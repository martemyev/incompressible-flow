#include "config.hpp"
#include "mfem.hpp"
#include "param.hpp"
#include "pressure_solver.hpp"
#include "saturation_solver.hpp"
#include "utilities.hpp"

#include <cmath>
#include <cfloat>

using namespace std;
using namespace mfem;

void run_serial(int argc, char **argv)
{
  if (argc == 1)
  {
    cout << "\nTo see the available options, run\n" << argv[0] << " -h\n" << endl;
    throw 1;
  }

  StopWatch total_time;
  total_time.Start();

  Param p;
  OptionsParser args(argc, argv);
  p.add_options(args);
  args.Parse();
  if (!args.Good())
  {
     args.PrintUsage(cout);
     throw 1;
  }
  if (p.info)
  {
    cout << p.get_info() << endl;
    throw 1;
  }
  args.PrintOptions(cout);

  p.init_arrays();

  string cmd = "mkdir -p " + string(p.outdir);
  system(cmd.c_str());

  const int gen_edges = 1;
  Mesh *mesh;
  if (p.spacedim == 2)
    mesh = new Mesh(p.nx, p.ny, Element::QUADRILATERAL, gen_edges, p.sx, p.sy);
  else if (p.spacedim == 3)
    mesh = new Mesh(p.nx, p.ny, p.nz, Element::HEXAHEDRON, gen_edges, p.sx, p.sy, p.sz);
  for (int el = 0; el < mesh->GetNE(); ++el)
    mesh->GetElement(el)->SetAttribute(el+1);

  const int dim = mesh->Dimension();
  MFEM_VERIFY(dim == p.spacedim, "Dimensions mismatch");

  FiniteElementCollection *hdiv_coll = new RT_FECollection(p.order_v, dim); // for velocity
  FiniteElementCollection *l2_coll   = new L2_FECollection(p.order_p, dim); // for pressure
  DG_FECollection         *dg_coll   = new DG_FECollection(p.order_s, dim); // for saturation

  FiniteElementSpace V_space(mesh, hdiv_coll);
  FiniteElementSpace P_space(mesh, l2_coll);
  FiniteElementSpace S_space(mesh, dg_coll);

  Array<int> block_offsets(3);
  block_offsets[0] = 0;
  block_offsets[1] = V_space.GetVSize();
  block_offsets[2] = P_space.GetVSize();
  block_offsets.PartialSum();

  cout << "dim(V)  = " << block_offsets[1] - block_offsets[0] << "\n";
  cout << "dim(P)  = " << block_offsets[2] - block_offsets[1] << "\n";
  cout << "dim(V+P)= " << block_offsets.Last() << "\n";
  cout << "dim(S)  = " << S_space.GetVSize() << "\n" << endl;

  GridFunction V, P;

  BlockVector x(block_offsets); // solution of the mixed problem
  V.Update(&V_space, x.GetBlock(0), 0);
  P.Update(&P_space, x.GetBlock(1), 0);

  GridFunction S(&S_space); // solution of the saturation problem
  S = 0.0;
  GridFunctionCoefficient saturation(&S);
  VectorGridFunctionCoefficient velocity(&V);

  const int n_cells = p.get_n_cells();
  Vector saturation_in_cells(n_cells);
  std::vector<int> saturation_flags(n_cells, 0);
  ValuesInCells vic(saturation, saturation_in_cells, saturation_flags, n_cells);

  if (p.seis_steps > 0)
  {
    S.ProjectCoefficient(vic);
    const std::string tstr = d2s(0, 0, 0, 0, 6); // ti = 0
    output_scalar_cells_serial(p, saturation_in_cells, tstr, "saturation");
  }

  VisItDataCollection visit("inc-flow-serial", mesh, p.outdir);
  visit.RegisterField("pressure", &P);
  visit.RegisterField("saturation", &S);

  StopWatch global_time_loop;
  global_time_loop.Start();

  int nt = ceil(p.t_final / p.dt_global);
  cout << "Number of global time steps: " << nt << endl;
  for (int ti = 1; ti <= nt; ++ti)
  {
    PressureSolver(block_offsets, *mesh, p, V_space, P_space, saturation, x);
    SaturationSolver(p, S, velocity, ti, p.dt_global);

    if (p.vis_steps_global > 0 && ti % p.vis_steps_global == 0)
    {
      visit.SetCycle(ti);
      visit.SetTime(ti*p.dt_global);
      visit.Save();
    }

    if (p.seis_steps > 0 && ti % p.seis_steps == 0) // update seismic properties
    {
//      Vector S_w(n_cells); // water saturation values in each cell
//      if (p.spacedim == 2)
//        compute_in_cells(p.sx, p.sy, p.nx, p.ny, *mesh, S, S_w);
//      else if (p.spacedim == 3)
//        compute_in_cells(p.sx, p.sy, p.sz, p.nx, p.ny, p.nz, *mesh, S, S_w);
//      else MFEM_ABORT("Not supported spacedim");
      S.ProjectCoefficient(vic);
      const std::string tstr = d2s(ti, 0, 0, 0, 6);
      output_scalar_cells_serial(p, saturation_in_cells, tstr, "saturation");
    }

    cout << "time step " << ti << " is done" << endl;

    if (S.CheckFinite() > 0)
    {
      cout << "\n\n\tSaturation has nonfinite values. The process stops.\n" << endl;
      break;
    }

    if (S.Max() > 1.1)
    {
      cout << "\n\n\tSaturation went up to more than 1.1. The process stops.\n" << endl;
      break;
    }
  }

  cout << "Time of the global time loop: " << global_time_loop.RealTime()
       << " sec" << endl;

  delete hdiv_coll;
  delete l2_coll;
  delete dg_coll;

  delete mesh;

  cout << "TOTAL TIME: " << total_time.RealTime() << " sec" << endl;
}
