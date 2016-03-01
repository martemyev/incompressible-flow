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
  Vector rho_array(n_cells);
  Vector vp_array(n_cells);
  Vector vs_array(n_cells);

  const double Kframe = K_frame(K_MINERAL_MATRIX, K_FLUID_COMPONENT,
                                F_MINERAL_MATRIX, F_FLUID_COMPONENT);
  cout << "Kframe = " << Kframe << endl;

  double minKsat = DBL_MAX, maxKsat = DBL_MIN;
  double rhomm[2], vpmm[2], vsmm[2];

  // compute seismic properties for S_w = 0
  for (int i = 0; i < n_cells; ++i)
  {
    const double Kfl   = K_func(VP_O, VS_O, RHO_O); // K_fluid = K_oil
    const double Ksat  = K_sat(Kframe, K_MINERAL_MATRIX, Kfl, p.phi_array[i]);
    rho_array[i] = rho_B(RHO_O, RHO_GRAIN, p.phi_array[i]); // rho_fluid = rho_oil
    vp_array[i]  = vp_func(Ksat, G_MINERAL_MATRIX, rho_array[i]);
    vs_array[i]  = vs_func(G_MINERAL_MATRIX, rho_array[i]);
    minKsat = min(minKsat, Ksat);
    maxKsat = max(maxKsat, Ksat);
    get_minmax(rho_array, n_cells, rhomm[0], rhomm[1]);
    get_minmax(vp_array, n_cells, vpmm[0], vpmm[1]);
    get_minmax(vs_array, n_cells, vsmm[0], vsmm[1]);
  }

  cout << "minKsat = " << minKsat << endl;
  cout << "maxKsat = " << maxKsat << endl;
  cout << "rhomin  = " << rhomm[0] << endl;
  cout << "rhomax  = " << rhomm[1] << endl;
  cout << "vpmin   = " << vpmm[0] << endl;
  cout << "vpmax   = " << vpmm[1] << endl;
  cout << "vsmin   = " << vsmm[0] << endl;
  cout << "vsmax   = " << vsmm[1] << endl;

  output_seismic_properties(p, 0, rho_array, vp_array, vs_array);

  VisItDataCollection visit("inc-flow-serial", mesh);
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
      Vector S_w(n_cells); // water saturation values in each cell
      if (p.spacedim == 2)
        compute_in_cells(p.sx, p.sy, p.nx, p.ny, *mesh, S, S_w);
      else if (p.spacedim == 3)
        compute_in_cells(p.sx, p.sy, p.sz, p.nx, p.ny, p.nz, *mesh, S, S_w);
      else MFEM_ABORT("Not supported spacedim");

      Gassmann(S_w, p, K_MINERAL_MATRIX, Kframe, RHO_GRAIN,
               p.phi_array, rho_array, vp_array, vs_array);

      get_minmax(rho_array, n_cells, rhomm[0], rhomm[1]);
      get_minmax(vp_array, n_cells, vpmm[0], vpmm[1]);
      get_minmax(vs_array, n_cells, vsmm[0], vsmm[1]);
      cout << "rhomin  = " << rhomm[0] << endl;
      cout << "rhomax  = " << rhomm[1] << endl;
      cout << "vpmin  = " << vpmm[0] << endl;
      cout << "vpmax  = " << vpmm[1] << endl;
      cout << "vsmin  = " << vsmm[0] << endl;
      cout << "vsmax  = " << vsmm[1] << endl;

      output_seismic_properties(p, ti, rho_array, vp_array, vs_array);
    }

    cout << "time step " << ti << " is done" << endl;

    if (S.Max() > 1.1)
    {
      cout << "Saturation went up to more than 1.1. The process stops." << endl;
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
