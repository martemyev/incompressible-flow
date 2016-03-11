#include "config.hpp"
#include "mfem.hpp"
#include "param.hpp"
#include "pressure_solver.hpp"
#include "saturation_solver.hpp"
#include "utilities.hpp"

#include <cmath>
#include <cfloat>
#include <fstream>

using namespace std;
using namespace mfem;

#if defined(MFEM_USE_MPI) // parallel mode



void run_parallel(int argc, char **argv)
{
  int num_procs, myid;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if (argc == 1)
  {
    if (myid == 0)
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
    if (myid == 0)
      args.PrintUsage(cout);
    throw 1;
  }
  if (p.info)
  {
    if (myid == 0)
      cout << p.get_info() << endl;
    throw 1;
  }
  if (myid == 0)
    args.PrintOptions(cout);

  p.init_arrays();

  if (myid == 0)
  {
    string cmd = "mkdir -p " + string(p.outdir);
    MFEM_VERIFY(system(cmd.c_str()) == 0, "Failed to create an output directory");
  }

  ParMesh *pmesh = nullptr;

  {
    // Create a serial mesh, and initialize the parallel one
    const bool gen_edges = true;
    Mesh *mesh;
    if (p.spacedim == 2)
      mesh = new Mesh(p.nx, p.ny, Element::QUADRILATERAL, gen_edges, p.sx, p.sy);
    else if (p.spacedim == 3)
      mesh = new Mesh(p.nx, p.ny, p.nz, Element::HEXAHEDRON, gen_edges, p.sx, p.sy, p.sz);
    for (int el = 0; el < mesh->GetNE(); ++el)
      mesh->GetElement(el)->SetAttribute(el+1);

    const int dim = mesh->Dimension();
    MFEM_VERIFY(dim == p.spacedim, "Dimensions mismatch");

    pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
  }

  FiniteElementCollection *hdiv_coll = new RT_FECollection(p.order_v, p.spacedim); // for velocity
  FiniteElementCollection *l2_coll   = new L2_FECollection(p.order_p, p.spacedim); // for pressure
  DG_FECollection         *dg_coll   = new DG_FECollection(p.order_s, p.spacedim); // for saturation

  ParFiniteElementSpace V_space(pmesh, hdiv_coll);
  ParFiniteElementSpace P_space(pmesh, l2_coll);
  ParFiniteElementSpace S_space(pmesh, dg_coll);

  HYPRE_Int dimV = V_space.GlobalTrueVSize();
  HYPRE_Int dimP = P_space.GlobalTrueVSize();
  HYPRE_Int dimS = S_space.GlobalTrueVSize();

  if (myid == 0)
  {
    cout << "dim(V)  = " << dimV << "\n";
    cout << "dim(P)  = " << dimP << "\n";
    cout << "dim(V+P)= " << dimV + dimP << "\n";
    cout << "dim(S)  = " << dimS << "\n" << endl;
  }

  Array<int> block_offsets(3); // number of variables + 1
  block_offsets[0] = 0;
  block_offsets[1] = V_space.GetVSize();
  block_offsets[2] = P_space.GetVSize();
  block_offsets.PartialSum();

  Array<int> block_trueOffsets(3); // number of variables + 1
  block_trueOffsets[0] = 0;
  block_trueOffsets[1] = V_space.TrueVSize();
  block_trueOffsets[2] = P_space.TrueVSize();
  block_trueOffsets.PartialSum();

  BlockVector x(block_offsets); // solution of the mixed problem
  BlockVector trueX(block_trueOffsets);

  ParGridFunction V, P;
  V.Update(&V_space, x.GetBlock(0), 0);
  P.Update(&P_space, x.GetBlock(1), 0);

  ParGridFunction S(&S_space); // solution of the saturation problem
  S = 0.0;
  GridFunctionCoefficient saturation(&S);
  VectorGridFunctionCoefficient velocity(&V);
  GridFunctionCoefficient pressure(&P);

  const int n_cells = p.get_n_cells();

  Vector saturation_in_cells(n_cells);
  std::vector<int> saturation_flags(n_cells, 0);
  ValuesInCells S_vic(saturation, saturation_in_cells, saturation_flags, n_cells);

  Vector pressure_in_cells(n_cells);
  std::vector<int> pressure_flags(n_cells, 0);
  ValuesInCells P_vic(pressure, pressure_in_cells, pressure_flags, n_cells);

  if (p.saturation_steps > 0)
  {
    S.ProjectCoefficient(S_vic);
    const std::string tstr = d2s(0, 0, 0, 0, 6); // ti = 0
    output_scalar_cells_parallel(p, saturation_in_cells, saturation_flags, tstr,
                                 "saturation");
    saturation_flags.clear();
    saturation_flags.resize(n_cells, 0);
  }

  VisItDataCollection visit_global("inc-flow-parallel-global", pmesh);
  visit_global.SetPrefixPath(p.outdir);
  visit_global.RegisterField("pressure", &P);
  visit_global.RegisterField("velocity", &V);
  visit_global.RegisterField("saturation", &S);

  VisItDataCollection visit_local("inc-flow-parallel-local", pmesh);
  visit_local.SetPrefixPath(p.outdir);
  visit_local.RegisterField("saturation", &S);

  StopWatch global_time_loop;
  global_time_loop.Start();

  int nt = ceil(p.t_final / p.dt_global);
  if (myid == 0)
    cout << "Number of global time steps: " << nt << endl;

  int valid_loop = 1;
  for (int ti = 1; ti <= nt && valid_loop; ++ti)
  {
    ParPressureSolver(block_offsets, block_trueOffsets, *pmesh, p,
                      V_space, P_space, saturation, x, trueX);

    V.Distribute(&trueX.GetBlock(0));
    ParSaturationSolver(p, S, velocity, ti, p.dt_global, visit_local);

    if (p.vis_steps_global > 0 && ti % p.vis_steps_global == 0)
    {
      P.Distribute(&trueX.GetBlock(1));

      visit_global.SetCycle(ti);
      visit_global.SetTime(ti*p.dt_global);
      visit_global.Save();

      MPI_Barrier(MPI_COMM_WORLD);
      P.ProjectCoefficient(P_vic);
      const std::string tstr = d2s(ti, 0, 0, 0, 6);
      output_scalar_cells_parallel(p, pressure_in_cells, pressure_flags, tstr,
                                   "pressure");
      pressure_flags.clear();
      pressure_flags.resize(n_cells, 0);
    }

    if (p.saturation_steps > 0 && ti % p.saturation_steps == 0)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      S.ProjectCoefficient(S_vic);
      const std::string tstr = d2s(ti, 0, 0, 0, 6);
      output_scalar_cells_parallel(p, saturation_in_cells, saturation_flags, tstr,
                                   "saturation");
      saturation_flags.clear();
      saturation_flags.resize(n_cells, 0);
    }

    if (myid == 0)
      cout << "time step " << ti << " is done" << endl;

    int nNonfinite = S.CheckFinite(); // number of nonfinite values (NaN, Inf)
    MPI_Allreduce(MPI_IN_PLACE, &nNonfinite, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
    if (nNonfinite > 0)
    {
      valid_loop = 0;
      if (myid == 0)
      {
        cout << "\n\n\tSaturation has nonfinite values. The process stops.\n"
             << endl;
      }
    }
//    else
//    {
//      double Smax = S.Max();
//      MPI_Allreduce(MPI_IN_PLACE, &Smax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//      if (Smax > 1.1)
//      {
//        valid_loop = 0;
//        if (myid == 0)
//        {
//          cout << "\n\n\tSaturation went up to more than 1.1. The process stops.\n"
//               << endl;
//        }
//      }
//    }
  }

  if (myid == 0)
  {
    cout << "Time of the global time loop: " << global_time_loop.RealTime()
         << " sec" << endl;
  }

  delete hdiv_coll;
  delete l2_coll;
  delete dg_coll;

  delete pmesh;

  if (myid == 0)
    cout << "TOTAL TIME: " << total_time.RealTime() << " sec" << endl;
}
#endif // MFEM_USE_MPI

