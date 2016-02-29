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

  Vector rho_array(p.n_cells);
  Vector vp_array(p.n_cells);
  Vector vs_array(p.n_cells);

  const double Kframe = K_frame(K_MINERAL_MATRIX, K_FLUID_COMPONENT,
                                F_MINERAL_MATRIX, F_FLUID_COMPONENT);
  if (myid == 0)
    cout << "Kframe = " << Kframe << endl;

  double minKsat = DBL_MAX, maxKsat = DBL_MIN;
  double rhomm[2], vpmm[2], vsmm[2];

  // compute seismic properties for S_w = 0
  for (int i = 0; i < p.n_cells; ++i)
  {
    const double Kfl   = K_func(VP_O, VS_O, RHO_O); // K_fluid = K_oil
    const double Ksat  = K_sat(Kframe, K_MINERAL_MATRIX, Kfl, p.phi_array[i]);
    rho_array[i] = rho_B(RHO_O, RHO_GRAIN, p.phi_array[i]); // rho_fluid = rho_oil
    vp_array[i]  = vp_func(Ksat, G_MINERAL_MATRIX, rho_array[i]);
    vs_array[i]  = vs_func(G_MINERAL_MATRIX, rho_array[i]);
    minKsat = min(minKsat, Ksat);
    maxKsat = max(maxKsat, Ksat);
    get_minmax(rho_array, p.n_cells, rhomm[0], rhomm[1]);
    get_minmax(vp_array, p.n_cells, vpmm[0], vpmm[1]);
    get_minmax(vs_array, p.n_cells, vsmm[0], vsmm[1]);
  }

  if (myid == 0)
  {
    cout << "minKsat = " << minKsat << endl;
    cout << "maxKsat = " << maxKsat << endl;
    cout << "rhomin  = " << rhomm[0] << endl;
    cout << "rhomax  = " << rhomm[1] << endl;
    cout << "vpmin   = " << vpmm[0] << endl;
    cout << "vpmax   = " << vpmm[1] << endl;
    cout << "vsmin   = " << vsmm[0] << endl;
    cout << "vsmax   = " << vsmm[1] << endl;
  }

  if (myid == 0)
    output_seismic_properties(p, 0, rho_array, vp_array, vs_array);

  VisItDataCollection visit("inc-flow-parallel", pmesh);
  visit.RegisterField("pressure", &P);
  visit.RegisterField("velocity", &V);
  visit.RegisterField("saturation", &S);

  VisItDataCollection visit_local("inc-flow-parallel-local", pmesh);
  visit_local.RegisterField("saturation", &S);

  StopWatch global_time_loop;
  global_time_loop.Start();

  int nt = ceil(p.t_final / p.dt);
  if (myid == 0)
    cout << "Number of global time steps: " << nt << endl;

  for (int ti = 1; ti <= nt; ++ti)
  {
    ParPressureSolver(block_offsets, block_trueOffsets, *pmesh, p,
                      V_space, P_space, saturation, x, trueX);

    V.Distribute(&trueX.GetBlock(0));
    ParSaturationSolver(p, S, velocity, ti, p.dt, visit_local);

    if (p.vis_steps_global > 0 && ti % p.vis_steps_global == 0)
    {
      V.Distribute(&trueX.GetBlock(0));
      P.Distribute(&trueX.GetBlock(1));

//      const std::string time_str =  d2s(ti, 0, 0, 0, 6);
//      const std::string rank_str =  d2s(myid, 0, 0, 0, 6);
//      std::string mesh_name = "output/mesh_t" + time_str + "." + rank_str;
//      std::string V_name = "output/solV_t" + time_str + "." + rank_str;
//      std::string P_name = "output/solP_t" + time_str + "." + rank_str;

//      ofstream mesh_ofs(mesh_name.c_str());
//      mesh_ofs.precision(8);
//      pmesh->Print(mesh_ofs);

//      ofstream V_ofs(V_name.c_str());
//      V_ofs.precision(8);
//      V.Save(V_ofs);

//      ofstream P_ofs(P_name.c_str());
//      P_ofs.precision(8);
//      P.Save(P_ofs);
      visit.SetCycle(ti);
      visit.SetTime(ti*p.dt);
      visit.Save();
    }




//    if (p.seis_steps > 0 && ti % p.seis_steps == 0) // update seismic properties
//    {
//      Vector S_w(p.n_cells); // water saturation values in each cell
//      if (p.spacedim == 2)
//        compute_in_cells(p.sx, p.sy, p.nx, p.ny, *mesh, S, S_w);
//      else if (p.spacedim == 3)
//        compute_in_cells(p.sx, p.sy, p.sz, p.nx, p.ny, p.nz, *mesh, S, S_w);
//      else MFEM_ABORT("Not supported spacedim");

//      Gassmann(S_w, p, K_MINERAL_MATRIX, Kframe, RHO_GRAIN,
//               p.phi_array, rho_array, vp_array, vs_array);

//      get_minmax(rho_array, p.n_cells, rhomm[0], rhomm[1]);
//      get_minmax(vp_array, p.n_cells, vpmm[0], vpmm[1]);
//      get_minmax(vs_array, p.n_cells, vsmm[0], vsmm[1]);
//      if (verbose)
//      {
//        cout << "rhomin  = " << rhomm[0] << endl;
//        cout << "rhomax  = " << rhomm[1] << endl;
//        cout << "vpmin  = " << vpmm[0] << endl;
//        cout << "vpmax  = " << vpmm[1] << endl;
//        cout << "vsmin  = " << vsmm[0] << endl;
//        cout << "vsmax  = " << vsmm[1] << endl;
//      }

//      output_seismic_properties(p, ti, rho_array, vp_array, vs_array);
//    }

//    cout << "time step " << ti << " is done" << endl;

//    if (S.Max() > 1.1)
//    {
//      if (verbose)
//        cout << "Saturation went up to more than 1.1. The process stops." << endl;
//      break;
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

