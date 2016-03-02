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



class ValuesInCells : public Coefficient
{
public:
  ValuesInCells(Coefficient &coef, Vector &in_cells, std::vector<int> &flags,
                int ncells)
    : coefficient(coef), values_in_cells(in_cells), flags_of_cells(flags),
      n_cells(ncells)
  {}

  virtual ~ValuesInCells() {}

  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
  {
    const int index = T.Attribute - 1; // use attribute as a cell number
    MFEM_ASSERT(index >= 0 && index < n_cells, "Element number (attribute) is "
                "out of range: " + d2s(index));

    const double val = coefficient.Eval(T, ip);
    values_in_cells(index) = val;
    flags_of_cells[index] = 1;
    return val;
  }

private:
  Coefficient &coefficient;
  Vector &values_in_cells;
  std::vector<int> &flags_of_cells;
  int n_cells;
};



//class GassmannRho: public Coefficient
//{
//public:
//  GassmannRho(Coefficient &S, double *rho, const double *phi)
//    : saturation(S), rho_array(rho), phi_array(phi)
//  {}

//  virtual ~GassmannRho() {}

//  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
//  {
//    const int index = T.Attribute - 1; // use attribute as a cell number
//    MFEM_ASSERT(index >= 0 && index < n_cells, "Element number (attribute) is "
//                "out of range: " + d2s(index));

//    const double S   = saturation.Eval(T, ip);
//    const double phi = phi_array[index];

//    const double rho_fl_mix = rho_fl(S, RHO_W, RHO_O);

//    const double rho = rho_B(rho_fl_mix, RHO_GRAIN, phi);
//    rho_array[index] = rho;

//    return rho;
//  }

//private:
//  Coefficient& saturation;
//  double *rho_array;
//  const double *phi_array;
//};



//class GassmannVp: public Coefficient
//{
//public:
//  GassmannVp(Coefficient &S, const double *rho, double *vp,
//             const double *vs, const double *phi)
//    : saturation(S), rho_array(rho), vp_array(vp), vs_array(vs), phi_array(phi)
//  {}

//  virtual ~GassmannVp() {}

//  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
//  {
//    const int index = T.Attribute - 1; // use attribute as a cell number
//    MFEM_ASSERT(index >= 0 && index < n_cells, "Element number (attribute) is "
//                "out of range: " + d2s(index));

//    const double K_w = K_func(VP_W, VS_W, RHO_W); // bulk modulus of water
//    const double K_o = K_func(VP_O, VS_O, RHO_O); // bulk modulus of oil

//    const double S   = saturation.Eval(T, ip);
//    const double rho = rho_array[index];
//    const double vs  = vs_array[index];
//    const double phi = phi_array[index];

//    const double G = G_func(vs, rho);

//    const double K_fl_mix   = K_fl(S, K_w, K_o);

//    const double Kframe = K_frame(K_MINERAL_MATRIX, K_FLUID_COMPONENT,
//                                  F_MINERAL_MATRIX, F_FLUID_COMPONENT);
//    const double Ksat = K_sat(Kframe, K_MINERAL_MATRIX, K_fl_mix, phi);

//    const double vp = vp_func(Ksat, G, rho);
//    vp_array[index] = vp;

//    return vp;
//  }

//private:
//  Coefficient& saturation;
//  const double *rho_array;
//  double *vp_array;
//  const double *vs_array;
//  const double *phi_array;
//};



//class GassmannVs: public Coefficient
//{
//public:
//  GassmannVs(const double *rho, double *vs)
//    : rho_array(rho), vs_array(vs)
//  {}

//  virtual ~GassmannVs() {}

//  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
//  {
//    const int index = T.Attribute - 1; // use attribute as a cell number
//    MFEM_ASSERT(index >= 0 && index < n_cells, "Element number (attribute) is "
//                "out of range: " + d2s(index));

//    const double rho = rho_array[index];
//    const double vs_old = vs_array[index];
//    const double G = G_func(vs_old, rho);
//    const double vs_new = vs_func(G, rho);
//    vs_array[index] = vs_new;

//    return vs_new;
//  }

//private:
//  const double *rho_array;
//  double *vs_array;
//};




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
    system(cmd.c_str());
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

  const int n_cells = p.get_n_cells();
  Vector saturation_in_cells(n_cells);
  std::vector<int> saturation_flags(n_cells, 0);
  ValuesInCells vic(saturation, saturation_in_cells, saturation_flags, n_cells);

  const double Kframe = K_frame(K_MINERAL_MATRIX, K_FLUID_COMPONENT,
                                F_MINERAL_MATRIX, F_FLUID_COMPONENT);
  if (myid == 0)
    cout << "Kframe = " << Kframe << endl;

  if (p.seis_steps > 0)
  {
    S.ProjectCoefficient(vic);
    const std::string tstr = d2s(0, 0, 0, 0, 6); // ti = 0
    output_scalar_cells(p, saturation_in_cells, saturation_flags, tstr,
                        "saturation");
  }

  VisItDataCollection visit_global("inc-flow-parallel-global", pmesh, p.outdir);
  visit_global.RegisterField("pressure", &P);
  visit_global.RegisterField("velocity", &V);
  visit_global.RegisterField("saturation", &S);

  VisItDataCollection visit_local("inc-flow-parallel-local", pmesh, p.outdir);
  visit_local.RegisterField("saturation", &S);

  StopWatch global_time_loop;
  global_time_loop.Start();

  int nt = ceil(p.t_final / p.dt_global);
  if (myid == 0)
    cout << "Number of global time steps: " << nt << endl;

  for (int ti = 1; ti <= nt; ++ti)
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
    }

    // update seismic properties
    if (p.seis_steps > 0 && ti % p.seis_steps == 0)
    {
      S.ProjectCoefficient(vic);
      const std::string tstr = d2s(ti, 0, 0, 0, 6);
      output_scalar_cells(p, saturation_in_cells, saturation_flags, tstr,
                          "saturation");
    }

    cout << "time step " << ti << " is done" << endl;

    if (S.Max() > 1.1)
    {
      cout << "Saturation went up to more than 1.1. The process stops." << endl;
      break;
    }
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

