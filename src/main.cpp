#include "config.hpp"
#include "mfem.hpp"
#include "param.hpp"
#include "pressure_solver.hpp"
#include "saturation_solver.hpp"
#include "utilities.hpp"

#include <cmath>

using namespace std;
using namespace mfem;



int main(int argc, char **argv)
{
  int spacedim = 2;
  int nx = 64, ny = 64, nz = 64;
  double sx = 1.0, sy = 1.0, sz = 1.0;
  int order_v = 1, order_p = 1, order_s = 0;
  double t_final = 3000;
  double dt = 100;
  int nt = ceil(t_final / dt);
  bool visualization = true;
  int vis_steps = 100;

  OptionsParser args(argc, argv);
  args.AddOption(&spacedim, "-d", "--dim", "Space dimension of the problem");
  args.AddOption(&nx, "-nx", "--nx", "Number of cells in x-direction");
  args.AddOption(&ny, "-ny", "--ny", "Number of cells in y-direction");
  args.AddOption(&nz, "-nz", "--nz", "Number of cells in z-direction");
  args.AddOption(&sx, "-sx", "--sx", "Size of domain in x-direction");
  args.AddOption(&sy, "-sy", "--sy", "Size of domain in y-direction");
  args.AddOption(&sz, "-sz", "--sz", "Size of domain in z-direction");
  args.AddOption(&order_v, "-ov", "--orderv", "Order (degree) of the finite elements for velocity");
  args.AddOption(&order_p, "-op", "--orderp", "Order (degree) of the finite elements for pressure");
  args.AddOption(&order_s, "-os", "--orders", "Order (degree) of the finite elements for saturation");
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step", "Time step.");
  args.AddOption(&visualization, "-vis", "--vis", "-no-vis", "--no-vis", "Enable or disable GLVis visualization.");
  args.AddOption(&vis_steps, "-vs", "--vis-steps", "Visualize every n-th timestep.");
  args.Parse();
  if (!args.Good())
  {
     args.PrintUsage(cout);
     return 1;
  }
  args.PrintOptions(cout);

  Param *param;
  if (spacedim == 2)
    param = new Param(nx, ny, sx, sz);
  else if (spacedim == 3)
    param = new Param(nx, ny, nz, sx, sy, sz);
  else MFEM_ABORT("Unsupported spacedim: " + d2s(spacedim));

  param->init_arrays();

  const bool generate_edges = true;
  Mesh *mesh;
  if (spacedim == 2)
    mesh = new Mesh(param->nx, param->ny, Element::QUADRILATERAL, generate_edges,
                    param->sx, param->sy);
  else if (spacedim == 3)
    mesh = new Mesh(param->nx, param->ny, param->nz, Element::QUADRILATERAL,
                    generate_edges, param->sx, param->sy, param->sz);

  const int dim = mesh->Dimension();

  FiniteElementCollection *hdiv_coll = new RT_FECollection(order_v, dim); // for velocity
  FiniteElementCollection *l2_coll   = new L2_FECollection(order_p, dim); // for pressure
  DG_FECollection         *dg_coll   = new DG_FECollection(order_s, dim); // for saturation

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

  Vector P_nodal, S_nodal, Vx_nodal, Vy_nodal, Vz_nodal;

  string fname;
  for (int ti = 1; ti <= nt; ++ti)
  {
    const string tstr = d2s(ti, 0, 0, 0, 6);

    PressureSolver(block_offsets, *mesh, *param, V_space, P_space, saturation, x);

    {
#if defined(TWO_PHASE_FLOW)
      fname = "pressure_2phase_" + tstr + ".vts";
#else
      fname = "pressure_1phase_" + tstr + ".vts";
#endif
      P.GetNodalValues(P_nodal);
      if (spacedim == 2)
        write_vts_scalar(fname, "pressure", param->sx, param->sy,
                         param->nx, param->ny, P_nodal);
      else if (spacedim == 3)
        write_vts_scalar(fname, "pressure", param->sx, param->sy, param->sz,
                         param->nx, param->ny, param->nz, P_nodal);
    }

    {
#if defined(TWO_PHASE_FLOW)
      fname = "velocity_2phase_" + tstr + ".vts";
#else
      fname = "velocity_1phase_" + tstr + ".vts";
#endif
      V.GetNodalValues(Vx_nodal, 1);
      V.GetNodalValues(Vy_nodal, 2);
      if (spacedim == 2)
        write_vts_vector(fname, "velocity", param->sx, param->sy,
                         param->nx, param->ny, Vx_nodal, Vy_nodal);
      else if (spacedim == 3)
      {
        V.GetNodalValues(Vz_nodal, 3);
        write_vts_vector(fname, "velocity", param->sx, param->sy, param->sz,
                         param->nx, param->ny, param->nz,
                         Vx_nodal, Vy_nodal, Vz_nodal);
      }
    }

    SaturationSolver(*param, S, velocity, ti, dt);

    {
#if defined(TWO_PHASE_FLOW)
      fname = "saturation_2phase_" + tstr + ".vts";
#else
      fname = "saturation_1phase_" + tstr + ".vts";
#endif
      S.GetNodalValues(S_nodal);
      if (spacedim == 2)
        write_vts_scalar(fname, "saturation", param->sx, param->sy,
                         param->nx, param->ny, S_nodal);
      else if (spacedim == 3)
        write_vts_scalar(fname, "saturation", param->sx, param->sy, param->sz,
                         param->nx, param->ny, param->nz, S_nodal);
    }

    cout << "time step " << ti << " is done" << endl;

    if (S.Max() > 1.1)
    {
      cout << "Saturation went up to more than 1.1. The process stops." << endl;
      break;
    }

//    Gassmann(S,
  }


  delete mesh;
  delete param;

  delete hdiv_coll;
  delete l2_coll;
  delete dg_coll;

  return 0;
}
