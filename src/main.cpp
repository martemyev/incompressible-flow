#include "config.hpp"
#include "mfem.hpp"
#include "param.hpp"
#include "pressure_solver.hpp"
#include "saturation_solver.hpp"
#include "utilities.hpp"

#include <cmath>

using namespace std;
using namespace mfem;

void output_scalar(const Param& p, const GridFunction& P, const string& tstr,
                   const string& name);

void output_vector(const Param& p, const GridFunction& V, const string& tstr,
                   const string& name);

void output_seismic_properties(const Param& p, int ti, const Vector& rho_array,
                               const Vector& vp_array, const Vector& vs_array);



int main(int argc, char **argv)
{
  StopWatch total_time;
  total_time.Start();

  Param p;
  OptionsParser args(argc, argv);
  p.add_options(args);
  args.Parse();
  if (!args.Good())
  {
     args.PrintUsage(cout);
     return 1;
  }
  args.PrintOptions(cout);

  p.init_arrays();

  const bool gen_edges = true;
  Mesh *mesh;
  if (p.spacedim == 2)
    mesh = new Mesh(p.nx, p.ny, Element::QUADRILATERAL, gen_edges, p.sx, p.sy);
  else if (p.spacedim == 3)
    mesh = new Mesh(p.nx, p.ny, p.nz, Element::QUADRILATERAL, gen_edges, p.sx, p.sy, p.sz);

  const int dim = mesh->Dimension();

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

  Vector rho_array(p.n_cells);
  Vector vp_array(p.n_cells);
  Vector vs_array(p.n_cells);
  rho_array = 2100.;
  vp_array  = 2400.;
  vs_array  = 1600.;
  output_seismic_properties(p, 0, rho_array, vp_array, vs_array);
  double K_saturated; // bulk modulus of the saturated media
  double Kframe = K_frame(K_MINERAL_MATRIX, K_FLUID_COMPONENT,
                          F_MINERAL_MATRIX, F_FLUID_COMPONENT);

  StopWatch global_time_loop;
  global_time_loop.Start();

  int nt = ceil(p.t_final / p.dt);
  cout << "Number of global time steps: " << nt << endl;
  for (int ti = 1; ti <= nt; ++ti)
  {
    const string tstr = d2s(ti, 0, 0, 0, 6);

    PressureSolver(block_offsets, *mesh, p, V_space, P_space, saturation, x);
    SaturationSolver(p, S, velocity, ti, p.dt);

    if (ti % p.vis_steps_global == 0) // visualize the solution
    {
      output_scalar(p, P, tstr, "pressure");
      output_vector(p, V, tstr, "velocity");
      output_scalar(p, S, tstr, "saturation");
    }

    if (ti % p.seis_steps == 0) // update seismic properties
    {
      Vector S_w(p.n_cells); // water saturation values in each cell
      if (p.spacedim == 2)
        compute_in_cells(p.sx, p.sy, p.nx, p.ny, *mesh, S, S_w);
      else if (p.spacedim == 3)
        compute_in_cells(p.sx, p.sy, p.sz, p.nx, p.ny, p.nz, *mesh, S, S_w);
      else MFEM_ABORT("Not supported spacedim");

      Gassmann(S, p, K_MINERAL_MATRIX, Kframe, RHO_GRAIN,
               p.phi_array, rho_array, vp_array, vs_array, K_saturated);

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

  delete mesh;

  delete hdiv_coll;
  delete l2_coll;
  delete dg_coll;

  cout << "TOTAL TIME: " << total_time.RealTime() << " sec" << endl;

  return 0;
}



void output_scalar(const Param& p, const GridFunction& P, const string& tstr,
                   const string& name)
{
#if defined(TWO_PHASE_FLOW)
  string fname = name + "_" + p.extra + "_" + d2s(p.spacedim) + "D_2phase_" +
                 tstr + ".vts";
#else
  string fname = name + "_" + p.extra + "_" + d2s(p.spacedim) + "D_1phase_" +
                 tstr + ".vts";
#endif
  Vector P_nodal;
  P.GetNodalValues(P_nodal);
  if (p.spacedim == 2)
    write_vts_scalar(fname, name, p.sx, p.sy, p.nx, p.ny, P_nodal);
  else if (p.spacedim == 3)
    write_vts_scalar(fname, name, p.sx, p.sy, p.sz, p.nx, p.ny, p.nz, P_nodal);
  else MFEM_ABORT("Not supported spacedim");
}



void output_vector(const Param& p, const GridFunction& V, const string& tstr,
                   const string& name)
{
#if defined(TWO_PHASE_FLOW)
  string fname = name + "_" + p.extra + "_" + d2s(p.spacedim) + "D_2phase_" +
                 tstr + ".vts";
#else
  string fname = name + "_" + p.extra + "_" + d2s(p.spacedim) + "D_1phase_" +
                 tstr + ".vts";
#endif
  Vector Vx_nodal, Vy_nodal;
  V.GetNodalValues(Vx_nodal, 1);
  V.GetNodalValues(Vy_nodal, 2);
  if (p.spacedim == 2)
    write_vts_vector(fname, name, p.sx, p.sy, p.nx, p.ny, Vx_nodal, Vy_nodal);
  else if (p.spacedim == 3)
  {
    Vector Vz_nodal;
    V.GetNodalValues(Vz_nodal, 3);
    write_vts_vector(fname, name, p.sx, p.sy, p.sz, p.nx, p.ny, p.nz,
                     Vx_nodal, Vy_nodal, Vz_nodal);
  }
  else MFEM_ABORT("Not supported spacedim");
}



void output_seismic_properties(const Param& p, int ti,
                               const Vector& rho_array, const Vector& vp_array,
                               const Vector& vs_array)
{
  string fname_rho_bin, fname_vp_bin, fname_vs_bin;
  string fname_rho_vts, fname_vp_vts, fname_vs_vts;
  string extra = (string)p.extra + "_";

#if defined(TWO_PHASE_FLOW)
  fname_rho_bin = extra + "2phase_" + d2s(p.n_cells) + "_t" + d2s(ti) + ".rho";
  fname_vp_bin  = extra + "2phase_" + d2s(p.n_cells) + "_t" + d2s(ti) + ".vp";
  fname_vs_bin  = extra + "2phase_" + d2s(p.n_cells) + "_t" + d2s(ti) + ".vs";
  fname_rho_vts = "rho_" + extra + d2s(p.spacedim) + "D_2phase_" + d2s(ti) + ".vts";
  fname_vp_vts  = "vp_" + extra + d2s(p.spacedim) + "D_2phase_"  + d2s(ti) + ".vts";
  fname_vs_vts  = "vs_" + extra + d2s(p.spacedim) + "D_2phase_"  + d2s(ti) + ".vts";
#else
  fname_rho_bin = extra + "1phase_" + d2s(p.n_cells) + "_t" + d2s(ti) + ".rho";
  fname_vp_bin  = extra + "1phase_" + d2s(p.n_cells) + "_t" + d2s(ti) + ".vp";
  fname_vs_bin  = extra + "1phase_" + d2s(p.n_cells) + "_t" + d2s(ti) + ".vs";
  fname_rho_vts = "rho_" + extra + d2s(p.spacedim) + "D_1phase_" + d2s(ti) + ".vts";
  fname_vp_vts  = "vp_" + extra + d2s(p.spacedim) + "D_1phase_"  + d2s(ti) + ".vts";
  fname_vs_vts  = "vs_" + extra + d2s(p.spacedim) + "D_1phase_"  + d2s(ti) + ".vts";
#endif

  write_binary(fname_rho_bin.c_str(), p.n_cells, rho_array.GetData());
  write_binary(fname_vp_bin.c_str(),  p.n_cells, vp_array.GetData());
  write_binary(fname_vs_bin.c_str(),  p.n_cells, vs_array.GetData());

  if (p.spacedim == 2)
  {
    write_vts_scalar_cells(fname_rho_vts, "density", p.sx, p.sy,
                           p.nx, p.ny, rho_array);
    write_vts_scalar_cells(fname_vp_vts, "vp", p.sx, p.sy,
                           p.nx, p.ny, vp_array);
    write_vts_scalar_cells(fname_vs_vts, "vs", p.sx, p.sy,
                           p.nx, p.ny, vs_array);
  }
  else if (p.spacedim == 3)
  {
    write_vts_scalar_cells(fname_rho_vts, "density", p.sx, p.sy, p.sz,
                           p.nx, p.ny, p.nz, rho_array);
    write_vts_scalar_cells(fname_vp_vts, "vp", p.sx, p.sy, p.sz,
                           p.nx, p.ny, p.nz, vp_array);
    write_vts_scalar_cells(fname_vs_vts, "vs", p.sx, p.sy, p.sz,
                           p.nx, p.ny, p.nz, vs_array);
  }
  else MFEM_ABORT("Not supported spacedim");
}
