#include "param.hpp"
#include "saturation_solver.hpp"
#include "utilities.hpp"

using namespace std;
using namespace mfem;

void FS(const Param& param, const GridFunction& S, Vector& fs)
{
  MFEM_VERIFY(S.Size() == fs.Size(), "Dimensions mismatch");
  for (int i = 0; i < S.Size(); ++i)
  {
    const double mw = param.Krw(S(i)) / MU_W;
    const double mo = param.Kro(S(i)) / MU_O;
    fs(i) = mw / (mw + mo);
  }
}

void SaturationSolver(const Param &param, GridFunction &S,
                      VectorCoefficient &velocity, int global_ti,
                      double global_dt)
{
  FiniteElementSpace &S_space = *S.FESpace();

  const int n_cells = param.get_n_cells();

  const bool own_array = false;
  CWConstCoefficient Phi(param.phi_array, n_cells, own_array);

  BilinearForm m(&S_space);
  m.AddDomainIntegrator(new MassIntegrator(Phi));
  BilinearForm k(&S_space);
  k.AddDomainIntegrator(new ConvectionIntegrator(velocity, -1.0));
  k.AddInteriorFaceIntegrator(
     new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, 0.5)));
  // the term for the boundary faces must be 0, because v \cdot n = 0 by
  // construction
//  k.AddBdrFaceIntegrator(
//     new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, -0.5)));

  SaturationSourceCoefficient R(param.injection, param.saturation_source,
                                param.spacedim);
  LinearForm b(&S_space);
  b.AddDomainIntegrator(new DomainLFIntegrator(R));

  m.Assemble();
  m.Finalize();
  int skip_zeros = 0;
  k.Assemble(skip_zeros);
  k.Finalize(skip_zeros);
  b.Assemble();

  SparseMatrix& K = k.SpMat();
  SparseMatrix *Kt = Transpose(K);
  (*Kt) *= -1.0;

  const SparseMatrix& M = m.SpMat();
  const int N = S.Size();

  DSmoother M_prec;
  CGSolver M_solver;
  M_solver.SetPreconditioner(M_prec);
  M_solver.SetOperator(M);
  M_solver.iterative_mode = false;
  M_solver.SetRelTol(1e-9);
  M_solver.SetAbsTol(0.0);
  M_solver.SetMaxIter(100);
  M_solver.SetPrintLevel(0);

  Vector y(N), z(N), fs(N);

  string extra = (string)param.extra + "_";

  double dt = param.dt_local;
  int nt = global_dt / dt;
  for (int ti = 1; ti <= nt; ++ti)
  {
     if (param.two_phase_flow)
     {
       FS(param, S, fs);
       (*Kt).Mult(fs, z); // S = S + dt M^-1( K F(S) + b)
     }
     else
     {
       (*Kt).Mult(S, z); // S = S + dt M^-1( K S + b)
     }

     z += b;
     M_solver.Mult(z, y);
     y *= dt;
     S += y;

     if (param.vis_steps_local > 0 && ti % param.vis_steps_local == 0)
     {
        Vector S_nodal;
        S.GetNodalValues(S_nodal);
        string tstr = d2s(ti, 0, 0, 0, 6);
        string phase = (param.two_phase_flow ? "2phase" : "1phase");
        string fname = "saturation_" + extra + d2s(param.spacedim) + "D_" +
                       phase + "_glob" + d2s(global_ti) + "_loc_" + tstr + ".vts";
        if (param.spacedim == 2)
          write_vts_scalar(fname, "saturation", param.sx, param.sy,
                           param.nx, param.ny, S_nodal);
        else if (param.spacedim == 3)
          write_vts_scalar(fname, "saturation", param.sx, param.sy, param.sz,
                           param.nx, param.ny, param.nz, S_nodal);
        else MFEM_ABORT("Not supported spacedim");
     }
  }

  delete Kt;
}

