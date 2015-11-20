#include "param.hpp"
#include "saturation_solver.hpp"
#include "utilities.hpp"

using namespace std;
using namespace mfem;

void FS(const GridFunction& S, Vector& fs)
{
  MFEM_VERIFY(S.Size() == fs.Size(), "Dimensions mismatch");
  for (int i = 0; i < S.Size(); ++i)
  {
    const double mw = Krw(S(i)) / MU_W;
    const double mo = Kro(S(i)) / MU_O;
    fs(i) = mw / (mw + mo);
  }
}

void SaturationSolver(const Param &param, GridFunction &S,
                      VectorCoefficient &velocity, int global_ti,
                      double global_dt)
{
  FiniteElementSpace &S_space = *S.FESpace();

  const bool own_array = false;
  CWConstCoefficient Por(param.por_array, param.n_cells, own_array);

  BilinearForm m(&S_space);
  m.AddDomainIntegrator(new MassIntegrator(Por));
  BilinearForm k(&S_space);
  k.AddDomainIntegrator(new ConvectionIntegrator(velocity, -1.0));
  k.AddInteriorFaceIntegrator(
     new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, 0.5)));
  // the term for the boundary faces must be 0, because v \cdot n = 0 by
  // construction
//  k.AddBdrFaceIntegrator(
//     new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, -0.5)));

  CWConstCoefficient R(param.R_array, param.n_cells, own_array);
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

  Vector y(N), z(N);
#if defined(TWO_PHASE_FLOW)
  Vector fs(N);
#endif

//  int nt = 400;
  double dt = 1e-2; //global_dt/nt;
  int nt = global_dt / dt;
  for (int ti = 1; ti <= nt; ++ti)
  {
#if defined(TWO_PHASE_FLOW)
     FS(S, fs);
     (*Kt).Mult(fs, z); // S = S + dt M^-1( K F(S) + b)
#else
     (*Kt).Mult(S, z); // S = S + dt M^-1( K S + b)
#endif
     z += b;
     M_solver.Mult(z, y);
     y *= dt;
     S += y;

     if (ti % 200 == 0)
     {
        Vector S_nodal;
        S.GetNodalValues(S_nodal);
#if defined(TWO_PHASE_FLOW)
        string fname = "saturation_2phase_" + d2s(global_ti) + "local_" +
                       d2s(ti, 0, 0, 0, 6) + ".vts";
#else
        string fname = "saturation_1phase_" + d2s(global_ti) + "local_" +
                       d2s(ti, 0, 0, 0, 6) + ".vts";
#endif
        if (param.n_cells == param.nx*param.ny)
          write_vts_scalar(fname, "saturation", param.sx, param.sy,
                           param.nx, param.ny, S_nodal);
        else if (param.n_cells == param.nx*param.ny*param.nz)
          write_vts_scalar(fname, "saturation", param.sx, param.sy, param.sz,
                           param.nx, param.ny, param.nz, S_nodal);
     }
  }

  delete Kt;
}

