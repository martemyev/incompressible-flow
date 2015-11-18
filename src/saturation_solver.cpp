#include "grid.hpp"
#include "saturation_solver.hpp"
#include "utilities.hpp"

using namespace std;
using namespace mfem;

void FS(const GridFunction& S, Vector& fs)
{
  MFEM_VERIFY(S.Size() == fs.Size(), "Sizes mismatch");
  for (int i = 0; i < S.Size(); ++i)
  {
    double s = S(i);
    double mw = Krw(s) / MU_W;
    double mo = Kro(s) / MU_O;
    fs(i) = mw / (mw+mo);
  }
}



void SaturationSolver(const Grid &grid, GridFunction &S,
                      VectorCoefficient &velocity, double global_dt)
{
  FiniteElementSpace &S_space = *S.FESpace();

  BilinearForm m(&S_space);
  m.AddDomainIntegrator(new MassIntegrator);
  BilinearForm k(&S_space);
//  k.AddDomainIntegrator(new ConvectionIntegrator(velocity, -1.0));
  k.AddInteriorFaceIntegrator(
     new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, -0.5)));
//  k.AddBdrFaceIntegrator(
//     new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, -0.5)));

  const bool own_array = false;
  CWConstCoefficient r_coef(grid.r_array, grid.n_cells, own_array);
  LinearForm b(&S_space);
  b.AddDomainIntegrator(new DomainLFIntegrator(r_coef));

  m.Assemble();
  m.Finalize();
  int skip_zeros = 0;
  k.Assemble(skip_zeros);
  k.Finalize(skip_zeros);
  b.Assemble();

  SparseMatrix& K = k.SpMat();  K *= -1.0;
  SparseMatrix *Kt = Transpose(K);
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

  double t_final = 100.;
  double dt = 1e-2; //global_dt/nt;
  const int nt = t_final/dt;
  double t = 0.0;
  for (int ti = 0; ti < nt; ++ti, t += dt)
  {
     // y = M^{-1} (K x + b)
//     FS(S, fs);
//     K.Mult(fs, z);
     (*Kt).Mult(S, z); // S = S + dt M^-1( K S + b)
     z += b;
     M_solver.Mult(z, y);
     y *= dt;
     S += y;

     if (ti % 200 == 0)
     {
        cout << "time step " << ti << " / " << nt << ", time: " << t << endl;
        Vector S_nodal;
        S.GetNodalValues(S_nodal);
        string fname = "saturation_" + d2s(ti, 0, 0, 0, 6) + ".vts";
        write_vts_scalar(fname, "saturation", grid.sx, grid.sy, grid.nx,
                         grid.ny, S_nodal);
     }
  }
}

