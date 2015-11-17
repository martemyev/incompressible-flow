#include "grid.hpp"
#include "saturation_solver.hpp"
#include "utilities.hpp"

using namespace std;
using namespace mfem;

void FS(const Vector& S, double mu_w, double mu_o, Vector& fs)
{
  for (int i = 0; i < S.Size(); ++i)
  {
    double s = S(i);
    double mw = CWCoefficient::Krw(s) / mu_w;
    double mo = CWCoefficient::Kro(s) / mu_o;
    fs(i) = mw / (mw+mo);
  }
}

void SaturationSolver(const Grid &grid, GridFunction &S,
                      VectorCoefficient &velocity, double global_dt)
{
  FiniteElementSpace &S_space = *S.FESpace();
  const int N = S.Size();
//  Vector S_nodal;

  BilinearForm m(&S_space);
  m.AddDomainIntegrator(new MassIntegrator);
  m.Assemble();
  m.Finalize();
  const SparseMatrix& M = m.SpMat();
  DSmoother M_prec(M);

  BilinearForm k(&S_space);
  k.AddDomainIntegrator(new ConvectionIntegrator(velocity, -1.0));
  k.AddInteriorFaceIntegrator(
        new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, -0.5)));
  k.AddBdrFaceIntegrator(
        new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, -0.5)));
  const int skip_zeros = 0;
  k.Assemble(skip_zeros);
  k.Finalize(skip_zeros);
  SparseMatrix& K = k.SpMat();
//  K *= -1.;

  const bool own_array = false;
  CWConstCoefficient r_coef(grid.r_array, own_array);

  LinearForm b(&S_space);
  b.AddDomainIntegrator(new DomainLFIntegrator(r_coef));
  b.Assemble();

  const int nt = 350;
  double dt = global_dt/nt;
  for (int t = 0; t < nt; ++t)
  {
    cout << "time step " << t+1 << " / " << nt << endl;

     Vector fs(N);
     FS(S, MU_W, MU_O, fs);

     Vector y(N);
     K.MultTranspose(fs, y); // y = K S_n

     Vector F = b;
     F -= y; // F = r_n - K S_n

     // z = M^-1 F
     Vector z(N);
     PCG(M, M_prec, F, z, 0, 100, 1e-9, 0.0);

     z *= dt;
     S += z; // u_n+1 = u_n + dt M^-1 (r_n - K S_n)

//     if (t % 200 == 0)
//     {
//       const string tstr = d2s(t, 0, 0, 0, 6);
//       string fname = "local_saturation_" + tstr + ".vts";
//       S.GetNodalValues(S_nodal);
//       write_vts_scalar(fname, "saturation", grid.sx, grid.sy, grid.nx, grid.ny,
//                        S_nodal);
//     }
  }
}

