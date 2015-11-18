#include "grid.hpp"
#include "pressure_solver.hpp"

using namespace std;
using namespace mfem;

double q_func(Vector& x)
{
  const double rad = 0.02;
  if (x(0) < rad && x(1) < rad) return 1.0; // injection
  if (x(0) > (1.-rad) && x(1) > (1.-rad)) return -1.0; // production
  return 0.0;
}

void PressureSolver(const Array<int> &block_offsets,
                    const Mesh& mesh,
                    const Grid &grid,
                    FiniteElementSpace &V_space,
                    FiniteElementSpace &P_space,
                    GridFunctionCoefficient &saturation,
                    BlockVector &x)
{
  const bool own_array = false;
//  CWCoefficient K(saturation, MU_W, MU_O, grid.K_array, grid.n_cells, own_array);
//  CWConstCoefficient Q(grid.Q_array, grid.n_cells, own_array);
  FunctionCoefficient Q(q_func);

  BlockVector rhs(block_offsets);

  rhs = 0.0;

  // LinearForm corresponding to the rhs.GetBlock(0) is zero

  LinearForm gform;
  gform.Update(&P_space, rhs.GetBlock(1), 0);
  gform.AddDomainIntegrator(new DomainLFIntegrator(Q));
  gform.Assemble();
  gform *= -1.0;

  BilinearForm mVarf(&V_space);
  MixedBilinearForm bVarf(&V_space, &P_space);

  mVarf.AddDomainIntegrator(new VectorFEMassIntegrator); // (K));
  mVarf.Assemble();

  bVarf.AddDomainIntegrator(new VectorFEDivergenceIntegrator()); //(K));
  bVarf.Assemble();

  Array<int> bdr_attr_is_ess(mesh.bdr_attributes.Max());
  bdr_attr_is_ess = 1;
  Array<int> ess_dofs;
  V_space.GetEssentialVDofs(bdr_attr_is_ess, ess_dofs);
  mVarf.EliminateEssentialBCFromDofs(ess_dofs, x.GetBlock(0), rhs.GetBlock(0));
  bVarf.EliminateEssentialBCFromTrialDofs(ess_dofs, x.GetBlock(0), rhs.GetBlock(0));

  mVarf.Finalize();
  bVarf.Finalize();

  SparseMatrix& M = mVarf.SpMat();
  SparseMatrix& B = bVarf.SpMat();
  B *= -1.0;
  SparseMatrix *BT = Transpose(B);

  BlockMatrix darcyMatrix(block_offsets);
  darcyMatrix.SetBlock(0, 0, &M);
  darcyMatrix.SetBlock(0, 1, BT);
  darcyMatrix.SetBlock(1, 0, &B);

  SparseMatrix *MinvBt = Transpose(B);
  Vector Md(M.Height());
  M.GetDiag(Md);
  for (int i = 0; i < Md.Size(); ++i)
    MinvBt->ScaleRow(i, 1.0/Md(i));
  SparseMatrix *S = Mult(B, *MinvBt);

  Solver *invM = new DSmoother(M);
  Solver *invS = new GSSmoother(*S);

  invM->iterative_mode = false;
  invS->iterative_mode = false;

  BlockDiagonalPreconditioner darcyPrec(block_offsets);
  darcyPrec.SetDiagonalBlock(0, invM);
  darcyPrec.SetDiagonalBlock(1, invS);

  const int maxiter = 1000;
  const double rtol = 1e-6;
  const double atol = 1e-10;

  MINRESSolver solver;
  solver.SetMaxIter(maxiter);
  solver.SetRelTol(rtol);
  solver.SetAbsTol(atol);
  solver.SetOperator(darcyMatrix);
  solver.SetPreconditioner(darcyPrec);
  solver.SetPrintLevel(0);

  x = 0.0;
  solver.Mult(rhs, x); cout << flush;
  MFEM_VERIFY(solver.GetConverged(), "The MINRES solver didn't converge");

  delete BT;
}
