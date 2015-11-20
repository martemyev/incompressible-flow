#include "param.hpp"
#include "pressure_solver.hpp"

using namespace std;
using namespace mfem;

void PressureSolver(const Array<int> &block_offsets,
                    const Mesh& mesh,
                    const Param& param,
                    FiniteElementSpace &V_space,
                    FiniteElementSpace &P_space,
                    GridFunctionCoefficient &saturation,
                    BlockVector &x)
{
  StopWatch timer;
  timer.Start();

  const bool own_array = false;
  CWCoefficient K(saturation, MU_W, MU_O, param.K_array, param.n_cells, own_array);
  CWConstCoefficient Q(param.Q_array, param.n_cells, own_array);

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

  mVarf.AddDomainIntegrator(new VectorFEMassIntegrator(K));
  mVarf.Assemble();

  bVarf.AddDomainIntegrator(new VectorFEDivergenceIntegrator); //(K));
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

  MINRESSolver solver;
  solver.SetMaxIter(param.darcy.maxiter);
  solver.SetRelTol(param.darcy.rtol);
  solver.SetAbsTol(param.darcy.atol);
  solver.SetOperator(darcyMatrix);
  solver.SetPreconditioner(darcyPrec);
  solver.SetPrintLevel(param.darcy.print_level);

  x = 0.0;
  solver.Mult(rhs, x);
  if (solver.GetConverged())
  {
    cout << "Darcy solver converged. Iter = " << solver.GetNumIterations()
         << ". Final norm = " << solver.GetFinalNorm()
         << ". Time = " << timer.RealTime() << " sec"
         << endl;
  }
  else
  {
    cerr << "Darcy solver DIDN'T converge. Iter = " << solver.GetNumIterations()
         << ". Final norm = " << solver.GetFinalNorm()
         << ". Time = " << timer.RealTime() << " sec"
         << endl;
    MFEM_ABORT("The solver didn't converge => terminating");
  }

  delete invS;
  delete invM;
  delete S;
  delete MinvBt;
  delete BT;
}
