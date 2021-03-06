#include "param.hpp"
#include "pressure_solver.hpp"

using namespace std;
using namespace mfem;

#if defined(MFEM_USE_MPI) // parallel mode
void ParPressureSolver(const Array<int>& block_offsets,
                       const Array<int>& block_trueOffsets,
                       const ParMesh& mesh,
                       const Param& param,
                       ParFiniteElementSpace& V_space,
                       ParFiniteElementSpace& P_space,
                       GridFunctionCoefficient &saturation,
                       BlockVector& x, BlockVector& trueX)
{
  StopWatch timer;
  timer.Start();

  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  const int n_cells = param.get_n_cells();

  const bool own_array = false;
  CWVectorCoefficient *K;
  if (param.spacedim == 2)
    K = new CWVectorCoefficient(saturation, MU_W, MU_O, param.K_array_x,
                                param.K_array_y, n_cells, param, own_array);
  else if (param.spacedim == 3)
    K = new CWVectorCoefficient(saturation, MU_W, MU_O, param.K_array_x,
                                param.K_array_y, param.K_array_z, n_cells,
                                param, own_array);
  else MFEM_ABORT("Unknown spacedim");

  WellFunctionCoefficient Q(param.injection, param.production,
                            param.inflow, param.outflow, param.spacedim);

  BlockVector rhs(block_offsets), trueRhs(block_trueOffsets);

  rhs = 0.0;
  trueRhs = 0.0;

  // ParLinearForm corresponding to the rhs.GetBlock(0) is zero
  // ParLinearForm corresponding to the trueRhs.GetBlock(0) is zero

  ParLinearForm gform;
  gform.Update(&P_space, rhs.GetBlock(1), 0);
  gform.AddDomainIntegrator(new DomainLFIntegrator(Q));
  gform.Assemble();
  gform *= -1.0;
  gform.ParallelAssemble(trueRhs.GetBlock(1));

  ParBilinearForm mVarf(&V_space);
  ParMixedBilinearForm bVarf(&V_space, &P_space);

  mVarf.AddDomainIntegrator(new VectorFEMassIntegrator(*K));
  mVarf.Assemble();

  bVarf.AddDomainIntegrator(new VectorFEDivergenceIntegrator);
  bVarf.Assemble();

  Array<int> bdr_attr_is_ess(mesh.bdr_attributes.Max());
  bdr_attr_is_ess = 1;
  Array<int> ess_dofs;
  V_space.GetEssentialVDofs(bdr_attr_is_ess, ess_dofs);
  mVarf.EliminateEssentialBCFromDofs(ess_dofs, x.GetBlock(0), rhs.GetBlock(0));
  bVarf.EliminateEssentialBCFromTrialDofs(ess_dofs, x.GetBlock(0), rhs.GetBlock(0));

  mVarf.Finalize();
  bVarf.Finalize();

  HypreParMatrix *M = mVarf.ParallelAssemble();
  HypreParMatrix *B = bVarf.ParallelAssemble();
  (*B) *= -1.0;
  HypreParMatrix *BT = B->Transpose();

  BlockOperator darcyOp(block_trueOffsets);
  darcyOp.SetBlock(0, 0, M);
  darcyOp.SetBlock(0, 1, BT);
  darcyOp.SetBlock(1, 0, B);

  HypreParMatrix *MinvBt = B->Transpose();
  HypreParVector Md(MPI_COMM_WORLD, M->GetGlobalNumRows(), M->GetRowStarts());
  M->GetDiag(Md);
  MinvBt->InvScaleRows(Md);
  HypreParMatrix *S = ParMult(B, MinvBt);

  HypreSolver *invM = new HypreDiagScale(*M);
  HypreSolver *invS = new HypreBoomerAMG(*S);

  invM->iterative_mode = false;
  invS->iterative_mode = false;

  BlockDiagonalPreconditioner darcyPrec(block_trueOffsets);
  darcyPrec.SetDiagonalBlock(0, invM);
  darcyPrec.SetDiagonalBlock(1, invS);

  MINRESSolver solver(MPI_COMM_WORLD);
  solver.SetMaxIter(param.darcy.maxiter);
  solver.SetRelTol(param.darcy.rtol);
  solver.SetAbsTol(param.darcy.atol);
  solver.SetOperator(darcyOp);
  solver.SetPreconditioner(darcyPrec);
  solver.SetPrintLevel(param.darcy.print_level);

  trueX = 0.0;
  solver.Mult(trueRhs, trueX);
  if (solver.GetConverged())
  {
    if (myid == 0)
      cout << "Darcy solver converged. Iter = " << solver.GetNumIterations()
           << ". Final norm = " << solver.GetFinalNorm()
           << ". Time = " << timer.RealTime() << " sec"
           << endl;
  }
  else
  {
    if (myid == 0)
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
  delete B;
  delete M;
  delete K;
}
#endif // MFEM_USE_MPI

