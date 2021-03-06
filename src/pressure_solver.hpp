#ifndef PRESSURE_SOLVER_HPP
#define PRESSURE_SOLVER_HPP

#include "config.hpp"
#include "mfem.hpp"

class Param;

void PressureSolver(const mfem::Array<int> &block_offsets,
                    const mfem::Mesh &mesh,
                    const Param &param,
                    mfem::FiniteElementSpace &V_space,
                    mfem::FiniteElementSpace &P_space,
                    mfem::GridFunctionCoefficient &saturation,
                    mfem::BlockVector &x);

#if defined(MFEM_USE_MPI) // parallel mode
void ParPressureSolver(const mfem::Array<int>& block_offsets,
                       const mfem::Array<int>& block_trueOffsets,
                       const mfem::ParMesh& mesh,
                       const Param& param,
                       mfem::ParFiniteElementSpace& V_space,
                       mfem::ParFiniteElementSpace& P_space,
                       mfem::GridFunctionCoefficient &saturation,
                       mfem::BlockVector& x, mfem::BlockVector& trueX);
#endif // MFEM_USE_MPI

#endif // PRESSURE_SOLVER_HPP
