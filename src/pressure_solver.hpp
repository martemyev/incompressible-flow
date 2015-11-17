#ifndef PRESSURE_SOLVER_HPP
#define PRESSURE_SOLVER_HPP

#include "config.hpp"
#include "mfem.hpp"

class Grid;

void PressureSolver(const mfem::Array<int> &block_offsets,
                    const mfem::Mesh &mesh,
                    const Grid &grid,
                    mfem::FiniteElementSpace &V_space,
                    mfem::FiniteElementSpace &P_space,
                    mfem::GridFunctionCoefficient &S,
                    mfem::BlockVector &x);

#endif // PRESSURE_SOLVER_HPP
