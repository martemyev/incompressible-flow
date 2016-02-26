#ifndef SATURATION_SOLVER_HPP
#define SATURATION_SOLVER_HPP

#include "config.hpp"
#include "mfem.hpp"

class Param;

void SaturationSolver(const Param &grid, mfem::GridFunction &S,
                      mfem::VectorCoefficient &velocity, int global_ti,
                      double global_dt);

#if defined(MFEM_USE_MPI) // parallel mode
void ParSaturationSolver(const Param &grid, mfem::GridFunction &S,
                         mfem::VectorCoefficient &velocity, int global_ti,
                         double global_dt);
#endif // MFEM_USE_MPI

#endif // SATURATION_SOLVER_HPP
