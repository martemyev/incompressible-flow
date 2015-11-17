#ifndef SATURATION_SOLVER_HPP
#define SATURATION_SOLVER_HPP

#include "config.hpp"
#include "mfem.hpp"

class Grid;

void SaturationSolver(const Grid &grid, mfem::GridFunction &S,
                      mfem::VectorCoefficient &velocity, double global_dt);

#endif // SATURATION_SOLVER_HPP
