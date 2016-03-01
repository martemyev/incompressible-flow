#ifndef PARAM_HPP
#define PARAM_HPP

#include "config.hpp"
#include "mfem.hpp"
#include "utilities.hpp"

#include <iostream>

double Krw(double S, bool two_phase_flow);
double Kro(double S, bool two_phase_flow);


struct DarcySolverParam
{
  int maxiter;
  double rtol;
  double atol;
  int print_level;

  DarcySolverParam();
  void add_options(mfem::OptionsParser& args);
};


struct Param
{
  int spacedim; ///< Spatial dimension of the problem (2 for 2D, 3 for 3D)
  int nx, ny, nz; ///< Number of cells of a Cartesian grid along each direction
  double sx, sy, sz; ///< Size of computational domain in each direction
  double *K_array; ///< Permeability array (cell-wise constant - one value per cell)
  double *Q_array; ///< Pressure source defined as an array (injection and
                   ///< production wells - inflow and outflow)
  double *phi_array; ///< Array of porosity values (cell-wise constant)
  double *R_array; ///< Saturation source defined as an array (cell-wise constant)

  int order_v, order_p, order_s; ///< Orders of finite elements for approximation
                                 ///< of (v)elocity, (p)ressure and (s)aturation
  double t_final; ///< Final time of simulation
  double dt_global; ///< Time step for pressure (global time loop)
  double dt_local; ///< Time step for saturation (local time loop)
  int vis_steps_global; ///< Visualization step for pressure solver (global)
  int vis_steps_local; ///< Visualization step for saturation solver (local)
  int seis_steps; ///< Output step for seismic media properties (based on Gassmann)

  double K; ///< Constant permeability
  double phi; ///< Constant porosity

  const char *K_file; ///< Binary file for heterogeneous permeability
  const char *phi_file; ///< Binary file for heterogeneous porosity

  const char *outdir; ///< Name of an output directory
  const char *extra; ///< Extra string to distinguish the output files

  bool two_phase_flow; ///< Simulate two phase flow (otherwise it's single phase)

  DarcySolverParam darcy; ///< Parameters of a Darcy solver

  bool info; ///< Whether or not print an info about the program
  std::string get_info() const;

  Param();
  ~Param();
  void init_arrays();
  void add_options(mfem::OptionsParser& args);

  /**
   * Get total number of cells
   */
  int get_n_cells() const
  {
    int n_cells = nx*ny;
    if (spacedim == 3) n_cells *= nz;
    return n_cells;
  }

  /**
   * Volume (area) of a one cell (all cells are of the same size)
   */
  double get_cell_volume() const
  {
    double V_domain = sx*sy;
    if (spacedim == 3) V_domain *= sz;
    return V_domain / get_n_cells();
  }
};



/**
 * Cell-wise constant coefficient.
 */
class CWConstCoefficient : public mfem::Coefficient
{
public:
  CWConstCoefficient(double *array, int ncells, bool own = 1)
    : val_array(array), n_cells(ncells), own_array(own)
  { }

  virtual ~CWConstCoefficient() { if (own_array) delete[] val_array; }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint&)
  {
    const int index = T.Attribute - 1; // use attribute as a cell number
    MFEM_ASSERT(index >= 0 && index < n_cells, "Element number (attribute) is "
                "out of range: " + d2s(index));
    return val_array[index];
  }

protected:
  double *val_array;
  int n_cells;
  bool own_array;
};



/**
 * Coefficient based on a product of a cell-wise coefficient (done via array),
 * and a grid function.
 */
class CWCoefficient : public mfem::Coefficient
{
public:
  CWCoefficient(mfem::GridFunctionCoefficient& f, double muw, double muo,
                double *array, int ncells, bool two_phase, bool own)
    : func(f)
    , mu_w(muw)
    , mu_o(muo)
    , val_array(array)
    , n_cells(ncells)
    , two_phase_flow(two_phase)
    , own_array(own)
  { }

  virtual ~CWCoefficient() { if (own_array) delete[] val_array; }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip)
  {
    const int index = T.Attribute - 1; // use attribute as a cell number
    MFEM_ASSERT(index >= 0 && index < n_cells, "Element number (attribute) is "
                "out of range: " + d2s(index));
    const double S = func.Eval(T, ip);
    const double val = Krw(S, two_phase_flow) / mu_w +
                       Kro(S, two_phase_flow) / mu_o;
    return 1./(val*val_array[index]);
  }

protected:
  mfem::GridFunctionCoefficient& func;
  double mu_w, mu_o;
  double *val_array;
  int n_cells;
  bool two_phase_flow;
  bool own_array;
};

#endif // PARAM_HPP
