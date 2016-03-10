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



/**
 * Wells are respresented by cylinders (circles in 2D).
 */
struct Well
{
  mfem::Vertex center;
  double radius;
  double height;
  std::string option_prefix; // to distinguish different wells via options

  Well();
};


struct Param
{
  int spacedim; ///< Spatial dimension of the problem (2 for 2D, 3 for 3D)
  int nx, ny, nz; ///< Number of cells of a Cartesian grid along each direction
  double sx, sy, sz; ///< Size of computational domain in each direction
  double *K_array_x; ///< Permeability array (cell-wise constant - one value per cell) (x-component of a diagonal tensor)
  double *K_array_y; ///< Permeability array (cell-wise constant - one value per cell) (y-component of a diagonal tensor)
  double *K_array_z; ///< Permeability array (cell-wise constant - one value per cell) (z-component of a diagonal tensor)
  double *phi_array; ///< Array of porosity values (cell-wise constant)

  int order_v, order_p, order_s; ///< Orders of finite elements for approximation
                                 ///< of (v)elocity, (p)ressure and (s)aturation
  double t_final; ///< Final time of simulation
  double dt_global; ///< Time step for pressure (global time loop)
  double dt_local; ///< Time step for saturation (local time loop)
  int vis_steps_global; ///< Visualization step for pressure solver (global)
  int vis_steps_local; ///< Visualization step for saturation solver (local)
  int saturation_steps; ///< Output step for saturation as a binary file

  double K_x; ///< Constant permeability (x-component of a tensor)
  double K_y; ///< Constant permeability (y-component of a tensor)
  double K_z; ///< Constant permeability (z-component of a tensor)
  double phi; ///< Constant porosity

  const char *K_file_x; ///< Binary file for heterogeneous permeability (x-component of a diagonal tensor)
  const char *K_file_y; ///< Binary file for heterogeneous permeability (y-component of a diagonal tensor)
  const char *K_file_z; ///< Binary file for heterogeneous permeability (z-component of a diagonal tensor)
  const char *phi_file; ///< Binary file for heterogeneous porosity

  const char *outdir; ///< Name of an output directory
  const char *extra; ///< Extra string to distinguish the output files

  bool two_phase_flow; ///< Simulate two phase flow (otherwise it's single phase)

  int ode_solver_type; ///< Type of the ODE solver for the saturation

  Well injection, production;
  double inflow, outflow, saturation_source;

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



/**
 * Coefficient based on a product of a cell-wise vector coefficient (done via
 * arrays - in (x and y) or (x, y and z0 directions depending on dimension of
 * a problem), and a grid function.
 */
class CWVectorCoefficient : public mfem::VectorCoefficient
{
public:
  CWVectorCoefficient(mfem::GridFunctionCoefficient& f, double muw, double muo,
                      double *arrayX, double *arrayY, int ncells,
                      bool two_phase, bool own)
    : VectorCoefficient(2)
    , func(f)
    , mu_w(muw)
    , mu_o(muo)
    , val_array_X(arrayX)
    , val_array_Y(arrayY)
    , val_array_Z(nullptr)
    , n_cells(ncells)
    , two_phase_flow(two_phase)
    , own_arrays(own)
  {}

  CWVectorCoefficient(mfem::GridFunctionCoefficient& f, double muw, double muo,
                      double *arrayX, double *arrayY, double *arrayZ, int ncells,
                      bool two_phase, bool own)
    : VectorCoefficient(3)
    , func(f)
    , mu_w(muw)
    , mu_o(muo)
    , val_array_X(arrayX)
    , val_array_Y(arrayY)
    , val_array_Z(arrayZ)
    , n_cells(ncells)
    , two_phase_flow(two_phase)
    , own_arrays(own)
  {}

  virtual ~CWVectorCoefficient()
  {
    if (own_arrays)
    {
      delete[] val_array_X;
      delete[] val_array_Y;
      delete[] val_array_Z;
    }
  }

  virtual void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
                    const mfem::IntegrationPoint &ip)
  {
    const int index = T.Attribute - 1; // use attribute as a cell number
    MFEM_ASSERT(index >= 0 && index < n_cells, "Element number (attribute) is "
                "out of range: " + d2s(index));
    const double S = func.Eval(T, ip);
    const double val = Krw(S, two_phase_flow) / mu_w +
                       Kro(S, two_phase_flow) / mu_o;

    V.SetSize(vdim);
    V(0) = 1./(val*val_array_X[index]);
    V(1) = 1./(val*val_array_Y[index]);
    if (vdim == 3)
      V(2) = 1./(val*val_array_Z[index]);
  }

protected:
  mfem::GridFunctionCoefficient& func;
  double mu_w, mu_o;
  double *val_array_X, *val_array_Y, *val_array_Z;
  int n_cells;
  bool two_phase_flow;
  bool own_arrays;
};



class WellFunctionCoefficient : public mfem::Coefficient
{
public:
  WellFunctionCoefficient(const Well& inject, const Well& product, double in,
                          double out, int dim)
    : injection(inject), production(product), inflow(in), outflow(out),
      spacedim(dim)
  {}

  virtual ~WellFunctionCoefficient() {}

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip)
  {
    double x[3];
    mfem::Vector transip(x, 3);
    T.Transform(ip, transip);
    return eval(transip);
  }

private:
  const Well& injection;
  const Well& production;
  double inflow;
  double outflow;
  int spacedim;

  double eval(const mfem::Vector& point)
  {
    if (spacedim == 2)
      return eval2D(point);
    else if (spacedim == 3)
      return eval3D(point);
    else throw std::runtime_error("Unknown spacedim");
  }

  double eval2D(const mfem::Vector& point)
  {
    const double x = point(0);
    const double y = point(1);

    if (fabs(x-injection.center(0)) <= injection.radius &&
        fabs(y-injection.center(1)) <= injection.radius)
      return inflow;

    if (fabs(x-production.center(0)) <= production.radius &&
        fabs(y-production.center(1)) <= production.radius)
      return outflow;

    return 0.;
  }

  double eval3D(const mfem::Vector& point)
  {
    const double x = point(0);
    const double y = point(1);
    const double z = point(2);

    if (fabs(x-injection.center(0)) <= injection.radius &&
        fabs(z-injection.center(2)) <= injection.radius &&
        y >= injection.center(1) &&
        y <= injection.center(1)+injection.height)
      return inflow;

    if (fabs(x-production.center(0)) <= production.radius &&
        fabs(z-production.center(2)) <= production.radius &&
        y >= production.center(1) &&
        y <= production.center(1)+production.height)
      return outflow;

    return 0.;
  }
};



class SaturationSourceCoefficient : public mfem::Coefficient
{
public:
  SaturationSourceCoefficient(const Well& inject, double source, int dim)
    : injection(inject), saturation_source(source), spacedim(dim)
  {}

  virtual ~SaturationSourceCoefficient() {}

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip)
  {
    double x[3];
    mfem::Vector transip(x, 3);
    T.Transform(ip, transip);
    return eval(transip);
  }

private:
  const Well& injection;
  double saturation_source;
  int spacedim;

  double eval(const mfem::Vector& point)
  {
    if (spacedim == 2)
      return eval2D(point);
    else if (spacedim == 3)
      return eval3D(point);
    else throw std::runtime_error("Unknown spacedim");
  }

  double eval2D(const mfem::Vector& point)
  {
    const double x = point(0);
    const double y = point(1);

    if (fabs(x-injection.center(0)) <= injection.radius &&
        fabs(y-injection.center(1)) <= injection.radius)
      return saturation_source;

    return 0.;
  }

  double eval3D(const mfem::Vector& point)
  {
    const double x = point(0);
    const double y = point(1);
    const double z = point(2);

    if (fabs(x-injection.center(0)) <= injection.radius &&
        fabs(z-injection.center(2)) <= injection.radius &&
        y >= injection.center(1) &&
        y <= injection.center(1)+injection.height)
      return saturation_source;

    return 0.;
  }
};




class ValuesInCells : public mfem::Coefficient
{
public:
  ValuesInCells(mfem::Coefficient &coef, mfem::Vector &in_cells,
                std::vector<int> &flags, int ncells)
    : coefficient(coef), values_in_cells(in_cells), flags_of_cells(flags),
      n_cells(ncells)
  {}

  virtual ~ValuesInCells() {}

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip)
  {
    const int index = T.Attribute - 1; // use attribute as a cell number
    MFEM_ASSERT(index >= 0 && index < n_cells, "Element number (attribute) is "
                "out of range: " + d2s(index));

    const double val = coefficient.Eval(T, ip);
    values_in_cells(index) = val;
    flags_of_cells[index] = 1;
    return val;
  }

private:
  mfem::Coefficient &coefficient;
  mfem::Vector &values_in_cells;
  std::vector<int> &flags_of_cells;
  int n_cells;
};

#endif // PARAM_HPP
