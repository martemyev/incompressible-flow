#ifndef PARAM_HPP
#define PARAM_HPP

#include "config.hpp"
#include "mfem.hpp"
#include "utilities.hpp"

#include <iostream>

#define TWO_PHASE_FLOW

double Krw(double S);
double Kro(double S);



struct Param
{
  int spacedim;
  int nx, ny, nz;
  int n_cells;
  double sx, sy, sz;
  double V; // volume (area) of a cell
  double *K_array; // permeability array (cell-wise constant - one value per cell)
  double *Q_array; // pressure source (injection and production wells - inflow and outflow)
  double *por_array; // array of porosity values (cell-wise constant)
  double *R_array; // saturation source (cell-wise constant)

  int order_v, order_p, order_s;
  double t_final, dt;
  int vis_steps;

  Param();
  ~Param();
  void init_arrays();
  void add_options(mfem::OptionsParser& args);

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
                      const mfem::IntegrationPoint &ip)
  {
    MFEM_ASSERT(T.ElementNo >= 0 && T.ElementNo < n_cells, "Element number is "
                "out of range: " + d2s(T.ElementNo));
    return val_array[T.ElementNo];
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
                double *array, int ncells, bool own = 1)
    : func(f)
    , mu_w(muw)
    , mu_o(muo)
    , val_array(array)
    , n_cells(ncells)
    , own_array(own)
  { }

  virtual ~CWCoefficient() { if (own_array) delete[] val_array; }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip)
  {
    MFEM_ASSERT(T.ElementNo >= 0 && T.ElementNo < n_cells, "Element number is "
                "out of range: " + d2s(T.ElementNo));
    const double S = func.Eval(T, ip);
    const double val = Krw(S) / mu_w + Kro(S) / mu_o;
    return 1./(val*val_array[T.ElementNo]);
  }

protected:
  mfem::GridFunctionCoefficient& func;
  double mu_w, mu_o;
  double *val_array;
  int n_cells;
  bool own_array;
};

#endif // PARAM_HPP
