#ifndef GRID_HPP
#define GRID_HPP

#include "config.hpp"
#include "mfem.hpp"
#include "utilities.hpp"

#include <iostream>

//#define TWO_PHASE_FLOW

double Krw(double S);
double Kro(double S);



struct Grid
{
  int nx, ny;
  double sx, sy;
  double V;
  double *K_array;
  double *Q_array;
  double *por_array;
  double *r_array;
  int n_cells;

  Grid(int _nx, int _ny, double _sx, double _sy)
    : nx(_nx), ny(_ny), sx(_sx), sy(_sy), V(0.0),
      K_array(nullptr), Q_array(nullptr), por_array(nullptr), r_array(nullptr)
  {
    double hx = sx / nx;
    double hy = sy / ny;
    V = hx*hy;
    n_cells = nx*ny;
  }

  ~Grid()
  {
    delete[] r_array;
    delete[] por_array;
    delete[] Q_array;
    delete[] K_array;
  }

  void init_arrays()
  {
    K_array   = new double[n_cells];
    Q_array   = new double[n_cells];
    por_array = new double[n_cells];
    r_array   = new double[n_cells];
    for (int i = 0; i < n_cells; ++i)
    {
      K_array[i]   = 1.0;
      Q_array[i]   = 0.0;
      por_array[i] = 1.0;
      r_array[i]   = 0.0;
    }
    Q_array[0]         =  1.0; // injection well
    Q_array[n_cells-1] = -1.0; // production well
//    r_array[n_cells-1] =  1.0;
    r_array[0] =  1.0;
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
    return val*val_array[T.ElementNo];
  }

protected:
  mfem::GridFunctionCoefficient& func;
  double mu_w, mu_o;
  double *val_array;
  int n_cells;
  bool own_array;
};

#endif // GRID_HPP
