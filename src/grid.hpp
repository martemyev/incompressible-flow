#ifndef GRID_HPP
#define GRID_HPP

#include "config.hpp"
#include "mfem.hpp"

#include <iostream>



struct Grid
{
  int nx, ny;
  double sx, sy;
  double V;
  double *K_array;
  double *Q_array;
  double *por_array;
  double *r_array;

  Grid(int _nx, int _ny, double _sx, double _sy)
    : nx(_nx), ny(_ny), sx(_sx), sy(_sy), V(0.0),
      K_array(nullptr), Q_array(nullptr), por_array(nullptr), r_array(nullptr)
  {
    double hx = sx / nx;
    double hy = sy / ny;
    V = hx*hy;
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
    K_array   = new double[nx*ny];
    Q_array   = new double[nx*ny];
    por_array = new double[nx*ny];
    r_array   = new double[nx*ny];
    for (int i = 0; i < nx*ny; ++i)
    {
      K_array[i]   = 1.0;
      Q_array[i]   = 0.0;
      por_array[i] = 1.0;
      r_array[i]   = 0.0;
    }
    Q_array[0]       = -1.0; // injection well
    Q_array[nx*ny-1] =  1.0; // production well
    r_array[0]       =  1.0;
  }

};



/**
 * Cell-wise constant coefficient.
 */
class CWConstCoefficient : public mfem::Coefficient
{
public:
  CWConstCoefficient(double *array, bool own = 1)
    : val_array(array), own_array(own)
  { }

  virtual ~CWConstCoefficient() { if (own_array) delete[] val_array; }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip)
  {
    return val_array[T.ElementNo];
  }

protected:
  double *val_array;
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
                double *array, bool own = 1)
    : func(f)
    , mu_w(muw)
    , mu_o(muo)
    , val_array(array)
    , own_array(own)
  { }

  virtual ~CWCoefficient() { if (own_array) delete[] val_array; }

  virtual double Eval(mfem::ElementTransformation &T,
                      const mfem::IntegrationPoint &ip)
  {
    const double S = func.Eval(T, ip);
    const double val = Krw(S) / mu_w + Kro(S) / mu_o;
    return val*val_array[T.ElementNo];
  }

  static double Krw(double S) { return S; } //*S; }
  static double Kro(double S) { return (1.0-S); } //*(1.0-S); }

protected:
  mfem::GridFunctionCoefficient& func;
  double mu_w, mu_o;
  double *val_array;
  bool own_array;
};

#endif // GRID_HPP
