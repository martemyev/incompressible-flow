#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include "config.hpp"

#include <sstream>
#include <stdexcept>

namespace mfem
{
  class GridFunction;
  class Mesh;
  class Vector;
}
class Param;

/**
 * Convert the data of any type which has oveloaded operator '<<' to string
 * @param data - the data
 * @param scientific - use scientific (exponential) format or not
 * @param precision - if scientific format is used, we can change the precision
 * @param noperiod - if true and the floating point number is being converted,
 * the period separating the number is removed
 * @param width - the width of a resulting string (filled with 0 if the number
 * is shorter)
 * @return data in string format
 */
template <typename T>
inline std::string d2s(T data,
                       bool scientific = false,
                       int precision = 0,
                       bool noperiod = false,
                       int width = 0)
{
  const char filler = '0';
  std::ostringstream o;
  if (scientific)    o.setf(std::ios::scientific);
  if (precision > 0) o.precision(precision);
  if (width > 0)     o << std::setfill(filler) << std::setw(width);

  if (!(o << data))
    throw std::runtime_error("Bad conversion of data to string!");

  // eliminate a period in case of floating-point numbers in non-scientific
  // format
  if (!scientific && noperiod)
  {
    std::string res = o.str(); // resulting string
    std::string::size_type pos = res.find('.');
    if (pos != std::string::npos)
      res.erase(pos, 1);
    if (width > 0 && static_cast<int>(res.size()) < width)
      res.insert(0, width-res.size(), filler);
    return res;
  }

  return o.str();
}

/**
 * Convert an angle in degrees to radians.
 */
double to_radians(double x);

/**
 * Read a binary file
 */
void read_binary(const char *filename, int n_values, double *values);

/**
 * Write a binary file
 */
void write_binary(const char *filename, int n_values, double *values);

/**
 * Find the max value of the given array (vector) a
 */
double get_max(double *a, int n_elements);

/**
 * Find min and max values of the given array (vector) a
 */
void get_minmax(double *a, int n_elements, double &min_val, double &max_val);

/**
 * Get an array of values of array 'a' that are less or equal than LE_value
 */
void get_LE(double *a, int n_elements, double LE_value, double *a_LE);

/**
 * Get an array of values of array 'a' that are greater or equal than GE_value
 */
void get_GE(double *a, int n_elements, double GE_value, double *a_GE);

/**
 * Write a snapshot of a vector wavefield in a VTS format
 * @param filename - output file name
 * @param solname - name of the wavefield
 * @param sx - size of domain in x-direction
 * @param sy - size of domain in y-direction
 * @param sz - size of domain in z-direction
 * @param nx - number of cells in x-direction
 * @param ny - number of cells in y-direction
 * @param nz - number of cells in z-direction
 * @param sol_x - x-component of the vector wavefield
 * @param sol_y - y-component of the vector wavefield
 */
void write_vts_vector(const std::string& filename, const std::string& solname,
                      double sx, double sy, double sz, int nx, int ny, int nz,
                      const mfem::Vector& sol_x, const mfem::Vector& sol_y,
                      const mfem::Vector& sol_z);
void write_vts_vector(const std::string& filename, const std::string& solname,
                      double sx, double sy, int nx, int ny,
                      const mfem::Vector& sol_x, const mfem::Vector& sol_y);

/**
 * Write scalar values in a VTS format
 * @param filename - output file name
 * @param solname - name of the wavefield
 * @param sx - size of domain in x-direction
 * @param sy - size of domain in y-direction
 * @param sz - size of domain in z-direction
 * @param nx - number of cells in x-direction
 * @param ny - number of cells in y-direction
 * @param nz - number of cells in z-direction
 * @param sol - scalar values
 */
void write_vts_scalar(const std::string& filename, const std::string& solname,
                      double sx, double sy, double sz, int nx, int ny, int nz,
                      const mfem::Vector& sol);
void write_vts_scalar(const std::string& filename, const std::string& solname,
                      double sx, double sy, int nx, int ny,
                      const mfem::Vector& sol);
void write_vts_scalar_cells(const std::string& filename, const std::string& solname,
                            double sx, double sy, int nx, int ny,
                            const mfem::Vector& sol);
void write_vts_scalar_cells(const std::string& filename, const std::string& solname,
                            double sx, double sy, double sz, int nx, int ny, int nz,
                            const mfem::Vector& sol);

void compute_in_cells(double sx, double sy, int nx, int ny,
                      const mfem::Mesh& mesh, const mfem::GridFunction& U,
                      mfem::Vector& values);
void compute_in_cells(double sx, double sy, double sz, int nx, int ny, int nz,
                      const mfem::Mesh& mesh, const mfem::GridFunction& U,
                      mfem::Vector& values);

void output_scalar(const Param& p, const mfem::Vector& P,
                   const std::string& tstr, const std::string& name);

void output_vector(const Param& p, const mfem::GridFunction& V,
                   const std::string& tstr, const std::string& name);

void output_seismic_properties(const Param& p, int ti,
                               const mfem::Vector& rho_array,
                               const mfem::Vector& vp_array,
                               const mfem::Vector& vs_array);

double K_func(double vp, double vs, double rho);
double G_func(double vs, double rho);
double vp_func(double K, double G, double rho);
double vs_func(double G, double rho);

double rho_B(double rho_fl, double rho_g, double phi);
double K_fl(double S_w, double K_w, double K_o);
double rho_fl(double S_w, double rho_w, double rho_o);
//double K_frame(double K_sat, double K_m, double K_fl, double phi);
double K_frame(double K1, double K2, double F1, double F2);
double K_sat(double K_frame, double K_m, double K_fl, double phi);

void Gassmann(const mfem::Vector& S, const Param& param, double K_m,
              double K_frame,double rho_gr, double *phi, double *rho,
              double *vp, double *vs);

#endif // UTILITIES_HPP
