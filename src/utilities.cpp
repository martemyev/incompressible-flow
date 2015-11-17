#include "mfem.hpp"
#include "utilities.hpp"

#include <cmath>
#include <fstream>

using namespace std;
using namespace mfem;



double to_radians(double x)
{
  return x*M_PI/180.0;
}



void read_binary(const char *filename, int n_values, double *values)
{
  ifstream in(filename, ios::binary);
  MFEM_VERIFY(in, "File '" + string(filename) + "' can't be opened");

  in.seekg(0, in.end); // jump to the end of the file
  int length = in.tellg(); // total length of the file in bytes
  int size_value = length / n_values; // size (in bytes) of one value

  MFEM_VERIFY(length % n_values == 0, "The number of bytes in the file '" +
              string(filename) + "' is not divisible by the number of elements "
              + d2s(n_values));

  in.seekg(0, in.beg); // jump to the beginning of the file

  if (size_value == sizeof(double))
  {
    in.read((char*)values, n_values*size_value); // read all at once

    MFEM_VERIFY(n_values == static_cast<int>(in.gcount()), "The number of "
                "successfully read elements is different from the expected one");
  }
  else if (size_value == sizeof(float))
  {
    float val = 0;
    for (int i = 0; i < n_values; ++i)  // read element-by-element
    {
      in.read((char*)&val, size_value); // read a 'float' value
      values[i] = val;                  // convert it to a 'double' value
    }
  }
  else MFEM_VERIFY(0, "Unknown size of an element (" + d2s(size_value) + ") in "
                   "bytes. Expected one is either sizeof(float) = " +
                   d2s(sizeof(float)) + ", or sizeof(double) = " +
                   d2s(sizeof(double)));

  in.close();
}



void write_binary(const char *filename, int n_values, double *values)
{
  ofstream out(filename, ios::binary);
  MFEM_VERIFY(out, "File '" + string(filename) + "' can't be opened");

  for (int i = 0; i < n_values; ++i)
  {
    float val = values[i];
    out.write(reinterpret_cast<char*>(&val), sizeof(float));
  }

  out.close();
}



double get_max(double *a, int n_elements)
{
  double max_val = a[0];
  for (int i = 1; i < n_elements; ++i)
    max_val = max(max_val, a[i]);
  return max_val;
}



void get_minmax(double *a, int n_elements, double &min_val, double &max_val)
{
  min_val = max_val = a[0];
  for (int i = 1; i < n_elements; ++i)
  {
    min_val = min(min_val, a[i]);
    max_val = max(max_val, a[i]);
  }
}



void get_LE(double *a, int n_elements, double LE_value, double *a_LE)
{
  for (int i = 0; i < n_elements; ++i)
  {
    if (a[i] <= LE_value)
      a_LE[i] = a[i];
    else
      a_LE[i] = LE_value;
  }
}



void get_GE(double *a, int n_elements, double GE_value, double *a_GE)
{
  for (int i = 0; i < n_elements; ++i)
  {
    if (a[i] >= GE_value)
      a_GE[i] = a[i];
    else
      a_GE[i] = GE_value;
  }
}



void write_vts_vector(const std::string& filename, const std::string& solname,
                      double sx, double sy, double sz, int nx, int ny, int nz,
                      const Vector& sol_x, const Vector& sol_y,
                      const Vector& sol_z)
{
  ofstream out(filename.c_str());
  MFEM_VERIFY(out, "File '" + filename + "' can't be opened");

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
  out << "  <StructuredGrid WholeExtent=\"1 " << nx+1 << " 1 " << ny+1 << " 1 " << nz+1 << "\">\n";
  out << "    <Piece Extent=\"1 " << nx+1 << " 1 " << ny+1 << " 1 " << nz+1 << "\">\n";
  out << "      <PointData Vectors=\"" << solname << "\" Scalars=\"" << solname << "_scalar_mag\">\n";
  out << "        <DataArray type=\"Float64\" Name=\"" << solname << "\" format=\"ascii\" NumberOfComponents=\"3\">\n";
  for (int iz = 0; iz < nz+1; ++iz)
    for (int iy = 0; iy < ny+1; ++iy)
      for (int ix = 0; ix < nx+1; ++ix)
      {
        const int glob_vert_index = iz*(nx+1)*(ny+1) + iy*(nx+1) + ix;
        out << sol_x(glob_vert_index) << " "
            << sol_y(glob_vert_index) << " "
            << sol_z(glob_vert_index) << " ";
      }
  out << "\n";
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Float64\" Name=\"" << solname << "_scalar_mag\" format=\"ascii\" NumberOfComponents=\"1\">\n";
  for (int iz = 0; iz < nz+1; ++iz)
    for (int iy = 0; iy < ny+1; ++iy)
      for (int ix = 0; ix < nx+1; ++ix)
      {
        const int glob_vert_index = iz*(nx+1)*(ny+1) + iy*(nx+1) + ix;
        double mag = pow(sol_x(glob_vert_index), 2) +
                     pow(sol_y(glob_vert_index), 2) +
                     pow(sol_z(glob_vert_index), 2);
        out << sqrt(mag) << " ";
      }
  out << "\n";
  out << "        </DataArray>\n";
  out << "      </PointData>\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  const double hx = sx / nx;
  const double hy = sy / ny;
  const double hz = sz / nz;

  for (int iz = 0; iz < nz+1; ++iz)
  {
    const double z = (iz == nz ? sz : iz*hz);
    for (int iy = 0; iy < ny+1; ++iy)
    {
      const double y = (iy == ny ? sy : iy*hy);
      for (int ix = 0; ix < nx+1; ++ix)
      {
        const double x = (ix == nx ? sx : ix*hx);
        out << x << " " << y << " " << z << " ";
      }
    }
  }

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "    </Piece>\n";
  out << "  </StructuredGrid>\n";
  out << "</VTKFile>\n";

  out.close();
}



void write_vts_vector(const std::string& filename, const std::string& solname,
                      double sx, double sy, int nx, int ny,
                      const Vector& sol_x, const Vector& sol_y)
{
  ofstream out(filename.c_str());
  MFEM_VERIFY(out, "File '" + filename + "' can't be opened");

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
  out << "  <StructuredGrid WholeExtent=\"1 " << nx+1 << " 1 1 1 " << ny+1 << "\">\n";
  out << "    <Piece Extent=\"1 " << nx+1 << " 1 1 1 " << ny+1 << "\">\n";
  out << "      <PointData Vectors=\"" << solname << "\">\n";
  out << "        <DataArray type=\"Float64\" Name=\"" << solname
      << "\" format=\"ascii\" NumberOfComponents=\"3\">\n";

  for (int i = 0; i < ny+1; ++i)
  {
    for (int j = 0; j < nx+1; ++j)
    {
      const int glob_vert_index = i*(nx+1) + j;
      out << sol_x(glob_vert_index) << " "
          << sol_y(glob_vert_index) << " 0.0 ";
    }
  }

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </PointData>\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         "format=\"ascii\">\n";

  const double hx = sx / nx;
  const double hy = sy / ny;

  for (int i = 0; i < ny+1; ++i)
  {
    const double y = (i == ny ? sy : i*hy);
    for (int j = 0; j < nx+1; ++j)
    {
      const double x = (j == nx ? sx : j*hx);
      out << x << " " << y << " 0.0 ";
    }
  }

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "    </Piece>\n";
  out << "  </StructuredGrid>\n";
  out << "</VTKFile>\n";

  out.close();
}



void write_vts_scalar(const std::string& filename, const std::string& solname,
                      double sx, double sy, double sz, int nx, int ny, int nz,
                      const Vector& sol)
{
  ofstream out(filename.c_str());
  MFEM_VERIFY(out, "File '" + filename + "' can't be opened");

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
  out << "  <StructuredGrid WholeExtent=\"1 " << nx+1 << " 1 " << ny+1 << " 1 " << nz+1 << "\">\n";
  out << "    <Piece Extent=\"1 " << nx+1 << " 1 " << ny+1 << " 1 " << nz+1 << "\">\n";
  out << "      <PointData Scalars=\"" << solname << "\">\n";
  out << "        <DataArray type=\"Float64\" Name=\"" << solname << "\" format=\"ascii\" NumberOfComponents=\"1\">\n";

  for (int iz = 0; iz < nz+1; ++iz)
    for (int iy = 0; iy < ny+1; ++iy)
      for (int ix = 0; ix < nx+1; ++ix)
      {
        const int glob_vert_index = iz*(nx+1)*(ny+1) + iy*(nx+1) + ix;
        out << sol(glob_vert_index) << " ";
      }

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </PointData>\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

  const double hx = sx / nx;
  const double hy = sy / ny;
  const double hz = sz / nz;

  for (int iz = 0; iz < nz+1; ++iz)
  {
    const double z = (iz == nz ? sz : iz*hz);
    for (int iy = 0; iy < ny+1; ++iy)
    {
      const double y = (iy == ny ? sy : iy*hy);
      for (int ix = 0; ix < nx+1; ++ix)
      {
        const double x = (ix == nx ? sx : ix*hx);
        out << x << " " << y << " " << z << " ";
      }
    }
  }

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "    </Piece>\n";
  out << "  </StructuredGrid>\n";
  out << "</VTKFile>\n";

  out.close();
}



void write_vts_scalar(const std::string& filename, const std::string& solname,
                      double sx, double sy, int nx, int ny,
                      const Vector& sol)
{
  ofstream out(filename.c_str());
  MFEM_VERIFY(out, "File '" + filename + "' can't be opened");

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
  out << "  <StructuredGrid WholeExtent=\"1 " << nx+1 << " 1 1 1 " << ny+1 << "\">\n";
  out << "    <Piece Extent=\"1 " << nx+1 << " 1 1 1 " << ny+1 << "\">\n";
  out << "      <PointData Scalars=\"" << solname << "\">\n";
  out << "        <DataArray type=\"Float64\" Name=\"" << solname
      << "\" format=\"ascii\" NumberOfComponents=\"1\">\n";

  for (int i = 0; i < ny+1; ++i)
  {
    for (int j = 0; j < nx+1; ++j)
    {
      const int glob_vert_index = i*(nx+1) + j;
      out << sol(glob_vert_index) << " ";
    }
  }

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </PointData>\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         "format=\"ascii\">\n";

  const double hx = sx / nx;
  const double hy = sy / ny;

  for (int i = 0; i < ny+1; ++i)
  {
    const double y = (i == ny ? sy : i*hy);
    for (int j = 0; j < nx+1; ++j)
    {
      const double x = (j == nx ? sx : j*hx);
      out << x << " " << y << " 0.0 ";
    }
  }

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  out << "    </Piece>\n";
  out << "  </StructuredGrid>\n";
  out << "</VTKFile>\n";

  out.close();
}
