#include "mfem.hpp"
#include "utilities.hpp"
#include "param.hpp"

#include <cfloat>
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



void write_vts_scalar_cells(const std::string& filename, const std::string& solname,
                            double sx, double sy, int nx, int ny, const Vector& sol)
{
  MFEM_VERIFY(sol.Size() == nx*ny, "Dimension mismatch");

  ofstream out(filename.c_str());
  MFEM_VERIFY(out, "File '" + filename + "' can't be opened");

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
  out << "  <StructuredGrid WholeExtent=\"1 " << nx+1 << " 1 " << ny+1 << " 1 1\">\n";
  out << "    <Piece Extent=\"1 " << nx+1 << " 1 " << ny+1 << " 1 1\">\n";
  out << "      <CellData Scalars=\"" << solname << "\">\n";
  out << "        <DataArray type=\"Float64\" Name=\"" << solname
      << "\" format=\"ascii\" NumberOfComponents=\"1\">\n";

  for (int i = 0; i < ny; ++i)
    for (int j = 0; j < nx; ++j)
      out << sol(i*nx + j) << " ";

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </CellData>\n";
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

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

void write_vts_scalar_cells(const std::string& filename, const std::string& solname,
                            double sx, double sy, double sz, int nx, int ny, int nz,
                            const Vector& sol)
{
  MFEM_VERIFY(sol.Size() == nx*ny*nz, "Dimension mismatch");

  ofstream out(filename.c_str());
  MFEM_VERIFY(out, "File '" + filename + "' can't be opened");

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
  out << "  <StructuredGrid WholeExtent=\"1 " << nx+1 << " 1 " << ny+1 << " 1 " << nz+1 << "\">\n";
  out << "    <Piece Extent=\"1 " << nx+1 << " 1 " << ny+1 << " 1 " << nz+1 << "\">\n";
  out << "      <CellData Scalars=\"" << solname << "\">\n";
  out << "        <DataArray type=\"Float64\" Name=\"" << solname
      << "\" format=\"ascii\" NumberOfComponents=\"1\">\n";

  for (int iz = 0; iz < nz; ++iz)
    for (int iy = 0; iy < ny; ++iy)
      for (int ix = 0; ix < nx; ++ix)
        out << sol(iz*nx*ny + iy*nx + ix) << " ";

  out << "\n";
  out << "        </DataArray>\n";
  out << "      </CellData>\n";
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



void compute_in_cells(double sx, double sy, int nx, int ny,
                      const Mesh& mesh, const GridFunction& U,
                      Vector& values)
{
  MFEM_VERIFY(values.Size() == nx*ny, "Dimensions mismatch");

  const double hx = sx / nx;
  const double hy = sy / ny;

  for (int iy = 0; iy < ny; ++iy)
  {
    const double y = (iy+0.5)*hy;
    for (int ix = 0; ix < nx; ++ix)
    {
      const double x = (ix+0.5)*hx;
      const int cell = iy*nx + ix;

      const Element *element = mesh.GetElement(cell);
      Array<int> vert_indices;
      element->GetVertices(vert_indices);
      const double *vert0 = mesh.GetVertex(vert_indices[0]); // min coords

      IntegrationPoint ip;
      ip.x = (x - vert0[0]) / hx;
      ip.y = (y - vert0[1]) / hy;

      values(cell) = U.GetValue(cell, ip);
    }
  }
}



void compute_in_cells(double sx, double sy, double sz, int nx, int ny, int nz,
                      const Mesh& mesh, const GridFunction& U,
                      Vector& values)
{
  MFEM_VERIFY(values.Size() == nx*ny*nz, "Dimensions mismatch");

  const double hx = sx / nx;
  const double hy = sy / ny;
  const double hz = sz / nz;

  for (int iz = 0; iz < nz; ++iz)
  {
    const double z = (iz+0.5)*hz;
    for (int iy = 0; iy < ny; ++iy)
    {
      const double y = (iy+0.5)*hy;
      for (int ix = 0; ix < nx; ++ix)
      {
        const double x = (ix+0.5)*hx;
        const int cell = iz*nx*ny + iy*nx + ix;

        const Element *element = mesh.GetElement(cell);
        Array<int> vert_indices;
        element->GetVertices(vert_indices);
        const double *vert0 = mesh.GetVertex(vert_indices[0]); // min coords

        IntegrationPoint ip;
        ip.x = (x - vert0[0]) / hx;
        ip.y = (y - vert0[1]) / hy;
        ip.z = (z - vert0[2]) / hz;

        values(cell) = U.GetValue(cell, ip);
      }
    }
  }
}



void output_scalar(const Param& p, const GridFunction& x, const string& tstr,
                   const string& name)
{
  string phase = (p.two_phase_flow ? "2phase" : "1phase");
  string fname = string(p.outdir) + "/" + name + "_" + p.extra + "_" +
                 d2s(p.spacedim) + "D_" + phase + "_" + tstr + ".vts";

  Vector x_nodal;
  x.GetNodalValues(x_nodal);
  if (p.spacedim == 2)
    write_vts_scalar(fname, name, p.sx, p.sy, p.nx, p.ny, x_nodal);
  else if (p.spacedim == 3)
    write_vts_scalar(fname, name, p.sx, p.sy, p.sz, p.nx, p.ny, p.nz, x_nodal);
  else MFEM_ABORT("Not supported spacedim");
}



void output_scalar_cells(const Param& p, Vector x, std::vector<int> flags,
                         const string& tstr, const string& name)
{
  string phase = (p.two_phase_flow ? "2phase" : "1phase");
  string fname = string(p.outdir) + "/" + name + "_" + p.extra + "_" +
                 d2s(p.spacedim) + "D_" + phase + "_" + tstr + ".vts";

#if defined(MFEM_USE_MPI)
  int num_procs, myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  const int n_cells = p.get_n_cells();
  const int tag_values = 11;
  const int tag_flags = 12;
  if (myid == 0)
  {
    double *Xvalues = new double[n_cells];
    int *Xflags = new int[n_cells];
    MPI_Status status;
    for (int rank = 1; rank < num_procs; ++rank)
    {
      MPI_Recv(Xvalues, n_cells, MPI_DOUBLE, MPI_ANY_SOURCE, tag_values,
               MPI_COMM_WORLD, &status);
      MPI_Recv(Xflags, n_cells, MPI_INTEGER, MPI_ANY_SOURCE, tag_flags,
               MPI_COMM_WORLD, &status);
      for (int el = 0; el < n_cells; ++el)
      {
        if (Xflags[el])
          x(el) = Xvalues[el];
      }
    }
    delete[] Xflags;
    delete[] Xvalues;

    if (p.spacedim == 2)
      write_vts_scalar_cells(fname, name, p.sx, p.sy, p.nx, p.ny, x);
    else if (p.spacedim == 3)
      write_vts_scalar_cells(fname, name, p.sx, p.sy, p.sz, p.nx, p.ny, p.nz, x);
    else MFEM_ABORT("Unknown spacedim");
  }
  else
  {
    MPI_Send(x.GetData(), n_cells, MPI_DOUBLE, 0, tag_values, MPI_COMM_WORLD);
    MPI_Send(&flags[0], n_cells, MPI_INTEGER, 0, tag_flags, MPI_COMM_WORLD);
  }
#else // MFEM_USE_MPI
  if (p.spacedim == 2)
    write_vts_scalar_cells(fname, name, p.sx, p.sy, p.nx, p.ny, x);
  else if (p.spacedim == 3)
    write_vts_scalar_cells(fname, name, p.sx, p.sy, p.sz, p.nx, p.ny, p.nz, x);
#endif // MFEM_USE_MPI
}



void output_vector(const Param& p, const GridFunction& V, const string& tstr,
                   const string& name)
{
  string phase = (p.two_phase_flow ? "2phase" : "1phase");
  string fname = name + "_" + p.extra + "_" + d2s(p.spacedim) + "D_" + phase +
                 "_" + tstr + ".vts";

  Vector Vx_nodal, Vy_nodal;
  V.GetNodalValues(Vx_nodal, 1);
  V.GetNodalValues(Vy_nodal, 2);
  if (p.spacedim == 2)
    write_vts_vector(fname, name, p.sx, p.sy, p.nx, p.ny, Vx_nodal, Vy_nodal);
  else if (p.spacedim == 3)
  {
    Vector Vz_nodal;
    V.GetNodalValues(Vz_nodal, 3);
    write_vts_vector(fname, name, p.sx, p.sy, p.sz, p.nx, p.ny, p.nz,
                     Vx_nodal, Vy_nodal, Vz_nodal);
  }
  else MFEM_ABORT("Not supported spacedim");
}



void output_seismic_properties(const Param& p, int ti,
                               const Vector& rho_array, const Vector& vp_array,
                               const Vector& vs_array)
{
  string fname_rho_bin, fname_vp_bin, fname_vs_bin;
  string fname_rho_vts, fname_vp_vts, fname_vs_vts;
  string extra = (string)p.extra + "_";

  const int n_cells = p.get_n_cells();

  string phase = (p.two_phase_flow ? "2phase" : "1phase");
  fname_rho_bin = extra + "1phase_" + d2s(n_cells) + "_t" + d2s(ti) + ".rho";
  fname_vp_bin  = extra + "1phase_" + d2s(n_cells) + "_t" + d2s(ti) + ".vp";
  fname_vs_bin  = extra + "1phase_" + d2s(n_cells) + "_t" + d2s(ti) + ".vs";
  fname_rho_vts = "rho_" + extra + d2s(p.spacedim) + "D_" + phase + "_" + d2s(ti) + ".vts";
  fname_vp_vts  = "vp_" + extra + d2s(p.spacedim) + "D_" + phase + "_"  + d2s(ti) + ".vts";
  fname_vs_vts  = "vs_" + extra + d2s(p.spacedim) + "D_" + phase + "_"  + d2s(ti) + ".vts";

  write_binary(fname_rho_bin.c_str(), n_cells, rho_array.GetData());
  write_binary(fname_vp_bin.c_str(),  n_cells, vp_array.GetData());
  write_binary(fname_vs_bin.c_str(),  n_cells, vs_array.GetData());

  if (p.spacedim == 2)
  {
    write_vts_scalar_cells(fname_rho_vts, "density", p.sx, p.sy,
                           p.nx, p.ny, rho_array);
    write_vts_scalar_cells(fname_vp_vts, "vp", p.sx, p.sy,
                           p.nx, p.ny, vp_array);
    write_vts_scalar_cells(fname_vs_vts, "vs", p.sx, p.sy,
                           p.nx, p.ny, vs_array);
  }
  else if (p.spacedim == 3)
  {
    write_vts_scalar_cells(fname_rho_vts, "density", p.sx, p.sy, p.sz,
                           p.nx, p.ny, p.nz, rho_array);
    write_vts_scalar_cells(fname_vp_vts, "vp", p.sx, p.sy, p.sz,
                           p.nx, p.ny, p.nz, vp_array);
    write_vts_scalar_cells(fname_vs_vts, "vs", p.sx, p.sy, p.sz,
                           p.nx, p.ny, p.nz, vs_array);
  }
  else MFEM_ABORT("Not supported spacedim");
}



double K_func(double vp, double vs, double rho)
{
  return rho * (vp*vp - 4./3.*vs*vs);
}

double G_func(double vs, double rho)
{
  return rho*vs*vs;
}

double rho_B(double rho_fl, double rho_g, double phi)
{
  return rho_g*(1.0-phi) + rho_fl*phi;
}

double K_fl(double S_w, double K_w, double K_o)
{
  return 1.0 / (S_w/K_w + (1.0-S_w)/K_o);
}

double rho_fl(double S_w, double rho_w, double rho_o)
{
  return S_w*rho_w + (1.0-S_w)*rho_o;
}

//double K_frame(double K_sat, double K_m, double K_fl, double phi)
//{
//  double a = phi*K_m/K_fl;
//  double numer = K_sat*(a+1.0-phi) - K_m;
//  double denom = a + K_sat/K_m - 1.0 - phi;
//  return numer/denom;
//}

double K_frame(double K1, double K2, double F1, double F2)
{
  MFEM_VERIFY(F1 >= 0. && F1 <= 1. && F2 >= 0. && F2 <= 1. && (F1+F2-1.) < 1e-12,
              "F1 or F2 is out of range");
  double K_Reuss = 1.0 / (F1/K1 + F2/K2);
  double K_Voigt = F1*K1 + F2*K2;
  return 0.5 * (K_Reuss + K_Voigt);
}

double K_sat(double K_frame, double K_m, double K_fl, double phi)
{
  double numer = (1.0-K_frame/K_m)*(1.0-K_frame/K_m);
  double denom = phi/K_fl + (1.0-phi)/K_m - K_frame/K_m/K_m;
  return K_frame + numer/denom;
}

double vp_func(double K, double G, double rho)
{
  return sqrt((K + 4./3.*G) / rho);
}

double vs_func(double G, double rho)
{
  return sqrt(G/rho);
}

void Gassmann(const Vector& S, const Param& param, double K_m, double Kframe,
              double rho_gr, double *phi, double *rho, double *vp, double *vs)
{
  const int n_cells = param.get_n_cells();
  MFEM_VERIFY(S.Size() == n_cells, "Sizes mismatch");

  const double K_w = K_func(VP_W, VS_W, RHO_W); // bulk modulus of water
  const double K_o = K_func(VP_O, VS_O, RHO_O); // bulk modulus of oil

  double K, G, K_fl_mix, rho_fl_mix, Ksat;
  double minKsat = DBL_MAX, maxKsat = DBL_MIN;

  for (int i = 0; i < n_cells; ++i)
  {
    K = K_func(vp[i], vs[i], rho[i]);
    G = G_func(vs[i], rho[i]);

    K_fl_mix   = K_fl(S(i), K_w, K_o);
    rho_fl_mix = rho_fl(S(i), RHO_W, RHO_O);

    Ksat = K_sat(Kframe, K_m, K_fl_mix, phi[i]);
    minKsat = min(minKsat, Ksat);
    maxKsat = max(maxKsat, Ksat);

    rho[i] = rho_B(rho_fl_mix, rho_gr, phi[i]);
    vp[i]  = vp_func(Ksat, G, rho[i]);
    vs[i]  = vs_func(G, rho[i]);
  }

  cout << "minKsat = " << minKsat << endl;
  cout << "maxKsat = " << maxKsat << endl;
}

double compute_rho(const Vector &x)
{
  return 0.;
}

double compute_vp(const Vector &x)
{
  return 0.;
}

double compute_vs(const Vector &x)
{
  return 0.;
}


