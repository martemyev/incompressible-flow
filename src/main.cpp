#include "config.hpp"
#include "grid.hpp"
#include "mfem.hpp"
#include "pressure_solver.hpp"
#include "saturation_solver.hpp"
#include "utilities.hpp"

#include <cmath>

using namespace std;
using namespace mfem;



int main(int argc, char **argv)
{
  Grid grid(64, 64, 1.0, 1.0);
  grid.init_arrays();

  const bool generate_edges = true;
  Mesh mesh(grid.nx, grid.ny, Element::QUADRILATERAL, generate_edges,
            grid.sx, grid.sy);
  const int dim = mesh.Dimension();

//  const int order = 1;
  FiniteElementCollection *hdiv_coll = new RT_FECollection(0, dim); // for velocity
  FiniteElementCollection *l2_coll   = new L2_FECollection(0, dim); // for pressure
  DG_FECollection         *dg_coll   = new DG_FECollection(0, dim); // for saturation

  FiniteElementSpace V_space(&mesh, hdiv_coll);
  FiniteElementSpace P_space(&mesh, l2_coll);
  FiniteElementSpace S_space(&mesh, dg_coll);

  Array<int> block_offsets(3);
  block_offsets[0] = 0;
  block_offsets[1] = V_space.GetVSize();
  block_offsets[2] = P_space.GetVSize();
  block_offsets.PartialSum();

  cout << "***********************************************************\n";
  cout << "dim(V)  = " << block_offsets[1] - block_offsets[0] << "\n";
  cout << "dim(P)  = " << block_offsets[2] - block_offsets[1] << "\n";
  cout << "dim(V+P)= " << block_offsets.Last() << "\n";
  cout << "dim(S)  = " << S_space.GetVSize() << "\n";
  cout << "***********************************************************\n";
  cout << endl;

  const double T = 0.7;
  const int nt = 25;
  const double dt = T / nt;

  GridFunction V, P;

  BlockVector x(block_offsets);
  V.Update(&V_space, x.GetBlock(0), 0);
  P.Update(&P_space, x.GetBlock(1), 0);

  GridFunction S(&S_space);
  S = 0.0;
  GridFunctionCoefficient saturation(&S);
  VectorGridFunctionCoefficient velocity(&V);

  Vector P_nodal, S_nodal, Vx_nodal, Vy_nodal;

  string fname;
  for (int t = 0; t < nt; ++t)
  {
    PressureSolver(block_offsets, mesh, grid, V_space, P_space, saturation, x);
    SaturationSolver(grid, S, velocity, dt);

    const string tstr = d2s(t, 0, 0, 0, 6);
    fname = "pressure_" + tstr + ".vts";
    P.GetNodalValues(P_nodal);
    write_vts_scalar(fname, "pressure", grid.sx, grid.sy, grid.nx, grid.ny, P_nodal);
    fname = "velocity_" + tstr + ".vts";
    V.GetNodalValues(Vx_nodal, 1);
    V.GetNodalValues(Vy_nodal, 2);
    write_vts_vector(fname, "velocity", grid.sx, grid.sy, grid.nx, grid.ny, Vx_nodal, Vy_nodal);
    fname = "saturation_" + tstr + ".vts";
    S.GetNodalValues(S_nodal);
    write_vts_scalar(fname, "saturation", grid.sx, grid.sy, grid.nx, grid.ny, S_nodal);

    cout << "time step " << t+1 << " is done" << endl;
  }

  return 0;
}
