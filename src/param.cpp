#include "param.hpp"

using namespace std;
using namespace mfem;

#if defined(TWO_PHASE_FLOW)
double Krw(double S) { return S*S; }
double Kro(double S) { return (1.0-S)*(1.0-S); }
#else
double Krw(double S) { return S; }
double Kro(double S) { return (1.0-S); }
#endif // TWO_PHASE_FLOW

Param::Param()
  : spacedim(2)
  , nx(0), ny(0), nz(0), n_cells(0)
  , sx(0.0), sy(0.0), sz(0.0), V(0.0)
  , K_array(nullptr)
  , Q_array(nullptr)
  , por_array(nullptr)
  , R_array(nullptr)
  , order_v(1), order_p(1), order_s(0)
  , t_final(200), dt(100)
  , vis_steps(100)
{ }

Param::~Param()
{
  delete[] R_array;
  delete[] por_array;
  delete[] Q_array;
  delete[] K_array;
}

void Param::init_arrays()
{
  n_cells = nx*ny;
  if (spacedim == 3) n_cells *= nz;
  double V_domain = sx*sy;
  if (spacedim == 3) V_domain *= sz;
  V = V_domain / n_cells;

  K_array   = new double[n_cells];
  Q_array   = new double[n_cells];
  por_array = new double[n_cells];
  R_array   = new double[n_cells];
  for (int i = 0; i < n_cells; ++i)
  {
    K_array[i]   = 1.0;
    Q_array[i]   = 0.0;
    por_array[i] = 1.0;
    R_array[i]   = 0.0;
  }
  Q_array[0]         =  1.0; // injection well
  Q_array[n_cells-1] = -1.0; // production well
  R_array[0]         =  1.0; // source of saturation
}

void Param::add_options(OptionsParser &args)
{
  args.AddOption(&spacedim, "-d", "--dim", "Space dimension of the problem");
  args.AddOption(&nx, "-nx", "--nx", "Number of cells in x-direction");
  args.AddOption(&ny, "-ny", "--ny", "Number of cells in y-direction");
  args.AddOption(&nz, "-nz", "--nz", "Number of cells in z-direction");
  args.AddOption(&sx, "-sx", "--sx", "Size of domain in x-direction");
  args.AddOption(&sy, "-sy", "--sy", "Size of domain in y-direction");
  args.AddOption(&sz, "-sz", "--sz", "Size of domain in z-direction");
  args.AddOption(&order_v, "-ov", "--orderv", "Order (degree) of the finite elements for velocity");
  args.AddOption(&order_p, "-op", "--orderp", "Order (degree) of the finite elements for pressure");
  args.AddOption(&order_s, "-os", "--orders", "Order (degree) of the finite elements for saturation");
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0.");
  args.AddOption(&dt, "-dt", "--time-step", "Time step.");
  args.AddOption(&vis_steps, "-vs", "--vis-steps", "Visualize every n-th timestep.");
}
