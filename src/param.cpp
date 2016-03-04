#include "param.hpp"

using namespace std;
using namespace mfem;


double Krw(double S, bool two_phase_flow)
{
  if (two_phase_flow)
    return S*S;
  else
    return S;
}

double Kro(double S, bool two_phase_flow)
{
  if (two_phase_flow)
    return (1.0-S)*(1.0-S);
  else
    return (1.0-S);
}



DarcySolverParam::DarcySolverParam()
  : maxiter(1000)
  , rtol(1e-6)
  , atol(1e-10)
  , print_level(0)
{ }

void DarcySolverParam::add_options(OptionsParser &args)
{
  args.AddOption(&maxiter, "-darcy-maxiter", "--darcy-maxiter", "Max number of iterations for the Darcy solver");
  args.AddOption(&rtol, "-darcy-rtol", "--darcy-rtol", "Rel tolerance for the Darcy solver");
  args.AddOption(&atol, "-darcy-atol", "--darcy-atol", "Abs tolerance for the Darcy solver");
  args.AddOption(&print_level, "-darcy-print-level", "--darcy-print-level", "Print level for the Darcy solver");
}





Param::Param()
  : spacedim(2)
  , nx(100), ny(100), nz(100)
  , sx(100.), sy(100.), sz(100.)
  , K_array(nullptr)
  , Q_array(nullptr)
  , phi_array(nullptr)
  , R_array(nullptr)
  , order_v(1), order_p(1), order_s(0)
  , t_final(200), dt_global(100), dt_local(1)
  , vis_steps_global(1)
  , vis_steps_local(-1) // negative means no output
  , seis_steps(1)
  , K(1.0)              // permeability
  , phi(1.0)            // porosity
  , K_file("no-file")
  , phi_file("no-file")
  , outdir("output")
  , extra("")
  , two_phase_flow(false)
  , ode_solver_type(1) // Forward Euler
  , info(false)
{ }

Param::~Param()
{
  delete[] R_array;
  delete[] phi_array;
  delete[] Q_array;
  delete[] K_array;
}

void Param::init_arrays()
{
  const int n_cells = get_n_cells();

  K_array   = new double[n_cells];
  Q_array   = new double[n_cells];
  phi_array = new double[n_cells];
  R_array   = new double[n_cells];

  if (!strcmp(K_file, "no-file"))
    for (int i = 0; i < n_cells; ++i) K_array[i] = K;
  else
    read_binary(K_file, n_cells, K_array);

  if (!strcmp(phi_file, "no-file"))
    for (int i = 0; i < n_cells; ++i) phi_array[i] = phi;
  else
    read_binary(phi_file, n_cells, phi_array);

  // sources and sinks
  for (int i = 0; i < n_cells; ++i)
  {
    Q_array[i]   = 0.0;
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
  args.AddOption(&t_final, "-tf", "--t-final", "Final time; start time is 0");
  args.AddOption(&dt_global, "-dtg", "--time-step-global", "Time step for pressure solver");
  args.AddOption(&dt_local, "-dtl", "--time-step-local", "Time step for saturation solver");
  args.AddOption(&vis_steps_global, "-vsg", "--vis-steps-global", "Visualize every n-th timestep in the pressure time loop (<=0 no output)");
  args.AddOption(&vis_steps_local, "-vsl", "--vis-steps-local", "Visualize every n-th timestep in the saturation time loops (<=0 no output)");
  args.AddOption(&seis_steps, "-ss", "--seis-steps", "Compute seismic properties with Gassmann and output them every n-th global timestep (<=0 no output)");
  args.AddOption(&K, "-K", "--K", "Constant permeability");
  args.AddOption(&phi, "-phi", "--phi", "Constant porosity");
  args.AddOption(&K_file, "-K-file", "--K-file", "Name of a binary file (single precision) for heterogeneous permeability");
  args.AddOption(&phi_file, "-phi-file", "--phi-file", "Name of a binary file (single precision) for heterogeneous porosity");
  args.AddOption(&outdir, "-outdir", "--output-directory", "Directory for output of the results");
  args.AddOption(&extra, "-extra", "--extra", "Extra string to distinguish output files");
  args.AddOption(&two_phase_flow, "-two-phase", "--two-phase-flow",
                 "-single-phase", "--single-phase-flow",
                 "Simulate two phase flow (otherwise it's single phase)");
  args.AddOption(&ode_solver_type, "-ode", "--ode-solver-type", "ODE solver for saturation (1=Forward Euler, 2=RK2, 3=RK3, 4=RK4, 6=RK6 ");

  darcy.add_options(args);

  args.AddOption(&info, "-info", "--info", "-no-info", "--no-info",
                 "Show info about the program and exit");
}

string Param::get_info() const
{
  const string str =
      "\nBuild type                 : " + BUILD_TYPE +
      "\nMFEM directory             : " + MFEM_DIR +
      "\nMFEM library               : " + MFEM_LIB +
      "\nHypre directory            : " + HYPRE_DIR +
      "\nHypre library              : " + HYPRE_LIB +
      "\nMetis directory            : " + METIS_DIR +
      "\nMetis library              : " + METIS_LIB +
      "\nConfig date and time (UTC) : " + CONFIG_TIME +
      "\nUser name                  : " + USER_NAME +
      "\nHost name                  : " + HOST_NAME +
      "\nSystem name                : " + SYSTEM_NAME +
      "\nSystem processor           : " + SYSTEM_PROCESSOR +
      "\nGit branch                 : " + GIT_BRANCH +
      "\nGit commit hash            : " + GIT_COMMIT_HASH +
      "\n";
  return str;
}
