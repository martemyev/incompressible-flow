#include "param.hpp"

using namespace std;
using namespace mfem;



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




Well::Well()
  : center(0., 0., 0.)
  , radius(1.)
  , height(1.)
{}




Param::Param()
  : spacedim(2)
  , nx(100), ny(100), nz(100)
  , sx(100.), sy(100.), sz(100.)
  , K_array_x(nullptr)
  , K_array_y(nullptr)
  , K_array_z(nullptr)
  , phi_array(nullptr)
  , order_v(1), order_p(1), order_s(0)
  , t_final(200), dt_global(100), dt_local(1)
  , vis_steps_global(1)
  , vis_steps_local(-1) // negative means no output
  , saturation_steps(1)
  , K_x(1.0)            // permeability (x-component of a tensor)
  , K_y(1.0)            // permeability (y-component of a tensor)
  , K_z(1.0)            // permeability (z-component of a tensor)
  , phi(1.0)            // porosity
  , K_file_x("no-file")
  , K_file_y("no-file")
  , K_file_z("no-file")
  , phi_file("no-file")
  , outdir("output")
  , extra("")
  , two_phase_flow(false)
  , ode_solver_type(1) // Forward Euler
  , injection()
  , production()
  , inflow(1.)
  , outflow(-1.)
  , saturation_source(1.)
  , Sor(0.)
  , Swc(0.)
  , info(false)
{ }

Param::~Param()
{
  delete[] phi_array;
  delete[] K_array_z;
  delete[] K_array_y;
  delete[] K_array_x;
}

void Param::init_arrays()
{
  const int n_cells = get_n_cells();

  K_array_x = new double[n_cells];
  K_array_y = new double[n_cells];
  if (spacedim == 3)
    K_array_z = new double[n_cells];
  phi_array = new double[n_cells];

  if (!strcmp(K_file_x, "no-file"))
    for (int i = 0; i < n_cells; ++i) K_array_x[i] = K_x;
  else
    read_binary(K_file_x, n_cells, K_array_x);

  if (!strcmp(K_file_y, "no-file"))
    for (int i = 0; i < n_cells; ++i) K_array_y[i] = K_y;
  else
    read_binary(K_file_y, n_cells, K_array_y);

  if (spacedim == 3)
  {
    if (!strcmp(K_file_z, "no-file"))
      for (int i = 0; i < n_cells; ++i) K_array_z[i] = K_z;
    else
      read_binary(K_file_z, n_cells, K_array_z);
  }

  if (!strcmp(phi_file, "no-file"))
    for (int i = 0; i < n_cells; ++i) phi_array[i] = phi;
  else
    read_binary(phi_file, n_cells, phi_array);
  for (int i = 0; i < n_cells; ++i)
    if (phi_array[i] < 1e-3)
      phi_array[i] = 1e-3;
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
  args.AddOption(&saturation_steps, "-ss", "--saturation-steps", "Output saturation solution as a binary file every n-th global timestep (<=0 no output)");
  args.AddOption(&K_x, "-Kx", "--Kx", "Constant permeability (x-component of a diagonal tensor)");
  args.AddOption(&K_y, "-Ky", "--Ky", "Constant permeability (y-component of a diagonal tensor)");
  args.AddOption(&K_z, "-Kz", "--Kz", "Constant permeability (z-component of a diagonal tensor)");
  args.AddOption(&phi, "-phi", "--phi", "Constant porosity");
  args.AddOption(&K_file_x, "-K-file-x", "--K-file-x", "Name of a binary file (single precision) for heterogeneous permeability (x-component of a tensor)");
  args.AddOption(&K_file_y, "-K-file-y", "--K-file-y", "Name of a binary file (single precision) for heterogeneous permeability (y-component of a tensor)");
  args.AddOption(&K_file_z, "-K-file-z", "--K-file-z", "Name of a binary file (single precision) for heterogeneous permeability (z-component of a tensor)");
  args.AddOption(&phi_file, "-phi-file", "--phi-file", "Name of a binary file (single precision) for heterogeneous porosity");
  args.AddOption(&outdir, "-outdir", "--output-directory", "Directory for output of the results");
  args.AddOption(&extra, "-extra", "--extra", "Extra string to distinguish output files");
  args.AddOption(&two_phase_flow, "-two-phase", "--two-phase-flow",
                 "-single-phase", "--single-phase-flow",
                 "Simulate two phase flow (otherwise it's single phase)");
  args.AddOption(&ode_solver_type, "-ode", "--ode-solver-type", "ODE solver for saturation (1=Forward Euler, 2=RK2, 3=RK3, 4=RK4, 6=RK6)");

  args.AddOption(&injection.center(0), "-inj-xcen", "--inj-center-x", "x-coordinate of center of well");
  args.AddOption(&injection.center(1), "-inj-ycen", "--inj-center-y", "y-coordinate of center of well");
  args.AddOption(&injection.center(2), "-inj-zcen", "--inj-center-z", "z-coordinate of center of well");
  args.AddOption(&injection.radius, "-inj-r", "--inj-radius", "radius of well");
  args.AddOption(&injection.height, "-inj-h", "--inj-height", "height of well");

  args.AddOption(&production.center(0), "-pro-xcen", "--pro-center-x", "x-coordinate of center of well");
  args.AddOption(&production.center(1), "-pro-ycen", "--pro-center-y", "y-coordinate of center of well");
  args.AddOption(&production.center(2), "-pro-zcen", "--pro-center-z", "z-coordinate of center of well");
  args.AddOption(&production.radius, "-pro-r", "--pro-radius", "radius of well");
  args.AddOption(&production.height, "-pro-h", "--pro-height", "height of well");

  args.AddOption(&inflow, "-inflow", "--inflow", "Inflow value (in injection well)");
  args.AddOption(&outflow, "-outflow", "--outflow", "Outflow value (in production well)");
  args.AddOption(&saturation_source, "-sat-source", "--saturation-source", "Value of saturation source");
  args.AddOption(&Sor, "-Sor", "--Sor", "Irreducible oil saturation");
  args.AddOption(&Swc, "-Swc", "--Swc", "Connate water saturation");

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



double Param::Krw(double S) const
{
  if (two_phase_flow)
  {
    double S_ = (S - Swc) / (1. - Sor - Swc);
    return S_*S_;
  }
  else
    return S;
}

double Param::Kro(double S) const
{
  if (two_phase_flow)
  {
    double S_ = (S - Swc) / (1. - Sor - Swc);
    return (1.-S_)*(1.-S_);
  }
  else
    return (1.-S);
}
