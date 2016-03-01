#include "param.hpp"
#include "saturation_solver.hpp"
#include "utilities.hpp"

using namespace std;
using namespace mfem;

#if defined(MFEM_USE_MPI) // parallel mode

void FS(const Vector& S, Vector& fs)
{
  MFEM_VERIFY(S.Size() == fs.Size(), "Dimensions mismatch");
  const bool two_phase_flow = true; // this is only valid for two phase flow
  for (int i = 0; i < S.Size(); ++i)
  {
    const double mw = Krw(S(i), two_phase_flow) / MU_W;
    const double mo = Kro(S(i), two_phase_flow) / MU_O;
    fs(i) = mw / (mw + mo);
  }
}

/** A time-dependent operator for the right-hand side of the ODE. The DG weak
    form of du/dt = v.grad(u) is M du/dt = K u + b, where M and K are the mass
    and advection matrices, and b describes the flow on the boundary. This can
    be written as a general ODE, du/dt = M^{-1} (K u + b), and this class is
    used to evaluate the right-hand side. */
class FE_Evolution : public TimeDependentOperator
{
public:
  FE_Evolution(HypreParMatrix &_M, HypreParMatrix &_K, const Vector &_b,
               bool two_phase)
    : TimeDependentOperator(_M.Height()),
      M(_M), K(_K), b(_b), M_solver(M.GetComm()), z(_M.Height()),
      two_phase_flow(two_phase)
  {
    M_prec.SetType(HypreSmoother::Jacobi);
    M_solver.SetPreconditioner(M_prec);
    M_solver.SetOperator(M);

    M_solver.iterative_mode = false;
    M_solver.SetRelTol(1e-9);
    M_solver.SetAbsTol(0.0);
    M_solver.SetMaxIter(100);
    M_solver.SetPrintLevel(0);
  }

  void Mult(const Vector &x, Vector &y) const
  {
    if (two_phase_flow)
    {
      Vector fs(x.Size());
      FS(x, fs);
      K.Mult(fs, z); // S = S + dt M^-1( K F(S) + b)
    }
    else
    {
      // S = S + dt M^-1( K S + b)
      // y = M^{-1} (K x + b)
      K.Mult(x, z);
    }
    z += b;
    M_solver.Mult(z, y);
  }

  virtual ~FE_Evolution() { }

private:
   HypreParMatrix &M, &K;
   const Vector &b;
   HypreSmoother M_prec;
   CGSolver M_solver;
   bool two_phase_flow;

   mutable Vector z;
};



void ParSaturationSolver(const Param &param, ParGridFunction &S,
                         VectorCoefficient &velocity, int global_ti,
                         double global_dt, VisItDataCollection &visit)
{
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  ParFiniteElementSpace &S_space = *(S.ParFESpace());

  const bool own_array = false;
  CWConstCoefficient Phi(param.phi_array, param.get_n_cells(), own_array);

  ParBilinearForm m(&S_space);
  m.AddDomainIntegrator(new MassIntegrator(Phi));
  ParBilinearForm k(&S_space);
  k.AddDomainIntegrator(new TransposeIntegrator(new ConvectionIntegrator(velocity, 1.0)));
  k.AddInteriorFaceIntegrator(new DGTraceIntegrator(velocity, -1.0, -0.5));

//  CWConstCoefficient R(param.R_array, param.n_cells, own_array);

  Vector R_array(param.nx*param.ny);
  R_array = 0.;
  R_array(0) = 1.;
  PWConstCoefficient R(R_array);

  ParLinearForm b(&S_space);
  b.AddDomainIntegrator(new DomainLFIntegrator(R));

  m.Assemble();
  m.Finalize();
  int skip_zeros = 0;
  k.Assemble(skip_zeros);
  k.Finalize(skip_zeros);
  b.Assemble();

  HypreParMatrix *M = m.ParallelAssemble();
  HypreParMatrix *K = k.ParallelAssemble();
  HypreParVector *B = b.ParallelAssemble();

  ODESolver *ode_solver = NULL;
  int ode_solver_type = 1;
  switch (ode_solver_type)
  {
    case 1: ode_solver = new ForwardEulerSolver; break;
    case 2: ode_solver = new RK2Solver(1.0); break;
    case 3: ode_solver = new RK3SSPSolver; break;
    case 4: ode_solver = new RK4Solver; break;
    case 6: ode_solver = new RK6Solver; break;
    default: throw runtime_error("AAA!");
  }

  HypreParVector *SV = S.GetTrueDofs();

  FE_Evolution adv(*M, *K, *B, param.two_phase_flow);
  ode_solver->Init(adv);

  double t = 0.0;
  double dt = param.dt_local;

  if (myid == 0)
  {
    std::cout << "t_final = " << global_dt << std::endl;
    std::cout << "dt      = " << dt << std::endl;
    std::cout << "SV norm = " << SV->Norml2() << std::endl;
    std::cout << "t       = " << t << std::endl;
  }


  for (int ti = 0; true; )
  {
    if (t >= global_dt - dt/2)
       break;

    ode_solver->Step(*SV, t, dt);
    ti++;

    if (param.vis_steps_local > 0 && ti % param.vis_steps_local == 0)
    {
       if (myid == 0)
         cout << "time step: " << ti << ", time: " << t << endl;

       S = *SV;

       visit.SetCycle(global_ti + ti);
       visit.SetTime(global_ti*global_dt + t);
       visit.Save();
    }
  }

  S = *SV;

  delete SV;
  delete ode_solver;
  delete B;
  delete M;
  delete K;
}

#endif // MFEM_USE_MPI
