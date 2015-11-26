#include "config.hpp"
#include "mfem.hpp"
#include <stdexcept>

using namespace std;

void run_parallel(int argc, char **argv);
void run_serial(int argc, char **argv);

int main(int argc, char **argv)
{
#if defined(MFEM_USE_MPI) // parallel mode
  MPI_Init(&argc, &argv);
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  try
  {
    run_parallel(argc, argv);
  }
  catch (int i)
  {
    MPI_Finalize();
    return i;
  }
  catch (const exception& e)
  {
    cout << "ID " << myid << ": Exception\n" << e.what() << endl;
    MPI_Finalize();
    return 1;
  }
  catch (...)
  {
    cout << "ID " << myid << ": Unknown exception" << endl;
    MPI_Finalize();
    return 1;
  }
  MPI_Finalize();
  return 0;
#else // serial mode
  try
  {
    run_serial(argc, argv);
  }
  catch (int i)
  {
    return i;
  }
  catch (const exception& e)
  {
    cout << "Exception\n" << e.what() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "Unknown exception" << endl;
    return 1;
  }
  return 0;
#endif // MFEM_USE_MPI
}
