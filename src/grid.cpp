#include "grid.hpp"

#if defined(TWO_PHASE_FLOW)
double Krw(double S) { return S*S; }
double Kro(double S) { return (1.0-S)*(1.0-S); }
#else
double Krw(double S) { return S; }
double Kro(double S) { return (1.0-S); }
#endif // TWO_PHASE_FLOW
