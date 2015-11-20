#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>

// MFEM config file
#include "/u/artemyev/projects/mfem/serial/config.hpp"

#ifndef nullptr
  #define nullptr NULL
#endif

const std::string SEISMOGRAMS_DIR = "seismograms/";
const std::string SNAPSHOTS_DIR   = "snapshots/";

const double VERY_SMALL_NUMBER = 1e-32;
const double FLOAT_NUMBERS_EQUALITY_TOLERANCE = 1e-12;
const double FLOAT_NUMBERS_EQUALITY_REDUCED_TOLERANCE = 1e-6;
const double FIND_CELL_TOLERANCE = 1e-12;

// Water properties
const double MU_W  = 1.0;    // Viscosity
const double VP_W  = 1500.0; // P-wave velocity
const double VS_W  = 0.0;    // S-wave velocity
const double RHO_W = 1000.0; // Density

// Oil properties
const double MU_O  = 1.0;    // Viscosity
const double VP_O  = 1200.0; // P-wave velocity
const double VS_O  = 0.0;    // S-wave velocity
const double RHO_O = 800.0;  // Density

const double K_QUARTZ = 37e+9;
const double G_QUARTZ = 44e+9;
const double RHO_QUARTZ = 2650;

const double K_WET_CLAY = 15e+9;

const double K_MINERAL_MATRIX = K_QUARTZ; // bulk modulus of the mineral matrix
const double G_MINERAL_MATRIX = G_QUARTZ; // shear modulus of the mineral matrix
const double RHO_GRAIN = RHO_QUARTZ; // grain density of the rock matrix

const double K_FLUID_COMPONENT = K_WET_CLAY;

const double F_MINERAL_MATRIX  = 0.9;
const double F_FLUID_COMPONENT = 1.0 - F_MINERAL_MATRIX;

#endif // CONFIG_HPP
