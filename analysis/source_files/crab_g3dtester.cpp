/* +++++++++++ Complex routines ++++++++++++++ */
/* This should work on any machine */
//#include "source_files/volya_complex.h"
//#define double_complex Complex

/* For g++, you need only */
// #include <complex.h>

/* For Borland C++ or for cxx on alpha vms, you need */
#define double_complex complex
#include <complex.h>
/* +++++++++++++++++++++++++++++++++++++++++++ */

/* If you wish to calc. the correlation function as a function of |q_2-q_1|/2
   rather than the relative momentum (the default), then you should
   define "REDUCED_MOM" */
//#define REDUCED_MOM

/* For the Random # generator */
long IDUM=-1234;

/* The max. number of impact parameters */
#define NBMAX 1
/* The max. # of phasespace pt.s saved for one impact parameter */
#define NPHASEMAX 5000

/* #define MIXED_PAIRS_FOR_DENOM */
#define NMAX_FOR_MIXING 10000

/* To smear momentum to model experimental accuracies */
//#define SMEAR_MOMENTA

#include "source_files/crab.h"
#include "interactions/crab_interaction_pipi.cpp" /* To define interactions,
						 choose/edit this file.
						 You also set the IDs of
						 particles whose phase-space
						 pt.s you wish to use here. */
#include "binnings/crab_bindefs_g3dtester.cpp"     /* Choose and edit to modify
						 binning. */
#include "source_files/crab_main.cpp"
#include "source_files/crab_prinput.cpp"      /* Edit to change input format */

#ifdef STRONG_INTERACTION
#include "source_files/crab_partwaveinit.cpp"
#endif

#ifdef COULOMB
#include "source_files/crab_coulomb.cpp"
#endif

#include "filters/crab_filter_acceptall.cpp" /* Choose/edit to set
						 exp. acceptance */

#ifdef SMEAR_MOMENTA
#include "smearings/crab_smearing_10mev.cpp" /* Choose/edit to set exp.
						 momentum smearing */
#endif

#include "source_files/crab_corrcalc.cpp"
#include "source_files/crab_misc.cpp"
#include "source_files/crab_random.cpp"
