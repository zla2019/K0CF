/* These first two numbers describe the number of ident numbers used to
   screen phase space pt.s from files that are to be used in the generation
   of correlation functions. These codes may be different than those of the
    particles you are investigating. For instance, you may wish to do pi-pi-
    correlations and use the phase space pt.s for all 3 species of pions */
#define N1TYPES 1
#define N2TYPES 1
/* These are the identifications. See pdg.lbl.gov/rpp/mcdata/all.mc for a
   list of identification numbers */
const int IDENT1[N1TYPES]={2212};
const int IDENT2[N2TYPES]={211};

#define MASS1 938.3
#define MASS2 139.58

/* Define if particles are identical */
//#define IDENTICAL
/* Turn off and on the Coulomb Interaction */
#define COULOMB
/* Turn off and on the Strong Interaction */
/* #define STRONG_INTERACTION */

#define Q1Q2 1

#define INTERACTION_WSYM 0.0
#define  INTERACTION_WANTI 0.0
#define INTERACTION_WNOSYM 1.0
/* fractions of symmetric and antisym weights of the various spin channels */

#define INTERACTION_DELK  0.5
/* spacing of mom. mesh for calc. of strong/coul. int. given in MeV/c 
   and mom. is reduced momentum (1/2 Qinv for m1=m2) */
/* these are used to define both the strong and coulomb meshes */
#define INTERACTION_NKMAX 100
/* number of momentum points in mesh for strong/coul. interaction */

/* For estimating wave-function squared with Breit-Wigner form, define the 
   parameters below */
#define BREIT_WIGNER */

#ifdef BREIT_WIGNER

#define N_BWRESONANCES 1
/* This is a "coalescence radius */
#define BWINFO_R 1.0
/* Make it a p-wave resonance with L=1 */
const int BWINFO_L[N_BWRESONANCES]={1};
const double BWINFO_J[N_BWRESONANCES]={1.5};
const double BWINFO_J1[N_BWRESONANCES]={0.5};
const double BWINFO_J2[N_BWRESONANCES]={0};
const double BWINFO_MINV[N_BWRESONANCES]={1232.0};
const double BWINFO_GAMMA[N_BWRESONANCES]={115.0};
const double BWINFO_BRANCHING[N_BWRESONANCES]={1.0};

#endif
