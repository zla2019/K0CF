/* These first two numbers describe the number of ident numbers used to
   screen phase space pt.s from files that are to be used in the generation
   of correlation functions. These codes may be different than those of the
    particles you are investigating. For instance, you may wish to do pi-pi-
    correlations and use the phase space pt.s for all 3 species of pions */
#define N1TYPES 3
#define N2TYPES 3
/* These are the identifications. See pdg.lbl.gov/rpp/mcdata/all.mc for a
   list of identification numbers */
const int IDENT1[N1TYPES]={211,111,-211};
const int IDENT2[N2TYPES]={211,111,-211};

#define MASS1 139.58
#define MASS2 139.58

/* Define if particles are identical */
#define IDENTICAL
/* Turn off and on the Coulomb Interaction */
//#define COULOMB
/* Turn off and on the Strong Interaction */
//#define STRONG_INTERACTION

#define Q1Q2 1

#define INTERACTION_WSYM 1.0
#define INTERACTION_WANTI 0.0
#define INTERACTION_WNOSYM 0.0
/* fractions of symmetric and antisym weights of the various spin channels */

#define INTERACTION_DELK  1.0
/* spacing of mom. mesh for calc. of strong/coul. int. given in MeV/c
   and mom. is reduced momentum (1/2 Qinv for m1=m2) */
/* these are used to define both the strong and coulomb meshes */
#define INTERACTION_NKMAX 100
/* number of momentum points in mesh for strong/coul. interaction */

