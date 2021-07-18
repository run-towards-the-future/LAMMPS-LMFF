#ifndef PAIR_LMFF_SW64_H_
#define PAIR_LMFF_SW64_H_
#ifdef __cplusplus
extern "C"{
#endif

#define NEIGHMASK 0x3FFFFFFF

#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG)
//#define LAMMPS_SMALLBIG
#define LAMMPS_BIGBIG
#endif

#ifdef LAMMPS_SMALLBIG
typedef int tagint;
typedef int64_t bigint;
#endif

#ifdef LAMMPS_BIGBIG
typedef int64_t tagint;
typedef int64_t bigint;
#endif

#ifdef LAMMPS_SMALLSMALL
typedef int tagint;
typedef int bigint;
#endif

#include "pair_lmff_parameters.h"

/* Buffer neighbor info for Tersoff*/
typedef struct ParamBf
{//buffer atom info for neighshort;
	double r, rsq, rinv, delr[3];
	int /*id,*/ tp;
	double padding;
	int ters_flag;
}ParamBf;


typedef struct atom_in
{
  double x[3];
  int type;
	/* padding into 64Bytes */
  int padding;
  double pad20, pad21;

	/* tagint */
	tagint molecule; 
  tagint tag;
}atom_in;


typedef struct pair_lmff_param_t
{
	double (*frc)[4];//wcache;
	long *mylock;

	double (*x)[3], (*f)[3];
	double cutshortsq;
	double eng_vdwl, eng_coul, virial[6], pvector[2];
  double (*vatom)[6], *eatom, *pack64;
	double *ilp_cutILPsq;//nelements^2;
	double *ilp_cutsq;//(ntypes+1)^2;
	
	Param ters_params[8];
	Param_ilp_str ilp_params[8];

	atom_in *my_atoms;
	
	tagint *tag, *molecule;
	bigint ntimestep, lastcall;
	int *type, *map;
	int *elem2param;//nelements^3;
	int nlocal, nghost, newton_pair;
	int ntypes, nelements, nparams, maxshort;
	int *iskip, *ijskip;//!!!
  int inum, *ilist,*numneigh,**firstneigh, gnum;
	int *nummol;
	
	int ilp_nelements, ilp_nparams, ilp_tap_flag;
	int *ilp_elem2param;//nelements^2;
	int *ilp_map;

	int eflag, vflag, evflag, eflag_either, vflag_either;
	int eflag_global, vflag_global;
	int eflag_atom, vflag_atom;
	
	int myrank;

	/* filter shortneighbors */
	int *flt_shortneigh;
	
	/* recalculate */
	int recalc;
}pair_lmff_param_t;


#ifdef __cplusplus
}
#endif
#endif
