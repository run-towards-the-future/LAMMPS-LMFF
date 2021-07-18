/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(ilp/graphene/hbn,PairILPGrapheneHBN)

#else

#ifndef LMP_PAIR_ILP_GRAPHENE_HBN_H
#define LMP_PAIR_ILP_GRAPHENE_HBN_H

#include "pair.h"

namespace LAMMPS_NS {

class PairILPGrapheneHBN : public Pair {
 public:
  PairILPGrapheneHBN(class LAMMPS *);
  virtual ~PairILPGrapheneHBN();

  virtual void compute(int, int);
  void compute_ilp(int, int);///
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void ILP_neigh();
  void calc_normal();
  void calc_FRep(int, int);
  void calc_FvdW(int, int);
  double single(int, int, int, int, double, double, double, double &);

 protected:
  class NeighList *list_ilp;//added by ping;

  int me;
  int maxlocal;                    // size of numneigh, firstneigh arrays
  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom
  MyPage<int> *ipage;              // neighbor list pages
  int *ILP_numneigh;                // # of pair neighbors for each atom
  int **ILP_firstneigh;             // ptr to 1st neighbor of each atom
  int tap_flag;                    // flag to turn on/off taper function

  struct Param {
    double z0,alpha,epsilon,C,delta,d,sR,reff,C6,S;
    double delta2inv,seff,lambda,rcut;
    int ielement,jelement;
  };
  Param *params;       // parameter set for I-J interactions
  char **elements;     // names of unique elements
  int **elem2param;    // mapping from element pairs to parameters
  int *map;            // mapping from atom types to elements
  int nelements;       // # of unique elements
  int nparams;         // # of stored parameter sets
  int maxparam;        // max # of parameter sets
  int nmax;            // max # of atoms

  double cut_global;
  double cut_normal;
  double **cut;
  double **cutILPsq;    // mapping the cutoff from element pairs to parameters
  double **offset;
  double **normal;
  double ***dnormdri;
  double ****dnormal;

  void read_file( char * );
  void allocate();

	inline void ev_tally_buffer(int i, int j, double *ei, double *vi, double *ej, double *vj, 
														int nlocal, int newton_pair,
                            double evdwl, double ecoul,
                            double fx, double fy, double fz,
                            double delx, double dely, double delz)
	{
		double evdwlhalf,ecoulhalf,epairhalf,v[6];
  	//printf("%p %p %p %p\n", ei, ej, vi, vj);
  	if (eflag_either) {
  	  if (eflag_global) {
  	    if (newton_pair) {
  	      eng_vdwl += evdwl;
  	      eng_coul += ecoul;
  	    } else {
  	      evdwlhalf = 0.5*evdwl;
  	      ecoulhalf = 0.5*ecoul;
  	      if (i < nlocal) {
  	        eng_vdwl += evdwlhalf;
  	        eng_coul += ecoulhalf;
  	      }
  	      if (j < nlocal) {
  	        eng_vdwl += evdwlhalf;
  	        eng_coul += ecoulhalf;
  	      }
  	    }
  	  }
  	  if (eflag_atom) {
  	    epairhalf = 0.5 * (evdwl + ecoul);
  	    if (newton_pair || i < nlocal) *ei += epairhalf;
  	    if (newton_pair || j < nlocal) *ej += epairhalf;
  	  }
  	}

  	if (vflag_either) {
  	  v[0] = delx*fx;
  	  v[1] = dely*fy;
  	  v[2] = delz*fz;
  	  v[3] = delx*fy;
  	  v[4] = delx*fz;
  	  v[5] = dely*fz;

  	  if (vflag_global) {
  	    if (newton_pair) {
  	      virial[0] += v[0];
  	      virial[1] += v[1];
  	      virial[2] += v[2];
  	      virial[3] += v[3];
  	      virial[4] += v[4];
  	      virial[5] += v[5];
  	    } else {
  	      if (i < nlocal) {
  	        virial[0] += 0.5*v[0];
  	        virial[1] += 0.5*v[1];
  	        virial[2] += 0.5*v[2];
  	        virial[3] += 0.5*v[3];
  	        virial[4] += 0.5*v[4];
  	        virial[5] += 0.5*v[5];
  	      }
  	      if (j < nlocal) {
  	        virial[0] += 0.5*v[0];
  	        virial[1] += 0.5*v[1];
  	        virial[2] += 0.5*v[2];
  	        virial[3] += 0.5*v[3];
  	        virial[4] += 0.5*v[4];
  	        virial[5] += 0.5*v[5];
  	      }
  	    }
  	  }

  	  if (vflag_atom) {
  	    if (newton_pair || i < nlocal) {
  	      vi[0] += 0.5*v[0];
  	      vi[1] += 0.5*v[1];
  	      vi[2] += 0.5*v[2];
  	      vi[3] += 0.5*v[3];
  	      vi[4] += 0.5*v[4];
  	      vi[5] += 0.5*v[5];
  	    }
  	    if (newton_pair || j < nlocal) {
  	      vj[0] += 0.5*v[0];
  	      vj[1] += 0.5*v[1];
  	      vj[2] += 0.5*v[2];
  	      vj[3] += 0.5*v[3];
  	      vj[4] += 0.5*v[4];
  	      vj[5] += 0.5*v[5];
  	    }
  	  }
  	}
	}

	inline void cross_deriv(double *pv, double (*dpv)[3][3], double (*vet)[3], int j, int k, int l) 
	{
	  pv[0] = vet[j][1]*vet[k][2] - vet[k][1]*vet[j][2];
	  pv[1] = vet[j][2]*vet[k][0] - vet[k][2]*vet[j][0];
	  pv[2] = vet[j][0]*vet[k][1] - vet[k][0]*vet[j][1];
	
	  dpv[j][0][0] =  0.0;
	  dpv[j][0][1] =  vet[k][2];
	  dpv[j][0][2] = -vet[k][1];
	  dpv[j][1][0] = -vet[k][2];
	  dpv[j][1][1] =  0.0;
	  dpv[j][1][2] =  vet[k][0];
	  dpv[j][2][0] =  vet[k][1];
	  dpv[j][2][1] = -vet[k][0];
	  dpv[j][2][2] =  0.0;
	
	  dpv[k][0][0] =  0.0;
	  dpv[k][0][1] = -vet[j][2];
	  dpv[k][0][2] =  vet[j][1];
	  dpv[k][1][0] =  vet[j][2];
	  dpv[k][1][1] =  0.0;
	  dpv[k][1][2] = -vet[j][0];
	  dpv[k][2][0] = -vet[j][1];
	  dpv[k][2][1] =  vet[j][0];
	  dpv[k][2][2] =  0.0;
	
	  dpv[l][0][0] =  0.0;
	  dpv[l][0][1] =  0.0;
	  dpv[l][0][2] =  0.0;
	  dpv[l][1][0] =  0.0;
	  dpv[l][1][1] =  0.0;
	  dpv[l][1][2] =  0.0;
	  dpv[l][2][0] =  0.0;
	  dpv[l][2][1] =  0.0;
	  dpv[l][2][2] =  0.0;
	}
	
	inline void calc_dnormal(double (*dnormal)[3][3], double (*dn1)[3][3], double *n1, double nn, double nn2)
	{
	  double dnn[3][3];
	  int m, id, ip;
	  // derivatives of nn, dnn:3x3 vector
	  // dnn[id][m]: the derivative of nn respect to r[id][m], id,m=0,1,2
	  // r[id][m]: the id's component of atom m
	  for (m = 0; m < 3; m++){
	    for (id = 0; id < 3; id++){
	      dnn[id][m] = (n1[0]*dn1[m][0][id] + n1[1]*dn1[m][1][id] + n1[2]*dn1[m][2][id])/nn;
	    }
	  }
	  // dnormal[m][id][ip]: the derivative of normal[id] respect to r[ip][m], id,ip=0,1,2
	  // for atom m, which is a neighbor atom of atom i, m=0,jnum-1
	  for (m = 0; m < 3; m++){
	    for (id = 0; id < 3; id++){
	      for (ip = 0; ip < 3; ip++){
	        dnormal[m][id][ip] = dn1[m][id][ip]/nn - n1[id]*dnn[ip][m]/nn2;
	      }
	    }
	  }
	}
	/************************************************************************************/

  /* ----Calculate the long-range cutoff term */
  inline double calc_Tap(double r_ij, double Rcut) {
    double Tap,r;
    double Tap_coeff[8] = {1.0,0.0,0.0,0.0,-35.0,84.0,-70.0,20.0};

    r = r_ij/Rcut;
    if(r >= 1.0) {Tap = 0.0;}
    else {
      Tap = Tap_coeff[7] * r + Tap_coeff[6];
      Tap = Tap * r  + Tap_coeff[5];
      Tap = Tap * r  + Tap_coeff[4];
      Tap = Tap * r  + Tap_coeff[3];
      Tap = Tap * r  + Tap_coeff[2];
      Tap = Tap * r  + Tap_coeff[1];
      Tap = Tap * r  + Tap_coeff[0];
    }
    return(Tap);
	}

  /* ----Calculate the derivatives of long-range cutoff term */
  inline double calc_dTap(double r_ij, double Rcut) {
    double dTap,r;
    double Tap_coeff[8] = {1.0,0.0,0.0,0.0,-35.0,84.0,-70.0,20.0};

    r = r_ij/Rcut;
    if(r >= 1.0) {dTap = 0.0;}
    else {
      dTap = 7.0*Tap_coeff[7] * r + 6.0*Tap_coeff[6];
      dTap = dTap * r  + 5.0*Tap_coeff[5];
      dTap = dTap * r  + 4.0*Tap_coeff[4];
      dTap = dTap * r  + 3.0*Tap_coeff[3];
      dTap = dTap * r  + 2.0*Tap_coeff[2];
      dTap = dTap * r  + Tap_coeff[1];
      dTap = dTap/Rcut;
    }
    return(dTap);
  }
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/

