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

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "pair_lmff.h"
#include <mpi.h>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairLMFF::PairLMFF(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
  map = NULL;

  maxshort = 10;
  neighshort = NULL;
  
	pvector = new double[nextra];

	/* filter shortneighbors */
	flt_shortneigh	= NULL;
	recalc = 0;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairLMFF::~PairLMFF()
{
  if (copymode) return;

  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(neighshort);
    delete [] map;

		/* filter shortneighbors */
    memory->destroy(flt_shortneigh);  //added by ping;
  }
  
	delete [] pvector;
}

/* ---------------------------------------------------------------------- */

void PairLMFF::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);
	/* Method(1):  not filtering shortneighbors  */
	//computeLMFFgeneral(eflag, vflag);
	
	/* Method(2): filtering shortneighbors  */
	if(flt_shortneigh == NULL)
	{
		memory->create(flt_shortneigh,	atom->nmax*4,	"pair:filter_shortneigh");
		memset(flt_shortneigh,	0, sizeof(int)*atom->nmax*4);
	}


	if(update->ntimestep == neighbor->lastcall)
	{
		computeLMFFTwice(eflag, vflag);
	}
	else
	{
		recalc = 0;
		computeLMFFOnce(eflag, vflag);
		if(recalc)
			computeLMFFTwice(eflag, vflag);
	}


	//computeTersoff(eflag, vflag);
	//computeILP(eflag, vflag);
}

void PairLMFF::computeLMFFgeneral(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,gnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3]/*,fk[3]*/;
  double zeta_ij,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;
  const double cutshortsq = cutmax*cutmax;

	int nt = atom->ntypes+1;
	int iskip[nt], ijskip[nt][nt];// When all is NULL;
	if(list->iskip == NULL && list->ijskip == NULL)
	{
		for(i = 0; i < nt; i++)
		{
			iskip[i] = 0;
			for(j = 0; j < nt; j++) ijskip[i][j] = 0;
		}
	}
	else
	{
		for(i = 0; i < nt; i++)
		{
			iskip[i] = neighbor->lists[0]->iskip[i];
			for(j = 0; j < nt; j++)
				ijskip[i][j] = neighbor->lists[0]->ijskip[i][j];
		}
	}

	int nlist = neighbor->nlist;
	int find_list = -1;
	for(i = 0; i < nlist; i++)
	{
		char *pairbuild = neighbor->pairnames[neighbor->lists[i]->pair_method-1];
		if(strcmp("full/bin/ghost", pairbuild) == 0)
		{
			inum				= neighbor->lists[i]->inum;
			gnum				= neighbor->lists[i]->gnum;
  		ilist				= neighbor->lists[i]->ilist;
  		numneigh		= neighbor->lists[i]->numneigh;
  		firstneigh	= neighbor->lists[i]->firstneigh;
			find_list = i;
			break;
		}
	}
	
	if(find_list == -1)
	{
		error->one(FLERR, "Can not find the right neighbor list");
	}


	/*** variables for ILP ***/
	pvector[0] = pvector[1] = 0.0;
  int l, ll;
  double erep, r,Rcut,r2inv,r6inv,r8inv;
	double Tap,dTap,Vilp,TSvdw,TSvdw2inv,fsum;
  erep = 0.0;
	//Vars for calc_normal
  int id,ip,m;
  double nn, nn2;
  double pv12[3],pv31[3],pv23[3],n1[3],dni[3],dnn[3][3],vet[3][3],dpvdri[3][3];
  double dn1[3][3][3],dpv12[3][3][3],dpv23[3][3][3],dpv31[3][3][3];
  double normal[3], dnormal[3][3][3], dnormdri[3][3];
  //Vars for calc_Frep
  double prodnorm1,fkcx,fkcy,fkcz;
  double rhosq1,exp0,exp1;
  double frho1,Erep,rdsq1,fpair1;
  double dprodnorm1[3] = {0.0, 0.0, 0.0};
  double fp1[3] = {0.0, 0.0, 0.0};
  double fprod1[3] = {0.0, 0.0, 0.0};
  double delkj[3] = {0.0, 0.0, 0.0};
  double fk[3] = {0.0, 0.0, 0.0};
  double fkk[3][3], ekk[3], vkk[3][6];
  //more for overflow
  int ilp_neigh[6];
	//transfer from force cls;
	int ilp_nelements					= force->ilp_nelements;       
  int ilp_nparams						= force->ilp_nparams;         
	int ilp_tap_flag					= force->ilp_tap_flag;
  double **ilp_cutILPsq			= force->ilp_cutILPsq;    
  double **ilp_cutsq				= force->ilp_cutsq;    
	Param_ilp_str *ilp_params = force->ilp_params;       
  int **ilp_elem2param			= force->ilp_elem2param;    
  int *ilp_map							= force->ilp_map;            
	/*** End variables for ILP ***/

	if(newton_pair == 0)
		error->one(FLERR, "The LMFF potential needs that newton_pair = on");
	
	double fxtmp,fytmp,fztmp;

  // loop over full neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) 
	{//ILP+Tersoff
		i = ilist[ii];
		xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    int itypem = ilp_map[type[i]];
    itag = tag[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int nilp = 0;
		
		//TERSOFF
    itype = map[type[i]];
    int numshort = 0;
    fxtmp = fytmp = fztmp = 0.0;

    for (jj = 0; jj < jnum; jj++) 
		{
      j = jlist[jj];
      j &= NEIGHMASK;
      int jtypem = ilp_map[type[j]];
      jtag = tag[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq != 0 && rsq < ilp_cutILPsq[itypem][jtypem] && atom->molecule[i] == atom->molecule[j])
			{
        ilp_neigh[nilp] = j;
        vet[nilp][0] = -delx;
        vet[nilp][1] = -dely;
        vet[nilp][2] = -delz;
        nilp ++;
      }

			//TERSOFF
      if(!iskip[type[i]] && !ijskip[type[i]][type[j]])
			{
				if (rsq < cutshortsq) 
				{
      	  neighshort[numshort++] = j;
      	  if (numshort >= maxshort) 
					{
      	    maxshort += maxshort/2;
      	    memory->grow(neighshort,maxshort,"pair:neighshort");
      	  }
      	}
			}//if-notSkip
    }//for-jj

		calc_normal(nilp, jnum, vet, normal, dnormal, dnormdri);
    
		for (kk = 0; kk < nilp; kk ++) 
		{
      fkk[kk][0] = 0.0;
      fkk[kk][1] = 0.0;
      fkk[kk][2] = 0.0;
      vkk[kk][0] = 0.0;
      vkk[kk][1] = 0.0;
      vkk[kk][2] = 0.0;
      vkk[kk][3] = 0.0;
      vkk[kk][4] = 0.0;
      vkk[kk][5] = 0.0;
      ekk[kk] = 0.0;
    }
    
		jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) 
		{
      j = jlist[jj];
      j &= NEIGHMASK;

      jtag = tag[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
		

			//TERSOFF
			int ters_flag = 1;
      int vdwflag = 1;
			if (itag > jtag) 
			{
    	  if ((itag+jtag) % 2 == 0) {ters_flag = 0; vdwflag = 0;}
    	} 
			else if (itag < jtag) 
			{
    	  if ((itag+jtag) % 2 == 1) {ters_flag = 0; vdwflag = 0;}
    	} 
			else 
			{
    	  if (x[j][2] < x[i][2]) {ters_flag = 0; vdwflag = 0;}
    	  if (x[j][2] == ztmp && x[j][1] < ytmp) {ters_flag = 0; vdwflag = 0;}
    	  if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) {ters_flag = 0; vdwflag = 0;}
    	}

      jtype = map[type[j]];
			r = sqrt(rsq);
      double rinv = 1.0 / r;
			
			if(ters_flag && !iskip[type[i]] && !ijskip[type[i]][type[j]])
			{
				iparam_ij = elem2param[itype][jtype][jtype];
				if (rsq < params[iparam_ij].cutsq)
				{
					ters_FRep(&params[iparam_ij],r,rinv,fpair,eflag,evdwl);

      		fxtmp += delx*fpair;
      		fytmp += dely*fpair;
      		fztmp += delz*fpair;
      		f[j][0] -= delx*fpair;
      		f[j][1] -= dely*fpair;
      		f[j][2] -= delz*fpair;

      		if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
				}
			}//if-notSkip

      if (rsq < ilp_cutsq[type[i]][type[j]] && atom->molecule[i] != atom->molecule[j]) 
			{
        int iparam_ij = ilp_elem2param[ilp_map[type[i]]][ilp_map[type[j]]];

        Param_ilp_str *p = &ilp_params[iparam_ij];

        if (ilp_tap_flag) 
				{
          Rcut = sqrt(ilp_cutsq[type[i]][type[j]]);
          Tap = calc_Tap(r,Rcut);
          dTap = calc_dTap(r,Rcut);
        } 
				else {Tap = 1.0; dTap = 0.0;}

        double fvdw = 0, evdw = 0;
        
        if (vdwflag) 
				{
					r2inv = rinv * rinv;
					r6inv = r2inv*r2inv*r2inv;
					r8inv = r6inv*r2inv;
					
					TSvdw = 1.0 + exp(-p->d*(r/p->seff - 1.0));
					TSvdw2inv = 1 / (TSvdw * TSvdw);//pow(TSvdw,-2.0);
					double vvdw = -p->C6*r6inv/TSvdw;
					
					// derivatives
					fpair = -6.0*p->C6*r8inv/TSvdw + p->C6*p->d/p->seff*(TSvdw-1.0)*TSvdw2inv*r8inv*r;
					fsum = fpair*Tap - vvdw*dTap*rinv;
					fvdw = fsum;
					evdw = vvdw * Tap;
					if(eflag)
          {
            pvector[0] += evdw;
          }
        }
				
				// Calculate the transverse distance
        prodnorm1 = normal[0]*delx + normal[1]*dely + normal[2]*delz;
        rhosq1 = rsq - prodnorm1*prodnorm1;  // rho_ij
        rdsq1 = rhosq1*p->delta2inv; // (rho_ij/delta)^2

        // store exponents
        exp0 = exp(-p->lambda*(r-p->z0));
        exp1 = exp(-rdsq1);

        frho1 = exp1*p->C;
        Erep = 0.5*p->epsilon + frho1;
        Vilp = exp0*Erep;

        // derivatives
        fpair  = p->lambda*exp0*rinv*Erep;
        fpair1 = 2.0*exp0*frho1*p->delta2inv;
        fsum = fpair + fpair1;

        // derivatives of the product of rij and ni, the result is a vector
        dprodnorm1[0] = dnormdri[0][0]*delx + dnormdri[1][0]*dely + dnormdri[2][0]*delz;
        dprodnorm1[1] = dnormdri[0][1]*delx + dnormdri[1][1]*dely + dnormdri[2][1]*delz;
        dprodnorm1[2] = dnormdri[0][2]*delx + dnormdri[1][2]*dely + dnormdri[2][2]*delz;
        fp1[0] = prodnorm1*normal[0]*fpair1;
        fp1[1] = prodnorm1*normal[1]*fpair1;
        fp1[2] = prodnorm1*normal[2]*fpair1;
        fprod1[0] = prodnorm1*dprodnorm1[0]*fpair1;
        fprod1[1] = prodnorm1*dprodnorm1[1]*fpair1;
        fprod1[2] = prodnorm1*dprodnorm1[2]*fpair1;

        fkcx = (delx*fsum - fp1[0])*Tap - Vilp*dTap*delx*rinv;
        fkcy = (dely*fsum - fp1[1])*Tap - Vilp*dTap*dely*rinv;
        fkcz = (delz*fsum - fp1[2])*Tap - Vilp*dTap*delz*rinv;

        //This should be no use because fkcx need a lot of variables
        double ftotx = fvdw * delx + fkcx;
        double ftoty = fvdw * dely + fkcy;
        double ftotz = fvdw * delz + fkcz;
        f[i][0] += ftotx - fprod1[0]*Tap;
        f[i][1] += ftoty - fprod1[1]*Tap;
        f[i][2] += ftotz - fprod1[2]*Tap;
        f[j][0] -= ftotx;
        f[j][1] -= ftoty;
        f[j][2] -= ftotz;
				
        for (kk = 0; kk < nilp; kk++) 
				{
          k = ilp_neigh[kk];
          dprodnorm1[0] = dnormal[kk][0][0]*delx + dnormal[kk][1][0]*dely + dnormal[kk][2][0]*delz;
          dprodnorm1[1] = dnormal[kk][0][1]*delx + dnormal[kk][1][1]*dely + dnormal[kk][2][1]*delz;
          dprodnorm1[2] = dnormal[kk][0][2]*delx + dnormal[kk][1][2]*dely + dnormal[kk][2][2]*delz;
          fk[0] = (-prodnorm1*dprodnorm1[0]*fpair1)*Tap;
          fk[1] = (-prodnorm1*dprodnorm1[1]*fpair1)*Tap;
          fk[2] = (-prodnorm1*dprodnorm1[2]*fpair1)*Tap;
          fkk[kk][0] += fk[0];
          fkk[kk][1] += fk[1];
          fkk[kk][2] += fk[2];
          delkj[0] = x[i][0] + vet[kk][0] - x[j][0];
          delkj[1] = x[i][1] + vet[kk][1] - x[j][1];
          delkj[2] = x[i][2] + vet[kk][2] - x[j][2];

          if (vflag_atom) 
          {
            if (evflag) 
							ev_tally_buffer(k, j, ekk + kk, vkk[kk], eatom + j, vatom[j], nlocal,newton_pair, 
															0.0,0.0,fk[0],fk[1],fk[2],delkj[0],delkj[1],delkj[2]);
          } 
          else 
          {
            if (evflag) 
							ev_tally_buffer(k, j, ekk + kk, vkk[kk], eatom + j, NULL, nlocal,newton_pair, 
															0.0,0.0,fk[0],fk[1],fk[2], delkj[0],delkj[1],delkj[2]);
          }
        }
        //erep = Tap*Vilp;
        if (eflag) pvector[1] += erep = Tap*Vilp;
        if (evflag) 
					ev_tally_xyz(i,j,nlocal,newton_pair,erep + evdw,0.0, ftotx,ftoty,ftotz,delx,dely,delz);
      }
    } // loop over jj
    for (kk = 0; kk < nilp; kk ++) 
		{
      int k = ilp_neigh[kk];
      f[k][0] += fkk[kk][0];
      f[k][1] += fkk[kk][1];
      f[k][2] += fkk[kk][2];
      if (eflag_atom) eatom[k] += ekk[kk];
      if (vflag_atom) 
			{
        vatom[k][0] += vkk[kk][0];
        vatom[k][1] += vkk[kk][1];
        vatom[k][2] += vkk[kk][2];
        vatom[k][3] += vkk[kk][3];
        vatom[k][4] += vkk[kk][4];
        vatom[k][5] += vkk[kk][5];
      }
    }//for-kk

		
		/*** Tersoff ***/
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // three-body interactions
    // skip immediately if I-J is not within cutoff
    double fjxtmp,fjytmp,fjztmp;
		double *dlr1, *dlr2, r1, r2, r1inv, r2inv;	
		//When iskip[type[i]] is ture, numshort==0;
		//			fxtmp == 0; fytmp == 0; fztmp == 0;
    for (jj = 0; jj < numshort; jj++) 
		{
      j = neighshort[jj];
      jtype = map[type[j]];
      iparam_ij = elem2param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 >= params[iparam_ij].cutsq) continue;
      
      // accumulate bondorder zeta for each i-j interaction via loop over k

      fjxtmp = fjytmp = fjztmp = 0.0;
      zeta_ij = 0.0;
      for (kk = 0; kk < numshort; kk++) 
			{
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];
        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[iparam_ijk].cutsq) continue;
        zeta_ij += zeta(&params[iparam_ijk],rsq1,rsq2,delr1,delr2);
      }

      // pairwise force due to zeta

      force_zeta(&params[iparam_ij],rsq1,zeta_ij,fpair,prefactor,eflag,evdwl);

      fxtmp += delr1[0]*fpair;
      fytmp += delr1[1]*fpair;
      fztmp += delr1[2]*fpair;
      fjxtmp -= delr1[0]*fpair;
      fjytmp -= delr1[1]*fpair;
      fjztmp -= delr1[2]*fpair;
      if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);

      // attractive term via loop over k
      for (kk = 0; kk < numshort; kk++) 
			{
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];
        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[iparam_ijk].cutsq) continue;
        attractive(&params[iparam_ijk],prefactor,rsq1,rsq2,delr1,delr2,fi,fj,fk);

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];
				
        if (vflag_atom) v_tally3(i,j,k,fj,fk,delr1,delr2);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
			
    }
		
		f[i][0] += fxtmp;
		f[i][1] += fytmp;
    f[i][2] += fztmp;
		
  }//for-ii

  if (vflag_fdotr) virial_fdotr_compute();
}


void PairLMFF::computeLMFFOnce(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,gnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3]/*,fk[3]*/;
  double zeta_ij,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;
  const double cutshortsq = cutmax*cutmax;

	int nt = atom->ntypes+1;
	int iskip[nt], ijskip[nt][nt];// When all is NULL;
	if(list->iskip == NULL && list->ijskip == NULL)
	{
		for(i = 0; i < nt; i++)
		{
			iskip[i] = 0;
			for(j = 0; j < nt; j++) ijskip[i][j] = 0;
		}
	}
	else
	{
		for(i = 0; i < nt; i++)
		{
			iskip[i] = neighbor->lists[0]->iskip[i];
			for(j = 0; j < nt; j++)
				ijskip[i][j] = neighbor->lists[0]->ijskip[i][j];
		}
	}

	int nlist = neighbor->nlist;
	int find_list = -1;
	for(i = 0; i < nlist; i++)
	{
		char *pairbuild = neighbor->pairnames[neighbor->lists[i]->pair_method-1];
		if(strcmp("full/bin/ghost", pairbuild) == 0)
		{
			inum				= neighbor->lists[i]->inum;
			gnum				= neighbor->lists[i]->gnum;
  		ilist				= neighbor->lists[i]->ilist;
  		numneigh		= neighbor->lists[i]->numneigh;
  		firstneigh	= neighbor->lists[i]->firstneigh;
			find_list = i;
			break;
		}
	}
	if(find_list == -1)
	{
		error->one(FLERR, "Can not find the right neighbor list");
	}


	/*** variables for ILP ***/
	pvector[0] = pvector[1] = 0.0;
  int l, ll;
  double erep, r,Rcut,r2inv,r6inv,r8inv;
	double Tap,dTap,Vilp,TSvdw,TSvdw2inv,fsum;
  erep = 0.0;
	//Vars for calc_normal
  int id,ip,m;
  double nn, nn2;
  double pv12[3],pv31[3],pv23[3],n1[3],dni[3],dnn[3][3],vet[3][3],dpvdri[3][3];
  double dn1[3][3][3],dpv12[3][3][3],dpv23[3][3][3],dpv31[3][3][3];
  double normal[3], dnormal[3][3][3], dnormdri[3][3];
  //Vars for calc_Frep
  double prodnorm1,fkcx,fkcy,fkcz;
  double rhosq1,exp0,exp1;
  double frho1,Erep,rdsq1,fpair1;
  double dprodnorm1[3] = {0.0, 0.0, 0.0};
  double fp1[3] = {0.0, 0.0, 0.0};
  double fprod1[3] = {0.0, 0.0, 0.0};
  double delkj[3] = {0.0, 0.0, 0.0};
  double fk[3] = {0.0, 0.0, 0.0};
  double fkk[3][3], ekk[3], vkk[3][6];
  //more for overflow
  int ilp_neigh[6];
	//transfer from force cls;
	int ilp_nelements					= force->ilp_nelements;       
  int ilp_nparams						= force->ilp_nparams;         
	int ilp_tap_flag					= force->ilp_tap_flag;
  double **ilp_cutILPsq			= force->ilp_cutILPsq;    
  double **ilp_cutsq				= force->ilp_cutsq;    
	Param_ilp_str *ilp_params = force->ilp_params;       
  int **ilp_elem2param			= force->ilp_elem2param;    
  int *ilp_map							= force->ilp_map;            
	/*** End variables for ILP ***/
	
	//filter short nbrs;
	int jshort, tot_short, jorg;
	int cur_ilpneigh[6], cur_nilp;
	int cur_neighshort[20], cur_numshort;
	
	double fxtmp,fytmp,fztmp;

  // loop over full neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) 
	{//ILP+Tersoff
		i = ilist[ii];
		xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    int itypem = ilp_map[type[i]];
    itag = tag[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int nilp = 0;
		
		//TERSOFF
    itype = map[type[i]];
    int numshort = 0;
    fxtmp = fytmp = fztmp = 0.0;
		
		/* filtering shortneighbors instead of tranversing long jlist */
		tot_short = flt_shortneigh[i*4+3];
		for(jshort = 0; jshort < tot_short; jshort++)
		{
			jorg = flt_shortneigh[i*4+jshort];
			j = jorg & NEIGHMASK;
			int jtypem = ilp_map[type[j]];
      jtag = tag[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
  		rsq = delx*delx + dely*dely + delz*delz;
			
			if((jorg >> 30) & 1)
			{
				ilp_neigh[nilp] = j;
  			vet[nilp][0] = -delx;
  			vet[nilp][1] = -dely;
  			vet[nilp][2] = -delz;
				nilp++;
			}
			
      if(!iskip[type[i]] && !ijskip[type[i]][type[j]])
			{
  			neighshort[numshort] = j;
				numshort++;
			}
		}//for

		calc_normal(nilp, jnum, vet, normal, dnormal, dnormdri);
    
		for (kk = 0; kk < nilp; kk ++) 
		{
      fkk[kk][0] = 0.0;
      fkk[kk][1] = 0.0;
      fkk[kk][2] = 0.0;
      vkk[kk][0] = 0.0;
      vkk[kk][1] = 0.0;
      vkk[kk][2] = 0.0;
      vkk[kk][3] = 0.0;
      vkk[kk][4] = 0.0;
      vkk[kk][5] = 0.0;
      ekk[kk] = 0.0;
    }
    
		jlist = firstneigh[i];
    jnum = numneigh[i];
			
		cur_numshort = 0;//recount;
		cur_nilp = 0;		 //recount;
    for (jj = 0; jj < jnum; jj++) 
		{
      j = jlist[jj];
      j &= NEIGHMASK;

      jtag = tag[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

			/* whether to recalculate */
			if (rsq < cutshortsq) 
			{
				if(!iskip[type[i]] && !ijskip[type[i]][type[j]])
				{
  	    	cur_neighshort[cur_numshort++] = j;
				}	
				
				int jtypem = ilp_map[type[j]];
				if(rsq != 0 && rsq < ilp_cutILPsq[itypem][jtypem] && atom->molecule[i] == atom->molecule[j])
				{
  	  	  cur_ilpneigh[cur_nilp++] = j;
  	  	}
  	  }

			//TERSOFF
			int ters_flag = 1;
      int vdwflag = 1;
			if (itag > jtag) 
			{
    	  if ((itag+jtag) % 2 == 0) {ters_flag = 0; vdwflag = 0;}
    	} 
			else if (itag < jtag) 
			{
    	  if ((itag+jtag) % 2 == 1) {ters_flag = 0; vdwflag = 0;}
    	} 
			else 
			{
    	  if (x[j][2] < x[i][2]) {ters_flag = 0; vdwflag = 0;}
    	  if (x[j][2] == ztmp && x[j][1] < ytmp) {ters_flag = 0; vdwflag = 0;}
    	  if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) {ters_flag = 0; vdwflag = 0;}
    	}

      jtype = map[type[j]];
			r = sqrt(rsq);
      double rinv = 1.0 / r;
			
			if(ters_flag && !iskip[type[i]] && !ijskip[type[i]][type[j]])
			{
				iparam_ij = elem2param[itype][jtype][jtype];
				if (rsq < params[iparam_ij].cutsq)
				{
					ters_FRep(&params[iparam_ij],r,rinv,fpair,eflag,evdwl);

      		fxtmp += delx*fpair;
      		fytmp += dely*fpair;
      		fztmp += delz*fpair;
      		f[j][0] -= delx*fpair;
      		f[j][1] -= dely*fpair;
      		f[j][2] -= delz*fpair;

      		if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
				}
			}//if-notSkip

      if (rsq < ilp_cutsq[type[i]][type[j]] && atom->molecule[i] != atom->molecule[j]) 
			{
        int iparam_ij = ilp_elem2param[ilp_map[type[i]]][ilp_map[type[j]]];

        Param_ilp_str *p = &ilp_params[iparam_ij];

        if (ilp_tap_flag) 
				{
          Rcut = sqrt(ilp_cutsq[type[i]][type[j]]);
          Tap = calc_Tap(r,Rcut);
          dTap = calc_dTap(r,Rcut);
        } 
				else {Tap = 1.0; dTap = 0.0;}

        double fvdw = 0, evdw = 0;
        
        if (vdwflag) 
				{
					r2inv = rinv * rinv;
					r6inv = r2inv*r2inv*r2inv;
					r8inv = r6inv*r2inv;
					
					TSvdw = 1.0 + exp(-p->d*(r/p->seff - 1.0));
					TSvdw2inv = 1 / (TSvdw * TSvdw);//pow(TSvdw,-2.0);
					double vvdw = -p->C6*r6inv/TSvdw;
					
					// derivatives
					fpair = -6.0*p->C6*r8inv/TSvdw + p->C6*p->d/p->seff*(TSvdw-1.0)*TSvdw2inv*r8inv*r;
					fsum = fpair*Tap - vvdw*dTap*rinv;
					fvdw = fsum;
					evdw = vvdw * Tap;
					if(eflag)
          {
            pvector[0] += evdw;
          }
        }
				
				// Calculate the transverse distance
        prodnorm1 = normal[0]*delx + normal[1]*dely + normal[2]*delz;
        rhosq1 = rsq - prodnorm1*prodnorm1;  // rho_ij
        rdsq1 = rhosq1*p->delta2inv; // (rho_ij/delta)^2

        // store exponents
        exp0 = exp(-p->lambda*(r-p->z0));
        exp1 = exp(-rdsq1);

        frho1 = exp1*p->C;
        Erep = 0.5*p->epsilon + frho1;
        Vilp = exp0*Erep;

        // derivatives
        fpair  = p->lambda*exp0*rinv*Erep;
        fpair1 = 2.0*exp0*frho1*p->delta2inv;
        fsum = fpair + fpair1;

        // derivatives of the product of rij and ni, the result is a vector
        dprodnorm1[0] = dnormdri[0][0]*delx + dnormdri[1][0]*dely + dnormdri[2][0]*delz;
        dprodnorm1[1] = dnormdri[0][1]*delx + dnormdri[1][1]*dely + dnormdri[2][1]*delz;
        dprodnorm1[2] = dnormdri[0][2]*delx + dnormdri[1][2]*dely + dnormdri[2][2]*delz;
        fp1[0] = prodnorm1*normal[0]*fpair1;
        fp1[1] = prodnorm1*normal[1]*fpair1;
        fp1[2] = prodnorm1*normal[2]*fpair1;
        fprod1[0] = prodnorm1*dprodnorm1[0]*fpair1;
        fprod1[1] = prodnorm1*dprodnorm1[1]*fpair1;
        fprod1[2] = prodnorm1*dprodnorm1[2]*fpair1;

        fkcx = (delx*fsum - fp1[0])*Tap - Vilp*dTap*delx*rinv;
        fkcy = (dely*fsum - fp1[1])*Tap - Vilp*dTap*dely*rinv;
        fkcz = (delz*fsum - fp1[2])*Tap - Vilp*dTap*delz*rinv;

        //This should be no use because fkcx need a lot of variables
        double ftotx = fvdw * delx + fkcx;
        double ftoty = fvdw * dely + fkcy;
        double ftotz = fvdw * delz + fkcz;
        f[i][0] += ftotx - fprod1[0]*Tap;
        f[i][1] += ftoty - fprod1[1]*Tap;
        f[i][2] += ftotz - fprod1[2]*Tap;
        f[j][0] -= ftotx;
        f[j][1] -= ftoty;
        f[j][2] -= ftotz;
				
        for (kk = 0; kk < nilp; kk++) 
				{
          k = ilp_neigh[kk];
          dprodnorm1[0] = dnormal[kk][0][0]*delx + dnormal[kk][1][0]*dely + dnormal[kk][2][0]*delz;
          dprodnorm1[1] = dnormal[kk][0][1]*delx + dnormal[kk][1][1]*dely + dnormal[kk][2][1]*delz;
          dprodnorm1[2] = dnormal[kk][0][2]*delx + dnormal[kk][1][2]*dely + dnormal[kk][2][2]*delz;
          fk[0] = (-prodnorm1*dprodnorm1[0]*fpair1)*Tap;
          fk[1] = (-prodnorm1*dprodnorm1[1]*fpair1)*Tap;
          fk[2] = (-prodnorm1*dprodnorm1[2]*fpair1)*Tap;
          fkk[kk][0] += fk[0];
          fkk[kk][1] += fk[1];
          fkk[kk][2] += fk[2];
          delkj[0] = x[i][0] + vet[kk][0] - x[j][0];
          delkj[1] = x[i][1] + vet[kk][1] - x[j][1];
          delkj[2] = x[i][2] + vet[kk][2] - x[j][2];


          if (vflag_atom) 
          {
            if (evflag) 
							ev_tally_buffer(k, j, ekk + kk, vkk[kk], eatom + j, vatom[j], nlocal,newton_pair,
															0.0,0.0,fk[0],fk[1],fk[2], delkj[0],delkj[1],delkj[2]);
          } 
          else 
          {
            if (evflag) 
							ev_tally_buffer(k, j, ekk + kk, vkk[kk], eatom + j, NULL, nlocal,newton_pair,
															0.0,0.0,fk[0],fk[1],fk[2],delkj[0],delkj[1],delkj[2]);
          }
        }
        //erep = Tap*Vilp;
        if (eflag) pvector[1] += erep = Tap*Vilp;
        if (evflag) 
					ev_tally_xyz(i,j,nlocal,newton_pair,erep + evdw,0.0, ftotx,ftoty,ftotz,delx,dely,delz);
      }
    } // loop over jj

		/*set recalc true,then return to recalculate */	
		if(cur_nilp != nilp || cur_numshort != numshort)
		{
			recalc = 1; return;
		}
		else
		{
			int re;
			if(cur_nilp == nilp)
			for(re = 0; re < nilp; re++)
			{
				if(ilp_neigh[re] != cur_ilpneigh[re])
				{
					recalc = 1; return;
				}
			}
			
			if(cur_numshort == numshort)
			for(re = 0; re < numshort; re++)
			{
				if(neighshort[re] != cur_neighshort[re])
				{
					recalc = 1; return;
				}
			}
		}

    for (kk = 0; kk < nilp; kk ++) 
		{
      int k = ilp_neigh[kk];
      f[k][0] += fkk[kk][0];
      f[k][1] += fkk[kk][1];
      f[k][2] += fkk[kk][2];
      if (eflag_atom) eatom[k] += ekk[kk];
      if (vflag_atom) 
			{
        vatom[k][0] += vkk[kk][0];
        vatom[k][1] += vkk[kk][1];
        vatom[k][2] += vkk[kk][2];
        vatom[k][3] += vkk[kk][3];
        vatom[k][4] += vkk[kk][4];
        vatom[k][5] += vkk[kk][5];
      }
    }//for-kk

		
		/*** Tersoff ***/
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // three-body interactions
    // skip immediately if I-J is not within cutoff
    double fjxtmp,fjytmp,fjztmp;
		double *dlr1, *dlr2, r1, r2, r1inv, r2inv;	
		//When iskip[type[i]] is ture, numshort==0;
		//			fxtmp == 0; fytmp == 0; fztmp == 0;
    for (jj = 0; jj < numshort; jj++) 
		{
      j = neighshort[jj];
      jtype = map[type[j]];
      iparam_ij = elem2param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 >= params[iparam_ij].cutsq) continue;
      
      // accumulate bondorder zeta for each i-j interaction via loop over k

      fjxtmp = fjytmp = fjztmp = 0.0;
      zeta_ij = 0.0;

      for (kk = 0; kk < numshort; kk++) 
			{
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];
        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[iparam_ijk].cutsq) continue;
        zeta_ij += zeta(&params[iparam_ijk],rsq1,rsq2,delr1,delr2);
      }

      // pairwise force due to zeta

      force_zeta(&params[iparam_ij],rsq1,zeta_ij,fpair,prefactor,eflag,evdwl);

      fxtmp += delr1[0]*fpair;
      fytmp += delr1[1]*fpair;
      fztmp += delr1[2]*fpair;
      fjxtmp -= delr1[0]*fpair;
      fjytmp -= delr1[1]*fpair;
      fjztmp -= delr1[2]*fpair;
      if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);
			

      // attractive term via loop over k
      for (kk = 0; kk < numshort; kk++) 
			{
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];
        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[iparam_ijk].cutsq) continue;
        attractive(&params[iparam_ijk],prefactor,rsq1,rsq2,delr1,delr2,fi,fj,fk);

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];
				
        if (vflag_atom) v_tally3(i,j,k,fj,fk,delr1,delr2);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
			
    }
		
		f[i][0] += fxtmp;
		f[i][1] += fytmp;
    f[i][2] += fztmp;
		
  }//for-ii

  if (vflag_fdotr) virial_fdotr_compute();
}

void PairLMFF::computeLMFFTwice(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,gnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3]/*,fk[3]*/;
  double zeta_ij,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;
  const double cutshortsq = cutmax*cutmax;

	int nt = atom->ntypes+1;
	int iskip[nt], ijskip[nt][nt];// When all is NULL;
	if(list->iskip == NULL && list->ijskip == NULL)
	{
		for(i = 0; i < nt; i++)
		{
			iskip[i] = 0;
			for(j = 0; j < nt; j++) ijskip[i][j] = 0;
		}
	}
	else
	{
		for(i = 0; i < nt; i++)
		{
			iskip[i] = neighbor->lists[0]->iskip[i];
			for(j = 0; j < nt; j++)
				ijskip[i][j] = neighbor->lists[0]->ijskip[i][j];
		}
	}

	int nlist = neighbor->nlist;
	int find_list = -1;
	for(i = 0; i < nlist; i++)
	{
		char *pairbuild = neighbor->pairnames[neighbor->lists[i]->pair_method-1];
		if(strcmp("full/bin/ghost", pairbuild) == 0)
		{
			inum				= neighbor->lists[i]->inum;
			gnum				= neighbor->lists[i]->gnum;
  		ilist				= neighbor->lists[i]->ilist;
  		numneigh		= neighbor->lists[i]->numneigh;
  		firstneigh	= neighbor->lists[i]->firstneigh;
			find_list = i;
			break;
		}
	}
	
	if(find_list == -1)
	{
		error->one(FLERR, "Can not find the right neighbor list");
	}


	/*** variables for ILP ***/
	pvector[0] = pvector[1] = 0.0;
  int l, ll;
  double erep, r,Rcut,r2inv,r6inv,r8inv;
	double Tap,dTap,Vilp,TSvdw,TSvdw2inv,fsum;
  erep = 0.0;
	//Vars for calc_normal
  int id,ip,m;
  double nn, nn2;
  double pv12[3],pv31[3],pv23[3],n1[3],dni[3],dnn[3][3],vet[3][3],dpvdri[3][3];
  double dn1[3][3][3],dpv12[3][3][3],dpv23[3][3][3],dpv31[3][3][3];
  double normal[3], dnormal[3][3][3], dnormdri[3][3];
  //Vars for calc_Frep
  double prodnorm1,fkcx,fkcy,fkcz;
  double rhosq1,exp0,exp1;
  double frho1,Erep,rdsq1,fpair1;
  double dprodnorm1[3] = {0.0, 0.0, 0.0};
  double fp1[3] = {0.0, 0.0, 0.0};
  double fprod1[3] = {0.0, 0.0, 0.0};
  double delkj[3] = {0.0, 0.0, 0.0};
  double fk[3] = {0.0, 0.0, 0.0};
  double fkk[3][3], ekk[3], vkk[3][6];
  //more for overflow
  int ilp_neigh[6];
	//transfer from force cls;
	int ilp_nelements					= force->ilp_nelements;       
  int ilp_nparams						= force->ilp_nparams;         
	int ilp_tap_flag					= force->ilp_tap_flag;
  double **ilp_cutILPsq			= force->ilp_cutILPsq;    
  double **ilp_cutsq				= force->ilp_cutsq;    
	Param_ilp_str *ilp_params = force->ilp_params;       
  int **ilp_elem2param			= force->ilp_elem2param;    
  int *ilp_map							= force->ilp_map;            
	/*** End variables for ILP ***/
	
	int jshort, tot_short, jorg;
	int cur_ilpneigh[6], cur_nilp;
	int cur_neighshort[20], cur_numshort;
	
	if(newton_pair == 0)
		error->one(FLERR, "The LMFF potential needs that newton_pair = on");
	
	double fxtmp,fytmp,fztmp;

  // loop over full neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) 
	{//ILP+Tersoff
		i = ilist[ii];
		xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    int itypem = ilp_map[type[i]];
    itag = tag[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int nilp = 0;
		
		//TERSOFF
    itype = map[type[i]];
    int numshort = 0;
    fxtmp = fytmp = fztmp = 0.0;
			
		tot_short = 0;
    for (jj = 0; jj < jnum; jj++) 
		{
      j = jlist[jj];
      j &= NEIGHMASK;
      int jtypem = ilp_map[type[j]];
      jtag = tag[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

			if (rsq < cutshortsq) 
			{
				flt_shortneigh[i*4+tot_short] = j;
				//TERSOFF
      	if(!iskip[type[i]] && !ijskip[type[i]][type[j]])
				{
      		neighshort[numshort++] = j;
      		if (numshort >= maxshort) 
					{
      		  maxshort += maxshort/2;
      		  memory->grow(neighshort,maxshort,"pair:neighshort");
      		}
				}//if-notSkip
				
				if (rsq != 0 && rsq < ilp_cutILPsq[itypem][jtypem] 
											&& atom->molecule[i] == atom->molecule[j])
				{
					flt_shortneigh[i*4+tot_short]	|= (1<<30);
					//ilp_index[nilp] = tot_short;
  	  	  ilp_neigh[nilp] = j;
      	  vet[nilp][0] = -delx;
      	  vet[nilp][1] = -dely;
      	  vet[nilp][2] = -delz;
      	  nilp ++;
      	}
				tot_short++;
			}//if-rsq
			flt_shortneigh[i*4+3] = tot_short;//num;
			
    }//for-jj
			
		#ifdef MY_DEBUG
		if(tot_short >= 4) printf("tot_short = %d\n", tot_short);
		#endif
		
		calc_normal(nilp, jnum, vet, normal, dnormal, dnormdri);
    
		for (kk = 0; kk < nilp; kk ++) 
		{
      fkk[kk][0] = 0.0;
      fkk[kk][1] = 0.0;
      fkk[kk][2] = 0.0;
      vkk[kk][0] = 0.0;
      vkk[kk][1] = 0.0;
      vkk[kk][2] = 0.0;
      vkk[kk][3] = 0.0;
      vkk[kk][4] = 0.0;
      vkk[kk][5] = 0.0;
      ekk[kk] = 0.0;
    }
    
		jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) 
		{
      j = jlist[jj];
      j &= NEIGHMASK;

      jtag = tag[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

			//TERSOFF
			int ters_flag = 1;
      int vdwflag = 1;
			if (itag > jtag) 
			{
    	  if ((itag+jtag) % 2 == 0) {ters_flag = 0; vdwflag = 0;}
    	} 
			else if (itag < jtag) 
			{
    	  if ((itag+jtag) % 2 == 1) {ters_flag = 0; vdwflag = 0;}
    	} 
			else 
			{
    	  if (x[j][2] < x[i][2]) {ters_flag = 0; vdwflag = 0;}
    	  if (x[j][2] == ztmp && x[j][1] < ytmp) {ters_flag = 0; vdwflag = 0;}
    	  if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) {ters_flag = 0; vdwflag = 0;}
    	}

      jtype = map[type[j]];
			r = sqrt(rsq);
      double rinv = 1.0 / r;
			
			if(ters_flag && !iskip[type[i]] && !ijskip[type[i]][type[j]])
			{
				iparam_ij = elem2param[itype][jtype][jtype];
				if (rsq < params[iparam_ij].cutsq)
				{
					ters_FRep(&params[iparam_ij],r,rinv,fpair,eflag,evdwl);

      		fxtmp += delx*fpair;
      		fytmp += dely*fpair;
      		fztmp += delz*fpair;
      		f[j][0] -= delx*fpair;
      		f[j][1] -= dely*fpair;
      		f[j][2] -= delz*fpair;

      		if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
				}
			}//if-notSkip

      if (rsq < ilp_cutsq[type[i]][type[j]] && atom->molecule[i] != atom->molecule[j]) 
			{
        int iparam_ij = ilp_elem2param[ilp_map[type[i]]][ilp_map[type[j]]];
        Param_ilp_str *p = &ilp_params[iparam_ij];

        if (ilp_tap_flag) 
				{
          Rcut = sqrt(ilp_cutsq[type[i]][type[j]]);
          Tap = calc_Tap(r,Rcut);
          dTap = calc_dTap(r,Rcut);
        } 
				else {Tap = 1.0; dTap = 0.0;}

        double fvdw = 0, evdw = 0;
        
        if (vdwflag) 
				{
					r2inv = rinv * rinv;
					r6inv = r2inv*r2inv*r2inv;
					r8inv = r6inv*r2inv;
					
					TSvdw = 1.0 + exp(-p->d*(r/p->seff - 1.0));
					TSvdw2inv = 1 / (TSvdw * TSvdw);//pow(TSvdw,-2.0);
					double vvdw = -p->C6*r6inv/TSvdw;
					
					// derivatives
					fpair = -6.0*p->C6*r8inv/TSvdw + p->C6*p->d/p->seff*(TSvdw-1.0)*TSvdw2inv*r8inv*r;
					fsum = fpair*Tap - vvdw*dTap*rinv;
					fvdw = fsum;
					evdw = vvdw * Tap;
					if(eflag)
          {
            pvector[0] += evdw;
          }
        }
				
				// Calculate the transverse distance
        prodnorm1 = normal[0]*delx + normal[1]*dely + normal[2]*delz;
        rhosq1 = rsq - prodnorm1*prodnorm1;  // rho_ij
        rdsq1 = rhosq1*p->delta2inv; // (rho_ij/delta)^2

        // store exponents
        exp0 = exp(-p->lambda*(r-p->z0));
        exp1 = exp(-rdsq1);

        frho1 = exp1*p->C;
        Erep = 0.5*p->epsilon + frho1;
        Vilp = exp0*Erep;

        // derivatives
        fpair  = p->lambda*exp0*rinv*Erep;
        fpair1 = 2.0*exp0*frho1*p->delta2inv;
        fsum = fpair + fpair1;

        // derivatives of the product of rij and ni, the result is a vector
        dprodnorm1[0] = dnormdri[0][0]*delx + dnormdri[1][0]*dely + dnormdri[2][0]*delz;
        dprodnorm1[1] = dnormdri[0][1]*delx + dnormdri[1][1]*dely + dnormdri[2][1]*delz;
        dprodnorm1[2] = dnormdri[0][2]*delx + dnormdri[1][2]*dely + dnormdri[2][2]*delz;
        fp1[0] = prodnorm1*normal[0]*fpair1;
        fp1[1] = prodnorm1*normal[1]*fpair1;
        fp1[2] = prodnorm1*normal[2]*fpair1;
        fprod1[0] = prodnorm1*dprodnorm1[0]*fpair1;
        fprod1[1] = prodnorm1*dprodnorm1[1]*fpair1;
        fprod1[2] = prodnorm1*dprodnorm1[2]*fpair1;

        fkcx = (delx*fsum - fp1[0])*Tap - Vilp*dTap*delx*rinv;
        fkcy = (dely*fsum - fp1[1])*Tap - Vilp*dTap*dely*rinv;
        fkcz = (delz*fsum - fp1[2])*Tap - Vilp*dTap*delz*rinv;

        //This should be no use because fkcx need a lot of variables
        double ftotx = fvdw * delx + fkcx;
        double ftoty = fvdw * dely + fkcy;
        double ftotz = fvdw * delz + fkcz;
        f[i][0] += ftotx - fprod1[0]*Tap;
        f[i][1] += ftoty - fprod1[1]*Tap;
        f[i][2] += ftotz - fprod1[2]*Tap;
        f[j][0] -= ftotx;
        f[j][1] -= ftoty;
        f[j][2] -= ftotz;

        for (kk = 0; kk < nilp; kk++) 
				{
          k = ilp_neigh[kk];
          dprodnorm1[0] = dnormal[kk][0][0]*delx + dnormal[kk][1][0]*dely + dnormal[kk][2][0]*delz;
          dprodnorm1[1] = dnormal[kk][0][1]*delx + dnormal[kk][1][1]*dely + dnormal[kk][2][1]*delz;
          dprodnorm1[2] = dnormal[kk][0][2]*delx + dnormal[kk][1][2]*dely + dnormal[kk][2][2]*delz;
          fk[0] = (-prodnorm1*dprodnorm1[0]*fpair1)*Tap;
          fk[1] = (-prodnorm1*dprodnorm1[1]*fpair1)*Tap;
          fk[2] = (-prodnorm1*dprodnorm1[2]*fpair1)*Tap;
          fkk[kk][0] += fk[0];
          fkk[kk][1] += fk[1];
          fkk[kk][2] += fk[2];
          delkj[0] = x[i][0] + vet[kk][0] - x[j][0];
          delkj[1] = x[i][1] + vet[kk][1] - x[j][1];
          delkj[2] = x[i][2] + vet[kk][2] - x[j][2];

          if (vflag_atom) 
          {
            if (evflag) 
							ev_tally_buffer(k, j, ekk + kk, vkk[kk], eatom + j, vatom[j], nlocal,newton_pair,
															0.0,0.0,fk[0],fk[1],fk[2], delkj[0],delkj[1],delkj[2]);
          } 
          else 
          {
            if (evflag) 
							ev_tally_buffer(k, j, ekk + kk, vkk[kk], eatom + j, NULL, nlocal,newton_pair,
															0.0,0.0,fk[0],fk[1],fk[2], delkj[0],delkj[1],delkj[2]);
          }
        }
        //erep = Tap*Vilp;
        if (eflag) pvector[1] += erep = Tap*Vilp;
        if (evflag) 
					ev_tally_xyz(i,j,nlocal,newton_pair,erep + evdw,0.0,ftotx,ftoty,ftotz,delx,dely,delz);
      }
    } // loop over jj
    for (kk = 0; kk < nilp; kk ++) 
		{
      int k = ilp_neigh[kk];
      f[k][0] += fkk[kk][0];
      f[k][1] += fkk[kk][1];
      f[k][2] += fkk[kk][2];
      if (eflag_atom) eatom[k] += ekk[kk];
      if (vflag_atom) 
			{
        vatom[k][0] += vkk[kk][0];
        vatom[k][1] += vkk[kk][1];
        vatom[k][2] += vkk[kk][2];
        vatom[k][3] += vkk[kk][3];
        vatom[k][4] += vkk[kk][4];
        vatom[k][5] += vkk[kk][5];
      }
    }//for-kk

		
		/*** Tersoff ***/
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // three-body interactions
    // skip immediately if I-J is not within cutoff
    double fjxtmp,fjytmp,fjztmp;
		double *dlr1, *dlr2, r1, r2, r1inv, r2inv;	
		//When iskip[type[i]] is ture, numshort==0;
		//			fxtmp == 0; fytmp == 0; fztmp == 0;
    for (jj = 0; jj < numshort; jj++) 
		{
      j = neighshort[jj];
      jtype = map[type[j]];
      iparam_ij = elem2param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 >= params[iparam_ij].cutsq) continue;
      
      // accumulate bondorder zeta for each i-j interaction via loop over k

      fjxtmp = fjytmp = fjztmp = 0.0;
      zeta_ij = 0.0;

      for (kk = 0; kk < numshort; kk++) 
			{
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];
        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[iparam_ijk].cutsq) continue;
        zeta_ij += zeta(&params[iparam_ijk],rsq1,rsq2,delr1,delr2);
      }

      // pairwise force due to zeta

      force_zeta(&params[iparam_ij],rsq1,zeta_ij,fpair,prefactor,eflag,evdwl);

      fxtmp += delr1[0]*fpair;
      fytmp += delr1[1]*fpair;
      fztmp += delr1[2]*fpair;
      fjxtmp -= delr1[0]*fpair;
      fjytmp -= delr1[1]*fpair;
      fjztmp -= delr1[2]*fpair;
      if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);

      // attractive term via loop over k
      for (kk = 0; kk < numshort; kk++) 
			{
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];
        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[iparam_ijk].cutsq) continue;
        attractive(&params[iparam_ijk],prefactor,rsq1,rsq2,delr1,delr2,fi,fj,fk);

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];
				
        if (vflag_atom) v_tally3(i,j,k,fj,fk,delr1,delr2);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
			
    }
		
		f[i][0] += fxtmp;
		f[i][1] += fytmp;
    f[i][2] += fztmp;
		
  }//for-ii

  if (vflag_fdotr) virial_fdotr_compute();
}

void PairLMFF::ters_FRep(Param *param, double r, double rinv, 
												double &fforce, int eflag, double &eng)
{
  double tmp_fc,tmp_fc_d,tmp_exp;

	if (r < param->bigr-param->bigd) {tmp_fc = 1.0; tmp_fc_d = 0.0;}
	else if (r > param->bigr+param->bigd) {tmp_fc = 0.0; tmp_fc_d = 0.0;}
	else 
	{
		tmp_fc = 0.5*(1.0 - sin(MY_PI2*(r - param->bigr)/param->bigd)); 
		tmp_fc_d =  -(MY_PI4/param->bigd) * cos(MY_PI2*(r - param->bigr)/param->bigd);
	}

  tmp_exp = exp(-param->lam1 * r);
  fforce = -param->biga * tmp_exp * (tmp_fc_d - tmp_fc*param->lam1) * rinv;
  if (eflag) eng = tmp_fc * param->biga * tmp_exp;
}


void PairLMFF::calc_normal(int nilp, int jnum, double vet[][3], 
													double normal[], double dnormal[][3][3], double dnormdri[][3])
{
	int id,ip,m;
  double nn, nn2;
  double pv12[3],pv31[3],pv23[3],n1[3],dni[3],dnn[3][3],/*vet[3][3],*/dpvdri[3][3];
  double dn1[3][3][3],dpv12[3][3][3],dpv23[3][3][3],dpv31[3][3][3];

	double jnuminv = 1.0 / jnum;

	//nilp;
  for (id = 0; id < 3; id++)
	{
    pv12[id] = 0.0;
    pv31[id] = 0.0;
    pv23[id] = 0.0;
    n1[id] = 0.0;
    dni[id] = 0.0;
    normal[id] = 0.0;
    for (ip = 0; ip < 3; ip++)
		{
      //vet[ip][id] = 0.0;
      dnn[ip][id] = 0.0;
      dpvdri[ip][id] = 0.0;
      dnormdri[ip][id] = 0.0;
      for (m = 0; m < 3; m++)
			{
        dpv12[m][ip][id] = 0.0;
        dpv31[m][ip][id] = 0.0;
        dpv23[m][ip][id] = 0.0;
        dn1[m][ip][id] = 0.0;
        dnormal[m][ip][id] = 0.0;
      }
    }
  }

  if (nilp <= 1) 
	{
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 1.0;
    for (id = 0; id < 3; id++)
		{
      for (ip = 0; ip < 3; ip++)
			{
        dnormdri[id][ip] = 0.0;
        for (m = 0; m < 3; m++)
				{
          dnormal[m][id][ip] = 0.0;
        }
      }
    }
  }
  else if (nilp == 2) 
	{
    cross_deriv(pv12, dpv12, vet, 0, 1, 2);
    // derivatives of pv12[0] to ri
    dpvdri[0][0] = 0.0;
    dpvdri[0][1] = vet[0][2]-vet[1][2];
    dpvdri[0][2] = vet[1][1]-vet[0][1];
    // derivatives of pv12[1] to ri
    dpvdri[1][0] = vet[1][2]-vet[0][2];
    dpvdri[1][1] = 0.0;
    dpvdri[1][2] = vet[0][0]-vet[1][0];
    // derivatives of pv12[2] to ri
    dpvdri[2][0] = vet[0][1]-vet[1][1];
    dpvdri[2][1] = vet[1][0]-vet[0][0];
    dpvdri[2][2] = 0.0;

    n1[0] = pv12[0];
    n1[1] = pv12[1];
    n1[2] = pv12[2];
    // the magnitude of the normal vector
    nn2 = n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2];
    nn = sqrt(nn2);
    if (nn == 0) error->one(FLERR,"The magnitude of the normal vector is zero");
    // the unit normal vector
    normal[0] = n1[0]/nn;
    normal[1] = n1[1]/nn;
    normal[2] = n1[2]/nn;
    // derivatives of nn, dnn:3x1 vector
    dni[0] = (n1[0]*dpvdri[0][0] + n1[1]*dpvdri[1][0] + n1[2]*dpvdri[2][0])/nn;
    dni[1] = (n1[0]*dpvdri[0][1] + n1[1]*dpvdri[1][1] + n1[2]*dpvdri[2][1])/nn;
    dni[2] = (n1[0]*dpvdri[0][2] + n1[1]*dpvdri[1][2] + n1[2]*dpvdri[2][2])/nn;
    for (id = 0; id < 3; id++)
		{
      for (ip = 0; ip < 3; ip++)
			{
        dnormdri[id][ip] = dpvdri[id][ip]/nn - n1[id]*dni[ip]/nn2;
      }
    }
    for (m = 0; m < 3; m++)
		{
      for (id = 0; id < 3; id++)
			{
        for (ip = 0; ip < 3; ip++)
				{
          dn1[m][id][ip] = dpv12[m][id][ip];
        }
      }
    }
    calc_dnormal(dnormal, dn1, n1, nn, nn2);

  }
  else if(nilp == 3) 
	{
    cross_deriv(pv12, dpv12, vet, 0, 1, 2);
    cross_deriv(pv31, dpv31, vet, 2, 0, 1);
    cross_deriv(pv23, dpv23, vet, 1, 2, 0);

    n1[0] = (pv12[0] + pv31[0] + pv23[0]) * jnuminv;
    n1[1] = (pv12[1] + pv31[1] + pv23[1]) * jnuminv;
    n1[2] = (pv12[2] + pv31[2] + pv23[2]) * jnuminv;
    // the magnitude of the normal vector
    nn2 = n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2];
    nn = sqrt(nn2);
    if (nn == 0) error->one(FLERR,"The magnitude of the normal vector is zero");
    // the unit normal vector
    normal[0] = n1[0]/nn;
    normal[1] = n1[1]/nn;
    normal[2] = n1[2]/nn;

    for (id = 0; id < 3; id++)
		{
      for (ip = 0; ip < 3; ip++)
			{
        dnormdri[id][ip] = 0.0;
      }
    }

    for (m = 0; m < 3; m++)
		{
      for (id = 0; id < 3; id++)
			{
        for (ip = 0; ip < 3; ip++)
				{
          dn1[m][id][ip] = (dpv12[m][id][ip] + dpv23[m][id][ip] + dpv31[m][id][ip])*jnuminv;///jnum;
        }
      }
    }
    calc_dnormal(dnormal, dn1, n1, nn, nn2);
  } 
	else 
	{
    error->one(FLERR,"There are too many neighbors for calculating normals");
  }
}


/*************************************************************************/
void PairLMFF::computeTersoff(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  //ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  const double cutshortsq = cutmax*cutmax;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
	
	bigint ntimestep  = update->ntimestep;
  bigint lastcall   = neighbor->lastcall;
  
	double fxtmp,fytmp,fztmp;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];
    int numshort = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutshortsq) {
        neighshort[numshort++] = j;
        if (numshort >= maxshort) 
					error->one(FLERR, "numshort exceeds maxshort");
      }

      jtag = tag[j];
      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < x[i][2]) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      jtype = map[type[j]];
      iparam_ij = elem2param[itype][jtype][jtype];
      if (rsq >= params[iparam_ij].cutsq) continue;

      repulsive(&params[iparam_ij],rsq,fpair,eflag,evdwl);

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }

    // three-body interactions
    // skip immediately if I-J is not within cutoff
    double fjxtmp,fjytmp,fjztmp;

    for (jj = 0; jj < numshort; jj++) {
      j = neighshort[jj];
      jtype = map[type[j]];
      iparam_ij = elem2param[itype][jtype][jtype];

      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 >= params[iparam_ij].cutsq) continue;

      // accumulate bondorder zeta for each i-j interaction via loop over k

      fjxtmp = fjytmp = fjztmp = 0.0;
      zeta_ij = 0.0;

      for (kk = 0; kk < numshort; kk++) {
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[iparam_ijk].cutsq) continue;

        zeta_ij += zeta(&params[iparam_ijk],rsq1,rsq2,delr1,delr2);
      }

      // pairwise force due to zeta

      force_zeta(&params[iparam_ij],rsq1,zeta_ij,fpair,prefactor,eflag,evdwl);

      fxtmp += delr1[0]*fpair;
      fytmp += delr1[1]*fpair;
      fztmp += delr1[2]*fpair;
      fjxtmp -= delr1[0]*fpair;
      fjytmp -= delr1[1]*fpair;
      fjztmp -= delr1[2]*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,-fpair,-delr1[0],-delr1[1],-delr1[2]);

      // attractive term via loop over k

      for (kk = 0; kk < numshort; kk++) {
        if (jj == kk) continue;
        k = neighshort[kk];
        ktype = map[type[k]];
        iparam_ijk = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 >= params[iparam_ijk].cutsq) continue;

        attractive(&params[iparam_ijk],prefactor,
                   rsq1,rsq2,delr1,delr2,fi,fj,fk);

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (vflag_atom) v_tally3(i,j,k,fj,fk,delr1,delr2);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();

	//computeILP(eflag, vflag);
}

void PairLMFF::computeILP(int eflag, int vflag)
{

	pvector[0] = pvector[1] = 0.0;
  int i,j,ii,jj,inum,jnum,itype,jtype,k,l,kk,ll;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,erep,fpair;
  double rsq,r,Rcut,r2inv,r6inv,r8inv,Tap,dTap,Vilp,TSvdw,TSvdw2inv,fsum;
  int *ilist,*jlist,*numneigh,**firstneigh;

  erep = 0.0;
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  //Vars for calc_normal
  int id,ip,m;
  double nn, nn2;
  double pv12[3],pv31[3],pv23[3],n1[3],dni[3],dnn[3][3],vet[3][3],dpvdri[3][3];
  double dn1[3][3][3],dpv12[3][3][3],dpv23[3][3][3],dpv31[3][3][3];
  double normal[3], dnormal[3][3][3], dnormdri[3][3];
  //Vars for calc_Frep
  double prodnorm1,fkcx,fkcy,fkcz;
  double rhosq1,exp0,exp1;
  double frho1,Erep,rdsq1,fpair1;
  double dprodnorm1[3] = {0.0, 0.0, 0.0};
  double fp1[3] = {0.0, 0.0, 0.0};
  double fprod1[3] = {0.0, 0.0, 0.0};
  double delkj[3] = {0.0, 0.0, 0.0};
  double fk[3] = {0.0, 0.0, 0.0};
  double fkk[3][3], ekk[3], vkk[3][6];// = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  //more for overflow
  int ilp_neigh[6];
  //int *ILP_neighs_i;

  //inum = list->inum;
  //ilist = list->ilist;
  //numneigh = list->numneigh;
  //firstneigh = list->firstneigh;

	list_ilp = neighbor->lists[1];//added by ping;
	inum				= list_ilp->inum;
  ilist				= list_ilp->ilist;
  numneigh		= list_ilp->numneigh;
  firstneigh	= list_ilp->firstneigh;
	
	bigint ntimestep  = update->ntimestep;
  bigint lastcall   = neighbor->lastcall;
	
	int ilp_nelements					= force->ilp_nelements;       
  int ilp_nparams						= force->ilp_nparams;         
	int ilp_tap_flag					= force->ilp_tap_flag;
  double **ilp_cutILPsq			= force->ilp_cutILPsq;    
  double **ilp_cutsq				= force->ilp_cutsq;    
	Param_ilp_str *ilp_params = force->ilp_params;       
  int **ilp_elem2param			= force->ilp_elem2param;    
  int *ilp_map							= force->ilp_map;            

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) 
	{
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    int itypem = ilp_map[itype];
    itag = tag[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int nilp = 0;
		
    for (jj = 0; jj < jnum; jj++) 
		{
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      int jtypem = ilp_map[jtype];
      jtag = tag[j];

      // two-body interactions from full neighbor list, skip half of them

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq != 0 && rsq < ilp_cutILPsq[itypem][jtypem] && atom->molecule[i] == atom->molecule[j])
			{
        ilp_neigh[nilp] = j;
        vet[nilp][0] = -delx;
        vet[nilp][1] = -dely;
        vet[nilp][2] = -delz;
        nilp ++;
      }
    }//for-jj
	
		//nilp;
    for (id = 0; id < 3; id++)
		{
      pv12[id] = 0.0;
      pv31[id] = 0.0;
      pv23[id] = 0.0;
      n1[id] = 0.0;
      dni[id] = 0.0;
      normal[id] = 0.0;
      for (ip = 0; ip < 3; ip++)
			{
        //vet[ip][id] = 0.0;
        dnn[ip][id] = 0.0;
        dpvdri[ip][id] = 0.0;
        dnormdri[ip][id] = 0.0;
        for (m = 0; m < 3; m++)
				{
          dpv12[m][ip][id] = 0.0;
          dpv31[m][ip][id] = 0.0;
          dpv23[m][ip][id] = 0.0;
          dn1[m][ip][id] = 0.0;
          dnormal[m][ip][id] = 0.0;
        }
      }
    }

    if (nilp <= 1) 
		{
      normal[0] = 0.0;
      normal[1] = 0.0;
      normal[2] = 1.0;
      for (id = 0; id < 3; id++)
			{
        for (ip = 0; ip < 3; ip++)
				{
          dnormdri[id][ip] = 0.0;
          for (m = 0; m < 3; m++)
					{
            dnormal[m][id][ip] = 0.0;
          }
        }
      }
    }
    else if (nilp == 2) 
		{
      cross_deriv(pv12, dpv12, vet, 0, 1, 2);
      // derivatives of pv12[0] to ri
      dpvdri[0][0] = 0.0;
      dpvdri[0][1] = vet[0][2]-vet[1][2];
      dpvdri[0][2] = vet[1][1]-vet[0][1];
      // derivatives of pv12[1] to ri
      dpvdri[1][0] = vet[1][2]-vet[0][2];
      dpvdri[1][1] = 0.0;
      dpvdri[1][2] = vet[0][0]-vet[1][0];
      // derivatives of pv12[2] to ri
      dpvdri[2][0] = vet[0][1]-vet[1][1];
      dpvdri[2][1] = vet[1][0]-vet[0][0];
      dpvdri[2][2] = 0.0;

      n1[0] = pv12[0];
      n1[1] = pv12[1];
      n1[2] = pv12[2];
      // the magnitude of the normal vector
      nn2 = n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2];
      nn = sqrt(nn2);
      if (nn == 0) error->one(FLERR,"The magnitude of the normal vector is zero");
      // the unit normal vector
      normal[0] = n1[0]/nn;
      normal[1] = n1[1]/nn;
      normal[2] = n1[2]/nn;
      // derivatives of nn, dnn:3x1 vector
      dni[0] = (n1[0]*dpvdri[0][0] + n1[1]*dpvdri[1][0] + n1[2]*dpvdri[2][0])/nn;
      dni[1] = (n1[0]*dpvdri[0][1] + n1[1]*dpvdri[1][1] + n1[2]*dpvdri[2][1])/nn;
      dni[2] = (n1[0]*dpvdri[0][2] + n1[1]*dpvdri[1][2] + n1[2]*dpvdri[2][2])/nn;
      // derivatives of unit vector ni respect to ri, the result is 3x3 matrix
      for (id = 0; id < 3; id++)
			{
        for (ip = 0; ip < 3; ip++)
				{
          dnormdri[id][ip] = dpvdri[id][ip]/nn - n1[id]*dni[ip]/nn2;
        }
      }
      // derivatives of non-normalized normal vector, dn1:3x3x3 array
      for (m = 0; m < 3; m++)
			{
        for (id = 0; id < 3; id++)
				{
          for (ip = 0; ip < 3; ip++)
					{
            dn1[m][id][ip] = dpv12[m][id][ip];
          }
        }
      }
      calc_dnormal(dnormal, dn1, n1, nn, nn2);

    }
    else if(nilp == 3) 
		{
      cross_deriv(pv12, dpv12, vet, 0, 1, 2);
      cross_deriv(pv31, dpv31, vet, 2, 0, 1);
      cross_deriv(pv23, dpv23, vet, 1, 2, 0);

      n1[0] = (pv12[0] + pv31[0] + pv23[0])/jnum;
      n1[1] = (pv12[1] + pv31[1] + pv23[1])/jnum;
      n1[2] = (pv12[2] + pv31[2] + pv23[2])/jnum;
      // the magnitude of the normal vector
      nn2 = n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2];
      nn = sqrt(nn2);
      if (nn == 0) error->one(FLERR,"The magnitude of the normal vector is zero");
      // the unit normal vector
      normal[0] = n1[0]/nn;
      normal[1] = n1[1]/nn;
      normal[2] = n1[2]/nn;

      // for the central atoms, dnormdri is always zero
      for (id = 0; id < 3; id++)
			{
        for (ip = 0; ip < 3; ip++)
				{
          dnormdri[id][ip] = 0.0;
        }
      }

      // derivatives of non-normalized normal vector, dn1:3x3x3 array
      for (m = 0; m < 3; m++)
			{
        for (id = 0; id < 3; id++)
				{
          for (ip = 0; ip < 3; ip++)
					{
            dn1[m][id][ip] = (dpv12[m][id][ip] + dpv23[m][id][ip] + dpv31[m][id][ip])/jnum;
          }
        }
      }
      calc_dnormal(dnormal, dn1, n1, nn, nn2);
    } 
		else 
		{
      error->one(FLERR,"There are too many neighbors for calculating normals");
    }
    
		for (kk = 0; kk < nilp; kk ++) 
		{
      fkk[kk][0] = 0.0;
      fkk[kk][1] = 0.0;
      fkk[kk][2] = 0.0;
      vkk[kk][0] = 0.0;
      vkk[kk][1] = 0.0;
      vkk[kk][2] = 0.0;
      vkk[kk][3] = 0.0;
      vkk[kk][4] = 0.0;
      vkk[kk][5] = 0.0;
      ekk[kk] = 0.0;
    }
    
		jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) 
		{
      j = jlist[jj];
      j &= NEIGHMASK;

      jtype = type[j];
      jtag = tag[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
			
      // only include the interation between different layers
      if (rsq < ilp_cutsq[itype][jtype] && atom->molecule[i] != atom->molecule[j]) 
			{

        int iparam_ij = ilp_elem2param[ilp_map[itype]][ilp_map[jtype]];
        Param_ilp_str *p = &ilp_params[iparam_ij];
        r = sqrt(rsq);
        double rinv = 1.0 / r;

        // turn on/off taper function
        if (ilp_tap_flag) 
				{
          Rcut = sqrt(ilp_cutsq[itype][jtype]);
          Tap = calc_Tap(r,Rcut);
          dTap = calc_dTap(r,Rcut);
        } 
				else {Tap = 1.0; dTap = 0.0;}

        int vdwflag = 1;
        if (itag > jtag) 
				{
          if ((itag+jtag) % 2 == 0) vdwflag = 0;
        } 
				else if (itag < jtag) 
				{
          if ((itag+jtag) % 2 == 1) vdwflag = 0;
        } 
				else 
				{
          if (x[j][2] < ztmp) vdwflag = 0;
          if (x[j][2] == ztmp && x[j][1] < ytmp) vdwflag = 0;
          if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) vdwflag = 0;
        }
        
        double fvdw = 0, evdw = 0;
        
        if (vdwflag) {
          r2inv = rinv * rinv;
          r6inv = r2inv*r2inv*r2inv;
          r8inv = r6inv*r2inv;

          TSvdw = 1.0 + exp(-p->d*(r/p->seff - 1.0));
          TSvdw2inv = 1 / (TSvdw * TSvdw);//pow(TSvdw,-2.0);
          double vvdw = -p->C6*r6inv/TSvdw;

          // derivatives
          fpair = -6.0*p->C6*r8inv/TSvdw + p->C6*p->d/p->seff*(TSvdw-1.0)*TSvdw2inv*r8inv*r;
          fsum = fpair*Tap - vvdw*dTap*rinv;

          fvdw = fsum;
          evdw = vvdw * Tap;
					if(eflag)
          {
            pvector[0] += evdw;
          }
          //pvector[0] += evdw;
        }
				
				// Calculate the transverse distance
        prodnorm1 = normal[0]*delx + normal[1]*dely + normal[2]*delz;
        rhosq1 = rsq - prodnorm1*prodnorm1;  // rho_ij
        rdsq1 = rhosq1*p->delta2inv; // (rho_ij/delta)^2

        // store exponents
        exp0 = exp(-p->lambda*(r-p->z0));
        exp1 = exp(-rdsq1);

        frho1 = exp1*p->C;
        Erep = 0.5*p->epsilon + frho1;
        Vilp = exp0*Erep;

        // derivatives
        fpair  = p->lambda*exp0*rinv*Erep;
        fpair1 = 2.0*exp0*frho1*p->delta2inv;
        fsum = fpair + fpair1;

        // derivatives of the product of rij and ni, the result is a vector
        dprodnorm1[0] = dnormdri[0][0]*delx + dnormdri[1][0]*dely + dnormdri[2][0]*delz;
        dprodnorm1[1] = dnormdri[0][1]*delx + dnormdri[1][1]*dely + dnormdri[2][1]*delz;
        dprodnorm1[2] = dnormdri[0][2]*delx + dnormdri[1][2]*dely + dnormdri[2][2]*delz;
        fp1[0] = prodnorm1*normal[0]*fpair1;
        fp1[1] = prodnorm1*normal[1]*fpair1;
        fp1[2] = prodnorm1*normal[2]*fpair1;
        fprod1[0] = prodnorm1*dprodnorm1[0]*fpair1;
        fprod1[1] = prodnorm1*dprodnorm1[1]*fpair1;
        fprod1[2] = prodnorm1*dprodnorm1[2]*fpair1;

        fkcx = (delx*fsum - fp1[0])*Tap - Vilp*dTap*delx*rinv;
        fkcy = (dely*fsum - fp1[1])*Tap - Vilp*dTap*dely*rinv;
        fkcz = (delz*fsum - fp1[2])*Tap - Vilp*dTap*delz*rinv;

        //This should be no use because fkcx need a lot of variables
        //fi + fj + sum(fk) = 0
        //-sum(fk) = fi + fj = -fprod*Tap
        //sum(fk) = fprod * Tap
        //fj = -fi - fprod*Tap
        double ftotx = fvdw * delx + fkcx;
        double ftoty = fvdw * dely + fkcy;
        double ftotz = fvdw * delz + fkcz;
        f[i][0] += ftotx - fprod1[0]*Tap;
        f[i][1] += ftoty - fprod1[1]*Tap;
        f[i][2] += ftotz - fprod1[2]*Tap;
        f[j][0] -= ftotx;
        f[j][1] -= ftoty;
        f[j][2] -= ftotz;
				
        // calculate the forces acted on the neighbors of atom i from atom j
        //ILP_neighs_i = ilp_neigh;
        for (kk = 0; kk < nilp; kk++) 
				{
          k = ilp_neigh[kk];
          //if (k == i) continue;
          // derivatives of the product of rij and ni respect to rk, k=0,1,2, where atom k is the neighbors of atom i
          dprodnorm1[0] = dnormal[kk][0][0]*delx + dnormal[kk][1][0]*dely + dnormal[kk][2][0]*delz;
          dprodnorm1[1] = dnormal[kk][0][1]*delx + dnormal[kk][1][1]*dely + dnormal[kk][2][1]*delz;
          dprodnorm1[2] = dnormal[kk][0][2]*delx + dnormal[kk][1][2]*dely + dnormal[kk][2][2]*delz;
          fk[0] = (-prodnorm1*dprodnorm1[0]*fpair1)*Tap;
          fk[1] = (-prodnorm1*dprodnorm1[1]*fpair1)*Tap;
          fk[2] = (-prodnorm1*dprodnorm1[2]*fpair1)*Tap;
          fkk[kk][0] += fk[0];
          fkk[kk][1] += fk[1];
          fkk[kk][2] += fk[2];
          delkj[0] = x[i][0] + vet[kk][0] - x[j][0];
          delkj[1] = x[i][1] + vet[kk][1] - x[j][1];
          delkj[2] = x[i][2] + vet[kk][2] - x[j][2];
          if (vflag_atom) 
          {
            if (evflag) 
							ev_tally_buffer(k, j, ekk + kk, vkk[kk], eatom + j, vatom[j], nlocal,newton_pair,
															0.0,0.0,fk[0],fk[1],fk[2],delkj[0],delkj[1],delkj[2]);
          } 
          else 
          {
            if (evflag) 
							ev_tally_buffer(k, j, ekk + kk, vkk[kk], eatom + j, NULL, nlocal,newton_pair,
															0.0,0.0,fk[0],fk[1],fk[2],delkj[0],delkj[1],delkj[2]);
          }
        }
        if (eflag) pvector[1] += erep = Tap*Vilp;
        if (evflag) 
					ev_tally_xyz(i,j,nlocal,newton_pair, erep + evdw,0.0,ftotx,ftoty,ftotz, delx,dely,delz);
      }
    } // loop over jj
    for (kk = 0; kk < nilp; kk ++) 
		{
      int k = ilp_neigh[kk];
      f[k][0] += fkk[kk][0];
      f[k][1] += fkk[kk][1];
      f[k][2] += fkk[kk][2];
      if (eflag_atom) eatom[k] += ekk[kk];
      if (vflag_atom) 
			{
        vatom[k][0] += vkk[kk][0];
        vatom[k][1] += vkk[kk][1];
        vatom[k][2] += vkk[kk][2];
        vatom[k][3] += vkk[kk][3];
        vatom[k][4] += vkk[kk][4];
        vatom[k][5] += vkk[kk][5];
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairLMFF::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(neighshort,maxshort,"pair:neighshort");
  map = new int[n+1];

	/*note: atom->nmax when using replicate command */
	
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLMFF::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLMFF::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params(); 

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLMFF::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style LMFF requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style LMFF requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLMFF::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairLMFF::read_file(char *file)
{
  int params_per_line = 17;
  char **words = new char*[params_per_line+1];

  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open LMFF potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int n,nwords,ielement,jelement,kelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in LMFF potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next line

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2],elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].powerm = atof(words[3]);
    params[nparams].gamma = atof(words[4]);
    params[nparams].lam3 = atof(words[5]);
    params[nparams].c = atof(words[6]);
    params[nparams].d = atof(words[7]);
    params[nparams].h = atof(words[8]);
    params[nparams].powern = atof(words[9]);
    params[nparams].beta = atof(words[10]);
    params[nparams].lam2 = atof(words[11]);
    params[nparams].bigb = atof(words[12]);
    params[nparams].bigr = atof(words[13]);
    params[nparams].bigd = atof(words[14]);
    params[nparams].lam1 = atof(words[15]);
    params[nparams].biga = atof(words[16]);
  	

    // currently only allow m exponent of 1 or 3

    params[nparams].powermint = int(params[nparams].powerm);

    if (params[nparams].c < 0.0 ||
        params[nparams].d < 0.0 ||
        params[nparams].powern < 0.0 ||
        params[nparams].beta < 0.0 ||
        params[nparams].lam2 < 0.0 ||
        params[nparams].bigb < 0.0 ||
        params[nparams].bigr < 0.0 ||
        params[nparams].bigd < 0.0 ||
        params[nparams].bigd > params[nparams].bigr ||
        params[nparams].lam1 < 0.0 ||
        params[nparams].biga < 0.0 ||
        params[nparams].powerm - params[nparams].powermint != 0.0 ||
        (params[nparams].powermint != 3 &&
         params[nparams].powermint != 1) ||
        params[nparams].gamma < 0.0)
      error->all(FLERR,"Illegal LMFF parameter");

    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairLMFF::setup_params()
{
  int i,j,k,m,n;

  // set elem2param for all element triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,nelements,"pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry");
        elem2param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].bigr + params[m].bigd;
    params[m].cutsq = params[m].cut*params[m].cut;

    params[m].c1 = pow(2.0*params[m].powern*1.0e-16,-1.0/params[m].powern);
    params[m].c2 = pow(2.0*params[m].powern*1.0e-8,-1.0/params[m].powern);
    params[m].c3 = 1.0/params[m].c2;
    params[m].c4 = 1.0/params[m].c1;	
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++)
    if (params[m].cut > cutmax) cutmax = params[m].cut;
}

/* ---------------------------------------------------------------------- */

void PairLMFF::repulsive(Param *param, double rsq, double &fforce,
                            int eflag, double &eng)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp;

  r = sqrt(rsq);
  tmp_fc = ters_fc(r,param);
  tmp_fc_d = ters_fc_d(r,param);
  tmp_exp = exp(-param->lam1 * r);
  fforce = -param->biga * tmp_exp * (tmp_fc_d - tmp_fc*param->lam1) / r;
  if (eflag) eng = tmp_fc * param->biga * tmp_exp;
}

/* ---------------------------------------------------------------------- */

double PairLMFF::zeta(Param *param, double rsqij, double rsqik,
                         double *delrij, double *delrik)
{
  double rij,rik,costheta,arg,ex_delr;

  rij = sqrt(rsqij);
  rik = sqrt(rsqik);
  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) / (rij*rik);

  if (param->powermint == 3) arg = cube(param->lam3 * (rij-rik));
  else arg = param->lam3 * (rij-rik);

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return ters_fc(rik,param) * ters_gijk(costheta,param) * ex_delr;
}

/* ---------------------------------------------------------------------- */

void PairLMFF::force_zeta(Param *param, double rsq, double zeta_ij,
                             double &fforce, double &prefactor,
                             int eflag, double &eng)
{
  double r,fa,fa_d,bij;

  r = sqrt(rsq);
  fa = ters_fa(r,param);
  fa_d = ters_fa_d(r,param);
  bij = ters_bij(zeta_ij,param);
  fforce = 0.5*bij*fa_d / r;
  prefactor = -0.5*fa * ters_bij_d(zeta_ij,param);
  if (eflag) eng = 0.5*bij*fa;
}

/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
------------------------------------------------------------------------- */

void PairLMFF::attractive(Param *param, double prefactor,
                             double rsqij, double rsqik,
                             double *delrij, double *delrik,
                             double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param);
}

/* ---------------------------------------------------------------------- */

double PairLMFF::ters_fc(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - sin(MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

double PairLMFF::ters_fc_d(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(MY_PI4/ters_D) * cos(MY_PI2*(r - ters_R)/ters_D);
}

/* ---------------------------------------------------------------------- */

double PairLMFF::ters_fa(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return -param->bigb * exp(-param->lam2 * r) * ters_fc(r,param);
}

/* ---------------------------------------------------------------------- */

double PairLMFF::ters_fa_d(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return param->bigb * exp(-param->lam2 * r) *
    (param->lam2 * ters_fc(r,param) - ters_fc_d(r,param));
}

/* ---------------------------------------------------------------------- */

double PairLMFF::ters_bij(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return 1.0/sqrt(tmp);
  if (tmp > param->c2)
    return (1.0 - pow(tmp,-param->powern) / (2.0*param->powern))/sqrt(tmp);
  if (tmp < param->c4) return 1.0;
  if (tmp < param->c3)
    return 1.0 - pow(tmp,param->powern)/(2.0*param->powern);
  return pow(1.0 + pow(tmp,param->powern), -1.0/(2.0*param->powern));
}

/* ---------------------------------------------------------------------- */

double PairLMFF::ters_bij_d(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return param->beta * -0.5*pow(tmp,-1.5);
  if (tmp > param->c2)
    return param->beta * (-0.5*pow(tmp,-1.5) *
                          // error in negligible 2nd term fixed 9/30/2015
                          // (1.0 - 0.5*(1.0 +  1.0/(2.0*param->powern)) *
                          (1.0 - (1.0 +  1.0/(2.0*param->powern)) *
                           pow(tmp,-param->powern)));
  if (tmp < param->c4) return 0.0;
  if (tmp < param->c3)
    return -0.5*param->beta * pow(tmp,param->powern-1.0);

  double tmp_n = pow(tmp,param->powern);
  return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*param->powern)))*tmp_n / zeta;
}

/* ---------------------------------------------------------------------- */

void PairLMFF::ters_zetaterm_d(double prefactor,
                                  double *rij_hat, double rij,
                                  double *rik_hat, double rik,
                                  double *dri, double *drj, double *drk,
                                  Param *param)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param);
  dfc = ters_fc_d(rik,param);
  if (param->powermint == 3) tmp = cube(param->lam3 * (rij-rik));
  else tmp = param->lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (param->powermint == 3)
    ex_delr_d = 3.0*cube(param->lam3) * square(rij-rik)*ex_delr;
  else ex_delr_d = param->lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(cos_theta,param);
  gijk_d = ters_gijk_d(cos_theta,param);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
}

/* ---------------------------------------------------------------------- */

void PairLMFF::costheta_d(double *rij_hat, double rij,
                             double *rik_hat, double rik,
                             double *dri, double *drj, double *drk)
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  double cos_theta = vec3_dot(rij_hat,rik_hat);

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(1.0/rij,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(1.0/rik,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}
