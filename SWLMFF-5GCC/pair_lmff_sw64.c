#include "sunway.h"
#include "pair_lmff_sw64.h"
//#include "pair_lmff_inline.h"
#include "pair_lmff_inline_sleef.h"
#include "sw_wcache.h"
#include "sleef_math_full.h"
//#define MY_DEBUG


#ifdef MPE
#define LWPF_UNITS U(LMFF)
#include "lwpf2.h"
#include "mpi.h"
long *mylock = NULL;
int flag_athread_init = 0;
extern SLAVE_FUN(pack_atoms_lmff_para)(pair_lmff_param_t *);
extern SLAVE_FUN(compute_lmff_para)(pair_lmff_param_t *);
extern SLAVE_FUN(copy_force_back_para)(pair_lmff_param_t *);
extern SLAVE_FUN(compute_lmff_onelonglist)(pair_lmff_param_t *);
extern SLAVE_FUN(compute_lmff_twolonglist)(pair_lmff_param_t *);

//int bubbleSort(int *jlist, int jnum, tagint imol, tagint *molecule)
//{
//	int tmp;
//	int i, j, iatom, jatom;
//	int samelayer = 0;
//	for(i = 0; i < jnum-1; i++)
//	{
//		if(i == 0)
//		{
//			if(imol == molecule[i])
//				samelayer++;
//		}
//		for(j = 0; j < jnum-1-i; j++)
//		{
//			if(molecule[jlist[j]] > molecule[jlist[j+1]])
//			{
//				tmp = jlist[j];
//				jlist[j] = jlist[j+1];
//				jlist[j+1] = tmp;
//			}
//		}
//	}
//	if(imol == molecule[jnum-1]) samelayer++;
//	return samelayer;
//}

void pack_atoms_lmff_serial(pair_lmff_param_t *pm)
{
	#ifdef SUNWAY
	//if(flag_athread_init == 0)
	//{
	//	CRTS_init();
	//	flag_athread_init = 1;
	//}
	
	
	
	if(athread_idle() == 0)
		athread_init();

	evt_conf_t conf;
	conf.pc_mask = 0xf;
	conf.evt[0] = PC0_CYCLE;
	conf.evt[1] = PC1_N_INST;
	conf.evt[2] = PC2_N_GLD;
	conf.evt[3] = PC3_CYCLE;
	lwpf_init(&conf);
	
	athread_spawn(pack_atoms_lmff_para, pm);
  athread_join();
  static reported = 0;
	if(myrank == 0 && !reported)
	{
	  lwpf_report_summary(stdout);
	  reported = 1;
	  //lwpf_report_detail(stdout);
	}
	
	#else
	int i, j;
	int inum					= pm->inum;
	int gnum 					= pm->gnum;
	int allnum				= inum + gnum;
	int *type					= pm->type;
	tagint *tag				= pm->tag;
	tagint *molecule	= pm->molecule;
	double (*x)[3] = pm->x;
	atom_in *my_atoms = pm->my_atoms;

	for(i = 0; i < allnum; i++)
	{
		my_atoms[i].x[0] = x[i][0];
		my_atoms[i].x[1] = x[i][1];
		my_atoms[i].x[2] = x[i][2];
		my_atoms[i].type = type[i];
		my_atoms[i].padding		= 0;
		my_atoms[i].molecule	= molecule[i];
		my_atoms[i].tag				= tag[i];
		my_atoms[i].pad20			= 0;
		my_atoms[i].pad21			= 0;
	}
	#endif
}

void copy_force_back_serial(pair_lmff_param_t *pm)
{
	#ifdef SUNWAY
	//if(flag_athread_init == 0)
	//{
	//	CRTS_init();
	//	flag_athread_init = 1;
	//}

	if(athread_idle() == 0)
		athread_init();
	athread_spawn(copy_force_back_para, pm);
  athread_join();

	#else
	double (*f)[3]		=  pm->f;
	double (*frc)[4]	=  pm->frc;
	int ntotal				= pm->inum + pm->gnum;
	int i;
	
	for(i = 0; i < ntotal; i++)
  {
    f[i][0] += frc[i][0];
    f[i][1] += frc[i][1];
    f[i][2] += frc[i][2];
  }
	#endif
}


void compute_lmff_serial(pair_lmff_param_t *pm)
{
  //if(flag_athread_init == 0)
  //{
  //	CRTS_init();
  //	flag_athread_init = 1;
  //}
	
  if(athread_idle() == 0)
    athread_init();
  int tt;
  /*** allocate for my_lock ***/	
  if(mylock == NULL)
    {
      int lbits = WRT_C_S;
      long lmask = (1 << lbits) - 1;
      int ntotal = pm->inum + pm->gnum + 1024;
      long nblk = (ntotal + lmask) & ~lmask;
      //mylock = (long *)libc_uncached_malloc(sizeof(long) * nblk * 2);
      mylock = (long *)malloc(sizeof(long) * nblk * 2);
      for(tt = 0; tt < nblk*2; tt++)
	mylock[tt] = 0;
    }
  pm->mylock = mylock;
	
	
  /* evt_conf_t conf; */
  /* conf.pc_mask = 0xf; */
  /* conf.evt[0] = PC0_CYCLE; */
  /* conf.evt[1] = PC1_INST; */
  /* conf.evt[2] = PC2_GLD; */
  /* conf.evt[3] = PC3_GST; */
  /* lwpf_init(&conf); */
	
  /* perf_config_t conf; */
  /* conf.pcrc = PCRC_ALL; */
  /* conf.pcr0 = PC0_CYCLE; */
  /* conf.pcr1 = PC1_CYCLE; */
  /* conf.pcr2 = PC2_N_GLD; */
  /* //conf.pcr2 = PC2_N_DMA_REQ; */
  /* lwpf_init(&conf); */
	
  pm->recalc = 0;

  if(pm->ntimestep == pm->lastcall)
    {
      athread_spawn(compute_lmff_twolonglist, pm);
      athread_join();
    }
  else if(pm->ntimestep != pm->lastcall)
    {
      athread_spawn(compute_lmff_onelonglist, pm);
      athread_join();
		
      if(pm->recalc)
	{
	  int allnum = pm->inum + pm->gnum;
	  memset(pm->frc, 0, sizeof(double)*4*(allnum+512));
	  athread_spawn(compute_lmff_twolonglist, pm);
	  athread_join();
	}
    }
  else
    {/* not using filtering */
      athread_spawn(compute_lmff_para, pm);
      athread_join();
    }
	
  if(pm->myrank == 0)
  {
    //lwpf_report_summary(stdout);
    //lwpf_report_summary(stdout, &conf);
    
  	//lwpf_report_detail(stdout);
  }
	
  return;


  int i,j,k,ii,jj,kk,inum,gnum,jnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3];
  double zeta_ij,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

  double (*x)[3]				= pm->x;
  double (*f)[3]				= pm->f;
  tagint *tag						= pm->tag;
  tagint *molecule			= pm->molecule;
  int *type							= pm->type;
  int ntypes						= pm->ntypes;
  int nlocal						= pm->nlocal;
  int newton_pair				= pm->newton_pair;
  double cutshortsq 		= pm->cutshortsq;
  int nelements					= pm->nelements;       
  int nparams						= pm->nparams;         
  Param *ters_params		= pm->ters_params;       
  int *elem2param				= pm->elem2param;//nelements^3;    
  int *map							= pm->map;            
  int maxshort					= pm->maxshort;
  atom_in *my_atoms			= pm->my_atoms;

  int eflag						= pm->eflag;
  int evflag					= pm->evflag;
  int eflag_either		= pm->eflag_either;
  int eflag_global		= pm->eflag_global;
  int eflag_atom			= pm->eflag_atom;
  int vflag						= pm->vflag;
  int vflag_either		= pm->vflag_either;
  int vflag_global		= pm->vflag_global;
  int vflag_atom			= pm->vflag_atom;
  double *p_eng_vdwl  = &(pm->eng_vdwl);
  double *p_eng_coul  = &(pm->eng_coul);
  double *p_virial    = pm->virial;
  double (*vatom)[6]  = pm->vatom;
  double *eatom       = pm->eatom;

  int nt = pm->ntypes+1;
  int iskip[nt], ijskip[nt][nt];
  for(i = 0; i < nt; i++)
    {
      iskip[i] = pm->iskip[i];
      for(j = 0; j < nt; j++)
	ijskip[i][j] = pm->ijskip[i*nt+j];
    }
  inum				= pm->inum;
  gnum				= pm->gnum;
  ilist				= pm->ilist;
  numneigh		= pm->numneigh;
  firstneigh	= pm->firstneigh;
	
  //buffer atom info for neighshort;
  ParamBf neighBf[20];
  int neighshort[20];

  /*** variables for ILP ***/
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
  int ilp_nelements					= pm->ilp_nelements;       
  int ilp_nparams						= pm->ilp_nparams;         
  int ilp_tap_flag					= pm->ilp_tap_flag;
  double *ilp_cutILPsq			= pm->ilp_cutILPsq;//nelements^2;        
  double *ilp_cutsq					= pm->ilp_cutsq;//(ntypes+1)^2;    
  Param_ilp_str *ilp_params = pm->ilp_params;       
  int *ilp_elem2param				= pm->ilp_elem2param;//nelements^2;        
  int *ilp_map							= pm->ilp_map;            
  double pvector[2];
  pvector[0] = pm->pvector[0];
  pvector[1] = pm->pvector[1];
  /*** End variables for ILP ***/

  bigint ntimestep  = pm->ntimestep;
  bigint lastcall   = pm->lastcall;
  int myrank				= pm->myrank;
  
  double fxtmp,fytmp,fztmp;
  atom_in *iatom, *jatom;

  int ist, ied, isz, ioff;
  int jst, jed, jsz, joff;

  // loop over full neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) 
    {//ILP+Tersoff
      i = ilist[ii];

#ifdef MY_DEBUG
      if(i != ii)
	printf("i = %d, ii = %d\n", i, ii);
#endif
      iatom = &my_atoms[i];
      xtmp = iatom->x[0];
      ytmp = iatom->x[1];
      ztmp = iatom->x[2];
      int itypem = ilp_map[iatom->type];
      itag = iatom->tag;
      jlist = firstneigh[i];
      jnum = numneigh[i];
      int nilp = 0;
		
      //TERSOFF
      itype = map[iatom->type];
      int numshort = 0;
      fxtmp = fytmp = fztmp = 0.0;

      for (jj = 0; jj < jnum; jj++) 
	{
	  j = jlist[jj];
	  jatom = &my_atoms[j];
	  int jtypem = ilp_map[jatom->type];
	  jtag = jatom->tag;

	  delx = xtmp - jatom->x[0];
	  dely = ytmp - jatom->x[1];
	  delz = ztmp - jatom->x[2];
	  rsq = delx*delx + dely*dely + delz*delz;

	  if(rsq != 0 && rsq < ilp_cutILPsq[itypem*nelements+jtypem] && iatom->molecule == jatom->molecule)
	    {
	      ilp_neigh[nilp] = j;
	      vet[nilp][0] = -delx;
	      vet[nilp][1] = -dely;
	      vet[nilp][2] = -delz;
	      nilp ++;
	    }


	  //TERSOFF
	  if(!iskip[iatom->type] && !ijskip[iatom->type][jatom->type])
	    {
	      if (rsq < cutshortsq) 
		{
		  neighBf[numshort].rsq = rsq;
		  neighBf[numshort].r = sqrt(rsq);
		  neighBf[numshort].rinv = 1/neighBf[numshort].r;
		  neighBf[numshort].delr[0] = -delx;
		  neighBf[numshort].delr[1] = -dely;
		  neighBf[numshort].delr[2] = -delz;
		  //neighBf[numshort].id = j;
		  neighBf[numshort].tp = map[jatom->type];
		  neighBf[numshort].padding = 0;

		  int ters_flag = 1;
		  if (itag > jtag) 
		    {
		      if ((itag+jtag) % 2 == 0) {ters_flag = 0; }
		    } 
		  else if (itag < jtag) 
		    {
		      if ((itag+jtag) % 2 == 1) {ters_flag = 0; }
		    } 
		  else 
		    {
		      if (jatom->x[2] < iatom->x[2]) {ters_flag = 0;}
		      if (jatom->x[2] == iatom->x[2] && jatom->x[1] <  iatom->x[1]) {ters_flag = 0; }
		      if (jatom->x[2] == iatom->x[2] && jatom->x[1] == iatom->x[1] && jatom->x[0] < iatom->x[0]) {ters_flag = 0;}
		    }
		  neighBf[numshort].ters_flag = ters_flag;


		  neighshort[numshort++] = j;
#ifdef MY_DEBUG
		  if (numshort >= maxshort && numshort > 4) 
		    {
		      printf("numshort > maxshot, numshort=%d\n", numshort);
		    }
#endif
		}
	    }//if-notSkip
	}//for-jj

      calc_normal(nilp, jnum, vet, normal, dnormal, dnormdri);
		
      /* Debugging */
#ifdef MY_DEBUG
      if(pm->myrank == 0)
	{
	  printf("i=%d, ntimestep=%d,lastcall=%d,", i,ntimestep,lastcall);
	  int tt;
	  for(tt = 0; tt < nilp; tt++) printf(" %d", ilp_neigh[tt]);
	  printf(";");
	  for(tt = 0; tt < numshort; tt++) printf(" %d", neighshort[tt]);
	  printf("\n");
	}
#endif
    
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
	  jatom = &my_atoms[j];

	  jtag = jatom->tag;
	  delx = xtmp - jatom->x[0];
	  dely = ytmp - jatom->x[1];
	  delz = ztmp - jatom->x[2];
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
	      if (jatom->x[2] < iatom->x[2]) {ters_flag = 0; vdwflag = 0;}
	      if (jatom->x[2] == ztmp && iatom->x[1] < ytmp) {ters_flag = 0; vdwflag = 0;}
	      if (jatom->x[2] == ztmp && iatom->x[1] == ytmp && iatom->x[0] < xtmp) {ters_flag = 0; vdwflag = 0;}
	    }

	  jtype = map[jatom->type];
	  r = sqrt(rsq);
	  double rinv = 1.0 / r;
			
	  if(ters_flag && !iskip[jatom->type] && !ijskip[iatom->type][jatom->type])
	    {
	      iparam_ij = elem2param[(itype*nelements+jtype)*nelements+jtype];
	      if (rsq < ters_params[iparam_ij].cutsq)
		{
		  ters_FRep(&ters_params[iparam_ij],r,rinv,&fpair,eflag,&evdwl);

		  fxtmp += delx*fpair;
		  fytmp += dely*fpair;
		  fztmp += delz*fpair;
		  f[j][0] -= delx*fpair;
		  f[j][1] -= dely*fpair;
		  f[j][2] -= delz*fpair;

					
		  if(evflag)
		    {
		      if (eflag_either && eflag_global) 
			*p_eng_vdwl += evdwl;
		      if (vflag_either && vflag_global) 
			{
			  p_virial[0] += delx*delx*fpair;
			  p_virial[1] += dely*dely*fpair;
			  p_virial[2] += delz*delz*fpair;
			  p_virial[3] += delx*dely*fpair;
			  p_virial[4] += delx*delz*fpair;
			  p_virial[5] += dely*delz*fpair;
			}
		    }//ev_tally_simple;

		}
	    }//if-notSkip

	  if (rsq < ilp_cutsq[iatom->type*nt+jatom->type] && iatom->molecule != jatom->molecule) 
	    {
	      int iparam_ij = ilp_elem2param[ilp_map[iatom->type] *nelements+ ilp_map[jatom->type]];
	      Param_ilp_str *p = &ilp_params[iparam_ij];

	      if (ilp_tap_flag) 
		{
		  Rcut = sqrt(ilp_cutsq[iatom->type*nt+jatom->type]);
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
		  delkj[0] = iatom->x[0] + vet[kk][0] - jatom->x[0];
		  delkj[1] = iatom->x[1] + vet[kk][1] - jatom->x[1];
		  delkj[2] = iatom->x[2] + vet[kk][2] - jatom->x[2];
					
		  if (evflag) 
		    {
		      if (vflag_either && vflag_global) 
			{
			  p_virial[0] += delkj[0]*fk[0];
			  p_virial[1] += delkj[1]*fk[1];
			  p_virial[2] += delkj[2]*fk[2];
			  p_virial[3] += delkj[0]*fk[1];
			  p_virial[4] += delkj[0]*fk[2];
			  p_virial[5] += delkj[1]*fk[2];
			}//ec_tally_buffer_simple;
		    }

		}
	      //erep = Tap*Vilp;
	      if (eflag) pvector[1] += erep = Tap*Vilp;
	      if (evflag) 
		{
		  if (eflag_either && eflag_global) 
		    *p_eng_vdwl += (erep+evdw);
		  if (vflag_either && vflag_global) 
		    {
		      p_virial[0] += delx*ftotx;
		      p_virial[1] += dely*ftoty;
		      p_virial[2] += delz*ftotz;
		      p_virial[3] += delx*ftoty;
		      p_virial[4] += delx*ftotz;
		      p_virial[5] += dely*ftotz;
		    }//ev_tally_xyz_simple

		}
	    }
	} // loop over jj
      for (kk = 0; kk < nilp; kk ++) 
	{
	  int k = ilp_neigh[kk];
	  f[k][0] += fkk[kk][0];
	  f[k][1] += fkk[kk][1];
	  f[k][2] += fkk[kk][2];
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
	  jtype = neighBf[jj].tp;
	  iparam_ij = elem2param[(itype*nelements+jtype)*nelements+jtype];
	  dlr1 = neighBf[jj].delr;
	  rsq1 = neighBf[jj].rsq;
	  r1 = neighBf[jj].r;
	  r1inv = neighBf[jj].rinv;
	  if (rsq1 >= ters_params[iparam_ij].cutsq) continue;

	  // accumulate bondorder zeta for each i-j interaction via loop over k

	  fjxtmp = fjytmp = fjztmp = 0.0;
	  zeta_ij = 0.0;

	  for (kk = 0; kk < numshort; kk++) 
	    {
	      if (jj == kk) continue;
	      k = neighshort[kk];
	      ktype = neighBf[kk].tp;
	      iparam_ijk = elem2param[(itype*nelements+jtype)*nelements+ktype];
	      dlr2 = neighBf[kk].delr;
	      rsq2 = neighBf[kk].rsq;
	      r2 = neighBf[kk].r;
	      r2inv = neighBf[kk].rinv;
	      if (rsq2 >= ters_params[iparam_ijk].cutsq) continue;
	      zeta_ij += ters_zeta(&ters_params[iparam_ijk],r1,r1inv,r2,r2inv,dlr1,dlr2);
	    }

	  // pairwise force due to zeta
	  ters_force_zeta(&ters_params[iparam_ij],r1,r1inv,zeta_ij,&fpair,&prefactor,eflag,&evdwl);
			
	  fxtmp += dlr1[0]*fpair;
	  fytmp += dlr1[1]*fpair;
	  fztmp += dlr1[2]*fpair;
	  fjxtmp -= dlr1[0]*fpair;
	  fjytmp -= dlr1[1]*fpair;
	  fjztmp -= dlr1[2]*fpair;

	  if(evflag)
	    {
	      if (eflag_either && eflag_global) 
		*p_eng_vdwl += evdwl;
	      if (vflag_either && vflag_global) 
		{
		  p_virial[0] += -dlr1[0]*dlr1[0]*fpair;
		  p_virial[1] += -dlr1[1]*dlr1[1]*fpair;
		  p_virial[2] += -dlr1[2]*dlr1[2]*fpair;
		  p_virial[3] += -dlr1[0]*dlr1[1]*fpair;
		  p_virial[4] += -dlr1[0]*dlr1[2]*fpair;
		  p_virial[5] += -dlr1[1]*dlr1[2]*fpair;
		}
	    }//ev_tally_simple;



	  // attractive term via loop over k
	  for (kk = 0; kk < numshort; kk++) 
	    {
	      if (jj == kk) continue;
	      k = neighshort[kk];
	      ktype = neighBf[kk].tp;
	      iparam_ijk = elem2param[(itype*nelements+jtype)*nelements+ktype];
	      dlr2 = neighBf[kk].delr;
	      rsq2 = neighBf[kk].rsq;
	      r2 = neighBf[kk].r;
	      r2inv = neighBf[kk].rinv;
	      if (rsq2 >= ters_params[iparam_ijk].cutsq) continue;
	      ters_Att(&ters_params[iparam_ijk],prefactor,r1,r1inv,r2,r2inv,dlr1,dlr2,fi,fj,fk);

	      fxtmp += fi[0];
	      fytmp += fi[1];
	      fztmp += fi[2];
	      fjxtmp += fj[0];
	      fjytmp += fj[1];
	      fjztmp += fj[2];
	      f[k][0] += fk[0];
	      f[k][1] += fk[1];
	      f[k][2] += fk[2];
	    }
	  f[j][0] += fjxtmp;
	  f[j][1] += fjytmp;
	  f[j][2] += fjztmp;
			
	}
		
      f[i][0] += fxtmp;
      f[i][1] += fytmp;
      f[i][2] += fztmp;
		
    }//for-ii

  pm->pvector[0] = pvector[0];
  pm->pvector[1] = pvector[1];
}

//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
void compute_lmff_layer_serial(pair_lmff_param_t *pm)
{
  int i,j,k,ii,jj,kk,inum,gnum,jnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3];
  double zeta_ij,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;

	double (*x)[3]				= pm->x;
  double (*f)[3]				= pm->f;
  tagint *tag						= pm->tag;
  tagint *molecule			= pm->molecule;
  int *type							= pm->type;
  int ntypes						= pm->ntypes;
  int nlocal						= pm->nlocal;
  int newton_pair				= pm->newton_pair;
  double cutshortsq 		= pm->cutshortsq;
	int nelements					= pm->nelements;       
  int nparams						= pm->nparams;         
	Param *ters_params		= pm->ters_params;       
  int *elem2param				= pm->elem2param;//nelements^3;    
  int *map							= pm->map;            
	int maxshort					= pm->maxshort;
	atom_in *my_atoms			= pm->my_atoms;

	int eflag						= pm->eflag;
	int evflag					= pm->evflag;
	int eflag_either		= pm->eflag_either;
	int eflag_global		= pm->eflag_global;
	int eflag_atom			= pm->eflag_atom;
	int vflag						= pm->vflag;
	int vflag_either		= pm->vflag_either;
	int vflag_global		= pm->vflag_global;
	int vflag_atom			= pm->vflag_atom;
	double *p_eng_vdwl  = &(pm->eng_vdwl);
  double *p_eng_coul  = &(pm->eng_coul);
  double *p_virial    = pm->virial;
	double (*vatom)[6]  = pm->vatom;
  double *eatom       = pm->eatom;

	int nt = pm->ntypes+1;
	int iskip[nt], ijskip[nt][nt];
	for(i = 0; i < nt; i++)
	{
		iskip[i] = pm->iskip[i];
		for(j = 0; j < nt; j++)
			ijskip[i][j] = pm->ijskip[i*nt+j];
	}
	inum				= pm->inum;
	gnum				= pm->gnum;
  ilist				= pm->ilist;
  numneigh		= pm->numneigh;
	int *nummol			= pm->nummol;
  firstneigh	= pm->firstneigh;
	
	//buffer atom info for neighshort;
	ParamBf neighBf[20];
	int neighshort[20];

	/*** variables for ILP ***/
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
	int ilp_nelements					= pm->ilp_nelements;       
  int ilp_nparams						= pm->ilp_nparams;         
	int ilp_tap_flag					= pm->ilp_tap_flag;
  double *ilp_cutILPsq			= pm->ilp_cutILPsq;//nelements^2;        
  double *ilp_cutsq					= pm->ilp_cutsq;//(ntypes+1)^2;    
	Param_ilp_str *ilp_params = pm->ilp_params;       
  int *ilp_elem2param				= pm->ilp_elem2param;//nelements^2;        
  int *ilp_map							= pm->ilp_map;            
	double pvector[2];
	pvector[0] = pm->pvector[0];
	pvector[1] = pm->pvector[1];
	/*** End variables for ILP ***/

	bigint ntimestep  = pm->ntimestep;
  bigint lastcall   = pm->lastcall;
	int myrank				= pm->myrank;
  
	double fxtmp,fytmp,fztmp;
	atom_in *iatom, *jatom;

	int ist, ied, isz, ioff;
	int jst, jed, jsz, joff;

  // loop over full neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) 
	{//ILP+Tersoff
		i = ilist[ii];

		#ifdef MY_DEBUG
		if(i != ii)
			printf("i = %d, ii = %d\n", i, ii);
		#endif
		iatom = &my_atoms[i];
		xtmp = iatom->x[0];
    ytmp = iatom->x[1];
    ztmp = iatom->x[2];
    int itypem = ilp_map[iatom->type];
    itag = iatom->tag;
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int nilp = 0;
	
		//sort jlist;
		//bubbleSort(jlist, jnum, iatom->molecule, molecule);

		//Gather neighbors in the same layer;
		//int pi, pj;
		//pi = 0; pj = 0;
		//int tmp, numlayer = 0;
		//for(jj = 1; jj < jnum; jj++)
		//{
		//	j = jlist[jj];
		//	if(molecule[j] == iatom->molecule)
		//	{
		//		tmp = jlist[pi];
		//		jlist[pi] = jlist[jj];
		//		jlist[jj] = tmp;
		//		pi++;
		//		numlayer++;
		//	}
		//}

		int numlayer = nummol[i];

		//TERSOFF
    itype = map[iatom->type];
    int numshort = 0;
    fxtmp = fytmp = fztmp = 0.0;
		//mytest:
		int samelayer = 0;
    //for (jj = 0; jj < jnum; jj++) 
    for (jj = 0; jj < numlayer; jj++) 
		{
      j = jlist[jj];
			jatom = &my_atoms[j];
      int jtypem = ilp_map[jatom->type];
      jtag = jatom->tag;

      delx = xtmp - jatom->x[0];
      dely = ytmp - jatom->x[1];
      delz = ztmp - jatom->x[2];
      rsq = delx*delx + dely*dely + delz*delz;

			//if(pm->myrank == 0)
			//{
			////	if(iatom->molecule ==  jatom->molecule)
			////		samelayer++;

			//	printf("i=%d, j=%d, imol=%d, jmol=%d\n", i, j, iatom->molecule, jatom->molecule);
			//}

      if(rsq != 0 && rsq < ilp_cutILPsq[itypem*nelements+jtypem] && iatom->molecule == jatom->molecule)
			{
        ilp_neigh[nilp] = j;
        vet[nilp][0] = -delx;
        vet[nilp][1] = -dely;
        vet[nilp][2] = -delz;
        nilp ++;
      }


			//TERSOFF
      if(!iskip[iatom->type] && !ijskip[iatom->type][jatom->type])
			{
				if (rsq < cutshortsq) 
				{
					neighBf[numshort].rsq = rsq;
					neighBf[numshort].r = sqrt(rsq);
					neighBf[numshort].rinv = 1/neighBf[numshort].r;
					neighBf[numshort].delr[0] = -delx;
					neighBf[numshort].delr[1] = -dely;
					neighBf[numshort].delr[2] = -delz;
					neighBf[numshort].tp = map[jatom->type];
					neighBf[numshort].padding = 0;

					int ters_flag = 1;
					if (itag > jtag) 
					{
    			  if ((itag+jtag) % 2 == 0) {ters_flag = 0; }
    			} 
					else if (itag < jtag) 
					{
    			  if ((itag+jtag) % 2 == 1) {ters_flag = 0; }
    			} 
					else 
					{
    			  if (jatom->x[2] < iatom->x[2]) {ters_flag = 0;}
    			  if (jatom->x[2] == iatom->x[2] && jatom->x[1] <  iatom->x[1]) {ters_flag = 0; }
    			  if (jatom->x[2] == iatom->x[2] && jatom->x[1] == iatom->x[1] && jatom->x[0] < iatom->x[0]) {ters_flag = 0;}
    			}
					neighBf[numshort].ters_flag = ters_flag;


      	  neighshort[numshort++] = j;
					#ifdef MY_DEBUG
      	  if (numshort >= maxshort && numshort > 4) 
					{
						printf("numshort > maxshot, numshort=%d\n", numshort);
      	  }
					#endif
      	}
			}//if-notSkip
    }//for-jj

		//if(myrank == 0)
		//{
		//	printf("i=%d, jnum=%d, samelayer=%d\n", i, jnum, samelayer);
		//}

		calc_normal(nilp, jnum, vet, normal, dnormal, dnormdri);
		
		/* Debugging */
		#ifdef MY_DEBUG
		if(pm->myrank == 0)
		{
			printf("i=%d, ntimestep=%d,lastcall=%d,", i,ntimestep,lastcall);
			int tt;
			for(tt = 0; tt < nilp; tt++) printf(" %d", ilp_neigh[tt]);
			printf(";");
			for(tt = 0; tt < numshort; tt++) printf(" %d", neighshort[tt]);
			printf("\n");
		}
		#endif
    
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

    //for (jj = 0; jj < jnum; jj++) 
    for (jj = numlayer; jj < jnum; jj++) 
		{
      j = jlist[jj];
      j &= NEIGHMASK;
			jatom = &my_atoms[j];

      jtag = jatom->tag;
      delx = xtmp - jatom->x[0];
      dely = ytmp - jatom->x[1];
      delz = ztmp - jatom->x[2];
      rsq = delx*delx + dely*dely + delz*delz;

      int vdwflag = 1;
			if (itag > jtag) 
			{
    	  if ((itag+jtag) % 2 == 0) {vdwflag = 0;}
    	} 
			else if (itag < jtag) 
			{
    	  if ((itag+jtag) % 2 == 1) {vdwflag = 0;}
    	} 
			else 
			{
    	  if (jatom->x[2] < iatom->x[2]) {vdwflag = 0;}
    	  if (jatom->x[2] == ztmp && iatom->x[1] < ytmp) {vdwflag = 0;}
    	  if (jatom->x[2] == ztmp && iatom->x[1] == ytmp && iatom->x[0] < xtmp) {vdwflag = 0;}
    	}

      jtype = map[jatom->type];
			r = sqrt(rsq);
      double rinv = 1.0 / r;
			
      if (rsq < ilp_cutsq[iatom->type*nt+jatom->type] && iatom->molecule != jatom->molecule) 
			{
        int iparam_ij = ilp_elem2param[ilp_map[iatom->type] *nelements+ ilp_map[jatom->type]];
        Param_ilp_str *p = &ilp_params[iparam_ij];

        if (ilp_tap_flag) 
				{
          Rcut = sqrt(ilp_cutsq[iatom->type*nt+jatom->type]);
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
          delkj[0] = iatom->x[0] + vet[kk][0] - jatom->x[0];
          delkj[1] = iatom->x[1] + vet[kk][1] - jatom->x[1];
          delkj[2] = iatom->x[2] + vet[kk][2] - jatom->x[2];
					
          if (evflag) 
					{
						if (vflag_either && vflag_global) 
						{
						  p_virial[0] += delkj[0]*fk[0];
						  p_virial[1] += delkj[1]*fk[1];
						  p_virial[2] += delkj[2]*fk[2];
						  p_virial[3] += delkj[0]*fk[1];
						  p_virial[4] += delkj[0]*fk[2];
						  p_virial[5] += delkj[1]*fk[2];
						}//ec_tally_buffer_simple;
					}

        }
        //erep = Tap*Vilp;
        if (eflag) pvector[1] += erep = Tap*Vilp;
        if (evflag) 
				{
					if (eflag_either && eflag_global) 
						*p_eng_vdwl += (erep+evdw);
  				if (vflag_either && vflag_global) 
					{
						p_virial[0] += delx*ftotx;
  				  p_virial[1] += dely*ftoty;
  				  p_virial[2] += delz*ftotz;
  				  p_virial[3] += delx*ftoty;
  				  p_virial[4] += delx*ftotz;
  				  p_virial[5] += dely*ftotz;
  				}//ev_tally_xyz_simple

				}
      }
    } // loop over jj
    for (kk = 0; kk < nilp; kk ++) 
		{
      int k = ilp_neigh[kk];
      f[k][0] += fkk[kk][0];
      f[k][1] += fkk[kk][1];
      f[k][2] += fkk[kk][2];
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
			jtype = neighBf[jj].tp;
      iparam_ij = elem2param[(itype*nelements+jtype)*nelements+jtype];
			dlr1 = neighBf[jj].delr;
			rsq1 = neighBf[jj].rsq;
			r1 = neighBf[jj].r;
			r1inv = neighBf[jj].rinv;
      if (rsq1 >= ters_params[iparam_ij].cutsq) continue;

      fjxtmp = fjytmp = fjztmp = 0.0;
			
			//Repulsive for Tersoff
			double Rep_evdwl = 0, Rep_fpair = 0;
			double fc, fc_d, tmp_exp;
			Param *param = &ters_params[iparam_ij];
			if (r1 < param->bigr-param->bigd) {fc = 1.0; fc_d = 0.0;}
			else if (r1 > param->bigr+param->bigd) {fc = 0.0; fc_d = 0.0;}
			else 
			{
				fc = 0.5*(1.0 - sin(MY_PI2*(r1 - param->bigr)*param->bigdinv )); 
				fc_d =  -(MY_PI4*param->bigdinv) * cos(MY_PI2*(r1 - param->bigr)*param->bigdinv);
			}

  		tmp_exp = exp(-param->lam1 * r1);
  		fpair = -param->biga * tmp_exp * (fc_d - fc*param->lam1) * r1inv;
  		if (eflag) evdwl = fc * param->biga * tmp_exp;
			
			if(neighBf[jj].ters_flag)
			{
				//ters_FRep(&ters_params[iparam_ij],r1,r1inv,&fpair,eflag,&evdwl);

      	fxtmp		+= -dlr1[0]*fpair;
      	fytmp 	+= -dlr1[1]*fpair;
      	fztmp 	+= -dlr1[2]*fpair;
				fjxtmp	-= -dlr1[0]*fpair;
      	fjytmp	-= -dlr1[1]*fpair;
      	fjztmp 	-= -dlr1[2]*fpair;

				Rep_evdwl = evdwl;
				Rep_fpair = fpair;
			}//if-ters_flag-notSkip


      // accumulate bondorder zeta for each i-j interaction via loop over k

      //fjxtmp = fjytmp = fjztmp = 0.0;
      zeta_ij = 0.0;

      for (kk = 0; kk < numshort; kk++) 
			{
        if (jj == kk) continue;
        k = neighshort[kk];
				ktype = neighBf[kk].tp;
				iparam_ijk = elem2param[(itype*nelements+jtype)*nelements+ktype];
				dlr2 = neighBf[kk].delr;
				rsq2 = neighBf[kk].rsq;
				r2 = neighBf[kk].r;
				r2inv = neighBf[kk].rinv;
      	if (rsq2 >= ters_params[iparam_ijk].cutsq) continue;
        zeta_ij += ters_zeta(&ters_params[iparam_ijk],r1,r1inv,r2,r2inv,dlr1,dlr2);
      }

      // pairwise force due to zeta
      //ters_force_zeta(&ters_params[iparam_ij],r1,r1inv,zeta_ij,&fpair,&prefactor,eflag,&evdwl);
      ters_force_zeta_fc(&ters_params[iparam_ij],r1,r1inv,zeta_ij,
													&fpair,&prefactor,eflag,&evdwl,fc,fc_d);

			fxtmp		+= dlr1[0]*fpair;
      fytmp 	+= dlr1[1]*fpair;
      fztmp 	+= dlr1[2]*fpair;
      fjxtmp	-= dlr1[0]*fpair;
      fjytmp 	-= dlr1[1]*fpair;
      fjztmp 	-= dlr1[2]*fpair;


      // attractive term via loop over k
      for (kk = 0; kk < numshort; kk++) 
			{
        if (jj == kk) continue;
        k = neighshort[kk];
				ktype = neighBf[kk].tp;
				iparam_ijk = elem2param[(itype*nelements+jtype)*nelements+ktype];
				dlr2 = neighBf[kk].delr;
				rsq2 = neighBf[kk].rsq;
				r2 = neighBf[kk].r;
				r2inv = neighBf[kk].rinv;
      	if (rsq2 >= ters_params[iparam_ijk].cutsq) continue;
        ters_Att(&ters_params[iparam_ijk],prefactor,r1,r1inv,r2,r2inv,dlr1,dlr2,fi,fj,fk);

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];
      }
			
			f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
			
			if(evflag)
			{
				if (eflag_either && eflag_global) 
				{
  			  *p_eng_vdwl += Rep_evdwl;
  			  *p_eng_vdwl += evdwl;
				}
  			if (vflag_either && vflag_global) 
				{
					p_virial[0] += -dlr1[0]*dlr1[0]*Rep_fpair;
  			  p_virial[1] += -dlr1[1]*dlr1[1]*Rep_fpair;
  			  p_virial[2] += -dlr1[2]*dlr1[2]*Rep_fpair;
  			  p_virial[3] += -dlr1[0]*dlr1[1]*Rep_fpair;
  			  p_virial[4] += -dlr1[0]*dlr1[2]*Rep_fpair;
  			  p_virial[5] += -dlr1[1]*dlr1[2]*Rep_fpair;

					p_virial[0] += -dlr1[0]*dlr1[0]*fpair;
  			  p_virial[1] += -dlr1[1]*dlr1[1]*fpair;
  			  p_virial[2] += -dlr1[2]*dlr1[2]*fpair;
  			  p_virial[3] += -dlr1[0]*dlr1[1]*fpair;
  			  p_virial[4] += -dlr1[0]*dlr1[2]*fpair;
  			  p_virial[5] += -dlr1[1]*dlr1[2]*fpair;
  			}
			}//ev_tally_simple;
			
    }//for-jj
		
		f[i][0] += fxtmp;
		f[i][1] += fytmp;
    f[i][2] += fztmp;
		
  }//for-ii

	pm->pvector[0] = pvector[0];
	pm->pvector[1] = pvector[1];
}

#endif






#ifdef CPE
#define LWPF_UNIT U(LMFF)
#define LWPF_KERNELS K(ALL) K(BEFORE) K(CMP) K(SHORT) K(ILP) K(UPDTJ) K(UPDTK)   K(TERSOFF) K(ZIJ) K(ZETA) K(ATT) K(UPDTKI) K(FLUSH) K(JLOOP) K(PART) K(PART1) K(PART2)
#include "lwpf2.h"

#define SAFER
#define ISTEP 64
#define JSTEP 64
#define PSTEP 1024
#define FSTEP 1024

void pack_atoms_lmff_para(pair_lmff_param_t *pm)
{
	dma_init();
  pair_lmff_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(pair_lmff_param_t));
  dma_syn();

  int i, ii, j;
  int inum              = l_pm.inum;
  int gnum              = l_pm.gnum;
  int *type             = l_pm.type;
  tagint *tag           = l_pm.tag;
  double (*x)[3]        = l_pm.x;
  tagint *molecule      = l_pm.molecule;
  atom_in *my_atoms     = l_pm.my_atoms;
  int allnum            = inum + gnum;
  int ist, ied, isz, ioff;
  double xi[PSTEP][3];
  int tpi[PSTEP];
	tagint mi[PSTEP];
  tagint tgi[PSTEP];
  atom_in ai[PSTEP];
  for(ist = _MYID * PSTEP; ist < allnum; ist += 64*PSTEP)
  {
    ied = ist + PSTEP;
    if(ied > allnum)
      ied = allnum;
    isz = ied - ist;
    pe_get(type     +ist, tpi,	sizeof(int)*isz);
    pe_get(molecule +ist, mi,   sizeof(tagint)*isz);
    pe_get(tag      +ist, tgi,  sizeof(tagint)*isz);
    pe_get(l_pm.x[ist],   xi,   sizeof(double)*isz*3);
    dma_syn();

    for(i = 0; i < isz; i++)
    {
      ioff = i - ist;
      ai[i].x[0]      = xi[i][0];
      ai[i].x[1]      = xi[i][1];
      ai[i].x[2]      = xi[i][2];
      ai[i].type      = tpi[i];//type[i];
      ai[i].padding   = 0;
      ai[i].tag       = tgi[i];//tag[i];
      ai[i].molecule  = mi[i];//molecule[i];
      ai[i].pad20			= 0.0;
      ai[i].pad21 		= 0.0;
    }//for-i
    pe_put(my_atoms+ist, ai, sizeof(atom_in)*isz);
    dma_syn();
  }//for-ist
}

void copy_force_back_para(pair_lmff_param_t *pm)
{
	dma_init();
  pair_lmff_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(pair_lmff_param_t));
  dma_syn();

  int nall  = l_pm.inum+l_pm.gnum;
  int i, ist, ied, isz, ioff;
  double fi[FSTEP][3], ifrc[FSTEP][4];
  
  for(ist = _MYID*FSTEP; ist < nall; ist+=64*FSTEP)
  {
    ied = ist + FSTEP;
    if(ied > nall)
      ied = nall;
    isz = ied - ist;
    pe_get(l_pm.f[ist],    fi,   sizeof(double)*isz*3);
    pe_get(l_pm.frc[ist],  ifrc, sizeof(double)*isz*4);
    dma_syn();
    for(i = ist; i < ied; i++)
    {
      ioff = i - ist;
      fi[ioff][0] += ifrc[ioff][0];
      fi[ioff][1] += ifrc[ioff][1];
      fi[ioff][2] += ifrc[ioff][2];
    }
    pe_put(l_pm.f[ist],    fi,   sizeof(double)*isz*3);
    dma_syn();
  }
}


/* filtering method */
#include "pair_lmff_cpe_onelonglist.h"
#include "pair_lmff_cpe_twolonglist.h"


/* not filtering shortneighbors */
void compute_lmff_para(pair_lmff_param_t *pm)
{
  lwpf_enter(LMFF);
  lwpf_start(ALL);

  lwpf_start(BEFORE);
  dma_init();
  pair_lmff_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(pair_lmff_param_t));
  dma_syn();

  int i,j,k,ii,jj,kk,inum,gnum,jnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3];
  double zeta_ij,prefactor;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;


  double (*x)[3]				= l_pm.x;
  double (*f)[3]				= l_pm.f;
  tagint *tag						= l_pm.tag;
  tagint *molecule			= l_pm.molecule;
  int *type							= l_pm.type;
  int ntypes						= l_pm.ntypes;
  int nlocal						= l_pm.nlocal;
  int newton_pair				= l_pm.newton_pair;
  double cutshortsq 		= l_pm.cutshortsq;
  int nelements					= l_pm.nelements;       
  int nparams						= l_pm.nparams;         
  Param *ters_params		= l_pm.ters_params;       
  int maxshort					= l_pm.maxshort;
  atom_in *my_atoms			= l_pm.my_atoms;

  int eflag						= l_pm.eflag;
  int evflag					= l_pm.evflag;
  int eflag_either		= l_pm.eflag_either;
  int eflag_global		= l_pm.eflag_global;
  int eflag_atom			= l_pm.eflag_atom;
  int vflag						= l_pm.vflag;
  int vflag_either		= l_pm.vflag_either;
  int vflag_global		= l_pm.vflag_global;
  int vflag_atom			= l_pm.vflag_atom;
	
  inum				= l_pm.inum;
  gnum				= l_pm.gnum;
  ilist				= l_pm.ilist;
  numneigh		= l_pm.numneigh;
  firstneigh	= l_pm.firstneigh;
	
  //buffer atom info for neighshort;
  ParamBf neighBf[20];
  int neighshort[20];

  /*** variables for ILP ***/
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
  int ilp_nelements					= l_pm.ilp_nelements;       
  int ilp_nparams						= l_pm.ilp_nparams;         
  int ilp_tap_flag					= l_pm.ilp_tap_flag;
  Param_ilp_str *ilp_params = l_pm.ilp_params;       
	
  double eng_virial[16];
  for(i = 0; i < 16; i++)
    eng_virial[i] = 0;
  double *p_eng_vdwl	= eng_virial;
  double *p_eng_coul	= p_eng_vdwl + 1; 
  double *p_virial		= p_eng_coul + 1;
  double *pvector			= p_virial	 + 6;
	
  bigint ntimestep  = l_pm.ntimestep;
  bigint lastcall   = l_pm.lastcall;
  int myrank				= l_pm.myrank;
  long *mylock			= l_pm.mylock;
  rvec4 *frc				= l_pm.frc;
  //write_cache;
  int wfctag[WRT_C_LCNT];
  rvec4 wfcache[WRT_C_LCNT][WRT_C_LSZ];
  for(i = 0; i < WRT_C_LCNT; i++)
    wfctag[i] = -1;
  rvec4 bwf[WRT_C_LSZ];
  rvec4 fwi, fwj;//, fwk;//ILP+Tersoff;
  fwi[0] = fwi[1] = fwi[2] = fwi[3] = 0;
  fwj[0] = fwj[1] = fwj[2] = fwj[3] = 0;

  rvec4 fwk[20];
  //read_cache;
  int read_ctag[READ_C_LCNT];
  atom_in atoms_cache[READ_C_LCNT][READ_C_LSZ];
  for(i = 0; i < READ_C_LCNT; i++)
    read_ctag[i] = -1;
	
  double fxtmp,fytmp,fztmp;
  atom_in *iatom, *jatom;

  int ist, ied, isz, ioff;
  int jst, jed, jsz, joff;
  int ni[ISTEP], *fni[ISTEP];
  int jlist_buf[JSTEP];
  atom_in ai[ISTEP];


  int nt = ntypes+1;
  int ne = nelements;
  int ne3 = ne * ne * ne;
  int map[nt], elem2param[ne*ne*ne];
  int iskip[nt], ijskip[nt][nt];
	
  int ilp_ne2 = ilp_nelements * ilp_nelements;
  int ilp_map[nt], ilp_elem2param[ilp_ne2];
  double ilp_cutsq[nt*nt], ilp_cutILPsq[ilp_ne2];
  pe_get(l_pm.map,									map,									sizeof(int)*nt);
  pe_get(&(l_pm.elem2param[0]),			&(elem2param[0]),			sizeof(int)*ne3);
  pe_get(&(l_pm.iskip[0]),					&(iskip[0]),					sizeof(int)*nt);
  pe_get(&(l_pm.ijskip[0]),					&(ijskip[0][0]),			sizeof(int)*nt*nt);

  pe_get(l_pm.ilp_map,        			ilp_map,							sizeof(int)*nt);
  pe_get(&(l_pm.ilp_cutsq[0]),      &(ilp_cutsq[0]),			sizeof(double)*nt*nt);
  pe_get(&(l_pm.ilp_cutILPsq[0]),		&(ilp_cutILPsq[0]),		sizeof(double)*ilp_ne2);
  pe_get(&(l_pm.ilp_elem2param[0]), &(ilp_elem2param[0]), sizeof(int)*ilp_ne2);
  dma_syn();
	
  lwpf_stop(BEFORE);

  lwpf_start(CMP);
  // loop over full neighbor list of my atoms
  for (ist = _MYID*ISTEP; ist < inum; ist += 64*ISTEP) 
    {
      ied = ist + ISTEP;
      if(ied > inum)
	ied = inum;
      isz = ied -ist;
	
      pe_get(firstneigh+ist,  fni,	sizeof(int*)*isz);
      pe_get(numneigh+ist,    ni,   sizeof(int)*isz);
      pe_get(my_atoms+ist,    ai,   sizeof(atom_in)*isz);
      dma_syn();
      for (ii = ist; ii < ied; ii++) 
	{//ILP+Tersoff
	  //i = ilist[ii];
	  i = ii;
	  ioff = ii - ist;
	  iatom = &ai[ioff];
	  xtmp = iatom->x[0];
  	  ytmp = iatom->x[1];
  	  ztmp = iatom->x[2];
  	  int itypem = ilp_map[iatom->type];
  	  itag = iatom->tag;
  	  
	  jlist = fni[ioff];
  	  jnum = ni[ioff];
  	  int nilp = 0;
			
	  //TERSOFF
  	  itype = map[iatom->type];
  	  int numshort = 0;
  	  fxtmp = fytmp = fztmp = 0.0;
		
	  fwi[0] = fwi[1] = fwi[2] = fwi[3] = 0;

	  lwpf_start(SHORT);
	  for (jst = 0; jst < jnum; jst+=JSTEP) 
	    {
	      jed = jst + JSTEP;
	      if(jed > jnum)
		jed = jnum;
	      jsz = jed - jst;
			
	      pe_get(jlist+jst, jlist_buf, sizeof(int)*jsz);
	      dma_syn();
			
	      for (jj = jst; jj < jed; jj++) 
		{
		  joff = jj - jst;
  	  	  //j = jlist[jj];
  	  	  j = jlist_buf[joff];
		  if (read_ctag[(j >> READ_C_S) & READ_C_LM] != j >> READ_C_S)
		    {
		      pe_get(my_atoms + (j & ~READ_C_MM), 
			     atoms_cache[(j >> READ_C_S) & READ_C_LM], 
			     sizeof(atom_in) * READ_C_LSZ);
		      dma_syn();
		      read_ctag[(j >> READ_C_S) & READ_C_LM] = j >> READ_C_S;
		    }
		  jatom = &atoms_cache[(j >> READ_C_S) & READ_C_LM][j & READ_C_MM];
		  //jatom = &my_atoms[j];

  	  	  int jtypem = ilp_map[jatom->type];
  	  	  jtag = jatom->tag;

  	  	  delx = xtmp - jatom->x[0];
  	  	  dely = ytmp - jatom->x[1];
  	  	  delz = ztmp - jatom->x[2];
  	  	  rsq = delx*delx + dely*dely + delz*delz;

  	  	  if(rsq != 0 && rsq < ilp_cutILPsq[itypem*nelements+jtypem] && iatom->molecule == jatom->molecule)
		    {
		      ilp_neigh[nilp] = j;
		      vet[nilp][0] = -delx;
		      vet[nilp][1] = -dely;
		      vet[nilp][2] = -delz;
		      nilp ++;
		    }
					
		  //TERSOFF
  	  	  if(!iskip[iatom->type] && !ijskip[iatom->type][jatom->type])
		    {
		      if (rsq < cutshortsq) 
			{
			  neighBf[numshort].rsq = rsq;
			  neighBf[numshort].r = sqrt(rsq);
			  neighBf[numshort].rinv = 1/neighBf[numshort].r;
			  neighBf[numshort].delr[0] = -delx;
			  neighBf[numshort].delr[1] = -dely;
			  neighBf[numshort].delr[2] = -delz;
			  //neighBf[numshort].id = j;
			  neighBf[numshort].tp = map[jatom->type];
			  neighBf[numshort].padding = 0;

			  int ters_flag = 1;
			  if (itag > jtag) 
			    {
			      if ((itag+jtag) % 2 == 0) {ters_flag = 0; }
			    } 
			  else if (itag < jtag) 
			    {
			      if ((itag+jtag) % 2 == 1) {ters_flag = 0; }
			    } 
			  else 
			    {
			      if (jatom->x[2] < iatom->x[2]) {ters_flag = 0;}
			      if (jatom->x[2] == iatom->x[2] && jatom->x[1] <  iatom->x[1]) {ters_flag = 0; }
			      if (jatom->x[2] == iatom->x[2] && jatom->x[1] == iatom->x[1] && jatom->x[0] < iatom->x[0]) {ters_flag = 0;}
			    }
			  neighBf[numshort].ters_flag = ters_flag;

  	  	  	  neighshort[numshort++] = j;

#ifdef MY_DEBUG
  	  	  	  if (numshort > 4) 
			    {
			      printf("numshort > 4, numshort=%d\n", numshort);
			    }
#endif 
  	  	  	}
		    }//if-notSkip
  	  	}//for-jj
	    }//for-jst
	  lwpf_stop(SHORT);

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
	  jlist = fni[ioff];
  	  jnum = ni[ioff];

	  lwpf_start(ILP);
	  for (jst = 0; jst < jnum; jst+=JSTEP) 
	    {
	      jed = jst + JSTEP;
	      if(jed > jnum)
		jed = jnum;
	      jsz = jed - jst;
				
	      pe_get(jlist+jst, jlist_buf, sizeof(int)*jsz);
	      dma_syn();
				
	      for (jj = jst; jj < jed; jj++) 
		{
		  joff = jj - jst;
  	  	  j = jlist_buf[joff];
  	  	  j &= NEIGHMASK;
		  if (read_ctag[(j >> READ_C_S) & READ_C_LM] != j >> READ_C_S)
		    {
		      pe_get(my_atoms + (j & ~READ_C_MM), 
			     atoms_cache[(j >> READ_C_S) & READ_C_LM], 
			     sizeof(atom_in) * READ_C_LSZ);
		      dma_syn();
		      read_ctag[(j >> READ_C_S) & READ_C_LM] = j >> READ_C_S;
		    }
		  jatom = &atoms_cache[(j >> READ_C_S) & READ_C_LM][j & READ_C_MM];
		
		  fwj[0] = fwj[1] = fwj[2] = fwj[3] = 0;

  	  	  jtag = jatom->tag;
  	  	  delx = xtmp - jatom->x[0];
  	  	  dely = ytmp - jatom->x[1];
  	  	  delz = ztmp - jatom->x[2];
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
		      if (jatom->x[2] < iatom->x[2]) {ters_flag = 0; vdwflag = 0;}
		      if (jatom->x[2] == ztmp && iatom->x[1] < ytmp) {ters_flag = 0; vdwflag = 0;}

		      if (jatom->x[2] == ztmp && iatom->x[1] == ytmp && iatom->x[0] < xtmp) {ters_flag = 0; vdwflag = 0;}
		    }

  	  	  jtype = map[jatom->type];
		  r = sqrt(rsq);
  	  	  double rinv = 1.0 / r;
					
		  if(ters_flag && !iskip[jatom->type] && !ijskip[iatom->type][jatom->type])
		    {
		      iparam_ij = elem2param[(itype*nelements+jtype)*nelements+jtype];
		      if (rsq < ters_params[iparam_ij].cutsq)
			{
			  ters_FRep(&ters_params[iparam_ij],r,rinv,&fpair,eflag,&evdwl);

			  fxtmp += delx*fpair;
			  fytmp += dely*fpair;
			  fztmp += delz*fpair;
			  fwj[0] -= delx*fpair;
			  fwj[1] -= dely*fpair;
			  fwj[2] -= delz*fpair;

							
			  if(evflag)
			    {
			      if (eflag_either && eflag_global) 
				*p_eng_vdwl += evdwl;
			      if (vflag_either && vflag_global) 
				{
				  p_virial[0] += delx*delx*fpair;
				  p_virial[1] += dely*dely*fpair;
				  p_virial[2] += delz*delz*fpair;
				  p_virial[3] += delx*dely*fpair;
				  p_virial[4] += delx*delz*fpair;
				  p_virial[5] += dely*delz*fpair;
				}
			    }//ev_tally_simple;

			}
		    }//if-notSkip



  	  	  if (rsq < ilp_cutsq[iatom->type*nt+jatom->type] && iatom->molecule != jatom->molecule) 
		    {

		      int iparam_ij = ilp_elem2param[ilp_map[iatom->type] *nelements+ ilp_map[jatom->type]];
		      Param_ilp_str *p = &ilp_params[iparam_ij];

		      if (ilp_tap_flag) 
			{
			  Rcut = sqrt(ilp_cutsq[iatom->type*nt+jatom->type]);
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
							
			  TSvdw = 1.0 + xexp(-p->d*(r/p->seff - 1.0));

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
		      exp0 = xexp(-p->lambda*(r-p->z0));
		      exp1 = xexp(-rdsq1);

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
		      fwi[0] += ftotx - fprod1[0]*Tap;
		      fwi[1] += ftoty - fprod1[1]*Tap;
		      fwi[2] += ftotz - fprod1[2]*Tap;
		      fwj[0] -= ftotx;
		      fwj[1] -= ftoty;
		      fwj[2] -= ftotz;
						

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
			  delkj[0] = iatom->x[0] + vet[kk][0] - jatom->x[0];
			  delkj[1] = iatom->x[1] + vet[kk][1] - jatom->x[1];
			  delkj[2] = iatom->x[2] + vet[kk][2] - jatom->x[2];
							
			  if (evflag) 
			    {

			      if (vflag_either && vflag_global) 
				{
				  p_virial[0] += delkj[0]*fk[0];
				  p_virial[1] += delkj[1]*fk[1];
				  p_virial[2] += delkj[2]*fk[2];
				  p_virial[3] += delkj[0]*fk[1];
				  p_virial[4] += delkj[0]*fk[2];
				  p_virial[5] += delkj[1]*fk[2];
				}//ec_tally_buffer_simple;
			    }

			}

		      //erep = Tap*Vilp;
		      if (eflag) pvector[1] += erep = Tap*Vilp;
		      if (evflag) 
			{
			  if (eflag_either && eflag_global) 
			    *p_eng_vdwl += (erep+evdw);
			  if (vflag_either && vflag_global) 
			    {
			      p_virial[0] += delx*ftotx;
			      p_virial[1] += dely*ftoty;
			      p_virial[2] += delz*ftotz;
			      p_virial[3] += delx*ftoty;
			      p_virial[4] += delx*ftotz;
			      p_virial[5] += dely*ftotz;
			    }//ev_tally_xyz_simple
			}
		    }//if-rsq-ilp
					
		  lwpf_start(UPDTJ);
		  update_cache_ters(j, fwj, wfcache, wfctag, frc, mylock);
		  lwpf_stop(UPDTJ);
  	  	} // loop over jj

	    }
	  lwpf_stop(ILP);

	  lwpf_start(UPDTK);
	  for (kk = 0; kk < nilp; kk ++) 
	    {
	      int k = ilp_neigh[kk];
	      update_cache_ters(k, fkk[kk], wfcache, wfctag, frc, mylock);
	    }//for-kk
	  lwpf_stop(UPDTK);

	  /* fix bug */
	  for(kk = 0; kk < numshort; kk++)
	    {
	      fwk[kk][0] = 0.0;
	      fwk[kk][1] = 0.0;
	      fwk[kk][2] = 0.0;
	      fwk[kk][3] = 0.0;
	    }

		
	  lwpf_start(TERSOFF);

	  /*** Tersoff ***/
  	  // three-body interactions
  	  // skip immediately if I-J is not within cutoff
  	  double fjxtmp,fjytmp,fjztmp;
	  double *dlr1, *dlr2, r1, r2, r1inv, r2inv;	
	  //When iskip[type[i]] is ture, numshort==0;
	  //			fxtmp == 0; fytmp == 0; fztmp == 0;
  	  for (jj = 0; jj < numshort; jj++) 
	    {
	      j = neighshort[jj];
	      jtype = neighBf[jj].tp;
	      iparam_ij = elem2param[(itype*nelements+jtype)*nelements+jtype];
	      dlr1 = neighBf[jj].delr;
	      rsq1 = neighBf[jj].rsq;
	      r1 = neighBf[jj].r;
	      r1inv = neighBf[jj].rinv;
	      if (rsq1 >= ters_params[iparam_ij].cutsq) continue;

	      // accumulate bondorder zeta for each i-j interaction via loop over k

	      fjxtmp = fjytmp = fjztmp = 0.0;
	      zeta_ij = 0.0;

	      for (kk = 0; kk < numshort; kk++) 
		{
		  if (jj == kk) continue;
		  k = neighshort[kk];
		  ktype = neighBf[kk].tp;
		  iparam_ijk = elem2param[(itype*nelements+jtype)*nelements+ktype];
		  dlr2 = neighBf[kk].delr;
		  rsq2 = neighBf[kk].rsq;
		  r2 = neighBf[kk].r;
		  r2inv = neighBf[kk].rinv;
		  if (rsq2 >= ters_params[iparam_ijk].cutsq) continue;
		  zeta_ij += ters_zeta(&ters_params[iparam_ijk],r1,r1inv,r2,r2inv,dlr1,dlr2);
		}

	      // pairwise force due to zeta
	      ters_force_zeta(&ters_params[iparam_ij],r1,r1inv,zeta_ij,&fpair,&prefactor,eflag,&evdwl);
				
	      fxtmp += dlr1[0]*fpair;
	      fytmp += dlr1[1]*fpair;
	      fztmp += dlr1[2]*fpair;
	      fjxtmp -= dlr1[0]*fpair;
	      fjytmp -= dlr1[1]*fpair;
	      fjztmp -= dlr1[2]*fpair;

	      if(evflag)
		{
		  if (eflag_either && eflag_global) 
		    *p_eng_vdwl += evdwl;
		  if (vflag_either && vflag_global) 
		    {
		      p_virial[0] += -dlr1[0]*dlr1[0]*fpair;
		      p_virial[1] += -dlr1[1]*dlr1[1]*fpair;
		      p_virial[2] += -dlr1[2]*dlr1[2]*fpair;
		      p_virial[3] += -dlr1[0]*dlr1[1]*fpair;
		      p_virial[4] += -dlr1[0]*dlr1[2]*fpair;
		      p_virial[5] += -dlr1[1]*dlr1[2]*fpair;
		    }
		}//ev_tally_simple;

	      // attractive term via loop over k
	      for (kk = 0; kk < numshort; kk++) 
		{
		  if (jj == kk) continue;
		  k = neighshort[kk];
		  ktype = neighBf[kk].tp;
		  iparam_ijk = elem2param[(itype*nelements+jtype)*nelements+ktype];
		  dlr2 = neighBf[kk].delr;
		  rsq2 = neighBf[kk].rsq;
		  r2 = neighBf[kk].r;
		  r2inv = neighBf[kk].rinv;
		  if (rsq2 >= ters_params[iparam_ijk].cutsq) continue;
		  ters_Att(&ters_params[iparam_ijk],prefactor,r1,r1inv,r2,r2inv,dlr1,dlr2,fi,fj,fk);

		  fxtmp += fi[0];
		  fytmp += fi[1];
		  fztmp += fi[2];
		  fjxtmp += fj[0];
		  fjytmp += fj[1];
		  fjztmp += fj[2];
		  fwk[kk][0] += fk[0];
		  fwk[kk][1] += fk[1];
		  fwk[kk][2] += fk[2];
		}
	      fwk[jj][0] += fjxtmp;
	      fwk[jj][1] += fjytmp;
	      fwk[jj][2] += fjztmp;
	    }
			
	  lwpf_stop(TERSOFF);
		
	  lwpf_start(UPDTKI);
	  for (kk = 0; kk < numshort; kk ++) 
	    {
	      int k = neighshort[kk];
	      update_cache_ters(k, fwk[kk], wfcache, wfctag, frc, mylock);
	    }//for-kk

			
	  fwi[0] += fxtmp;
	  fwi[1] += fytmp;
  	  fwi[2] += fztmp;
			
	  update_cache_ters(i, fwi, wfcache, wfctag, frc, mylock);
	  lwpf_stop(UPDTKI);
			
  	}//for-ii
    }//for-ist
  lwpf_stop(CMP);

  lwpf_start(FLUSH);
  flush_cache_ters(frc, wfcache, wfctag, mylock);
  lwpf_stop(FLUSH);

  pe_put(l_pm.pack64+_MYID*16, eng_virial, sizeof(double)*16);
  dma_syn();
	
  lwpf_stop(ALL);
  lwpf_exit(LMFF);
}
#endif
