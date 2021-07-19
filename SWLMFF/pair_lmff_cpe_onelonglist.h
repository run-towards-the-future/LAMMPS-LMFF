void compute_lmff_onelonglist(pair_lmff_param_t *pm)
{
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
	int neighshort[20], neighindex[20];

	/* whether to recalculate */
	int cur_neighshort[20], cur_numshort;
	int cur_ilpneigh[6], cur_nilp;
	int recalc = l_pm.recalc;


	/*** variables for ILP ***/
  int l, ll;
  double erep, r,Rcut,r2inv,r6inv,r8inv;
	double Tap,dTap,Vilp,TSvdw,TSvdwinv,TSvdw2inv,fsum;
  erep = 0.0;
	//Vars for calc_normal
  int id,ip,m;
  double nn, nn2;
  double pv12[3],pv31[3],pv23[3],n1[3],dni[3],dnn[3][3],vet[3][3],dpvdri[3][3];
  double dn1[3][3][3],dpv12[3][3][3],dpv23[3][3][3],dpv31[3][3][3];
  double normal[3], dnormal[3][3][3], dnormdri[3][3];
  //Vars for calc_Frep
  double prodnorm1,fkcx,fkcy,fkcz, pdf1[3];
  double rhosq1,exp0,exp1;
  double frho1,Erep,rdsq1,fpair1;
  double dprodnorm1[3] = {0.0, 0.0, 0.0};
  double fp1[3] = {0.0, 0.0, 0.0};
  double fprod1[3] = {0.0, 0.0, 0.0};
  double delkj[3] = {0.0, 0.0, 0.0};
  double fk[3] = {0.0, 0.0, 0.0};
  //more for overflow
  int ilp_neigh[6], ilp_index[6];
	//transfer from force cls;
	int ilp_nelements					= l_pm.ilp_nelements;       
  int ilp_nparams						= l_pm.ilp_nparams;         
	int ilp_tap_flag					= l_pm.ilp_tap_flag;
	Param_ilp_str *ilp_params = l_pm.ilp_params;       
	/*** End variables for ILP ***/

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
	rvec4 fwi, fwj;
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

	int jshort, tot_short, jorg;
	int flt_shortneigh[ISTEP][4];			

	int nt = ntypes+1;
	int ne = nelements;
	int ne3 = ne * ne * ne;
	int map[nt], elem2param[ne3];
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
	
	
	int my_start, my_end;
	int nipe		= inum / 64;
	int nextra = inum % 64;
	if(_MYID < nextra)
	{
		my_start = (nipe+1) * _MYID;
		my_end		= my_start + nipe + 1;
	}
	else
	{
		my_start = nipe * _MYID + nextra;
		my_end		= my_start + nipe;
	}
  
	// loop over full neighbor list of my atoms
  for (ist = my_start; ist < my_end; ist += ISTEP) 
	{
		ied = ist + ISTEP;
    if(ied > my_end)
      ied = my_end;
    isz = ied -ist;
	
	 	pe_get(l_pm.flt_shortneigh+ist*4,	&flt_shortneigh[0][0],sizeof(int)*isz*4);
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


			/* filtering shortneighbors instead of tranversing long jlist */
			tot_short = flt_shortneigh[ioff][3];
			for(jshort = 0; jshort < tot_short; jshort++)
			{
				jorg = flt_shortneigh[ioff][jshort];
				j = jorg & NEIGHMASK;

				if (read_ctag[(j >> READ_C_S) & READ_C_LM] != j >> READ_C_S)
  			{
  			  pe_get(my_atoms + (j & ~READ_C_MM), 
  			        atoms_cache[(j >> READ_C_S) & READ_C_LM], 
  			        sizeof(atom_in) * READ_C_LSZ);
  			  dma_syn();
  			  read_ctag[(j >> READ_C_S) & READ_C_LM] = j >> READ_C_S;
  			}
  			jatom = &atoms_cache[(j >> READ_C_S) & READ_C_LM][j & READ_C_MM];
				
				int jtypem = ilp_map[jatom->type];
  			jtag = jatom->tag;
  			delx = xtmp - jatom->x[0];
  			dely = ytmp - jatom->x[1];
  			delz = ztmp - jatom->x[2];
  			rsq = delx*delx + dely*dely + delz*delz;
				
				if((jorg >> 30) & 1)
				{
					ilp_neigh[nilp] = j;
					ilp_index[nilp] = jshort;
  				vet[nilp][0] = -delx;
  				vet[nilp][1] = -dely;
  				vet[nilp][2] = -delz;
					nilp++;
				}
							
  		  if(!iskip[iatom->type] && !ijskip[iatom->type][jatom->type])
				{
					neighBf[numshort].rsq = rsq;
					neighBf[numshort].r = sqrt(rsq);
					neighBf[numshort].rinv = 1/neighBf[numshort].r;
					neighBf[numshort].delr[0] = -delx;
					neighBf[numshort].delr[1] = -dely;
					neighBf[numshort].delr[2] = -delz;
					neighBf[numshort].id = j;
					neighBf[numshort].tp = map[jatom->type];
					neighBf[numshort].padding = 0;
  				neighshort[numshort] = j;
					neighindex[numshort] = jshort;
					numshort++;
				}
			}//for

			calc_normal(nilp, jnum, vet, normal, dnormal, dnormdri);
  	 
			for (kk = 0; kk < 10; kk ++) 
			{
  	    fwk[kk][0] = 0.0;
  	    fwk[kk][1] = 0.0;
  	    fwk[kk][2] = 0.0;
  	    fwk[kk][3] = 0.0;
  	  }
			jlist = fni[ioff];
  	  jnum = ni[ioff];
			
			cur_numshort = 0;
			cur_nilp = 0;
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
					
					/* whether to recalculate */
					if (rsq < cutshortsq) 
					{
						if(!iskip[iatom->type] && !ijskip[iatom->type][jatom->type])
						{
  	  	    	cur_neighshort[cur_numshort++] = j;
						}	
						
						int jtypem = ilp_map[jatom->type];
						if(rsq != 0 && rsq < ilp_cutILPsq[itypem*nelements+jtypem] 
												&& iatom->molecule == jatom->molecule)
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
  	  		  if (jatom->x[2] < iatom->x[2]) {ters_flag = 0; vdwflag = 0;}
  	  		  if (jatom->x[2] == ztmp && jatom->x[1] < ytmp) {ters_flag = 0; vdwflag = 0;}

  	  		  if (jatom->x[2] == ztmp && jatom->x[1] == ytmp && jatom->x[0] < xtmp) {ters_flag = 0; vdwflag = 0;}
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
							
							if (evflag && eflag_either && eflag_global) 
  						  *p_eng_vdwl += evdwl;
							if (evflag && vflag_either && vflag_global) 
							{
								p_virial[0] += delx*delx*fpair;
  						  p_virial[1] += dely*dely*fpair;
  						  p_virial[2] += delz*delz*fpair;
  						  p_virial[3] += delx*dely*fpair;
  						  p_virial[4] += delx*delz*fpair;
  						  p_virial[5] += dely*delz*fpair;
  						}

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
							
							TSvdw = 1.0 + exp(-p->d*(r*p->seffinv - 1.0));
							TSvdwinv = 1 / TSvdw;
							TSvdw2inv = TSvdwinv * TSvdwinv;
							double vvdw = -p->C6*r6inv*TSvdwinv;
							
							// derivatives
							fpair = -6.0*p->C6*r8inv*TSvdwinv + p->C6*p->d*p->seffinv*(TSvdw-1.0)*TSvdw2inv*r8inv*r;
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
  	  	    
						pdf1[0] = prodnorm1*fpair1*Tap;
  	  	    pdf1[1] = prodnorm1*fpair1*Tap;
  	  	    pdf1[2] = prodnorm1*fpair1*Tap;
						fkcx = (delx*fsum*Tap - normal[0]*pdf1[0]) - Vilp*dTap*delx*rinv;
  	  	    fkcy = (dely*fsum*Tap - normal[1]*pdf1[1]) - Vilp*dTap*dely*rinv;
  	  	    fkcz = (delz*fsum*Tap - normal[2]*pdf1[2]) - Vilp*dTap*delz*rinv;

  	  	    //This should be no use because fkcx need a lot of variables
  	  	    double ftotx = fvdw * delx + fkcx;
  	  	    double ftoty = fvdw * dely + fkcy;
  	  	    double ftotz = fvdw * delz + fkcz;
						
						fwi[0] += ftotx - dprodnorm1[0]*pdf1[0];
  	  	    fwi[1] += ftoty - dprodnorm1[1]*pdf1[1];
  	  	    fwi[2] += ftotz - dprodnorm1[2]*pdf1[2];
						fwj[0] -= ftotx;
  	  	    fwj[1] -= ftoty;
  	  	    fwj[2] -= ftotz;
						

  	  	    for (kk = 0; kk < nilp; kk++) 
						{
  	  	      k = ilp_neigh[kk];
  	  	      dprodnorm1[0] = dnormal[kk][0][0]*delx + dnormal[kk][1][0]*dely + dnormal[kk][2][0]*delz;
  	  	      dprodnorm1[1] = dnormal[kk][0][1]*delx + dnormal[kk][1][1]*dely + dnormal[kk][2][1]*delz;
  	  	      dprodnorm1[2] = dnormal[kk][0][2]*delx + dnormal[kk][1][2]*dely + dnormal[kk][2][2]*delz;
							fk[0] = -pdf1[0]*dprodnorm1[0];
  	  	      fk[1] = -pdf1[1]*dprodnorm1[1];
  	  	      fk[2] = -pdf1[2]*dprodnorm1[2];
							fwk[ilp_index[kk]][0] += fk[0];
  	  	      fwk[ilp_index[kk]][1] += fk[1];
  	  	      fwk[ilp_index[kk]][2] += fk[2];
							
  	  	      if (evflag && vflag_either && vflag_global) 
							{
								p_virial[0] += (delx + vet[kk][0])*fk[0];
								p_virial[1] += (dely + vet[kk][1])*fk[1];
								p_virial[2] += (delz + vet[kk][2])*fk[2];
								p_virial[3] += (delx + vet[kk][0])*fk[1];
								p_virial[4] += (delx + vet[kk][0])*fk[2];
								p_virial[5] += (dely + vet[kk][1])*fk[2];
							}//ev_tally_xyz_simple
  	  	    }

  	  	    if (eflag) pvector[1] += erep = Tap*Vilp;
						if (evflag && eflag_either && eflag_global) 
							*p_eng_vdwl += (erep+evdw);
						if (evflag && vflag_either && vflag_global) 
						{
							p_virial[0] += delx*ftotx;
  					  p_virial[1] += dely*ftoty;
  					  p_virial[2] += delz*ftotz;
  					  p_virial[3] += delx*ftoty;
  					  p_virial[4] += delx*ftotz;
  					  p_virial[5] += dely*ftotz;
  					}//ev_tally_xyz_simple
  	  	  }//if-rsq-ilp
					
					update_cache_ters(j, fwj, wfcache, wfctag, frc, mylock);
  	  	} // loop over jj

			}

			/*set recalc true,then return to recalculate */	
			if(cur_nilp != nilp || cur_numshort != numshort)
			{
				pm->recalc = 1; return;
			}
			else
			{
				int re;
				if(cur_nilp == nilp)
				for(re = 0; re < nilp; re++)
				{
					if(ilp_neigh[re] != cur_ilpneigh[re])
					{
						pm->recalc = 1; return;
					}
				}
				
				if(cur_numshort == numshort)
				for(re = 0; re < numshort; re++)
				{
					if(neighshort[re] != cur_neighshort[re])
					{
						pm->recalc = 1; return;
					}
				}
			}


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

				if (evflag && eflag_either && eflag_global) 
  			  *p_eng_vdwl += evdwl;
  			if (evflag && vflag_either && vflag_global) 
				{
					p_virial[0] += -dlr1[0]*dlr1[0]*fpair;
  			  p_virial[1] += -dlr1[1]*dlr1[1]*fpair;
  			  p_virial[2] += -dlr1[2]*dlr1[2]*fpair;
  			  p_virial[3] += -dlr1[0]*dlr1[1]*fpair;
  			  p_virial[4] += -dlr1[0]*dlr1[2]*fpair;
  			  p_virial[5] += -dlr1[1]*dlr1[2]*fpair;
  			}

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
					fwk[neighindex[kk]][0] += fk[0];
  	      fwk[neighindex[kk]][1] += fk[1];
  	      fwk[neighindex[kk]][2] += fk[2];
  	    }
  	    
				fwk[neighindex[jj]][0] += fjxtmp;
  	    fwk[neighindex[jj]][1] += fjytmp;
  	    fwk[neighindex[jj]][2] += fjztmp;
  	  }
			
			for (kk = 0; kk < flt_shortneigh[ioff][3]; kk ++) 
			{
  	    int k = flt_shortneigh[ioff][kk] & NEIGHMASK;
				update_cache_ters(k, fwk[kk], wfcache, wfctag, frc, mylock);
  	  }//for-kk
			
			fwi[0] += fxtmp;
			fwi[1] += fytmp;
  	  fwi[2] += fztmp;
			update_cache_ters(i, fwi, wfcache, wfctag, frc, mylock);
			
  	}//for-ii
	}//for-ist

	flush_cache_ters(frc, wfcache, wfctag, mylock);

	pe_put(l_pm.pack64+_MYID*16, eng_virial, sizeof(double)*16);
	dma_syn();
}
