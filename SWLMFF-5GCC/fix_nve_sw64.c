#include "sunway.h"
#include "fix_nve_sw64.h"

#ifdef MPE
extern SLAVE_FUN(fix_nve_initial_integrate_para)(fix_nve_param_t *);
extern SLAVE_FUN(fix_nve_final_integrate_para)(fix_nve_param_t *);
void fix_nve_initial_integrate_serial(fix_nve_param_t *pm)
{
  //if(flag_athread_init == 0)
	//{
	//	CRTS_init();
	//	flag_athread_init = 1;
	//}

	if(athread_idle() == 0)
		athread_init();
	athread_spawn(fix_nve_initial_integrate_para, pm);
  athread_join();
	return;
}

void fix_nve_final_integrate_serial(fix_nve_param_t *pm)
{
  //if(flag_athread_init == 0)
	//{
	//	CRTS_init();
	//	flag_athread_init = 1;
	//}

	if(athread_idle() == 0)
		athread_init();
	athread_spawn(fix_nve_final_integrate_para, pm);
  athread_join();
	return;
}
#endif




#ifdef CPE
#define ISTEP 512
void fix_nve_initial_integrate_para(fix_nve_param_t *pm)
{
	dma_init();
	fix_nve_param_t l_pm;
	pe_get(pm, &l_pm, sizeof(fix_nve_param_t));
	dma_syn();

  double dtfm;

  double (*x)[3]	= l_pm.x;
  double (*v)[3] 	= l_pm.v;
  double (*f)[3] 	= l_pm.f;
  double *rmass		= l_pm.rmass;
  double *mass		= l_pm.mass;
	double dtv			= l_pm.dtv;
	double dtf			= l_pm.dtf;
  int *type				= l_pm.type;
  int *mask 			= l_pm.mask;
  int nlocal			= l_pm.nlocal;
	int groupbit		= l_pm.groupbit;
	int ntypes			= l_pm.ntypes;

	int i, ist, ied, isz, ioff;
	double xi[ISTEP][3], vi[ISTEP][3], fi[ISTEP][3];
	int ti[ISTEP], mi[ISTEP];
	
	double massinv[ntypes+1];
	pe_get(l_pm.massinv, massinv, sizeof(double)*(ntypes+1));
	dma_syn();
  
  for (ist = _MYID*ISTEP; ist < nlocal; ist+=64*ISTEP)
	{
		ied = ist + ISTEP;
		if(ied > nlocal)
			ied = nlocal;
		isz = ied - ist;
		
		pe_get(mask+ist, mi, sizeof(int)*isz);
		pe_get(type+ist, ti, sizeof(int)*isz);
		pe_get(&x[ist][0], &xi[0][0], sizeof(double)*isz*3);
		pe_get(&v[ist][0], &vi[0][0], sizeof(double)*isz*3);
		pe_get(&f[ist][0], &fi[0][0], sizeof(double)*isz*3);
		dma_syn();
		for(i = ist; i < ied; i++)
		{
			ioff = i - ist;
  	  if (mi[ioff]& groupbit) 
			{
  	    //dtfm = dtf / mass[ti[ioff]];
				dtfm = dtf * massinv[ti[ioff]];
  	    vi[ioff][0] += dtfm * fi[ioff][0];
  	    vi[ioff][1] += dtfm * fi[ioff][1];
  	    vi[ioff][2] += dtfm * fi[ioff][2];
  	    xi[ioff][0] += dtv * vi[ioff][0];
  	    xi[ioff][1] += dtv * vi[ioff][1];
  	    xi[ioff][2] += dtv * vi[ioff][2];
  	  }
		}//for-i
		
		pe_put(&x[ist][0], &xi[0][0], sizeof(double)*isz*3);
		pe_put(&v[ist][0], &vi[0][0], sizeof(double)*isz*3);
		dma_syn();

	}//for-ist
}


void fix_nve_final_integrate_para(fix_nve_param_t *pm)
{
	dma_init();
	fix_nve_param_t l_pm;
	pe_get(pm, &l_pm, sizeof(fix_nve_param_t));
	dma_syn();

  double dtfm;

  double (*v)[3] 	= l_pm.v;
  double (*f)[3] 	= l_pm.f;
  double *mass		= l_pm.mass;
	double dtf			= l_pm.dtf;
  int *type				= l_pm.type;
  int *mask 			= l_pm.mask;
  int nlocal			= l_pm.nlocal;
	int groupbit		= l_pm.groupbit;
	int ntypes			= l_pm.ntypes;

	int i, ist, ied, isz, ioff;
	double vi[ISTEP][3], fi[ISTEP][3];
	int ti[ISTEP], mi[ISTEP];
	
	double massinv[ntypes+1];
	pe_get(l_pm.massinv, massinv, sizeof(double)*(ntypes+1));
	dma_syn();
  
  for (ist = _MYID*ISTEP; ist < nlocal; ist+=64*ISTEP)
	{
		ied = ist + ISTEP;
		if(ied > nlocal)
			ied = nlocal;
		isz = ied - ist;
		
		pe_get(mask+ist, mi, sizeof(int)*isz);
		pe_get(type+ist, ti, sizeof(int)*isz);
		pe_get(&v[ist][0], &vi[0][0], sizeof(double)*isz*3);
		pe_get(&f[ist][0], &fi[0][0], sizeof(double)*isz*3);
		dma_syn();
		for(i = ist; i < ied; i++)
		{
			ioff = i - ist;
  	  if (mi[ioff]& groupbit) 
			{
				dtfm = dtf * massinv[ti[ioff]];
        vi[ioff][0] += dtfm * fi[ioff][0];
        vi[ioff][1] += dtfm * fi[ioff][1];
        vi[ioff][2] += dtfm * fi[ioff][2];
  	  }
		}//for-i
		
		pe_put(&v[ist][0], &vi[0][0], sizeof(double)*isz*3);
		dma_syn();

	}//for-ist
}

#endif
