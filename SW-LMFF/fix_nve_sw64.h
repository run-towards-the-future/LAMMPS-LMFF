#ifndef FIX_NVE_SW64_H_
#define FIX_NVE_SW64_H_
#ifdef __cplusplus
extern "C"{
#endif

typedef struct fix_nve_param_t
{
	double (*x)[3], (*v)[3], (*f)[3];
	double *rmass, *mass, *massinv, dtv, dtf;
	int *type, *mask;
	int nlocal, groupbit, ntypes;

}fix_nve_param_t;

void fix_nve_initial_integrate_serial(fix_nve_param_t *pm);
void fix_nve_final_integrate_serial(fix_nve_param_t *pm);


#ifdef __cplusplus
}
#endif
#endif
