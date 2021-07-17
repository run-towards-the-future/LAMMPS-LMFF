#ifndef PAIR_LMFF_PARAMETERS_H_
#define PAIR_LMFF_PARAMETERS_H_
#ifdef __cplusplus
extern "C"{
#endif

/* ILP */
typedef struct Param_ilp_str 
{
	double z0,alpha,epsilon,C,delta,d,sR,reff,C6,S;
  double delta2inv,seff,lambda,rcut;
	double seffinv;
  int ielement,jelement;
}Param_ilp_str;

/* Tersoff  */
typedef struct Param 
{
	double lam1,lam2,lam3;
  double c,d,h;
  double gamma,powerm;
  double powern,beta;
  double biga,bigb,bigd,bigr;
  double cut,cutsq;
  double c1,c2,c3,c4;
  int ielement,jelement,kelement;
  int powermint;

	double dinv, bigdinv;
	double powerninvhalf;

  //double Z_i,Z_j;              // added for TersoffZBL
  //double ZBLcut,ZBLexpscale;
  //double c5,ca1,ca4;           // added for TersoffMOD
  //double powern_del;
  //double c0;                   // added for TersoffMODC
}Param;


#ifdef __cplusplus
}
#endif
#endif
