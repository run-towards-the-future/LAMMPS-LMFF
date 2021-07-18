#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_lmff_parameters.h"
#include "sleef_math_full.h"

#define always_inline __attribute__((always_inline)) inline

#define MY_PI2  (1.57079632679489661923)
#define MY_PI4  (0.78539816339744830962)
#define  cube(x)  (x*x*x)
#define  square(x) (x*x)

#define SLEEF
/************************************************************************************/
inline double ters_gijk(double costheta,Param *param) 
{
  const double ters_c = param->c * param->c;
  const double ters_d = param->d * param->d;
  const double ters_dinv = param->dinv * param->dinv;
  const double hcth = param->h - costheta;

  return param->gamma*(1.0 + ters_c*ters_dinv - ters_c / (ters_d + hcth*hcth));
}

inline double ters_gijk_d(double costheta, Param * param) 
{
  const double ters_c = param->c * param->c;
  const double ters_d = param->d * param->d;
  const double hcth = param->h - costheta;
  const double numerator = -2.0 * ters_c * hcth;
  const double denominator = 1.0/(ters_d + hcth*hcth);
  return param->gamma*numerator*denominator*denominator;
}

inline double vec3_dot(double x[3], double y[3]) 
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

inline void vec3_add(double x[3], double y[3],double * z) 
{
  z[0] = x[0]+y[0];  z[1] = x[1]+y[1];  z[2] = x[2]+y[2];
}

inline void vec3_scale(double k, double x[3], double y[3]) 
{
  y[0] = k*x[0];  y[1] = k*x[1];  y[2] = k*x[2];
}

inline void vec3_scaleadd(const double k, const double x[3],
        const double y[3], double * const z)  
{
  z[0] = k*x[0]+y[0];
  z[1] = k*x[1]+y[1];
  z[2] = k*x[2]+y[2];
}

void costheta_d(double *rij_hat, double rij, double rijinv,
                double *rik_hat, double rik, double rikinv,
                double *dri, double *drj, double *drk)
{
  double cos_theta = vec3_dot(rij_hat,rik_hat);

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(rijinv,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(rikinv,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}
void ters_FRep(Param *param, double r, double rinv, 
							double *fforce, int eflag, double *eng)
{
  double tmp_fc,tmp_fc_d,tmp_exp;

	if (r < param->bigr-param->bigd) {tmp_fc = 1.0; tmp_fc_d = 0.0;}
	else if (r > param->bigr+param->bigd) {tmp_fc = 0.0; tmp_fc_d = 0.0;}
	else 
	{
		tmp_fc = 0.5*(1.0 - xsin(MY_PI2*(r - param->bigr)*param->bigdinv )); 
		tmp_fc_d =  -(MY_PI4*param->bigdinv) * xcos(MY_PI2*(r - param->bigr)*param->bigdinv);
	}

  tmp_exp = xexp(-param->lam1 * r);

  *fforce = -param->biga * tmp_exp * (tmp_fc_d - tmp_fc*param->lam1) * rinv;
  if (eflag) *eng = tmp_fc * param->biga * tmp_exp;
}

void ters_Att(Param *param, double prefactor,
              double rij, double rijinv, double rik, double rikinv,
              double *delrij, double *delrik,
              double *dri, double *drj, double *drk)
{
  double rij_hat[3],rik_hat[3];

  vec3_scale(rijinv,delrij,rij_hat);

  vec3_scale(rikinv,delrik,rik_hat);

	double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

	if (rik < param->bigr-param->bigd) {fc = 1.0; dfc = 0.0;}
	else if (rik > param->bigr+param->bigd) {fc = 0.0; dfc = 0.0;}
	else 
	{
		fc = 0.5*(1.0 - xsin(MY_PI2*(rik - param->bigr)*param->bigdinv)); 
		dfc =  -(MY_PI4*param->bigdinv) * xcos(MY_PI2*(rik - param->bigr)*param->bigdinv);
	}
  
	if (param->powermint == 3) tmp = cube(param->lam3 * (rij-rik));
  else tmp = param->lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = xexp(tmp);

  if (param->powermint == 3)
    ex_delr_d = 3.0*cube(param->lam3) * square(rij-rik)*ex_delr;
  else ex_delr_d = param->lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(cos_theta,param);
  gijk_d = ters_gijk_d(cos_theta,param);
  costheta_d(rij_hat,rij,rijinv,rik_hat,rik,rikinv,dcosdri,dcosdrj,dcosdrk);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
	
}
double ters_fc(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;
  double ters_Dinv = param->bigdinv;

  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - xsin(MY_PI2*(r - ters_R)*ters_Dinv));
}

/* ---------------------------------------------------------------------- */

double ters_fc_d(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;
  double ters_Dinv = param->bigdinv;

  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(MY_PI4*ters_Dinv) * xcos(MY_PI2*(r - ters_R)*ters_Dinv);
}
void ters_force_zeta(Param *param, double r,double rinv, double zeta_ij,
                     double *fforce, double *prefactor,
                     int eflag, double *eng)
{
  double fa,fa_d,bij;
	double powerninvhalf = param->powerninvhalf;

	if (r > param->bigr + param->bigd) {fa = 0.0; fa_d = 0.0;}
	else
	{
		double fc, dfc, tmp_exp;
		if (r < param->bigr-param->bigd) {fc = 1.0; dfc = 0.0;}
		else if (r > param->bigr+param->bigd) {fc = 0.0; dfc = 0.0;}
		else 
		{
			fc = 0.5*(1.0 - xsin(MY_PI2*(r - param->bigr)*param->bigdinv)); 
			dfc =  -(MY_PI4*param->bigdinv) * xcos(MY_PI2*(r - param->bigr)*param->bigdinv);
		}
		tmp_exp = xexp(-param->lam2 * r);

		fa = -param->bigb * tmp_exp * fc; 
		fa_d = param->bigb * tmp_exp * (param->lam2 * fc - dfc); 
	}

	double dbij;
	double tmp = param->beta * zeta_ij;
  if (tmp > param->c1) 
	{
		bij = 1.0/sqrt(tmp);
		dbij = param->beta * -0.5*xpow(tmp,-1.5);
	}
  else if (tmp > param->c2)
	{
		double tmp_nn = xpow(tmp,-param->powern);
		bij = (1.0 - tmp_nn *powerninvhalf)/sqrt(tmp);
		dbij = param->beta * (-0.5*xpow(tmp,-1.5) * (1.0 - (1.0+powerninvhalf)*tmp_nn));

	}
  else if (tmp < param->c4) 
	{
		bij = 1.0; dbij = 0.0;
	}
  else if (tmp < param->c3)
	{
		double tmp_n1 = xpow(tmp,param->powern-1.0);
    bij = 1.0 - tmp_n1*tmp*powerninvhalf;
    dbij = -0.5*param->beta * tmp_n1;
	}
  else 
	{
		double tmp_n = xpow(tmp,param->powern);
		bij = xpow(1.0 + tmp_n, -powerninvhalf);
		dbij = -0.5 * xpow(1.0+tmp_n, -1.0-powerninvhalf)*tmp_n / zeta_ij;
	}


  *fforce = 0.5*bij*fa_d *rinv;
  *prefactor = -0.5*fa * dbij;
  if (eflag) *eng = 0.5*bij*fa;
}

void ters_force_zeta_fc(Param *param, double r,double rinv, double zeta_ij,
                     double *fforce, double *prefactor,
                     int eflag, double *eng, 
										 double fc, double dfc)
{
  double fa,fa_d,bij;
	double powerninvhalf = param->powerninvhalf;

	if (r > param->bigr + param->bigd) {fa = 0.0; fa_d = 0.0;}
	else
	{
		double tmp_exp;
		tmp_exp = xexp(-param->lam2 * r);

		fa = -param->bigb * tmp_exp * fc; 
		fa_d = param->bigb * tmp_exp * (param->lam2 * fc - dfc); 
	}

	double dbij;
	double tmp = param->beta * zeta_ij;
  if (tmp > param->c1) 
	{
		bij = 1.0/sqrt(tmp);
		dbij = param->beta * -0.5*xpow(tmp,-1.5);
	}
  else if (tmp > param->c2)
	{
		double tmp_nn = xpow(tmp,-param->powern);
		bij = (1.0 - tmp_nn *powerninvhalf)/sqrt(tmp);
		dbij = param->beta * (-0.5*xpow(tmp,-1.5) * (1.0 - (1.0+powerninvhalf)*tmp_nn));

	}
  else if (tmp < param->c4) 
	{
		bij = 1.0; dbij = 0.0;
	}
  else if (tmp < param->c3)
	{
		double tmp_n1 = xpow(tmp,param->powern-1.0);
    bij = 1.0 - tmp_n1*tmp*powerninvhalf;
    dbij = -0.5*param->beta * tmp_n1;
	}
  else 
	{
		double tmp_n = xpow(tmp,param->powern);
		bij = xpow(1.0 + tmp_n, -powerninvhalf);
		dbij = -0.5 * xpow(1.0+tmp_n, -1.0-powerninvhalf)*tmp_n / zeta_ij;
	}


  *fforce = 0.5*bij*fa_d *rinv;
  *prefactor = -0.5*fa * dbij;
  if (eflag) *eng = 0.5*bij*fa;
}



double ters_zeta(Param *param, double rij,double rijinv, 
				double rik, double rikinv,
        double *delrij, double *delrik)
{
  double costheta,arg,ex_delr;

  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) *rijinv*rikinv;

  if (param->powermint == 3) arg = cube(param->lam3 * (rij-rik));
  else arg = param->lam3 * (rij-rik);

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = xexp(arg);

  return ters_fc(rik,param) * ters_gijk(costheta,param) * ex_delr;
}

/************************************************************************************/
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
	
inline void calc_dnormal(double (*dnormal)[3][3], double (*dn1)[3][3], double *n1, 
												double nn, double nn2, double nninv, double nn2inv)
{
  double dnn[3][3];
  int m, id, ip;
  for (m = 0; m < 3; m++){
    for (id = 0; id < 3; id++){
      dnn[id][m] = (n1[0]*dn1[m][0][id] + n1[1]*dn1[m][1][id] + n1[2]*dn1[m][2][id])*nninv;
    }
  }
  for (m = 0; m < 3; m++){
    for (id = 0; id < 3; id++){
      for (ip = 0; ip < 3; ip++){
        dnormal[m][id][ip] = dn1[m][id][ip]*nninv - n1[id]*dnn[ip][m]*nn2inv;
      }
    }
  }
}

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


void calc_normal(int nilp, int jnum, double vet[][3], 
								double normal[], double dnormal[][3][3], double dnormdri[][3])
{
	int id,ip,m;
  double nn, nn2, nninv, nn2inv;
  double pv12[3],pv31[3],pv23[3],n1[3],dni[3],dnn[3][3],/*vet[3][3],*/dpvdri[3][3];
  double dn1[3][3][3],dpv12[3][3][3],dpv23[3][3][3],dpv31[3][3][3];
  //double normal[3], dnormal[3][3][3], dnormdri[3][3];

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
		nninv = 1.0 / nn;
		nn2inv = nninv * nninv;
    
		if (nn == 0) printf("The magnitude of the normal vector is zero\n");

    // the unit normal vector
    normal[0] = n1[0]*nninv;
    normal[1] = n1[1]*nninv;
    normal[2] = n1[2]*nninv;
    // derivatives of nn, dnn:3x1 vector
    dni[0] = (n1[0]*dpvdri[0][0] + n1[1]*dpvdri[1][0] + n1[2]*dpvdri[2][0])*nninv;
    dni[1] = (n1[0]*dpvdri[0][1] + n1[1]*dpvdri[1][1] + n1[2]*dpvdri[2][1])*nninv;
    dni[2] = (n1[0]*dpvdri[0][2] + n1[1]*dpvdri[1][2] + n1[2]*dpvdri[2][2])*nninv;
    for (id = 0; id < 3; id++)
		{
      for (ip = 0; ip < 3; ip++)
			{
        dnormdri[id][ip] = dpvdri[id][ip]*nninv - n1[id]*dni[ip]*nn2inv;
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
    calc_dnormal(dnormal, dn1, n1, nn, nn2, nninv, nn2inv);

  }
  else if(nilp == 3) 
	{
    cross_deriv(pv12, dpv12, vet, 0, 1, 2);
    cross_deriv(pv31, dpv31, vet, 2, 0, 1);
    cross_deriv(pv23, dpv23, vet, 1, 2, 0);

    n1[0] = (pv12[0] + pv31[0] + pv23[0]) * jnuminv;;
    n1[1] = (pv12[1] + pv31[1] + pv23[1]) * jnuminv;;
    n1[2] = (pv12[2] + pv31[2] + pv23[2]) * jnuminv;;
    // the magnitude of the normal vector
    nn2 = n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2];
    nn = sqrt(nn2);
		nninv = 1.0 / nn;
		nn2inv = nninv * nninv;

    if (nn == 0) printf("The magnitude of the normal vector is zero\n");
    // the unit normal vector
    normal[0] = n1[0]*nninv;
    normal[1] = n1[1]*nninv;
    normal[2] = n1[2]*nninv;

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
          dn1[m][id][ip] = (dpv12[m][id][ip] + dpv23[m][id][ip] + dpv31[m][id][ip])*jnuminv;
        }
      }
    }
    calc_dnormal(dnormal, dn1, n1, nn, nn2, nninv, nn2inv);
  } 
	else 
	{
    printf("There are too many neighbors for calculating normals\n");
  }
}
/************************************************************************************/
