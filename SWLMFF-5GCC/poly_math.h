extern inline double p_lnd(double in);
extern inline double p_expd(double in);
extern inline double p_powd(double base, double power);
extern inline double p_sind(double in);
extern inline double p_cosd(double in);
extern inline double p_sinnpi_pid(double in);
extern inline double p_cosnpi_pid(double in);

extern inline doublev4 simd_vfloord(doublev4 in);
extern inline doublev4 simd_vceild(doublev4 in);
extern inline doublev4 simd_vroundd(doublev4 in);
extern inline doublev4 simd_vltod(doublev4 in);
extern inline doublev4 simd_vlnd(doublev4 in);
extern inline doublev4 simd_vln1_2d(doublev4 in);
extern inline doublev4 simd_vexpd(doublev4 in);
extern inline doublev4 simd_vpowd(doublev4 base, doublev4 power);
extern inline doublev4 simd_vsind(doublev4 in);
extern inline doublev4 simd_vcosd(doublev4 in);
extern inline doublev4 simd_vsinnpi_pid(doublev4 in);
extern inline doublev4 simd_vcosnpi_pid(doublev4 in);
extern inline doublev4 simd_vinv_sqrtd(doublev4 in);

extern __thread_local long exp_ln_coef[24];
extern __thread_local long sin_cos_coef[24];
extern __thread_local long inv_sqrt_coef[4];

//////////sqrt////////////
#define inv_sqrt(x, r) {                        \
    double y = x;                               \
    double xhalf = 0.5 * y;                     \
    long i = *(long*)(&y);                      \
    i = 0x5fe6ec85e7de30daLL - (i >> 1);        \
    y = *(double*)(&i);                         \
    y = y * (1.5 - xhalf * y * y );             \
    y = y * (1.5 - xhalf * y * y );             \
    y = y * (1.5 - xhalf * y * y );             \
    r = y;                                      \
  }




#include "math_macro.h"
#define msimd_inv_sqrtd(out, in)                                        \
  {                                                                     \
    float t0, t1, t2, t3;                                               \
    asm(                                                                \
        asm_vinv_sqrtd("%0", "%1", "%3", "%4", "%5", "%6", "%2")        \
        : "=&r"(out)                                                    \
        : "r"(in), "r"(inv_sqrt_coef),                                  \
          "r"(t0), "r"(t1), "r"(t2), "r"(t3));                          \
  }                                                                     \

