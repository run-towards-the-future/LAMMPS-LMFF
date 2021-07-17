/* WTCache*/
#define WRT_C_H    11
#define WRT_C_S    4
#define WRT_C_LSZ  (1 << WRT_C_S)
#define WRT_C_LCNT (1 << (WRT_C_H - WRT_C_S))
#define WRT_C_MM   (WRT_C_LSZ - 1)
#define WRT_C_LM   (WRT_C_LCNT - 1)

/* RDCache */
#define READ_C_H    10 //10 //10
#define READ_C_S		3  //2  //3
#define READ_C_LSZ  (1 << READ_C_S)
#define READ_C_LCNT (1 << (READ_C_H - READ_C_S))
#define READ_C_MM   (READ_C_LSZ - 1)
#define READ_C_LM   (READ_C_LCNT - 1)

//#define READ_C_S_DBL		(READ_C_S+1) // two-way
//#define READ_C_LSZ_DBL  (1 << READ_C_S_DBL)
//#define READ_C_MM_DBL   (READ_C_LSZ_DBL - 1)

/* RDCache for marking atoms */
#define MREAD_C_H			10 //10
#define MREAD_C_S			3
#define MREAD_C_LSZ  (1 << MREAD_C_S)
#define MREAD_C_LCNT (1 << (MREAD_C_H - MREAD_C_S))
#define MREAD_C_MM   (MREAD_C_LSZ - 1)
#define MREAD_C_LM   (MREAD_C_LCNT - 1)


#ifdef CPE
#include "sunway.h"

/*** SWCACHE-Write ***/
#define SWCACHE_LOCK(lock) {                    \
    long t0, t1;                                 \
    asm volatile("faal %0, 0(%2)\n\t"           \
                 "1000:"                        \
                 "ldl %1, 8(%2)\n\t"            \
                 "subl %1, %0, %1\n\t"          \
                 "bne %1, 1000b\n\t"            \
                 : "=&r"(t0), "=&r"(t1)         \
                 : "r"(lock) : "memory");                  \
  }
#define SWCACHE_UNLOCK(lock) asm volatile("faal $31, 8(%0)\n\t" :: "r"(lock): "memory")



void update_cache_ters(int i, double *fi, 
									rvec4 wfcache[][WRT_C_LSZ], int *wfctag,
									rvec4 *frc, long *mylock)
{
  dma_init();
	int glb_line		= i >> WRT_C_S;
	int cache_line	= glb_line & WRT_C_LM;
	int cache_head	= i & ~WRT_C_MM;
	int cache_off		= i & WRT_C_MM;
	rvec4 bwf[WRT_C_LSZ];
	long *ticket;
	int j;
	if (wfctag[cache_line] == -1)
  {
    wfctag[cache_line] = cache_head;			/*load a new line */
		for (j = 0; j < WRT_C_LSZ; j ++)
		{
			wfcache[cache_line][j][0] = 0;
			wfcache[cache_line][j][1] = 0;
			wfcache[cache_line][j][2] = 0;
			wfcache[cache_line][j][3] = 0;
		}
  }
	else if (wfctag[cache_line] !=  cache_head)
  {
		//miss_wc++;
		/* write back the current line*/
		ticket = mylock + (wfctag[cache_line] >> WRT_C_S)*2;
		SWCACHE_LOCK(ticket);									/* LOCK */
		pe_get(frc+wfctag[cache_line], bwf, sizeof(double)*4*WRT_C_LSZ);
		dma_syn();
		for(j = 0; j < WRT_C_LSZ; j++)
		{
			bwf[j][0] += wfcache[cache_line][j][0];
			bwf[j][1] += wfcache[cache_line][j][1];
			bwf[j][2] += wfcache[cache_line][j][2];
			bwf[j][3] = 0;
		}
		pe_put(frc+wfctag[cache_line], bwf, sizeof(double)*4*WRT_C_LSZ);
		dma_syn();
		SWCACHE_UNLOCK(ticket);							/* UNLOCK */

    wfctag[cache_line] = cache_head;		/*load a new line */
		for (j = 0; j < WRT_C_LSZ; j ++)
		{
			wfcache[cache_line][j][0] = 0;
			wfcache[cache_line][j][1] = 0;
			wfcache[cache_line][j][2] = 0;
			wfcache[cache_line][j][3] = 0;
		}
  }
	//tot_wc++;
	wfcache[cache_line][cache_off][0] += fi[0];
	wfcache[cache_line][cache_off][1] += fi[1];
	wfcache[cache_line][cache_off][2] += fi[2];
	wfcache[cache_line][cache_off][3] = 0;
}

void flush_cache_ters(rvec4 *frc, rvec4 wfcache[][WRT_C_LSZ], int *wfctag, long *mylock)
{/* Flush whole write_cache */
	dma_init();
	int i, j;
	rvec4 bwf[WRT_C_LSZ];
	long *ticket;
	for(i = 0; i < WRT_C_LCNT; i++)
	{
		if(wfctag[i] != -1)
		{
			ticket = mylock + (wfctag[i] >> WRT_C_S) * 2;
			SWCACHE_LOCK(ticket);
			pe_get(frc+wfctag[i], bwf, sizeof(double)*4*WRT_C_LSZ);
			dma_syn();
			for(j = 0; j < WRT_C_LSZ; j++)
			{
				bwf[j][0] += wfcache[i][j][0];
				bwf[j][1] += wfcache[i][j][1];
				bwf[j][2] += wfcache[i][j][2];
				bwf[j][3] = 0;
			}
			pe_put(frc+wfctag[i], bwf, sizeof(double)*4*WRT_C_LSZ);
			dma_syn();
			SWCACHE_UNLOCK(ticket);
		}
	}
}

#endif
