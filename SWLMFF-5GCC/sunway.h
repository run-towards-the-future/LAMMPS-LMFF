#ifndef SUNWAY_H_
#define SUNWAY_H_
//extern long *ptr_lock;
extern long *mylock;

typedef double rvec4[4];

#ifdef MPE
#include <athread.h>
//#include <crts.h>
//extern int flag_athread_init = 0;
extern int flag_athread_init;

#endif




#ifdef CPE
#include <slave.h>
#include <simd.h>
#include "dma_macros_sw5.h"

/*** DMA for SW9***/
//#ifdef SAFER
//	#include "dma_macros_sw9.h"
//#else
//	#define dma_init() volatile int reply = 0; int pe_cnt = 0;
//	#define pe_get(mem, ldm, size) {if(size){__builtin_sw_slave_athread_dma(PE_MODE, DMA_GET, (mem), (ldm), (size), (void*)&reply, 0, 0, 0); pe_cnt ++;}}
//	#define pe_put(mem, ldm, size) {if(size){__builtin_sw_slave_athread_dma(PE_MODE, DMA_PUT, (mem), (ldm), (size), (void*)&reply, 0, 0, 0); pe_cnt ++;}}
//	#define dma_syn() {while (reply != pe_cnt); asm volatile("memb");}
//	#define dma_rpl(desc, mem, ldm, reply) asm("dma %0, %1, %2\n\t" : : "r"(desc), "r"(mem), "r"(ldm), "r"(&reply) : "memory");
//#endif


#endif


#endif

