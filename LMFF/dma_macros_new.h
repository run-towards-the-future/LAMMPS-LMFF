//#ifdef __sw_256bit_simd__
//typedef int desc_t __attribute__ ((__mode__(__V1OI__)));
//#elif defined(__sw_512bit_simd__)
typedef int desc_t __attribute__ ((__mode__(__V1XI__)));
//#endif
#define dma_desc_init(desc, dir, mode, rpl_addr)   \
  asm volatile("vinsw %[RPL], $31, 2, %[DESC]\n\t" \
               "sll %[MODE], 56, %[DESC]\n\t"      \
               : [ DESC ] "=&r"(desc)              \
               : [ MODE ] "r"((mode << 4) | (dir)), [ RPL ] "r"(rpl_addr))
#define dma_set_mask(desc, mask)                              \
  asm volatile("vinsw %[DESC_IN], %[DESC_IN], 4, %[DESC]\n\t" \
               "vextw %[DESC], 2, %[DESC]\n\t"                \
               "zap   %[DESC], 0xf8, %[DESC]\n\t"             \
               "sll   %[MASK], 24, %[MASK]\n\t"               \
               "or    %[MASK], %[DESC], %[DESC]\n\t"          \
               "srl   %[MASK], 24, %[MASK]\n\t"               \
               "vinsw %[DESC], %[DESC], 2, %[DESC]\n\t"       \
               "vextw %[DESC], 4, %[DESC]\n\t"                \
               : [ DESC ] "=r"(desc)                          \
               : [ MASK ] "r"(mask), [ DESC_IN ] "0"(desc))
#define dma_set_stride(desc, stride, bsize)                  \
  asm volatile("vinsw %[STRIDE], %[DESC_IN], 3, %[DESC]\n\t" \
               "vinsf %[DESC], 2, %[DESC]\n\t"               \
               "zap %[DESC], 0x70, %[DESC]\n\t"              \
               "sll %[BSIZE], 32, %[BSIZE]\n\t"              \
               "or %[BSIZE], %[DESC], %[DESC]\n\t"           \
               "srl %[BSIZE], 32, %[BSIZE]\n\t"              \
               : [ DESC ] "=r"(desc),                        \
               : [ STRIDE ] "r"(stride), [ BSIZE ] "r"(bsize), [ DESC_IN ] "0"(desc))
#define dma_set_size(desc, size)                        \
  asm volatile("vinsw %[SIZE], %[DESC_IN], %[DESC]\n\t" \
               : [ DESC ] "=r"(desc),                   \
               : [ SIZE ] "r"(size), [ DESC_IN ] "0"(desc));
#define dma_rpl(desc, mem, ldm, rpl)             \
  asm volatile("dma %[DESC], %[MEM], %[LDM]\n\t" \
               : "+m"(rpl)                       \
               : [ DESC ] "r"(desc), [ MEM ] "r"(mem), [ LDM ] "r"(ldm))
#define dma_size_rpl(desc, mem, ldm, size, rpl)                                             \
  asm volatile("vinsw %[SIZE], %[DESC], 0, %[DESC]\n\t"                                     \
               "dma %[DESC], %[MEM], %[LDM]\n\t"                                            \
               : "+m"(rpl), "=m"(*(char(*)[size])(mem)), "=m"(*(char(*)[size])(ldm))        \
               : [ DESC ] "r"(desc), [ MEM ] "r"(mem), [ LDM ] "r"(ldm), [ SIZE ] "r"(size) \
               : "memory")
#ifdef SAFER
#define reply_t volatile int
#else
#define reply_t int
#endif
#define dma_init()                                                      \
  desc_t pe_get_desc_shadow, pe_put_desc_shadow;                        \
  reply_t reply_shadow = 0, count_shadow = 0;                           \
  {                                                                     \
    dma_desc_init(pe_get_desc_shadow, DMA_GET, PE_MODE, &reply_shadow); \
    dma_desc_init(pe_put_desc_shadow, DMA_PUT, PE_MODE, &reply_shadow); \
  }
#define pe_get(mem, ldm, size)                                        \
  {                                                                   \
    if (size) {                                                       \
      dma_size_rpl(pe_get_desc_shadow, mem, ldm, size, reply_shadow); \
      count_shadow++;                                                 \
    }                                                                 \
  }
#define pe_put(mem, ldm, size)                                        \
  {                                                                   \
    if (size) {                                                       \
      dma_size_rpl(pe_put_desc_shadow, mem, ldm, size, reply_shadow); \
      count_shadow++;                                                 \
    }                                                                 \
  }
#define pe_getn(mem, ldm, size) pe_get(mem, ldm, ((size) * sizeof(*(ldm))));
#define pe_putn(mem, ldm, size) pe_put(mem, ldm, ((size) * sizeof(*(ldm))));
#define dma_syn()                                                  \
  {                                                                \
    int tmp;                                                       \
    asm volatile("1:\n\t"                                          \
                 "ldw %[RPL], %[MRPL]\n\t"                         \
                 "subw %[RPL], %[CNT], %[RPL]\n\t"                 \
                 "bne %[RPL], 1b\n\t"                              \
                 "memb\n\t"                                        \
                 : [ RPL ] "=&r"(tmp), [ MRPL ] "+m"(reply_shadow) \
                 : [ CNT ] "r"(count_shadow)                       \
                 : "memory");                                      \
  }

#ifdef __cplusplus
#define DMA_CLASS                                \
  dest_t pe_get_desc_shadow, pe_put_desc_shadow; \
  volatile int reply_shadow = 0, count_shadow = 0;
#define dma_init_class()                                                \
  reply_shadow = 0, count_shadow = 0;                                   \
  {                                                                     \
    dma_desc_init(pe_get_desc_shadow, DMA_GET, PE_MODE, &reply_shadow); \
    dma_desc_init(pe_put_desc_shadow, DMA_PUT, PE_MODE, &reply_shadow); \
  }
#endif
