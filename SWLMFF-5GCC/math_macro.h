#define asm_vroundd(out, in, mode, t0, t1, t2, t3)      \
  "vextf   " in ", 0, " t0 "\n\t"                       \
  "vextf   " in ", 1, " t1 "\n\t"                       \
  "vextf   " in ", 2, " t2 "\n\t"                       \
  "vextf   " in ", 3, " t3 "\n\t"                       \
                                                        \
  "fcvtdlr " t0 ", " mode ", " t0 "\n\t"                \
  "fcvtdlr " t1 ", " mode ", " t1 "\n\t"                \
  "fcvtdlr " t2 ", " mode ", " t2 "\n\t"                \
  "fcvtdlr " t3 ", " mode ", " t3 "\n\t"                \
                                                        \
  "fcvtld  " t0 ", " t0 "\n\t"                          \
  "fcvtld  " t1 ", " t1 "\n\t"                          \
  "fcvtld  " t2 ", " t2 "\n\t"                          \
  "fcvtld  " t3 ", " t3 "\n\t"                          \
                                                        \
  "vinsf   " t0 ", " out ", 0, " out "\n\t"             \
  "vinsf   " t1 ", " out ", 1, " out "\n\t"             \
  "vinsf   " t2 ", " out ", 2, " out "\n\t"             \
  "vinsf   " t3 ", " out ", 3, " out "\n\t"             \

#define asm_vmodinvd(out, divee, diver, diverinv, t0, t1, t2, t3)       \
  "vmuld " divee ", " diverinv "," out "\n\t"                           \
  asm_vroundd(out, out, "5", t0, t1, t2, t3)                            \
  "vnmad " out ", " diver ", " divee ", " out "\n\t"                    \
    
#define msimd_modinvd(out, divee, diver, diverinv)                      \
  {                                                                     \
    double t0, t1, t2, t3;                                              \
    asm (                                                               \
         asm_vmodinvd("%0", "%5", "%6", "%7", "%1", "%2", "%3", "%4")   \
         : "=r"(out)                                                    \
         : "r"(t0), "r"(t1), "r"(t2), "r"(t3),                          \
           "r"(divee), "r"(diver), "r"(diverinv));                      \
  }                                                                     \

#define asm_vltod(out, in, t0, t1, t2, t3)      \
  "vextf   " in ", 1, " t1 "\n\t"               \
  "vextf   " in ", 2, " t2 "\n\t"               \
  "vextf   " in ", 3, " t3 "\n\t"               \
                                                \
  "fcvtld  " in ", " t0 "\n\t"                  \
  "fcvtld  " t1 ", " t1 "\n\t"                  \
  "fcvtld  " t2 ", " t2 "\n\t"                  \
  "fcvtld  " t3 ", " t3 "\n\t"                  \
                                                \
  "vinsf   " t1 ", " t0 ", 1, " out "\n\t"      \
  "vinsf   " t2 ", " out ", 2, " out "\n\t"     \
  "vinsf   " t3 ", " out ", 3, " out "\n\t"     \

#define msimd_ltod(out, in)                             \
  {                                                     \
    double t0, t1, t2, t3;                              \
    asm(                                                \
        asm_vltod("%0", "%1", "%2", "%3", "%4", "%5")   \
        : "=r"(out)                                     \
        : "r"(in), "r"(t0), "r"(t1), "r"(t2), "r"(t3)); \
  }                                                     \

#define asm_vptol(out, in, tmp)                         \
  "vinsw  " "$31,  " in  ", 0, " out       "\n\t"       \
                                                        \
  "ldi    " tmp ", "  "769($31)"           "\n\t"       \
  "ldih   " tmp ", " "1797(" tmp ")"       "\n\t"       \
                                                        \
  "vshfw  " out ", " out ", " tmp ", " out "\n\t"       \
  "vsrlw  " out ", 20, " out               "\n\t"       \
  "ldi    " tmp ", "  "1023($31)"          "\n\t"       \
  "vshff  " tmp ", " tmp ", $31" ", " tmp  "\n\t"       \
  "vsubl  " out ", " tmp ", " out          "\n\t"       \

#define msimd_ptol(out, in)                     \
  {                                             \
    int tmp;                                    \
    asm(                                        \
        asm_vptol("%0", "%1", "%2")            \
        : "=r"(out)                             \
        : "r"(in), "r"(tmp));                   \
  }                                             \

#define asm_vln1_2d(out, in, y2, y, t0, t1, t2, t3, tab)        \
  "ldde  " t0 ", 184(" tab ")\n\t"                              \
  "vaddd " in ", " t0 ", " t1 "\n\t"                            \
  "vsubd " in ", " t0 ", " t2 "\n\t"                            \
  "vdivd " t2 ", " t1 ", " y  "\n\t"                            \
  "vmuld " y  ", " y  ", " y2 "\n\t"                            \
  "vmuld  " y2 ", " y2 ", " out "\n\t"                          \
  "ldde  " t0 ", 160(" tab ")\n\t"                              \
  "ldde  " t1 ", 152(" tab ")\n\t"                              \
  "ldde  " t2 ", 144(" tab ")\n\t"                              \
  "ldde  " t3 ", 136(" tab ")\n\t"                              \
  "vmad  " t0 ", " out ", " t2 ", " t2 "\n\t"                   \
  "vmad  " t1 ", " out ", " t3 ", " t3 "\n\t"                   \
  "ldde  " t0 ", 128(" tab ")\n\t"                              \
  "ldde  " t1 ", 120(" tab ")\n\t"                              \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
  "ldde  " t0 ", 112(" tab ")\n\t"                              \
  "ldde  " t1 ", 104(" tab ")\n\t"                              \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
  "ldde  " t0 ", 96(" tab ")\n\t"                               \
  "ldde  " t1 ", 88(" tab ")\n\t"                               \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
  "vmad  " t2 ", " y2 ", " t3 ", " out "\n\t"                   \
  "vmuld " out ", " y ", " out "\n\t"                           \

#define asm_ln1_2d(out, in, y2, y, t0, t1, t2, t3, tab)         \
  "fldd  " t0 ", 184(" tab ")\n\t"                              \
  "faddd " in ", " t0 ", " t1 "\n\t"                            \
  "fsubd " in ", " t0 ", " t2 "\n\t"                            \
  "fdivd " t2 ", " t1 ", " y  "\n\t"                            \
  "fmuld " y  ", " y  ", " y2 "\n\t"                            \
  "fmuld  " y2 ", " y2 ", " out "\n\t"                          \
  "fldd  " t0 ", 160(" tab ")\n\t"                              \
  "fldd  " t1 ", 152(" tab ")\n\t"                              \
  "fldd  " t2 ", 144(" tab ")\n\t"                              \
  "fldd  " t3 ", 136(" tab ")\n\t"                              \
  "fmad  " t0 ", " out ", " t2 ", " t2 "\n\t"                   \
  "fmad  " t1 ", " out ", " t3 ", " t3 "\n\t"                   \
  "fldd  " t0 ", 128(" tab ")\n\t"                              \
  "fldd  " t1 ", 120(" tab ")\n\t"                              \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
  "fldd  " t0 ", 112(" tab ")\n\t"                              \
  "fldd  " t1 ", 104(" tab ")\n\t"                              \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
  "fldd  " t0 ", 96(" tab ")\n\t"                               \
  "fldd  " t1 ", 88(" tab ")\n\t"                               \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
  "fmad  " t2 ", " y2 ", " t3 ", " out "\n\t"                   \
  "fmuld " out ", " y ", " out "\n\t"                           \

  //  asm_vptol(y, in, t1)                              \

#define asm_lnd(out, in, y, y2, t0, t1, t2, t3, tab)    \
  "fldd   " t0 ", 184(" tab ")\n\t"                     \
  "fcpyse " t0 ", " in ", " out "\n\t"                  \
  asm_ln1_2d(out, out, y, y2, t0, t1, t2, t3, tab)      \
  "srl    " in ", 52, " y "\n\t"                        \
  "ldi    " t1 ", "  "1023($31)"          "\n\t"        \
  "subl   " y  ", " t1 ", " y "\n\t"                    \
  "fcvtld " y ", " y2 "\n\t"                            \
  "fldd   " t3 ", 168(" tab ")\n\t"                     \
  "fmad " y2 ", " t3 ", " out ", " out "\n\t"           \

#define asm_vlnd(out, in, y, y2, t0, t1, t2, t3, tab)   \
  "ldde   " t0 ", 184(" tab ")\n\t"                     \
  "vcpyse " t0 ", " in ", " out "\n\t"                  \
  asm_vln1_2d(out, out, y, y2, t0, t1, t2, t3, tab)     \
  asm_vptol(y, in, t1)                                  \
  asm_vltod(y2, y, t0, t1, t2, t3)                      \
  "ldde   " t3 ", 168(" tab ")\n\t"                     \
  "vmad " y2 ", " t3 ", " out ", " out "\n\t"           \

#define asm_vdivmodinvd(outres, outmod, divee, diver, diverinv, t0, t1, t2, t3) \
  "vmuld " divee ", " diverinv "," outres "\n\t"                        \
  asm_vroundd(outres, outres, "5", t0, t1, t2, t3)                      \
  "vnmad " outres ", " diver ", " divee ", " outmod "\n\t"              \

/*
  asm_vlddi_ln2inv(t0)                                  \
  "vmuld " t0 ", " in ", " eexp  "\n\t"                 \
*/

#define asm_expd(out, in, eexp, mod, t0, t1, t2, t3, tab)       \
  "fldd  " t0 ", 176(" tab ") \n\t"                             \
  "fmuld " t0 ", " in ", " eexp  "\n\t"                         \
                                                                \
  "fcvtdlr " eexp ", 7, " eexp "\n\t"                           \
  "fldd  " t0 ", 168(" tab ") \n\t"                             \
  "fldd  " t1 ", 184(" tab ") \n\t"                             \
  "fcvtld  " eexp ", " mod "\n\t"                               \
                                                                \
  "fnmad " mod ", " t0 ", " in "," mod "\n\t"                   \
  "sll   " eexp ", 52, " eexp "\n\t"                            \
  "addl  " eexp ", " t1 ", " eexp "\n\t"                        \
  "fcpys      $31, " mod ", " mod "\n\t"                        \
                                                                \
  "fmuld " mod ", " mod ", " out "\n\t"                         \
  "fldd  " t0 ", 80(" tab ") \n\t"                              \
  "fldd  " t1 ", 72(" tab ") \n\t"                              \
  "fldd  " t2 ", 64(" tab ") \n\t"                              \
  "fldd  " t3 ", 56(" tab ") \n\t"                              \
  "fmad  " out ", " t0 ", " t2 ", " t2 "\n\t"                   \
  "fmad  " out ", " t1 ", " t3 ", " t3 "\n\t"                   \
  "fldd  " t0 ", 48(" tab ") \n\t"                              \
  "fldd  " t1 ", 40(" tab ") \n\t"                              \
  "fmad  " out ", " t2 ", " t0 ", " t2 "\n\t"                   \
  "fmad  " out ", " t3 ", " t1 ", " t3 "\n\t"                   \
  "fldd  " t0 ", 32(" tab ") \n\t"                              \
  "fldd  " t1 ", 24(" tab ") \n\t"                              \
  "fmad  " out ", " t2 ", " t0 ", " t2 "\n\t"                   \
  "fmad  " out ", " t3 ", " t1 ", " t3 "\n\t"                   \
  "fldd  " t0 ", 16(" tab ") \n\t"                              \
  "fldd  " t1  ", 8(" tab ") \n\t"                              \
  "fmad  " out ", " t2 ", " t0 ", " t2 "\n\t"                   \
  "fmad  " out ", " t3 ", " t1 ", " t3 "\n\t"                   \
  "fldd  " t0  ", 0(" tab ") \n\t"                              \
  "fmuld " mod ", " t3 ", " t3 "\n\t"                           \
  "fmad  " out ", " t2 ", " t0 ", " t2 "\n\t"                   \
  "faddd " t3 ", " t2 ", " out "\n\t"                           \
  "fmuld " eexp ", " out ", " out "\n\t"                        \

#define asm_vexpd(out, in, eexp, mod, t0, t1, t2, t3, tab)      \
  "ldde  " t0 ", 176(" tab ") \n\t"                             \
  "vmuld " t0 ", " in ", " eexp  "\n\t"                         \
  "vextf   " eexp ", 1, " t1 "\n\t"                             \
  "vextf   " eexp ", 2, " t2 "\n\t"                             \
  "vextf   " eexp ", 3, " t3 "\n\t"                             \
                                                                \
  "fcvtdlr " eexp ", 7, " t0 "\n\t"                             \
  "fcvtdlr " t1 ", 7, " t1 "\n\t"                               \
  "fcvtdlr " t2 ", 7, " t2 "\n\t"                               \
  "fcvtdlr " t3 ", 7, " t3 "\n\t"                               \
                                                                \
  "vinsw   " t0 ", $31, 1, " eexp "\n\t"                        \
  "vinsw   " t1 ", " eexp ", 3, " eexp "\n\t"                   \
  "vinsw   " t2 ", " eexp ", 5, " eexp "\n\t"                   \
  "vinsw   " t3 ", " eexp ", 7, " eexp "\n\t"                   \
                                                                \
  "fcvtld  " t0 ", " t0 "\n\t"                                  \
  "fcvtld  " t1 ", " t1 "\n\t"                                  \
  "fcvtld  " t2 ", " t2 "\n\t"                                  \
  "fcvtld  " t3 ", " t3 "\n\t"                                  \
                                                                \
  "vinsf   " t1 ", " t0  ", 1, " mod "\n\t"                     \
  "ldde  " t0 ", 168(" tab ") \n\t"                             \
  "ldde  " t1 ", 184(" tab ") \n\t"                             \
  "vinsf   " t2 ", " mod ", 2, " mod "\n\t"                     \
  "vinsf   " t3 ", " mod ", 3, " mod "\n\t"                     \
                                                                \
  "vnmad " mod ", " t0 ", " in "," mod "\n\t"                   \
  "vsllw    " eexp ", 20, " eexp "\n\t"                         \
  "vaddl " eexp ", " t1 ", " eexp "\n\t"                        \
  "vcpys      $31, " mod ", " mod "\n\t"                        \
                                                                \
  "vmuld " mod ", " mod ", " out "\n\t"                         \
  "ldde  " t0 ", 80(" tab ") \n\t"                              \
  "ldde  " t1 ", 72(" tab ") \n\t"                              \
  "ldde  " t2 ", 64(" tab ") \n\t"                              \
  "ldde  " t3 ", 56(" tab ") \n\t"                              \
  "vmad  " out ", " t0 ", " t2 ", " t2 "\n\t"                   \
  "vmad  " out ", " t1 ", " t3 ", " t3 "\n\t"                   \
  "ldde  " t0 ", 48(" tab ") \n\t"                              \
  "ldde  " t1 ", 40(" tab ") \n\t"                              \
  "vmad  " out ", " t2 ", " t0 ", " t2 "\n\t"                   \
  "vmad  " out ", " t3 ", " t1 ", " t3 "\n\t"                   \
  "ldde  " t0 ", 32(" tab ") \n\t"                              \
  "ldde  " t1 ", 24(" tab ") \n\t"                              \
  "vmad  " out ", " t2 ", " t0 ", " t2 "\n\t"                   \
  "vmad  " out ", " t3 ", " t1 ", " t3 "\n\t"                   \
  "ldde  " t0 ", 16(" tab ") \n\t"                              \
  "ldde  " t1  ", 8(" tab ") \n\t"                              \
  "vmad  " out ", " t2 ", " t0 ", " t2 "\n\t"                   \
  "vmad  " out ", " t3 ", " t1 ", " t3 "\n\t"                   \
  "ldde  " t0  ", 0(" tab ") \n\t"                              \
  "vmuld " mod ", " t3 ", " t3 "\n\t"                           \
  "vmad  " out ", " t2 ", " t0 ", " t2 "\n\t"                   \
  "vaddd " t3 ", " t2 ", " out "\n\t"                           \
  "vmuld " eexp ", " out ", " out "\n\t"                        \


#define msimd_expd(out, in)                                             \
  {                                                                     \
    doublev4 eexp, mod, t0, t1, t2, t3;                                 \
    asm(                                                                \
        asm_vexpd("%0", "%1", "%2", "%3", "%4", "%5", "%6", "%7", "%8") \
        : "=r"(out)                                                     \
        :"r"(in), "r"(eexp), "r"(mod),                                  \
         "r"(t0), "r"(t1), "r"(t2), "r"(t3), "r"(exp_ln_coef));         \
  }                                                                     \

#define asm_vsind(out, in, x, x2, t0, t1, t2, t3, tab)  \
  "ldde  " t0 ", 168(" tab ")\n\t"                      \
  "vmuld " in ", " t0 ", " x "\n\t"                     \
                                                        \
  "vextf " x ", 1, " t1 "\n\t"                          \
  "vextf " x ", 2, " t2 "\n\t"                          \
  "vextf " x ", 3, " t3 "\n\t"                          \
                                                        \
  "fcvtdlr "  x ", 5, " t0 "\n\t"                       \
  "fcvtdlr " t1 ", 5, " t1 "\n\t"                       \
  "fcvtdlr " t2 ", 5, " t2 "\n\t"                       \
  "fcvtdlr " t3 ", 5, " t3 "\n\t"                       \
                                                        \
  "fcvtld " t0 ", " t0 "\n\t"                           \
  "fcvtld " t1 ", " t1 "\n\t"                           \
  "fcvtld " t2 ", " t2 "\n\t"                           \
  "fcvtld " t3 ", " t3 "\n\t"                           \
                                                        \
  "ldde  "  x ", 160(" tab ")\n\t"                      \
  "vinsf " t1 ", " t0 ", 1, " x2 "\n\t"                 \
  "vinsf " t2 ", " x2 ", 2, " x2 "\n\t"                 \
  "vinsf " t3 ", " x2 ", 3, " x2 "\n\t"                 \
                                                        \
  "vaddd " x ", $31, " t3 "\n\t"                        \
  "vaddd " x ", " x ", " x "\n\t"                       \
                                                        \
  "vnmad " x ", " x2 ", " in ", " x "\n\t"              \
  "vsubd " x ", " t3 ", " x "\n\t"                      \
  "vmuld " x ", " x ", " x2 "\n\t"                      \
                                                        \
  "ldde  " t0 ", 72(" tab ")\n\t"                       \
  "ldde  " t1 ", 64(" tab ")\n\t"                       \
  "vmuld " x2 ", " x2 ", " out "\n\t"                   \
  "ldde  " t2 ", 56(" tab ")\n\t"                       \
  "ldde  " t3 ", 48(" tab ")\n\t"                       \
  "vmad  " t0 ", " out ", " t2 ", " t2 "\n\t"           \
  "vmad  " t1 ", " out ", " t3 ", " t3 "\n\t"           \
                                                        \
  "ldde  " t0 ", 40(" tab ")\n\t"                       \
  "ldde  " t1 ", 32(" tab ")\n\t"                       \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"           \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"           \
                                                        \
  "ldde  " t0 ", 24(" tab ")\n\t"                       \
  "ldde  " t1 ", 16(" tab ")\n\t"                       \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"           \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"           \
                                                        \
  "ldde  " t0 ", 8(" tab ")\n\t"                        \
  "ldde  " t1 ", 0(" tab ")\n\t"                        \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"           \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"           \
                                                        \
  "vmad  " t2 ", " x2 ", " t3 ", " out "\n\t"           \
  "vmuld " out ", " x ", " out "\n\t"                   \

#define asm_vcosd(out, in, x, x2, t0, t1, t2, t3, tab)  \
  "ldde  " t0 ", 168(" tab ")\n\t"                      \
  "vmuld " in ", " t0 ", " x "\n\t"                     \
                                                        \
  "vextf " x ", 1, " t1 "\n\t"                          \
  "vextf " x ", 2, " t2 "\n\t"                          \
  "vextf " x ", 3, " t3 "\n\t"                          \
                                                        \
  "fcvtdlr "  x ", 5, " t0 "\n\t"                       \
  "fcvtdlr " t1 ", 5, " t1 "\n\t"                       \
  "fcvtdlr " t2 ", 5, " t2 "\n\t"                       \
  "fcvtdlr " t3 ", 5, " t3 "\n\t"                       \
                                                        \
  "fcvtld " t0 ", " t0 "\n\t"                           \
  "fcvtld " t1 ", " t1 "\n\t"                           \
  "fcvtld " t2 ", " t2 "\n\t"                           \
  "fcvtld " t3 ", " t3 "\n\t"                           \
                                                        \
  "ldde  "  x ", 160(" tab ")\n\t"                      \
  "vinsf " t1 ", " t0 ", 1, " x2 "\n\t"                 \
  "vinsf " t2 ", " x2 ", 2, " x2 "\n\t"                 \
  "vinsf " t3 ", " x2 ", 3, " x2 "\n\t"                 \
                                                        \
  "vaddd " x ", $31, " t3 "\n\t"                        \
  "vaddd " x ", " x ", " x "\n\t"                       \
                                                        \
  "vnmad " x ", " x2 ", " in ", " x "\n\t"              \
  "vsubd " x ", " t3 ", " x "\n\t"                      \
  "vmuld " x ", " x ", " x2 "\n\t"                      \
                                                        \
  "ldde  " t0 ", 152(" tab ")\n\t"                      \
  "ldde  " t1 ", 144(" tab ")\n\t"                      \
  "vmuld " x2 ", " x2 ", " out "\n\t"                   \
  "ldde  " t2 ", 136(" tab ")\n\t"                      \
  "ldde  " t3 ", 128(" tab ")\n\t"                      \
  "vmad  " t0 ", " out ", " t2 ", " t2 "\n\t"           \
  "vmad  " t1 ", " out ", " t3 ", " t3 "\n\t"           \
                                                        \
  "ldde  " t0 ", 120(" tab ")\n\t"                      \
  "ldde  " t1 ", 112(" tab ")\n\t"                      \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"           \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"           \
                                                        \
  "ldde  " t0 ", 104(" tab ")\n\t"                      \
  "ldde  " t1 ", 96(" tab ")\n\t"                       \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"           \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"           \
                                                        \
  "ldde  " t0 ", 88(" tab ")\n\t"                       \
  "ldde  " t1 ", 80(" tab ")\n\t"                       \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"           \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"           \
                                                        \
  "vmad  " t2 ", " x2 ", " t3 ", " out "\n\t"           \

#define asm_vsinnpi_pid(out, in, x, x2, t0, t1, t2, t3, tab)    \
  "vmuld " in ", " in ", " x2 "\n\t"                            \
                                                                \
  "ldde  " t0 ", 72(" tab ")\n\t"                               \
  "ldde  " t1 ", 64(" tab ")\n\t"                               \
  "vmuld " x2 ", " x2 ", " out "\n\t"                           \
  "ldde  " t2 ", 56(" tab ")\n\t"                               \
  "ldde  " t3 ", 48(" tab ")\n\t"                               \
  "vmad  " t0 ", " out ", " t2 ", " t2 "\n\t"                   \
  "vmad  " t1 ", " out ", " t3 ", " t3 "\n\t"                   \
                                                                \
  "ldde  " t0 ", 40(" tab ")\n\t"                               \
  "ldde  " t1 ", 32(" tab ")\n\t"                               \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
                                                                \
  "ldde  " t0 ", 24(" tab ")\n\t"                               \
  "ldde  " t1 ", 16(" tab ")\n\t"                               \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
                                                                \
  "ldde  " t0 ", 8(" tab ")\n\t"                                \
  "ldde  " t1 ", 0(" tab ")\n\t"                                \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
                                                                \
  "vmad  " t2 ", " x2 ", " t3 ", " out "\n\t"                   \
  "vnmad " out ", " in ", $31," out "\n\t"                      \

#define asm_vcosnpi_pid(out, in, x2, t0, t1, t2, t3, tab)       \
  "vmuld " in ", " in ", " x2 "\n\t"                            \
                                                                \
  "ldde  " t0 ", 152(" tab ")\n\t"                              \
  "ldde  " t1 ", 144(" tab ")\n\t"                              \
  "vmuld " x2 ", " x2 ", " out "\n\t"                           \
  "ldde  " t2 ", 136(" tab ")\n\t"                              \
  "ldde  " t3 ", 128(" tab ")\n\t"                              \
  "vmad  " t0 ", " out ", " t2 ", " t2 "\n\t"                   \
  "vmad  " t1 ", " out ", " t3 ", " t3 "\n\t"                   \
                                                                \
  "ldde  " t0 ", 120(" tab ")\n\t"                              \
  "ldde  " t1 ", 112(" tab ")\n\t"                              \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
                                                                \
  "ldde  " t0 ", 104(" tab ")\n\t"                              \
  "ldde  " t1 ", 96(" tab ")\n\t"                               \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
                                                                \
  "ldde  " t0 ", 88(" tab ")\n\t"                               \
  "ldde  " t1 ", 80(" tab ")\n\t"                               \
  "vmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "vmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
                                                                \
  "vnmsd  " t2 ", " x2 ", " t3 ", " out "\n\t"                  \

#define asm_sind(out, in, x, x2, t0, t1, t2, t3, tab)  \
  "fldd  " t0 ", 168(" tab ")\n\t"                      \
  "fmuld " in ", " t0 ", " x "\n\t"                     \
                                                        \
  "fcvtdlr "  x ", 5, " t0 "\n\t"                       \
  "fcvtld " t0 ", " x2 "\n\t"                           \
                                                        \
  "fldd  "  x ", 160(" tab ")\n\t"                      \
                                                        \
  "faddd " x ", $31, " t3 "\n\t"                        \
  "faddd " x ", " x ", " x "\n\t"                       \
                                                        \
  "fnmad " x ", " x2 ", " in ", " x "\n\t"              \
  "fsubd " x ", " t3 ", " x "\n\t"                      \
  "fmuld " x ", " x ", " x2 "\n\t"                      \
                                                        \
  "fldd  " t0 ", 72(" tab ")\n\t"                       \
  "fldd  " t1 ", 64(" tab ")\n\t"                       \
  "fmuld " x2 ", " x2 ", " out "\n\t"                   \
  "fldd  " t2 ", 56(" tab ")\n\t"                       \
  "fldd  " t3 ", 48(" tab ")\n\t"                       \
  "fmad  " t0 ", " out ", " t2 ", " t2 "\n\t"           \
  "fmad  " t1 ", " out ", " t3 ", " t3 "\n\t"           \
                                                        \
  "fldd  " t0 ", 40(" tab ")\n\t"                       \
  "fldd  " t1 ", 32(" tab ")\n\t"                       \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"           \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"           \
                                                        \
  "fldd  " t0 ", 24(" tab ")\n\t"                       \
  "fldd  " t1 ", 16(" tab ")\n\t"                       \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"           \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"           \
                                                        \
  "fldd  " t0 ", 8(" tab ")\n\t"                        \
  "fldd  " t1 ", 0(" tab ")\n\t"                        \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"           \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"           \
                                                        \
  "fmad  " t2 ", " x2 ", " t3 ", " out "\n\t"           \
  "fmuld " out ", " x ", " out "\n\t"                   \

#define asm_cosd(out, in, x, x2, t0, t1, t2, t3, tab)  \
  "fldd  " t0 ", 168(" tab ")\n\t"                      \
  "fmuld " in ", " t0 ", " x "\n\t"                     \
                                                        \
  "fcvtdlr "  x ", 5, " t0 "\n\t"                       \
  "fcvtld " t0 ", " x2 "\n\t"                           \
                                                        \
  "fldd  "  x ", 160(" tab ")\n\t"                      \
                                                        \
  "faddd " x ", $31, " t3 "\n\t"                        \
  "faddd " x ", " x ", " x "\n\t"                       \
                                                        \
  "fnmad " x ", " x2 ", " in ", " x "\n\t"              \
  "fsubd " x ", " t3 ", " x "\n\t"                      \
  "fmuld " x ", " x ", " x2 "\n\t"                      \
                                                        \
  "fldd  " t0 ", 152(" tab ")\n\t"                      \
  "fldd  " t1 ", 144(" tab ")\n\t"                      \
  "fmuld " x2 ", " x2 ", " out "\n\t"                   \
  "fldd  " t2 ", 136(" tab ")\n\t"                      \
  "fldd  " t3 ", 128(" tab ")\n\t"                      \
  "fmad  " t0 ", " out ", " t2 ", " t2 "\n\t"           \
  "fmad  " t1 ", " out ", " t3 ", " t3 "\n\t"           \
                                                        \
  "fldd  " t0 ", 120(" tab ")\n\t"                      \
  "fldd  " t1 ", 112(" tab ")\n\t"                      \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"           \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"           \
                                                        \
  "fldd  " t0 ", 104(" tab ")\n\t"                      \
  "fldd  " t1 ", 96(" tab ")\n\t"                       \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"           \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"           \
                                                        \
  "fldd  " t0 ", 88(" tab ")\n\t"                       \
  "fldd  " t1 ", 80(" tab ")\n\t"                       \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"           \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"           \
                                                        \
  "fmad  " t2 ", " x2 ", " t3 ", " out "\n\t"           \

#define asm_sinnpi_pid(out, in, x, x2, t0, t1, t2, t3, tab)    \
  "fmuld " in ", " in ", " x2 "\n\t"                            \
                                                                \
  "fldd  " t0 ", 72(" tab ")\n\t"                               \
  "fldd  " t1 ", 64(" tab ")\n\t"                               \
  "fmuld " x2 ", " x2 ", " out "\n\t"                           \
  "fldd  " t2 ", 56(" tab ")\n\t"                               \
  "fldd  " t3 ", 48(" tab ")\n\t"                               \
  "fmad  " t0 ", " out ", " t2 ", " t2 "\n\t"                   \
  "fmad  " t1 ", " out ", " t3 ", " t3 "\n\t"                   \
                                                                \
  "fldd  " t0 ", 40(" tab ")\n\t"                               \
  "fldd  " t1 ", 32(" tab ")\n\t"                               \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
                                                                \
  "fldd  " t0 ", 24(" tab ")\n\t"                               \
  "fldd  " t1 ", 16(" tab ")\n\t"                               \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
                                                                \
  "fldd  " t0 ", 8(" tab ")\n\t"                                \
  "fldd  " t1 ", 0(" tab ")\n\t"                                \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
                                                                \
  "fmad  " t2 ", " x2 ", " t3 ", " out "\n\t"                   \
  "fnmad " out ", " in ", $31," out "\n\t"                      \

#define asm_cosnpi_pid(out, in, x2, t0, t1, t2, t3, tab)       \
  "fmuld " in ", " in ", " x2 "\n\t"                            \
                                                                \
  "fldd  " t0 ", 152(" tab ")\n\t"                              \
  "fldd  " t1 ", 144(" tab ")\n\t"                              \
  "fmuld " x2 ", " x2 ", " out "\n\t"                           \
  "fldd  " t2 ", 136(" tab ")\n\t"                              \
  "fldd  " t3 ", 128(" tab ")\n\t"                              \
  "fmad  " t0 ", " out ", " t2 ", " t2 "\n\t"                   \
  "fmad  " t1 ", " out ", " t3 ", " t3 "\n\t"                   \
                                                                \
  "fldd  " t0 ", 120(" tab ")\n\t"                              \
  "fldd  " t1 ", 112(" tab ")\n\t"                              \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
                                                                \
  "fldd  " t0 ", 104(" tab ")\n\t"                              \
  "fldd  " t1 ", 96(" tab ")\n\t"                               \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
                                                                \
  "fldd  " t0 ", 88(" tab ")\n\t"                               \
  "fldd  " t1 ", 80(" tab ")\n\t"                               \
  "fmad  " t2 ", " out ", " t0 ", " t2 "\n\t"                   \
  "fmad  " t3 ", " out ", " t1 ", " t3 "\n\t"                   \
                                                                \
  "fnmsd  " t2 ", " x2 ", " t3 ", " out "\n\t"                  \

#define asm_vinv_sqrtd(y, x, t0, t1, t2, t3, tab)       \
  "ldde    " t0 ", 16(" tab ")\n\t"                     \
  "ldde    " t1 ", 8(" tab ")\n\t"                      \
  "ldde    " t2 ", 24(" tab ")\n\t"                     \
  "ldde    " t3 ", 0(" tab ")\n\t"                      \
  "vsubl   " x ", " t1 ", " t1 "\n\t"                   \
                                                        \
  "#t1 = x / 2, t3 = 1.5\n\t"                           \
                                                        \
  "vlog2x4 " x ", " t0 ", " t0 "\n\t"                   \
  "srlow   " t0 ", 1, " t0 "\n\t"                       \
  "vsubl   " t2 ", " t0 ", " y "\n\t"                   \
                                                        \
  "vmuld   " t3 ", " y  ", " t2 "\n\t"                  \
  "vmuld   "  y ", " y  ", " t0 "\n\t"                  \
  "vmuld   " t1 ", " y  ", " y "\n\t"                   \
  "vnmad   " t0 ", " y ", " t2 ", " y"\n\t"             \
                                                        \
  "vmuld   " t3 ", " y  ", " t2 "\n\t"                  \
  "vmuld   "  y ", " y  ", " t0 "\n\t"                  \
  "vmuld   " t1 ", " y  ", " y "\n\t"                   \
  "vnmad   " t0 ", " y ", " t2 ", " y"\n\t"             \
                                                        \
  "vmuld   " t3 ", " y  ", " t2 "\n\t"                  \
  "vmuld   "  y ", " y  ", " t0 "\n\t"                  \
  "vmuld   " t1 ", " y  ", " y "\n\t"                   \
  "vnmad   " t0 ", " y ", " t2 ", " y "\n\t"            \


/* #define asm_vtrailing_aggresive(y, x, t0, t1, t2, t3)   \ */
/*   "vextf   " x  ", 1, " t1 "\n\t"                       \ */
/*   "vextf   " x  ", 2, " t2 "\n\t"                       \ */
/*   "vextf   " x  ", 3, " t3 "\n\t"                       \ */
/*                                                         \ */
/*   "fcvtdlr " x  ", 0, " t4 */
/*   "fcvtdlr " x  ", " y  */
