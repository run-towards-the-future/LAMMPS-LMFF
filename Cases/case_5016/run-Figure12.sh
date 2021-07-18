##Validation of results for a six-layer graphene system of 5016 atoms with a time length of 1 ns. ORG is the results of Tersoff and ILP potential. OPT is the results of optimized LMFF potential, SWLMFF

##SW39000(Sunway SW39000)
 #Baseline
  bsub -cache_size 0  -o out.ORG	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 2304 -b -cgsp 64  ../../Baseline/lmp_sunway_big -in in.replicate.ilp_tersoff.validation
 #SWLMFF
  bsub -cache_size 0  -o out.OPT	-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 2304 -b -cgsp 64  ../../LMFF/lmp_sunway_big -in in.replicate.validation		