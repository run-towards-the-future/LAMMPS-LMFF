##Table 6: Comparison with Previous Work on a Single Node of Sunway
##SW26010(Sunway SW26010)
 #Tersoff
 #Please refer to previous work [10] in paper.

##SW26010(Sunway SW26010)
 #SWLMFF
 #When run on SW26010, change the library function supported by SW39000 to the corresponding function of SW26010.
 #such as:
 #	CRTS_init() --> athread_init() ,	 	
 # 	libc_uncached_malloc() --> malloc() ,
 #  and math functions.
  rm out.SW5.1Node.Baseline
  bsub  -o out.SW5.1Node.Baseline	-q q_test_ss -host_stack 256 -share_size 6144 -priv_size 4 -n 4 -b -cgsp 64  ../../Baseline/lmp_sunway_big -in in.replicate.ilp_tersoff
  rm out.SW5.1Node.SWLMFF
  bsub  -o out.SW5.1Node.SWLMFF	    -q q_test_ss -host_stack 256 -share_size 6144 -priv_size 4 -n 4 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate

 #We also provide a version on SW26010, named SWLMFF-5GCC
  rm out.SW5.1Node.Baseline
  bsub -o out.SW5.1Node.Baseline -q q_test_ss -host_stack 256 -share_size 6144 -priv_size 4 -n 4 -b -cgsp 64  ../../Baseline/lmp_sunway_big -in in.replicate.ilp_tersoff
  rm out.SW5.1Node.SWLMFF
  bsub -o out.SW5.1Node.SWLMFF	    -q q_test_ss -host_stack 256 -share_size 6144 -priv_size 4 -n 4 -b -cgsp 64  ../../SWLMFF-5GCC/lmp_sunway_big -in in.replicate
  


##SW39000(Sunway SW39000)
 #SWLMFF	
  rm out.SW9.1Node.Baseline
  bsub -cache_size 0  -o out.SW9.1Node.Baseline -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 6 -b -cgsp 64  ../../Baseline/lmp_sunway_big -in in.replicate.ilp_tersoff	
  rm out.SW9.1Node.SWLMFF	
  bsub -cache_size 0  -o out.SW9.1Node.SWLMFF	    -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate	
