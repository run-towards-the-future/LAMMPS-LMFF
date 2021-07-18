##Figure 18: Strong Scaling with small processes
 #128M expr
 #133,760,000 atoms
 bsub -cache_size 0  -o out.13376wn48		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 48   -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate
 bsub -cache_size 0  -o out.13376wn96		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 96   -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate
 bsub -cache_size 0  -o out.13376wn192		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 192  -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate
 bsub -cache_size 0  -o out.13376wn384		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 384  -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate
 bsub -cache_size 0  -o out.13376wn768	    -q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 768  -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate 
 bsub -cache_size 0  -o out.13376wn1536		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 1536 -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate 
 bsub -cache_size 0  -o out.13376wn3072		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 3072 -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate 
 bsub -cache_size 0  -o out.13376wn6144		-q q_test_ss -host_stack 512 -share_size 15000 -priv_size 4 -n 6144 -np 6 -b -cgsp 64  ../../SWLMFF/lmp_sunway_big -in in.replicate 
