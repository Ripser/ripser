/usr/bin/time -l ~/Bitbucket/dipha/dipha --benchmark --upper_dim 2 --dual ~/Bitbucket/phat-paper/benchmark/dipha/sphere_3_192.complex /dev/null 2>&1 | tee dipha.dim_1.out.txt
/usr/bin/time -l ~/Bitbucket/dipha/dipha --benchmark --upper_dim 3 --dual ~/Bitbucket/phat-paper/benchmark/dipha/sphere_3_192.complex /dev/null 2>&1 | tee dipha.dim_2.out.txt

