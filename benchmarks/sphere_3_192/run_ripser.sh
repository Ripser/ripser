/usr/bin/time -l ~/Bitbucket/ripser/ripser ~/Bitbucket/ripser/examples/sphere_3_192.lower_distance_matrix 2>&1 | tee ripser.dim_1.out.txt
/usr/bin/time -l ~/Bitbucket/ripser/ripser ~/Bitbucket/ripser/examples/sphere_3_192.lower_distance_matrix --dim 2 2>&1 | tee ripser.dim_2.out.txt

