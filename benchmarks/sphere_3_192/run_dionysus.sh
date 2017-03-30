/usr/bin/time -l ~/Source/Dionysus/examples/cohomology/rips-pairwise-cohomology -s2 -p2 ~/Bitbucket/phat-paper/benchmark/point\ cloud/sphere_3_192_points.dat 2>&1 | tee dionysus.dim_1.out.txt
/usr/bin/time -l ~/Source/Dionysus/examples/cohomology/rips-pairwise-cohomology -s3 -p2 ~/Bitbucket/phat-paper/benchmark/point\ cloud/sphere_3_192_points.dat 2>&1 | tee dionysus.dim_2.out.txt

