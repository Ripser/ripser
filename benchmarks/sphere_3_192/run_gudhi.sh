/usr/bin/time -l ~/Source/Gudhi_library_1.3.1/example/Persistent_cohomology/rips_persistence -d2 -p2 ~/Bitbucket/phat-paper/benchmark/point\ cloud/sphere_3_192_points.dat 2>&1 | tee gudhi.dim_1.out.txt
/usr/bin/time -l ~/Source/Gudhi_library_1.3.1/example/Persistent_cohomology/rips_persistence -d3 -p2 ~/Bitbucket/phat-paper/benchmark/point\ cloud/sphere_3_192_points.dat 2>&1 | tee gudhi.dim_2.out.txt

