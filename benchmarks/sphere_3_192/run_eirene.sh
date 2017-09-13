/usr/bin/time -l julia --load /Users/uli/Source/Eirene0_3_5/Eirene0_3_5.jl --eval \
'p = readdlm("/Users/uli/Bitbucket/phat-paper/benchmark/point cloud/sphere_3_192_points.dat"); eirene(p, rowsare="points", bettimax=0); eirene(p, rowsare="points", bettimax=1);' \
2>&1 | tee eirene.dim_1.out.txt

/usr/bin/time -l julia --load /Users/uli/Source/Eirene0_3_5/Eirene0_3_5.jl --eval \
'p = readdlm("/Users/uli/Bitbucket/phat-paper/benchmark/point cloud/sphere_3_192_points.dat"); eirene(p, rowsare="points", bettimax=0); eirene(p, rowsare="points", bettimax=2);' \
2>&1 | tee eirene.dim_2.out.txt
