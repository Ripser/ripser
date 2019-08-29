# FROM ubuntu:latest

#FROM ubuntu:18.10 as benchmark-setup
# RUN echo "deb http://archive.ubuntu.com/ubuntu/ cosmic main restricted universe" > /etc/apt/sources.list

FROM ubuntu:19.04 as benchmark-setup

RUN apt-get update && \
    apt-get install -y \
    apt-utils \
    cmake \
    curl \
    g++ \
    git \
    gudhi-utils \
    hdf5-tools \
    libcgal-dev \
    libboost-dev \
    libboost-filesystem-dev \
    libboost-test-dev \
    libgmp3-dev \
    libmpfr-dev \
    libtbb-dev \
    libopenmpi-dev \
    make \
    python-minimal \
    python-pip \
    time \
    unzip \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /external

RUN curl -LO https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.0-linux-x86_64.tar.gz \
	&& tar xf julia-1.1.0-linux-x86_64.tar.gz
ENV PATH="${PATH}:/external/julia-1.1.0/bin"

RUN julia --eval 'using Pkg; Pkg.add("Plotly"); Pkg.add("Eirene");  using Eirene;'


WORKDIR /benchmark
COPY sphere_3_192_points.dat .
COPY sphere_3_192.complex .
COPY o3_1024.txt .
COPY o3_1024_1.8.dipha .
COPY o3_4096.txt .
# COPY o3_4096_1.4.dipha .
COPY clifford_torus_50000.points.txt .
# COPY dragon_2000.dipha .



RUN curl -LO https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_distmat/fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt \
&& curl -LO https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_distmat/dipha/fractal_9_5_2_random.bin \
\
&& curl -LO https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_point_cloud/dragon_vrip.ply.txt_2000_.txt \
\
# && curl -LO https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_distmat/human_gene_distmat.txt \
\
&& curl -LO https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_point_cloud/random_point_cloud_50_16_.txt \
&& curl -LO https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_distmat/dipha/random_50_16.bin

RUN (echo "OFF\n2000 0 0" && cat dragon_vrip.ply.txt_2000_.txt) >dragon_vrip.ply.txt_2000_.off \
&& (echo "OFF\n50 0 0" && cat random_point_cloud_50_16_.txt) >random_point_cloud_50_16_.off \
&& (echo "OFF\n1024 0 0" && cat sphere_3_192_points.dat) > sphere_3_192.off \
&& (echo "OFF\n1024 0 0" && cat o3_1024.txt) >o3_1024.off \
&& (echo "OFF\n4096 0 0" && cat o3_4096.txt) >o3_4096.off \
&& (echo "nOFF\n4 50000 0 0" && cat clifford_torus_50000.points.txt) >clifford_torus_50000.points.off


FROM benchmark-setup as benchmark-ripser

WORKDIR /ripser
RUN git clone https://github.com/Ripser/ripser.git \
&& cd ripser \
&& git checkout dev \
&& make

ENV PATH="${PATH}:/ripser/ripser"

WORKDIR /benchmark
RUN time -v -o sphere_3_192.ripser.txt    ripser sphere_3_192_points.dat --format point-cloud --dim 2
RUN time -v -o o3_1024.ripser.txt    ripser o3_1024.txt --format point-cloud --dim 3 --threshold 1.8 --ratio 2
RUN time -v -o o3_4096.ripser.txt    ripser o3_4096.txt --format point-cloud --dim 3 --threshold 1.4 --ratio 2
RUN time -v -o clifford_torus_50000.ripser.txt    ripser clifford_torus_50000.points.txt --format point-cloud --dim 2 --threshold .15 --ratio 2

RUN time -v -o random.ripser.txt    ripser random_point_cloud_50_16_.txt --format point-cloud --dim 7
RUN time -v -o fractal-r.ripser.txt    ripser fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt --dim 2
RUN time -v -o dragon-2.ripser.txt    ripser dragon_vrip.ply.txt_2000_.txt --format point-cloud --dim 1
# RUN time -v -o genome.ripser.txt    ripser human_gene_distmat.txt --dim 1


FROM benchmark-ripser as benchmark-ripser-no-emergent

WORKDIR /ripser/ripser
RUN git pull && git checkout benchmarks/disable-emergent-pairs \
&& make

WORKDIR /benchmark
RUN time -v -o sphere_3_192.ripser-no-emergent.txt    ripser sphere_3_192_points.dat --format point-cloud --dim 2
RUN time -v -o random.ripser-no-emergent.txt    ripser random_point_cloud_50_16_.txt --format point-cloud --dim 7
RUN time -v -o fractal-r.ripser-no-emergent.txt    ripser fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt --dim 2
RUN time -v -o dragon-2.ripser-no-emergent.txt    ripser dragon_vrip.ply.txt_2000_.txt --format point-cloud --dim 1
# RUN time -v -o genome.ripser-no-emergent.txt    ripser human_gene_distmat.txt --dim 1
RUN time -v -o o3_1024.ripser-no-emergent.txt    ripser o3_1024.txt --format point-cloud --dim 3 --threshold 1.8 --ratio 2
RUN time -v -o o3_4096.ripser-no-emergent.txt    ripser o3_4096.txt --format point-cloud --dim 3 --threshold 1.4 --ratio 2
RUN time -v -o clifford_torus_50000.ripser-no-emergent.txt    ripser clifford_torus_50000.points.txt --format point-cloud --dim 2 --threshold .15 --ratio 2


FROM benchmark-ripser as benchmark-ripser-store-reduced

WORKDIR /ripser/ripser
RUN git pull && git checkout benchmarks/store-reduced-matrix \
&& make

WORKDIR /benchmark
RUN time -v -o sphere_3_192.ripser-store-reduced.txt    ripser sphere_3_192_points.dat --format point-cloud --dim 2
RUN time -v -o random.ripser-store-reduced.txt    ripser random_point_cloud_50_16_.txt --format point-cloud --dim 7
# RUN time -v -o fractal-r.ripser-store-reduced.txt    ripser fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt --dim 2
# RUN time -v -o dragon-2.ripser-store-reduced.txt    ripser dragon_vrip.ply.txt_2000_.txt --format point-cloud --dim 1
# RUN time -v -o genome.ripser-store-reduced.txt    ripser human_gene_distmat.txt --dim 1
RUN time -v -o o3_1024.ripser-store-reduced.txt    ripser o3_1024.txt --format point-cloud --dim 3 --threshold 1.8 --ratio 2
# RUN time -v -o o3_4096.ripser-store-reduced.txt    ripser o3_4096.txt --format point-cloud --dim 3 --threshold 1.4 --ratio 2
# RUN time -v -o clifford_torus_50000.ripser-store-reduced.txt    ripser clifford_torus_50000.points.txt --format point-cloud --dim 2 --threshold .15 --ratio 2


FROM benchmark-ripser as benchmark-ripser-use-reduced

WORKDIR /ripser/ripser
RUN git pull && git checkout benchmarks/use-reduced-matrix \
&& make

WORKDIR /benchmark
RUN time -v -o sphere_3_192.ripser-use-reduced.txt    ripser sphere_3_192_points.dat --format point-cloud --dim 2
RUN time -v -o random.ripser-use-reduced.txt    ripser random_point_cloud_50_16_.txt --format point-cloud --dim 7
# RUN time -v -o fractal-r.ripser-use-reduced.txt    ripser fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt --dim 2
# RUN time -v -o dragon-2.ripser-use-reduced.txt    ripser dragon_vrip.ply.txt_2000_.txt --format point-cloud --dim 1
# RUN time -v -o genome.ripser-use-reduced.txt    ripser human_gene_distmat.txt --dim 1
RUN time -v -o o3_1024.ripser-use-reduced.txt    ripser o3_1024.txt --format point-cloud --dim 3 --threshold 1.8 --ratio 2
# RUN time -v -o o3_4096.ripser-use-reduced.txt    ripser o3_4096.txt --format point-cloud --dim 3 --threshold 1.4 --ratio 2
# RUN time -v -o clifford_torus_50000.ripser-use-reduced.txt    ripser clifford_torus_50000.points.txt --format point-cloud --dim 2 --threshold .15 --ratio 2


FROM benchmark-setup as benchmark-gudhi

# WORKDIR /gudhi
# 
# RUN curl -LO "https://gforge.inria.fr/frs/download.php/file/37696/2018-09-04-14-25-00_GUDHI_2.3.0.tar.gz" \
# && tar xf 2018-09-04-14-25-00_GUDHI_2.3.0.tar.gz \
# && cd 2018-09-04-14-25-00_GUDHI_2.3.0 \
# && mkdir build \
# && cd build \
# && cmake -DCMAKE_BUILD_TYPE="Release" .. \
# && make rips_distance_matrix_persistence 
# 
# ENV PATH="${PATH}:/gudhi/2018-09-04-14-25-00_GUDHI_2.3.0/build/utilities/Rips_complex"
#
# WORKDIR /benchmark
# RUN time -v -o sphere_3_192.gudhi.txt     rips_distance_matrix_persistence -d3 -p2 sphere_3_192.distance_matrix

WORKDIR /benchmark
RUN time -v -o sphere_3_192.gudhi.txt     gudhi-rips-persistence -d3 -p2 sphere_3_192.off >sphere_3_192.gudhi.out.txt
# RUN time -v -o fractal-r.gudhi.txt     gudhi-rips-distance-matrix-persistence -d3 -p2 fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt >fractal-r.gudhi.out.txt
RUN time -v -o o3_1024.gudhi.txt     gudhi-rips-persistence -d4 -p2 -r 1.8 o3_1024.off >o3_1024.gudhi.out.txt
RUN time -v -o o3_4096.gudhi.txt     gudhi-rips-persistence -d4 -p2 -r 1.4 o3_4096.off >o3_4096.gudhi.out.txt
RUN time -v -o clifford_torus_50000.gudhi.txt     gudhi-rips-persistence -d3 -p2 -r 0.15 clifford_torus_50000.points.off >clifford_torus_50000.gudhi.out.txt


FROM benchmark-setup as benchmark-dipha

WORKDIR /dipha

RUN git clone https://github.com/DIPHA/dipha.git \
&& cd dipha \
&& mkdir build \
&& cd build \
&& cmake .. \
&& make

ENV PATH="${PATH}:/dipha/dipha/build"

WORKDIR /benchmark
RUN time -v -o sphere_3_192.dipha.txt     dipha --dual --benchmark --upper_dim 3 sphere_3_192.complex sphere_3_192.dipha.diagram
RUN time -v -o o3_1024.dipha.txt     dipha --dual --benchmark --upper_dim 4 o3_1024_1.8.dipha o3_1024_1.8.dipha.diagram
# RUN time -v -o o3_4096.dipha.txt     dipha --dual --benchmark --upper_dim 4 o3_4096_1.4.dipha o3_4096_1.4.dipha.diagram
# RUN time -v -o dragon-2.dipha.txt     dipha --dual --benchmark --upper_dim 2 dragon_2000.dipha dragon_2000.dipha.diagram



# FROM benchmark-dipha as benchmark-dipha-multicore

RUN time -v -o sphere_3_192.dipha-multicore.txt     mpiexec --allow-run-as-root --mca btl_vader_single_copy_mechanism none \
dipha --dual --benchmark --upper_dim 3 sphere_3_192.complex sphere_3_192.dipha-multicore.diagram

RUN time -v -o o3_1024.dipha-multicore.txt     mpiexec --allow-run-as-root --mca btl_vader_single_copy_mechanism none \
dipha --dual --benchmark --upper_dim 4 o3_1024_1.8.dipha o3_1024_1.8.dipha-multicore.diagram


### FROM benchmark-setup as benchmark-eirene037
### 
### WORKDIR /eirene
### 
### RUN curl -LO https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz \
### && tar xf julia-0.6.4-linux-x86_64.tar.gz
### 
### RUN curl -LO "http://gregoryhenselman.org/eirene/ewExternalFiles/Eirene0_3_7.zip" \
### && unzip Eirene0_3_7.zip \
### && cd Eirene0_3_7
### 
### # ENV PATH="${PATH}:/eirene/julia-9d11f62bcb/bin"
### 
### # RUN julia --eval 'Pkg.clone("https://github.com/Eetion/Eirene.jl.git"); using Eirene'
### # RUN julia --eval ' \
### #	Pkg.clone("https://github.com/bnels/PlotlyJS.jl.git"); \
### #	Pkg.checkout("PlotlyJS","fix-jsstring"); \
### #	Pkg.add("Eirene");'
### RUN /eirene/julia-9d11f62bcb/bin/julia --eval ' \
### 	Pkg.clone("https://github.com/bnels/PlotlyJS.jl.git"); \
### 	Pkg.checkout("PlotlyJS","fix-jsstring");'
### RUN /eirene/julia-9d11f62bcb/bin/julia --eval 'Pkg.add("Distances"); Pkg.add("JLD"); Pkg.add("Blink"); Pkg.add("Plotly"); Pkg.add("MultivariateStats");'
### RUN /eirene/julia-9d11f62bcb/bin/julia --load /eirene/Eirene0_3_7/Eirene0_3_7_mutable.jl --eval 'eirene([0 1; 1 0]);'
### 
### WORKDIR /benchmark
### RUN time -v -o sphere_3_192.eirene037.txt \
###     /eirene/julia-9d11f62bcb/bin/julia --load /eirene/Eirene0_3_7/Eirene0_3_7_mutable.jl --eval 'eirene([0 1; 1 0]); p = readdlm("sphere_3_192.distance_matrix"); \
###     tic(); eirene(p, record="none", bettimax=2); toc();' >sphere_3_192.eirene037.out.txt
### 
### RUN time -v -o dragon-2.eirene037.txt \
###     /eirene/julia-9d11f62bcb/bin/julia --load /eirene/Eirene0_3_7/Eirene0_3_7_mutable.jl --eval 'eirene([0 1; 1 0]); p = readdlm("dragon_vrip.ply.txt_2000_.txt"); \
###     tic(); eirene(p, rowsare="points", record="none", bettimax=1); toc();' >dragon-2.eirene037.out.txt
### RUN time -v -o fractal-r.eirene037.txt \
###     /eirene/julia-9d11f62bcb/bin/julia --load /eirene/Eirene0_3_7/Eirene0_3_7_mutable.jl --eval 'eirene([0 1; 1 0]); p = readdlm("fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt"); \
###     tic(); eirene(p, record="none", bettimax=2); toc();' >fractal-r.eirene037.out.txt
### RUN time -v -o random.eirene037.txt  \
###     /eirene/julia-9d11f62bcb/bin/julia --load /eirene/Eirene0_3_7/Eirene0_3_7_mutable.jl --eval 'eirene([0 1; 1 0]); p = readdlm("random_point_cloud_50_16_.txt"); \
###     tic(); eirene(p, rowsare="points", record="none", bettimax=7); toc();' >random.eirene037.out.txt



FROM benchmark-setup as benchmark-eirene

WORKDIR /eirene

#RUN curl -LO https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.0-linux-x86_64.tar.gz \
#&& tar xf julia-1.1.0-linux-x86_64.tar.gz
#ENV PATH="${PATH}:/eirene/julia-1.1.0/bin"

#RUN julia --eval 'using Pkg; Pkg.add("Plotly"); Pkg.add("Eirene"); using Eirene'

WORKDIR /benchmark

RUN time -v -o sphere_3_192.eirene.txt \
	julia --eval 'maxdim = 2; using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed (d = transpose(readdlm("sphere_3_192_points.dat")); result = eirene(d, model="pc", record="none", maxdim=maxdim))); for i in 0:maxdim print("dim $(i): "); show(stdout, "text/plain", (barcode(result , dim = i))); println(); end' >sphere_3_192.eirene.out.txt

RUN time -v -o o3_1024.eirene.txt \
	julia --eval 'maxdim = 3; maxrad = 1.8; using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed (d = transpose(readdlm("o3_1024.txt")); result = eirene(d, model="pc", record="none", maxdim=maxdim, maxrad=maxrad))); for i in 0:maxdim print("dim $(i): "); show(stdout, "text/plain", (barcode(result , dim = i))); println(); end' > o3_1024.eirene.out.txt

RUN time -v -o o3_4096.eirene.txt \
	julia --eval 'maxdim = 3; maxrad = 1.4; using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed (d = transpose(readdlm("o3_4096.txt")); result = eirene(d, model="pc", record="none", maxdim=maxdim, maxrad=maxrad))); for i in 0:maxdim print("dim $(i): "); show(stdout, "text/plain", (barcode(result , dim = i))); println(); end' > o3_4096.eirene.out.txt

# RUN time -v -o clifford_torus_50000.points.eirene.txt \
#	julia --eval 'maxdim = 2; maxrad = .15; using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed (d = transpose(readdlm("clifford_torus_50000.points.txt")); result = eirene(d, model="pc", record="none", maxdim=maxdim, maxrad=maxrad))); for i in 0:maxdim print("dim $(i): "); show(stdout, "text/plain", (barcode(result , dim = i))); println(); end' > clifford_torus_50000.eirene.out.txt

RUN time -v -o dragon-2.eirene.txt \
	julia --eval 'maxdim = 1; using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed (d = transpose(readdlm("dragon_vrip.ply.txt_2000_.txt")); result = eirene(d, model="pc", record="none", maxdim=maxdim))); for i in 0:maxdim print("dim $(i): "); show(stdout, "text/plain", (barcode(result , dim = i))); println(); end' >dragon-2.eirene.out.txt

RUN time -v -o fractal-r.eirene.txt \
	julia --eval 'maxdim = 2; using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed (d = (m->(m+transpose(m))/2)(readdlm("fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt")); result = eirene(d, record="none", maxdim=maxdim))); for i in 0:maxdim print("dim $(i): "); show(stdout, "text/plain", (barcode(result , dim = i))); println(); end' >fractal-r.eirene.out.txt

RUN time -v -o random.eirene.txt  \
	julia --eval 'maxdim = 7; maxrad = 1.7753855751520569; using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed (d = transpose(readdlm("random_point_cloud_50_16_.txt")); result = eirene(d, model="pc", record="none", maxdim=maxdim, maxrad=maxrad))); for i in 0:maxdim print("dim $(i): "); show(stdout, "text/plain", (barcode(result , dim = i))); println(); end' >random.eirene.out.txt

#
# using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed result=eirene(transpose(readdlm("/Users/uli/Bitbucket/ripser/examples/o3_1024.txt")), model="pc", record="none", maxdim=3, maxrad=1.8))
# using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed result=eirene(transpose(readdlm("/Users/uli/Bitbucket/ripser/examples/o3_2048.txt")), model="pc", record="none", maxdim=3, maxrad=1.6))
# using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed result=eirene(transpose(readdlm("/Users/uli/Bitbucket/ripser/examples/o3_4096.txt")), model="pc", record="none", maxdim=3, maxrad=1.4))
#

FROM benchmark-setup as benchmark-dionysus2

WORKDIR /dionysus2

RUN pip install numpy scipy dionysus

# RUN git clone https://github.com/mrzv/dionysus.git \
# && cd dionysus \
# && mkdir build \
# && cd build \
# && cmake .. \
# && make rips-pairwise
# 
# ENV PATH="${PATH}:/dionysus2/dionysus/build/examples/rips"

WORKDIR /benchmark
# RUN time -v -o sphere_3_192.dionysus-cpp.txt  rips-pairwise sphere_3_192_points.dat -s3

#RUN time -v -o sphere.dionysus-py.txt python -c 'import dionysus as d; from numpy import loadtxt, inf; import math; from scipy.spatial.distance import squareform; \
#	d.cohomology_persistence(d.fill_rips(squareform(loadtxt("sphere_3_192.distance_matrix")), 3, inf));'

#RUN time -v -o dragon-2.dionysus-py.txt python -c 'import dionysus as d; from numpy import loadtxt, inf; import math; from scipy.spatial.distance import squareform; \
#	d.cohomology_persistence(d.fill_rips(loadtxt("dragon_vrip.ply.txt_2000_.txt"), 2, inf));'

#RUN time -v -o fractal.dionysus-py.txt python -c 'import dionysus as d; from numpy import loadtxt, inf; import math; from scipy.spatial.distance import squareform; \
#	d.cohomology_persistence(d.fill_rips(squareform(loadtxt("fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt")), 3, inf));'


FROM benchmark-setup as benchmark-output

COPY --from=benchmark-ripser /benchmark /benchmark
COPY --from=benchmark-gudhi /benchmark /benchmark
COPY --from=benchmark-dipha /benchmark /benchmark
COPY --from=benchmark-eirene /benchmark /benchmark
# COPY --from=benchmark-eirene037 /benchmark /benchmark
COPY --from=benchmark-dionysus2 /benchmark /benchmark

COPY --from=benchmark-ripser-no-emergent /benchmark /benchmark
COPY --from=benchmark-ripser-store-reduced /benchmark /benchmark
COPY --from=benchmark-ripser-use-reduced /benchmark /benchmark

RUN ls -l /benchmark

