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
RUN curl -LO https://raw.githubusercontent.com/Ripser/ripser/dev/benchmarks/sphere_3_192/sphere_3_192_points.dat  \ 
&& curl -LO https://raw.githubusercontent.com/Ripser/ripser/dev/benchmarks/sphere_3_192/sphere_3_192.distance_matrix \
&& curl -LO https://raw.githubusercontent.com/Ripser/ripser/dev/benchmarks/sphere_3_192/sphere_3_192.complex \
\
&& curl -LO https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_distmat/fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt \
&& curl -LO https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_distmat/dipha/fractal_9_5_2_random.bin \
\
&& curl -LO https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_point_cloud/dragon_vrip.ply.txt_2000_.txt \
\
# && curl -LO https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_distmat/human_gene_distmat.txt \
\
&& curl -LO https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_point_cloud/random_point_cloud_50_16_.txt \
&& curl -LO https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_distmat/dipha/random_50_16.bin

COPY o3_1024.txt .

FROM benchmark-setup as benchmark-ripser

WORKDIR /ripser
RUN git clone https://github.com/Ripser/ripser.git \
&& cd ripser \
&& make

ENV PATH="${PATH}:/ripser/ripser"

WORKDIR /benchmark
RUN time -v -o sphere_3_192.ripser.txt    ripser sphere_3_192.distance_matrix --dim 2
RUN time -v -o random.ripser.txt    ripser random_point_cloud_50_16_.txt --format point-cloud --dim 7
RUN time -v -o fractal-r.ripser.txt    ripser fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt --dim 2
RUN time -v -o dragon-2.ripser.txt    ripser dragon_vrip.ply.txt_2000_.txt --format point-cloud --dim 1
# RUN time -v -o genome.ripser.txt    ripser human_gene_distmat.txt --dim 1
RUN time -v -o o3_1024.ripser.txt    ripser o3_1024.txt --format point-cloud --dim 3 --threshold 1.8 --ratio 2



FROM benchmark-setup as benchmark-gudhi

WORKDIR /gudhi

RUN curl -LO "https://gforge.inria.fr/frs/download.php/file/37696/2018-09-04-14-25-00_GUDHI_2.3.0.tar.gz" \
&& tar xf 2018-09-04-14-25-00_GUDHI_2.3.0.tar.gz \
&& cd 2018-09-04-14-25-00_GUDHI_2.3.0 \
&& mkdir build \
&& cd build \
&& cmake -DCMAKE_BUILD_TYPE="Release" .. \
&& make rips_distance_matrix_persistence 

ENV PATH="${PATH}:/gudhi/2018-09-04-14-25-00_GUDHI_2.3.0/build/utilities/Rips_complex"

WORKDIR /benchmark
RUN time -v -o sphere_3_192.gudhi.txt     rips_distance_matrix_persistence -d3 -p2 sphere_3_192.distance_matrix



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



FROM benchmark-dipha as benchmark-dipha-multicore

RUN time -v -o sphere_3_192.dipha-multicore.txt     mpiexec --allow-run-as-root --mca btl_vader_single_copy_mechanism none \
dipha --dual --benchmark --upper_dim 3 sphere_3_192.complex sphere_3_192.dipha-multicore.diagram



FROM benchmark-setup as benchmark-eirene037

WORKDIR /eirene

RUN curl -LO https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz \
&& tar xf julia-0.6.4-linux-x86_64.tar.gz

RUN curl -LO "http://gregoryhenselman.org/eirene/ewExternalFiles/Eirene0_3_7.zip" \
&& unzip Eirene0_3_7.zip \
&& cd Eirene0_3_7

# ENV PATH="${PATH}:/eirene/julia-9d11f62bcb/bin"

# RUN julia --eval 'Pkg.clone("https://github.com/Eetion/Eirene.jl.git"); using Eirene'
# RUN julia --eval ' \
#	Pkg.clone("https://github.com/bnels/PlotlyJS.jl.git"); \
#	Pkg.checkout("PlotlyJS","fix-jsstring"); \
#	Pkg.add("Eirene");'
RUN /eirene/julia-9d11f62bcb/bin/julia --eval ' \
	Pkg.clone("https://github.com/bnels/PlotlyJS.jl.git"); \
	Pkg.checkout("PlotlyJS","fix-jsstring");'
RUN /eirene/julia-9d11f62bcb/bin/julia --eval 'Pkg.add("Distances"); Pkg.add("JLD"); Pkg.add("Blink"); Pkg.add("Plotly"); Pkg.add("MultivariateStats");'
RUN /eirene/julia-9d11f62bcb/bin/julia --load /eirene/Eirene0_3_7/Eirene0_3_7_mutable.jl --eval 'eirene([0 1; 1 0]);'

WORKDIR /benchmark
RUN time -v -o sphere_3_192.eirene.txt \
    /eirene/julia-9d11f62bcb/bin/julia --load /eirene/Eirene0_3_7/Eirene0_3_7_mutable.jl --eval 'eirene([0 1; 1 0]); p = readdlm("sphere_3_192.distance_matrix"); \
    tic(); eirene(p, record="none", bettimax=2); toc();' >sphere_3_192.eirene037.out.txt
RUN time -v -o dragon-2.eirene.txt \
    /eirene/julia-9d11f62bcb/bin/julia --load /eirene/Eirene0_3_7/Eirene0_3_7_mutable.jl --eval 'eirene([0 1; 1 0]); p = readdlm("dragon_vrip.ply.txt_2000_.txt"); \
    tic(); eirene(p, rowsare="points", record="none", bettimax=1); toc();' >dragon-2.eirene037.out.txt
RUN time -v -o fractal-r.eirene.txt \
    /eirene/julia-9d11f62bcb/bin/julia --load /eirene/Eirene0_3_7/Eirene0_3_7_mutable.jl --eval 'eirene([0 1; 1 0]); p = readdlm("fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt"); \
    tic(); eirene(p, record="none", bettimax=2); toc();' >fractal-r.eirene037.out.txt
RUN time -v -o random.eirene.txt  \
    /eirene/julia-9d11f62bcb/bin/julia --load /eirene/Eirene0_3_7/Eirene0_3_7_mutable.jl --eval 'eirene([0 1; 1 0]); p = readdlm("random_point_cloud_50_16_.txt"); \
    tic(); eirene(p, rowsare="points", record="none", bettimax=7); toc();' >random.eirene037.out.txt



FROM benchmark-setup as benchmark-eirene

WORKDIR /eirene

#RUN curl -LO https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.0-linux-x86_64.tar.gz \
#&& tar xf julia-1.1.0-linux-x86_64.tar.gz
#ENV PATH="${PATH}:/eirene/julia-1.1.0/bin"

#RUN julia --eval 'using Pkg; Pkg.add("Plotly"); Pkg.add("Eirene"); using Eirene'

WORKDIR /benchmark
RUN time -v -o sphere_3_192.eirene.txt \
	julia --eval 'using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed eirene(transpose(readdlm("sphere_3_192.distance_matrix")), record="none", maxdim=2))' >sphere_3_192.eirene.out.txt
RUN time -v -o dragon-2.eirene.txt \
	julia --eval 'using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed eirene(transpose(readdlm("dragon_vrip.ply.txt_2000_.txt")), model="pc", record="none", maxdim=1))' >dragon-2.eirene.out.txt
RUN time -v -o fractal-r.eirene.txt \
	julia --eval 'using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed eirene(transpose(readdlm("fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt")), record="none", maxdim=2))' >fractal-r.eirene.out.txt
#RUN time -v -o random.eirene.txt  \
#	julia --eval 'using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed eirene(transpose(readdlm("random_point_cloud_50_16_.txt")), model="pc", record="none", maxdim=7))' >random.eirene.out.txt
RUN time -v -o o3_1024.eirene.txt \
	julia --eval 'using Eirene; using DelimitedFiles; eirene([0 1; 1 0]); println(@elapsed eirene(transpose(readdlm("o3_1024.txt")), model="pc", record="none", maxdim=3, maxrad=1.8))' > o3_1024.eirene.out.txt
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

RUN time -v -o sphere.dionysus-py.txt python -c 'import dionysus as d; from numpy import loadtxt, inf; import math; from scipy.spatial.distance import squareform; \
	d.cohomology_persistence(d.fill_rips(squareform(loadtxt("sphere_3_192.distance_matrix")), 3, inf));'

RUN time -v -o dragon-2.dionysus-py.txt python -c 'import dionysus as d; from numpy import loadtxt, inf; import math; from scipy.spatial.distance import squareform; \
	d.cohomology_persistence(d.fill_rips(loadtxt("dragon_vrip.ply.txt_2000_.txt"), 2, inf));'

RUN time -v -o sphere.dionysus-py.txt python -c 'import dionysus as d; from numpy import loadtxt, inf; import math; from scipy.spatial.distance import squareform; \
	d.cohomology_persistence(d.fill_rips(squareform(loadtxt("fractal_9_5_2_random_edge_list.txt_0.19795_distmat.txt")), 3, inf));'
