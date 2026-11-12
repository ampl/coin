#!/bin/bash
set -ex
cd `dirname $0`
toexecute=do_build_highs_simple.sh
if [ "$1" = "x64" ]; then
    image="quay.io/pypa/manylinux_2_28_x86_64"
    builddir="build-x64"
    nthreads=2
    user=""
    preprocdef="-DCUPDLP_GPU=off  -DHIPO=ON -DBUILD_SHARED_EXTRAS_LIB=OFF -DBUILD_OPENBLAS=ON  -DTARGET=GENERIC -DDYNAMIC_ARCH=OFF  -DCMAKE_RELEASE_POSTFIX=-ampl"
elif [ "$1" = "x64-gpu" ]; then
    image="mapgccv/manylinux_2_28_cuda13:latest"
    builddir="build-x64-gpu"
    nthreads=2
    user=""
    preprocdef="-DCUPDLP_GPU=on  -DHIPO=OFF -DHIGHS_GPU_LIB=ON -DCUPDLP_FIND_CUDA=OFF  -DCMAKE_RELEASE_POSTFIX=-ampl-gpu"
elif [ "$1" = "arm64" ]; then
    image="quay.io/pypa/manylinux_2_28_aarch64"
    builddir="build-arm64"
    nthreads=4
    user="-u $(id -u):$(id -g)"
    preprocdef="-DCUPDLP_GPU=off -DHIPO=ON -DBUILD_SHARED_EXTRAS_LIB=OFF -DBUILD_OPENBLAS=ON -DTARGET=GENERIC -DDYNAMIC_ARCH=OFF  -DCMAKE_RELEASE_POSTFIX=-ampl"
elif [ "$1" = "arm64-gpu" ]; then
    image="mapgccv/manylinux_2_28_cuda13_aarch64:latest"
    builddir="build-arm64-gpu"
    nthreads=4
    user="-u $(id -u):$(id -g)"
    preprocdef="-DCUPDLP_GPU=on -DHIPO=OFF -DHIGHS_GPU_LIB=ON -DCUPDLP_FIND_CUDA=OFF  -DCMAKE_RELEASE_POSTFIX=-ampl-gpu"
    toexecute=do_build_highs_simple_arm.sh
else
    echo "Usage: $0 [x64|x64-gpu|arm64|arm64-gpu]"
    exit 1
fi

docker run --rm $image /bin/bash -c "ldd --version"
docker run -v `pwd`:/coin-build $image /bin/bash -c "chmod +x /coin-build/$toexecute; /coin-build/$toexecute $builddir $nthreads $preprocdef"
