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
    image="mapgccv/manylinux_2_28_cuda12:latest"
    builddir="build-x64-gpu"
    nthreads=2
    user=""
    preprocdef="-DCUPDLP_GPU=on  -DHIPO=OFF -DHIGHS_GPU_LIB=ON -DCUPDLP_FIND_CUDA=OFF  -DCMAKE_RELEASE_POSTFIX=-ampl-gpu"
elif [ "$1" = "arm64" ]; then
    image="quay.io/pypa/manylinux_2_28_aarch64"
    builddir="build-arm64"
    nthreads=4
    user="-u $(id -u):$(id -g)"
    # Since v1.16 HiGHS forces OpenBLAS's DYNAMIC_ARCH on. With GCC >= 14,
    # OpenBLAS 0.3.30 then adds the ARMV9SME core, whose CMake build calls
    # the SME sgemm_direct kernels without compiling them (undefined
    # sgemm_direct_sme1_* when linking). DYNAMIC_LIST is OpenBLAS's default
    # list (ARMV8 is always added) without ARMV9SME. The single quotes keep
    # the ;-separated list in one piece through the docker bash -c.
    preprocdef="-DCUPDLP_GPU=off -DHIPO=ON -DBUILD_SHARED_EXTRAS_LIB=OFF -DBUILD_OPENBLAS=ON -DTARGET=GENERIC -DDYNAMIC_ARCH=OFF '-DDYNAMIC_LIST=CORTEXA53;CORTEXA57;THUNDERX;THUNDERX2T99;TSV110;EMAG8180;NEOVERSEN1;THUNDERX3T110;NEOVERSEV1;NEOVERSEN2;ARMV8SVE;A64FX' -DCMAKE_RELEASE_POSTFIX=-ampl"
elif [ "$1" = "arm64-gpu" ]; then
    image="mapgccv/manylinux_2_28_cuda12_aarch64:latest"
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
