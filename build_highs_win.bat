set HIGHS_TAG=v1.14.0

git clone --depth 1 --branch %HIGHS_TAG% https://github.com/ERGO-Code/HiGHS.git

cd highs
md build
cd build
cmake ..  -DBUILD_SHARED_LIBS=ON -DBUILD_EXAMPLES=OFF -DZLIB=OFF -DHIPO=ON -DBUILD_OPENBLAS=ON  -G "Visual Studio 17 2022" -DCMAKE_RELEASE_POSTFIX="-ampl"
cmake --build . --config Release --target package --parallel 14

cd ..
md build-gpu
cd build-gpu
cmake ..  -DBUILD_SHARED_LIBS=ON -DBUILD_EXAMPLES=OFF -DZLIB=OFF -DCUPDLP_GPU=ON -DHIGHS_GPU_LIB=ON -DFIND_CUDA=ON -G "Visual Studio 17 2022" -DCMAKE_RELEASE_POSTFIX="-ampl-gpu"
cmake --build . --config Release --target package --parallel 14