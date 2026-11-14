REM Same configurations as the windows job in .github/workflows/build.yml
set HIGHS_TAG=v1.16.0

git clone --depth 1 --branch %HIGHS_TAG% https://github.com/ERGO-Code/HiGHS.git

cd highs
md build-x64
cd build-x64
REM ONLY_CBLAS: HiGHS means to set it (FetchContent ignores the CMAKE_ARGS it
REM uses); without it OpenBLAS also builds its utest, whose post-build step
REM copies openblas.dll, which the -ampl postfix renamed
REM REQUIRE_OPTIMIZED: fail rather than silently fall back to the slow
REM MSVC-built OpenBLAS if clang-cl or Ninja is missing (install the Visual
REM Studio component "C++ Clang Compiler for Windows")
cmake .. -G "Visual Studio 17 2022" -A x64 -DBUILD_SHARED_LIBS=ON -DBUILD_EXAMPLES=OFF -DZLIB=OFF -DHIPO=ON -DBUILD_SHARED_EXTRAS_LIB=OFF -DBUILD_OPENBLAS=ON -DONLY_CBLAS=ON -DHIGHS_OPENBLAS_REQUIRE_OPTIMIZED=ON -DCMAKE_RELEASE_POSTFIX=-ampl
cmake --build . --config Release --target package --parallel 14

cd ..
md build-x64-gpu
cd build-x64-gpu
cmake .. -G "Visual Studio 17 2022" -A x64 -DBUILD_SHARED_LIBS=ON -DBUILD_EXAMPLES=OFF -DZLIB=OFF -DCUPDLP_GPU=ON -DHIGHS_GPU_LIB=ON -DCMAKE_RELEASE_POSTFIX=-ampl-gpu
cmake --build . --config Release --target package --parallel 14
