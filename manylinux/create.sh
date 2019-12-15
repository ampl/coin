#!/bin/bash
cd `dirname $0`

docker pull quay.io/pypa/manylinux1_i686
docker pull quay.io/pypa/manylinux1_x86_64

docker rm manylinux64 manylinux64_tmp
docker run -v `pwd`:/manylinux --name manylinux64_tmp quay.io/pypa/manylinux1_x86_64 /manylinux/bootstrap.sh
docker commit manylinux64_tmp fdabrandao/manylinux64:coin_cmake-3.12.4
docker rm manylinux64_tmp

docker rm manylinux32 manylinux32_tmp
docker run -v `pwd`:/manylinux --name manylinux32_tmp quay.io/pypa/manylinux1_i686 /manylinux/bootstrap.sh
docker commit manylinux32_tmp fdabrandao/manylinux32:coin_cmake-3.12.4
docker rm manylinux32_tmp

# Upload images:
# $ docker login --username=fdabrandao
# $ docker push fdabrandao/manylinux32:coin_cmake-3.12.4
# $ docker push fdabrandao/manylinux64:coin_cmake-3.12.4
