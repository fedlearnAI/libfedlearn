# Threshold-Paillier-without-Trust-Dealer
基于阈值解密和分布式秘钥生成的Pailliar实现。
依赖 NTL，GMP和OpenMP。
由于分布式 RSA 模数生成模块基于 BGW 协议，因此该加密系统对于诚实的大多数是安全的。

## Install guide

1. 静态编译 libntl 
./configure PREFIX='../../libntl' NTL_THREADS=on NTL_THREAD_BOOST=on 
make CXXFLAGS=-fPIC -j4
make install
build libDistPaillier:
cd build

// release mode
cmake -DCMAKE_BUILD_TYPE=Release ..
// debug mode
cmake -DCMAKE_BUILD_TYPE=Debug   ..

make clean
make CXXFLAGS=-O3 -j4

2. 在JAVA Run/Debug Configuration 中
VM options 加入-ea -Djava.library.path=./src/main/cxx/build


参考文献: Nishide, T., & Sakurai, K. (2010, August). Distributed paillier cryptosystem without trusted dealer. In International Workshop on Information Security Applications (pp. 44-60). Springer, Berlin, Heidelberg.


