## Compile libntl statically 
    ./configure PREFIX='../../libntl' NTL_THREADS=on NTL_THREAD_BOOST=on 
    make CXXFLAGS=-fPIC -j4
    make install
## build libDistPaillier:
    cd build
    
    // release mode
    cmake -DCMAKE_BUILD_TYPE=Release ..
    // debug mode
    cmake -DCMAKE_BUILD_TYPE=Debug   ..
    
    make clean
    make CXXFLAGS=-O3 -j4
## 在JAVA Run/Debug Configuration 中
    VM options 加入-ea -Djava.library.path=./src/main/cxx/build