# Threshold-Paillier-Without-Trust-Dealer
基于阈值解密和分布式秘钥生成的Pailliar实现。
依赖 NTL，GMP和 (可选的)OpenMP。
由于分布式 RSA 模数生成模块基于 BGW 协议，因此该加密系统对于诚实的大多数是安全的。

## 安装指南

本项目依赖

#### 1.系统和编辑器要求

1. 本项目使用cmake作为编译工具，需要先安装cmake，详情参考cmake官网。
2. 依赖libm，pthread等包，大部分Linux和Mac平台已内置，如果出现包缺失，请手动安装。
3. 目前我们只对部分平台和编译期进行了验证，包括
   - Linux平台  Ubuntu 18，Ubuntu 20 ，编译器采用GCC 9.3
   - MacOS 平台 MacOS 10.14， 10.15 编译期 Clang 11.0
   - Windows平台   to be continued

#### 2. 编译与链接

- 切换到项目根目录并执行配置
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release  # 如需采用debug模式，将 Release替换为Debug
```
- 编译和链接
```bash
cmake --build build
```

如果想指定make参数，进行精细配置，也可以采用make方式编译

```bash
cd ./build
make clean
make CXXFLAGS=-O3 -j4
```

编译成功后可以在 build 目录下看到

```shell
libdistpaillier.so # macos系统中是libdistpaillier.dylib
```

文件，复制到算法包中即可作为JNI的依赖使用，操作详情参考fedlearn项目README



## 参考文献: 
1. Nishide, T., & Sakurai, K. (2010, August). Distributed paillier cryptosystem without trusted dealer. In International Workshop on Information Security Applications (pp. 44-60). Springer, Berlin, Heidelberg.

