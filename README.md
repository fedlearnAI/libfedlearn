# Threshold-Paillier-without-Trust-Dealer
基于阈值解密和分布式秘钥生成的Pailliar实现。
依赖 NTL，GMP和 (可选的)OpenMP。
由于分布式 RSA 模数生成模块基于 BGW 协议，因此该加密系统对于诚实的大多数是安全的。

## Install guide

#### 1. linux平台编译
- 需要安装Cmake,  以及依赖libm，pthread等包，如果系统中没有需要手动安装

- 并切换到项目根目录并执行配置
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release  # 如需采用debug模式，将 Release替换为Debug
```
- 编译和链接
```bash
cmake --build build
```

也可以采用make方式编译

```bash
cd ./build
make clean
make CXXFLAGS=-O3 -j4
```

编译成功后可以在 build 目录下看到

```shell
libdistpaillier.so
```

文件，复制到算法包中即可作为JNI的依赖使用，操作详情参考fedlearn项目README

注：目前测试过Ubuntu18、Ubuntu20  GCC 9.3，其他版本如有问题请反馈我们



#### 2. MACOS 平台编译

- 需要安装Cmake,  以及依赖libm，pthread等包，如果系统中没有需要手动安装

  在 macos平台默认编译器为 clang

  配置编译和链接方式与linux基本相同，因mac系统默认无gmp，需要手动安装

- 安装GMP
  a. 
  ```bash
  cd thirdparty/gmp-6.2.1/
  ./configure --prefix=/PROJECT_ROOT/thirdparty/libgmp/ #PROJECT_ROOT指用户根目录
  make clean
  make -j8 CXXFLAGS=-fPIC # -j 参数指定编译使用的CPU核数，根据实际情况自行修改
  make install
  ```
  
  安装完后

链接完成后的包名为
```shell
libdistpaillier.dylib
```


#### 3.Windows 平台编译

 TO DO

## 参考文献: 
1. Nishide, T., & Sakurai, K. (2010, August). Distributed paillier cryptosystem without trusted dealer. In International Workshop on Information Security Applications (pp. 44-60). Springer, Berlin, Heidelberg.

