LIBS:= ntl gmp 
LIBFLAGS:=$(addprefix -l, $(LIBS));
CCFLAGS:= -Wall -std=c++11  -g -fopenmp
JAVALIBPATH1:= /usr/lib/jvm/java-1.8.0-openjdk-amd64/include/
JAVALIBPATH2:= /usr/lib/jvm/java-1.8.0-openjdk-amd64/include/linux

objects = test_shamir_share test_distributed_key_gene test_paillier libDistpaillier.so test_c_java_cast test_ntl

CC = clang++

.PHONY : clean
clean :
	-rm -i -f $(objects) *.o *.so

test_shamir_share : test_shamir_share.cpp shamir_secret_sharing_integer.cpp shamir_secret_sharing.cpp
	$(CC) $(CCFLAGS) $^ -o $@ $(LIBFLAGS) 

test_distributed_key_gene : test_distributed_key_gene.cpp distributed_paillier.cpp shamir_secret_sharing_integer.cpp shamir_secret_sharing.cpp
	$(CC) $(CCFLAGS) $^ -o $@ $(LIBFLAGS) 

libDistpaillier.so : DistributedPaillier.cpp shamir_secret_sharing_integer.cpp shamir_secret_sharing.cpp distributed_paillier.cpp
	rm -f jnilib/*
	$(CC) -fPIC $(CCFLAGS) -I$(JAVALIBPATH1)  -I$(JAVALIBPATH2) -shared  $^ -o $@ $(LIBFLAGS) 
	mv $@ ../jnilib

test_paillier : test_paillier.c paillier.c
	$(CC) $(CCFLAGS) $^ -o $@ $(LIBFLAGS) 

test_c_java_cast : test_c_java_cast.cpp
	$(CC) $(CCFLAGS) $^ -o $@ $(LIBFLAGS) 

test_ntl : test_ntl.cpp
	$(CC) $(CCFLAGS) $^ -o $@ $(LIBFLAGS) 