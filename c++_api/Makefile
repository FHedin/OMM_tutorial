CC=gcc
CXX=g++
FC=gfortran

INC_PATH=/usr/local/openmm/include
LD_PATH=/usr/local/openmm/lib
LIB_NAME=-lOpenMM

all:
	$(CXX) -O2 -march=native -I$(INC_PATH) -L$(LD_PATH) -Wl,-rpath,$(LD_PATH) HelloSodiumChloride.cpp  -o HelloSodiumChloride    $(LIB_NAME)
	$(CC)  -O2 -march=native -I$(INC_PATH) -L$(LD_PATH) -Wl,-rpath,$(LD_PATH) HelloSodiumChlorideInC.c -o HelloSodiumChlorideInC $(LIB_NAME) -lstdc++
	$(FC)  -O2 -march=native -I$(INC_PATH) -L$(LD_PATH) -Wl,-rpath,$(LD_PATH) HelloSodiumChlorideInFortran.f90 -o HelloSodiumChlorideInFortran $(LIB_NAME) -lstdc++

clean:
	rm *.mod HelloSodiumChloride HelloSodiumChlorideInC HelloSodiumChlorideInFortran

benchmark:
	time HelloSodiumChloride > cpp.pdb
	time HelloSodiumChlorideInC > c.pdb
	time HelloSodiumChlorideInFortran > f90.pdb

