CC=g++
LIBS=
IDIR =-I/usr/include/boost -I/usr/include/eigen3
CFLAGS= -w -ggdb -fdiagnostics-color=auto
OBJ = main.o Parameters.o randomsNumb.o ModelData.o Smooth.o filterS.o
DEP = EM_Classes.h RN_Class.h
CHK_SOURCES = main.cpp Parameters.cpp randomsNumb.cpp ModelData.cpp Smooth.cpp filterS.cpp
all: MS

MS: $(OBJ) 
	$(CC) $^ -o $@ $(LIBS) $(IDIR) -std=c++1y

%.o: %.cpp $(DEP) 
	$(CC) $(CFLAGS) $(IDIR) $(LIBS) -c -o $@ $< -std=c++1y

clean:
	rm -rf *o MS

check-syntax:
	/usr/bin/g++ -o nul -S${CHK_SOURCES} -std=c++1y  
