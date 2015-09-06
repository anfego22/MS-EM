CC=g++
LIBS=
IDIR =-I/usr/include/boost -I/usr/include/eigen3
CFLAGS= -w -ggdb
OBJ = main.o Parameters.o randomsNumb.o ModelData.o
DEP = EM_Classes.h RN_Class.h
all: MS

MS: $(OBJ) 
	$(CC) $^ -o $@ $(LIBS) $(IDIR) -std=c++11

%.o: %.cpp $(DEP) 
	$(CC) $(CFLAGS) $(IDIR) $(LIBS) -c -o $@ $< -std=c++11

clean:
	rm -rf *o MS
