CC=g++
LIBS=
IDIR =-I/usr/include/boost -I/usr/include/eigen3
CFLAGS= -w -ggdb -fdiagnostics-color=auto
OBJ = main.o Parameters.o randomsNumb.o ModelData.o
DEP = EM_Classes.h RN_Class.h
all: MS

MS: $(OBJ) 
	$(CC) $^ -o $@ $(LIBS) $(IDIR) -std=c++1y

%.o: %.cpp $(DEP) 
	$(CC) $(CFLAGS) $(IDIR) $(LIBS) -c -o $@ $< -std=c++1y

clean:
	rm -rf *o MS
