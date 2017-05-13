IDIR=/Users/MariaBJonsdottir/seqan-library-2.2.0/include
CXX=clang++
DEBUGFLAG=-g
LIBS=-lbz2 -lz 
CFLAGS= -arch x86_64 -std=c++1y -stdlib=libc++ -I$(IDIR) $(DEBUGFLAG) $(LIBS)

all:
	$(CXX) $(CFLAGS) *.cpp

clean:
	rm -rv a.out a.out*
