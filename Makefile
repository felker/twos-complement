CXX := mpicc
CXXFLAGS := -acc -Minline -Minfo=accel -ta=tesla,cuda8.0,cc60,pinned

all: mandelbrot

mandelbrot: mandelbrot_mpi.o
	$(CXX) $(CXXFLAGS) -o $@

mandelbrot_mpi.o: mandelbrot_mpi.c
	$(CXX) -c $<

clean:
	rm *o mandelbrot
