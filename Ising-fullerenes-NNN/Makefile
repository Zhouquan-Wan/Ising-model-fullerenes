CC = /usr/bin/c++ -O3 -std=c++17 -m64 -fopenmp
# BLAS_LAPACK_LIBFLAGS=-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -lpthread -L/opt/intel/compilers_and_libraries_2019.4.243/linux/compiler/lib/intel64 -Wl,--no-as-needed -liomp5 -lirc -lsvml
# BLAS_LAPACK_INCLUDEFLAGS=-I/opt/intel/mkl/include -I/opt/intel/compilers_and_libraries_2019.4.243/linux/compiler/include/intel64
IDIR = inc
ODIR = obj
SDIR = src
BDIR = bin
CFLAGS = -I $(IDIR)
_OBJ = main.o gsf_tensor.o fullerenes.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))
_DEPS = gsf_tensor.h fullerenes.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
dir_guard=@mkdir -p $(@D)
$(ODIR)/%.o: src/%.cpp $(DEPS)
	$(dir_guard)
	$(CC) -c -o $@ $< $(CFLAGS) $(BLAS_LAPACK_INCLUDEFLAGS) $(BLAS_LAPACK_LIBFLAGS)

$(BDIR)/main: $(OBJ)
	$(dir_guard)
	$(CC) -o $@ $^ $(CFLAGS) $(BLAS_LAPACK_INCLUDEFLAGS) $(BLAS_LAPACK_LIBFLAGS)

.PHONY:clean
clean:
	-rm -rf $(ODIR)/*.o  
	
