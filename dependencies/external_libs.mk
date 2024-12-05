#
# external libraries compilation
#
libs: $(wildcard $(LIBS_DIR)/2decomp-fft/src/*.f90)
	cd $(LIBS_DIR)/2decomp-fft && make
libsclean: $(wildcard $(LIBS_DIR)/2decomp-fft/src/*.f90)
	cd $(LIBS_DIR)/2decomp-fft && make clean
