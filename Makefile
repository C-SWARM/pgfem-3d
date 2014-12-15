# Be sure to change the top_level_dir variable to point to this directory

include Makefile.path
all: LIB SRC

LIB:
	cd $(pfem_lib_dir); make;

SRC:
	cd $(pfem_src_dir); make depend tags all;

clean:
	cd $(pfem_lib_dir); make clean;
	cd $(pfem_src_dir); make realclean;
