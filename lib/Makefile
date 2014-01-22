include ../Makefile.path

all: renumbering vtk_io

renumbering:
	cd $(RENUMBERING_dir)/src; make clean; make;

vtk_io:
	cd $(VTK_IO_dir)/src; make clean; make;

clean:
	cd $(RENUMBERING_dir)/src; make clean;
	cd $(VTK_IO_dir)/src; make clean;
