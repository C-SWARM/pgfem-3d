include ../../../Makefile.path

CC = mpicxx
CFLAGS = -Wall -O3

AR		= ar
ARFLAGS		= -rcs

INCLUDES = $(VTK_IO_incl) $(VTK_incl) -I$(pfem_incl_dir)

SRC = PGFem3D_to_VTK.cxx
OBJ = $(SRC:.cxx=.o)
LIB = $(VTK_IO_dir)/lib/libPGFem3D_to_VTK.a

all: $(OBJ) $(LIB)

%.o: %.cxx
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

$(LIB): $(OBJ)
	$(AR) $(ARFLAGS) $(LIB) $(OBJ)

clean:
	rm -rf *.o *~ 
realclean: clean
	rm -rf $(LIB)
tags:
	etags $(SRC)
depend:
	makedepend -- $(INCLUDES) $(CFLAGS) -- $(SRC)


# DO NOT DELETE

PGFem3D_to_VTK.o: /afs/crc.nd.edu/user/m/mmosby1/SVN/research/trunk/pfem3d/PFem3d/branch/lib/VTK_IO/include/PGFem3D_to_VTK.hpp
PGFem3D_to_VTK.o: /afs/crc.nd.edu/user/m/mmosby1/SVN/research/trunk/pfem3d/PFem3d/branch/include/element.h
PGFem3D_to_VTK.o: /usr/include/stdio.h /usr/include/features.h
PGFem3D_to_VTK.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
PGFem3D_to_VTK.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
PGFem3D_to_VTK.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
PGFem3D_to_VTK.o: /usr/include/libio.h /usr/include/_G_config.h
PGFem3D_to_VTK.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
PGFem3D_to_VTK.o: /usr/include/bits/sys_errlist.h
PGFem3D_to_VTK.o: /afs/crc.nd.edu/user/m/mmosby1/SVN/research/trunk/pfem3d/PFem3d/branch/include/supp.h
PGFem3D_to_VTK.o: /afs/crc.nd.edu/user/m/mmosby1/SVN/research/trunk/pfem3d/PFem3d/branch/include/node.h
PGFem3D_to_VTK.o: /afs/crc.nd.edu/user/m/mmosby1/SVN/research/trunk/pfem3d/PFem3d/branch/include/sig.h
PGFem3D_to_VTK.o: /afs/crc.nd.edu/user/m/mmosby1/SVN/research/trunk/pfem3d/PFem3d/branch/include/eps.h
PGFem3D_to_VTK.o: /afs/crc.nd.edu/user/m/mmosby1/SVN/research/trunk/pfem3d/PFem3d/branch/include/volumetric_damage.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkVersion.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkCommonCoreModule.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkObject.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkObjectBase.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkIndent.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkSystemIncludes.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkWin32Header.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkConfigure.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkABI.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkIOStream.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkType.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkOStreamWrapper.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkOStrStreamWrapper.h
PGFem3D_to_VTK.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
PGFem3D_to_VTK.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
PGFem3D_to_VTK.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
PGFem3D_to_VTK.o: /usr/include/sys/types.h /usr/include/time.h
PGFem3D_to_VTK.o: /usr/include/sys/select.h /usr/include/bits/select.h
PGFem3D_to_VTK.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
PGFem3D_to_VTK.o: /usr/include/sys/sysmacros.h
PGFem3D_to_VTK.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
PGFem3D_to_VTK.o: /usr/include/string.h /usr/include/xlocale.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkSetGet.h
PGFem3D_to_VTK.o: /usr/include/math.h /usr/include/bits/huge_val.h
PGFem3D_to_VTK.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
PGFem3D_to_VTK.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
PGFem3D_to_VTK.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkTimeStamp.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkWeakPointerBase.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkVersionMacros.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkSmartPointer.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkSmartPointerBase.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkXMLUnstructuredGridReader.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkIOXMLModule.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkIOGeometryModule.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkXMLUnstructuredDataReader.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkXMLDataReader.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkXMLReader.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkAlgorithm.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkCommonExecutionModelModule.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkXMLWriter.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkXMLUnstructuredGridWriter.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkXMLUnstructuredDataWriter.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkUnstructuredGrid.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkCommonDataModelModule.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkPointSet.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkDataSet.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkDataObject.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkPoints.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkDataArray.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkAbstractArray.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkVariant.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkStdString.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkUnicodeString.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkVariantInlineOperators.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkDoubleArray.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkDataArrayTemplate.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkIdTypeArray.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkCellArray.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkCell.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkIdList.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkCellType.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkPointData.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkDataSetAttributes.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkFieldData.h
PGFem3D_to_VTK.o: /opt/crc/vtk/5.10.1/gcc/include/vtk-6.0/vtkCellData.h
PGFem3D_to_VTK.o: /afs/crc.nd.edu/user/m/mmosby1/SVN/research/trunk/pfem3d/PFem3d/branch/include/enumerations.h
PGFem3D_to_VTK.o: /afs/crc.nd.edu/user/m/mmosby1/SVN/research/trunk/pfem3d/PFem3d/branch/include/interface_macro.h
