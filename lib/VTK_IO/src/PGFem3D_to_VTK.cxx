#include "PGFem3D_to_VTK.hpp"
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <iostream>

/** Create a vtkUnstructuredGrid object from a PGFem3D mesh and
    associated information. ***CURRENTLY ONLY SUPPORT LINEAR TETRAS
    FOR DAMAGE*** */
void* PGFem3D_to_vtkUnstructuredGrid(const int nnode,
				     const int nelems,
				     const NODE *nodes,
				     const ELEMENT *elems,
				     const SUPP supports,
				     const SIG *stress,
				     const EPS *strain,
				     const double *dofs,
				     const int analysis_type)
{
  int have_pressure = 0;
  // Create Points and point data
  vtkSmartPointer<vtkPoints> points = 
    vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkDoubleArray> disp =
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> macro_disp =
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> pressure =
    vtkSmartPointer<vtkDoubleArray>::New();

  // Allocate points and point data
  points->SetNumberOfPoints(nnode);
  disp->SetNumberOfComponents(3);
  disp->SetNumberOfTuples(nnode);
  disp->SetName("Displacement");
  if(supports->multi_scale){
    macro_disp->SetNumberOfComponents(3);
    macro_disp->SetNumberOfTuples(nnode);
    macro_disp->SetName("MacroDisplacement");
    // compute macro gradient
    // double *jump = new double[3];
    // compute_interface_macro_jump_u(jump,supports);
    // compute_interface_macro_grad_u(supports->F0,supports->lc,
    // 				   jump,supports->N0);
    // delete[] jump;

    compute_macro_grad_u(supports->F0,supports,analysis_type);
  }
  if(analysis_type == STABILIZED
     || analysis_type == MINI
     || analysis_type == MINI_3F){
    have_pressure = 1;
    pressure->SetNumberOfComponents(1);
    pressure->SetNumberOfTuples(nnode);
    pressure->SetName("Pressure");
  }
  for(int i=0; i<nnode; i++){
    double displ[3] = {0.0,0.0,0.0};
    for(int j=0; j<3; j++){
      if(nodes[i].id[j] == 0){
	displ[j] = 0.0;
      }else if(nodes[i].id[j] > 0){
	displ[j] = dofs[nodes[i].id[j]-1];
      } else if(nodes[i].id[j] <  0){
	displ[j] = supports->defl[abs(nodes[i].id[j])-1];
      }
    }
    points->SetPoint(i,nodes[i].x1_fd,nodes[i].x2_fd,nodes[i].x3_fd);
    disp->SetTuple3(i,displ[0],displ[1],displ[2]);
    if(supports->multi_scale){
      compute_interface_macro_disp_at_node(displ,&nodes[i],supports->F0,analysis_type);
      macro_disp->SetTuple3(i,displ[0],displ[1],displ[2]);
    }
    if(have_pressure){
      pressure->SetValue(i,dofs[nodes[i].id[3]-1]);
    }
  }

  // Cell connectivity and data
  vtkSmartPointer<vtkCellArray> cells = 
    vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkIdTypeArray> conn = 
    vtkSmartPointer<vtkIdTypeArray>::New();
  int *type;
  vtkSmartPointer<vtkDoubleArray> vtkStress = 
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> vtkStrain = 
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> effStress = 
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> effStrain = 
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> damage = 
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> chi = 
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIdTypeArray> prop = 
    vtkSmartPointer<vtkIdTypeArray>::New();

  /* output F, P, and W at 1st integration point */
  vtkSmartPointer<vtkDoubleArray> vtkF = 
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> vtkP = 
    vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> vtkW = 
    vtkSmartPointer<vtkDoubleArray>::New();

  // pre-allocate assuming quadradic tets
  conn->SetNumberOfValues(11*nelems); 
  conn->Reset(); // Set insertion point to begining and "size" to zero

  type = new int[nelems];

  vtkStress->SetNumberOfComponents(6);
  vtkStress->SetNumberOfTuples(nelems);
  vtkStress->SetName("CauchyStress");

  vtkStrain->SetNumberOfComponents(6);
  vtkStrain->SetNumberOfTuples(nelems);
  vtkStrain->SetName("EulerStrain");

  effStress->SetNumberOfValues(nelems);
  effStress->SetName("EffectiveStress");

  effStrain->SetNumberOfValues(nelems);
  effStrain->SetName("EffectiveStrain");

  damage->SetNumberOfValues(nelems);
  damage->SetName("Damage");

  chi->SetNumberOfValues(nelems);
  chi->SetName("Chi");

  prop->SetNumberOfValues(nelems);
  prop->SetName("CellProperty");

  vtkF->SetNumberOfComponents(9);
  vtkF->SetNumberOfTuples(nelems);
  vtkF->SetName("F");

  vtkP->SetNumberOfComponents(9);
  vtkP->SetNumberOfTuples(nelems);
  vtkP->SetName("P");

  vtkW->SetNumberOfValues(nelems);
  vtkW->SetName("W");

  // build connectivity and populate fields
  unsigned int count = 0;
  for(int i=0; i<nelems; i++){
    switch(elems[i].toe){
    case 4: type[i] = VTK_TETRA; break;
    case 8: type[i] = VTK_HEXAHEDRON; break;
    case 10: type[i] = VTK_QUADRATIC_TETRA; break;
    default:
      std::cerr << "Unsupported element type!\n" << std::endl;
      MPI_Abort(MPI_COMM_WORLD,0);
      break;
    }
    // conn->InsertValue(count,elems[i].toe); count++;
    conn->InsertNextValue(elems[i].toe);
    for(int j=0;j<elems[i].toe; j++){
      // conn->InsertValue(count,elems[i].nod[j]);count++;
      conn->InsertNextValue(elems[i].nod[j]);
    }

    // swap stress and strain values
    double tmp_stress[9],tmp_strain[9];
    tmp_stress[0] = stress[i].el.o[0];
    tmp_stress[1] = stress[i].el.o[1];
    tmp_stress[2] = stress[i].el.o[2];
    tmp_stress[3] = stress[i].el.o[5];
    tmp_stress[4] = stress[i].el.o[3];
    tmp_stress[5] = stress[i].el.o[4];

    tmp_strain[0] = strain[i].el.o[0];
    tmp_strain[1] = strain[i].el.o[1];
    tmp_strain[2] = strain[i].el.o[2];
    tmp_strain[3] = strain[i].el.o[5]/2;
    tmp_strain[4] = strain[i].el.o[3]/2;
    tmp_strain[5] = strain[i].el.o[4]/2;

    vtkStress->SetTuple(i,tmp_stress);
    vtkStrain->SetTuple(i,tmp_strain);
    effStress->SetValue(i,stress[i].el.eq);
    effStrain->SetValue(i,strain[i].el.eq);
    damage->SetValue(i,strain[i].dam[0].wn);
    chi->SetValue(i,strain[i].dam[0].Xn);
    prop->SetValue(i,elems[i].pr);

    vtkF->SetTuple(i,strain[i].il[0].F);
    vtkW->SetValue(i,strain[i].il[0].Y);

    // compute P from FS
    double *pF = strain[i].il[0].F;
    tmp_stress[0] = stress[i].il[0].o[0];
    tmp_stress[1] = stress[i].il[0].o[5];
    tmp_stress[2] = stress[i].il[0].o[4];

    tmp_stress[3] = stress[i].il[0].o[5];
    tmp_stress[4] = stress[i].il[0].o[1];
    tmp_stress[5] = stress[i].il[0].o[3];

    tmp_stress[6] = stress[i].il[0].o[3];
    tmp_stress[7] = stress[i].il[0].o[4];
    tmp_stress[8] = stress[i].il[0].o[2];

    // compute tmp_strain (I know...) = FS
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
	int i_jk = j*3+k;
	tmp_strain[i_jk] = 0.0;
	for(int l=0; l<3; l++){
	  int i_jl = j*3+l;
	  int i_lk = l*3+k;
	  tmp_strain[i_jk] += pF[i_jl]*tmp_stress[i_lk];
	}
      }
    }
    vtkP->SetTuple(i,tmp_strain);
  }

  conn->Squeeze();
  cells->SetCells(nelems,conn);

  // Create Grid
  vtkSmartPointer<vtkUnstructuredGrid> grid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->SetPoints(points);
  grid->SetCells(type,cells);

  // POINT DATA
  grid->GetPointData()->AddArray(disp);
  if(supports->multi_scale){
    grid->GetPointData()->AddArray(macro_disp);
  }
  if(have_pressure){
    grid->GetPointData()->AddArray(pressure);
  }

  // CELL DATA
  grid->GetCellData()->AddArray(vtkStress);
  grid->GetCellData()->AddArray(vtkStrain);
  grid->GetCellData()->AddArray(effStrain);
  grid->GetCellData()->AddArray(effStress);
  grid->GetCellData()->AddArray(prop);
  grid->GetCellData()->AddArray(damage);
  grid->GetCellData()->AddArray(chi);

  grid->GetCellData()->AddArray(vtkF);
  grid->GetCellData()->AddArray(vtkP);
  grid->GetCellData()->AddArray(vtkW);

  grid->Register(grid); // increment ref counter

  delete[] type;
  return (void*) grid;
}

void PGFem3D_destroy_vtkUnstructuredGrid(void *grid)
{
  vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = 
    (vtkUnstructuredGrid*) grid;
  vtkGrid->Delete();
}

void PGFem3D_write_vtkUnstructuredGrid(const char* filename,
				       const void *grid,
				       const int ascii)
{
  vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = 
    (vtkUnstructuredGrid*) grid;
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = 
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(filename);

  if(ascii){
    writer->SetDataModeToAscii();
  }

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(vtkGrid);
#else
  writer->SetInputData(vtkGrid);
#endif

  writer->Write();
}
