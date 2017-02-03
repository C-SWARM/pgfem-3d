#include "PGFem3D_to_VTK.hpp"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#ifndef NO_VTK_LIB
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
     || analysis_type == MINI_3F
     || (analysis_type == TF && elems[0].toe!=10)){
    have_pressure = 1;
    pressure->SetNumberOfComponents(1);
    pressure->SetNumberOfTuples(nnode);
    pressure->SetName("Pressure");
  }
  for(int i=0; i<nnode; i++){
    double displ[3] = {0.0,0.0,0.0};
    for(int j=0; j<3; j++){
      if(nodes[i].id_map[0].id[j] == 0){
	displ[j] = 0.0;
      }else if(nodes[i].id_map[0].id[j] > 0){
	displ[j] = dofs[nodes[i].id_map[0].id[j]-1];
      } else if(nodes[i].id_map[0].id[j] <  0){
	displ[j] = supports->defl[abs(nodes[i].id_map[0].id[j])-1];
      }
    }
    points->SetPoint(i,nodes[i].x1_fd,nodes[i].x2_fd,nodes[i].x3_fd);
    disp->SetTuple3(i,displ[0],displ[1],displ[2]);
    if(supports->multi_scale){
      compute_interface_macro_disp_at_node(displ,&nodes[i],supports->F0,analysis_type);
      macro_disp->SetTuple3(i,displ[0],displ[1],displ[2]);
    }
    if(have_pressure){
      pressure->SetValue(i,dofs[nodes[i].id_map[0].id[3]-1]);
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
      
  vtkSmartPointer<vtkDoubleArray> tf_pressure = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> tf_volume = vtkSmartPointer<vtkDoubleArray>::New();              

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

  if(analysis_type==TF)
  { 
    if(elems[0].toe==10)
    {
      tf_pressure->SetNumberOfValues(nelems);
      tf_pressure->SetName("TF_Pressure"); 
    }
    
    tf_volume->SetNumberOfValues(nelems);
    tf_volume->SetName("TF_Volume"); 
  }  

  // build connectivity and populate fields
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
    conn->InsertNextValue(elems[i].toe);
    for(int j=0;j<elems[i].toe; j++){
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
    
    if(analysis_type==TF)
    {
      if(elems[i].toe==10)
        tf_pressure->SetValue(i, strain[i].d_T[0]);  
      
      tf_volume->SetValue(i, strain[i].T[0]);  
    }
      

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
  if(analysis_type==TF)
  {
    if(elems[0].toe==10)
      grid->GetCellData()->AddArray(tf_pressure);
      
    grid->GetCellData()->AddArray(tf_volume);    
  }
  
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

int read_VTK_file4TF(char fn[], double *r, double *P, double *V)
{
  //parse command line arguments

  std::string filename = fn;


  //************ read all the data from the file **************
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update(); //read in data from file.

  vtkSmartPointer<vtkUnstructuredGrid> grid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid = reader->GetOutput(); //data is now in unstructured grid.

  //*************************************************************

  double *temp = 0;
  int nelems = 0;

  vtkCellData *cData = grid->GetCellData(); // store the grid's cell data
  vtkPointData *pData = grid->GetPointData();

  //************* Example: If(data exists) {}
  /* vtkDataArray *CauchyStress;
  if(cData->HasArray("CauchyStress")) {
    CauchyStress = cData->GetArray("CauchyStress");
  }
 */
 
 
 
   //********************* Print out displacement data ***********
  //  vtkPointData *pData = grid->GetPointData(); --This is already initialized above
  if(pData->HasArray("Displacement")) 
  {
    vtkDataArray *Displacement = pData->GetArray("Displacement");
    nelems = Displacement->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) 
    {
      temp = Displacement->GetTuple(i);
      r[i*3+0] = temp[0];
      r[i*3+1] = temp[1];      
      r[i*3+2] = temp[2];      
    }
  }
  //**************************************************************

  if(cData->HasArray("TF_Pressure")) 
  {
    vtkDataArray *pressure = cData->GetArray("TF_Pressure");
    nelems = pressure->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++)
    {
      temp = pressure->GetTuple(i);
      P[i] = temp[0];
    } 
  }
  
  if(cData->HasArray("TF_Volume")) 
  {
    vtkDataArray *volume = cData->GetArray("TF_Volume");
    nelems = volume->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++)
    {
      temp = volume->GetTuple(i);
      V[i] = temp[0];
    } 
  }

  return EXIT_SUCCESS;
}


int read_VTK_file(char fn[], double *r)
{
  //parse command line arguments

  std::string filename = fn;


  //************ read all the data from the file **************
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update(); //read in data from file.

  vtkSmartPointer<vtkUnstructuredGrid> grid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid = reader->GetOutput(); //data is now in unstructured grid.

  //*************************************************************

  double *temp = 0;
  int nelems = 0;

  vtkCellData *cData = grid->GetCellData(); // store the grid's cell data
  vtkPointData *pData = grid->GetPointData();
  vtkDataArray *vtkData; // used to temporarily store vtk data arrays. CauchyStress, etc.

  //************* Example: If(data exists) {}
  /* vtkDataArray *CauchyStress;
  if(cData->HasArray("CauchyStress")) {
    CauchyStress = cData->GetArray("CauchyStress");
  }
 */
 
 
 
   //********************* Print out displacement data ***********
  //  vtkPointData *pData = grid->GetPointData(); --This is already initialized above
  if(pData->HasArray("Displacement")) 
  {
    vtkDataArray *Displacement = pData->GetArray("Displacement");
    nelems = Displacement->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) 
    {
      temp = Displacement->GetTuple(i);
      r[i*3+0] = temp[0];
      r[i*3+1] = temp[1];      
      r[i*3+2] = temp[2];      
    }
  }
  //**************************************************************
  
  return EXIT_SUCCESS;
 
 
 
 
  if(cData->HasArray("CauchyStress")) {
    vtkData = cData->GetArray("CauchyStress");
    nelems = vtkData->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) {
      //******************** Print out CauchyStress data ************
      temp = vtkData->GetTuple(i);
      std::cout << "----------CauchyStress: " << i << " ---------- " << std::endl;
      for(int j = 0; j < 6; j++) {
        std::cout <<  "Component " << j << ": " << temp[j] << std::endl;
      }
      std::cout << std::endl;
      //*************************************************************
     } 
  }

  if(cData->HasArray("EulerStrain")) {
    vtkData = cData->GetArray("EulerStrain");
    nelems = vtkData->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) {
      //******************** Print out EulerStrain data ************
      temp = vtkData->GetTuple(i);
      std::cout << "----------EulerStrain: " << i << " ---------- " << std::endl;
      for(int j = 0; j < 6; j++) {
        std::cout <<  "Component " << j << ": " << temp[j] << std::endl;
      }
      std::cout << std::endl;
      //*************************************************************
     } 
  }

  if(cData->HasArray("EffectiveStress")) {
    vtkData = cData->GetArray("EffectiveStress");
    nelems = vtkData->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) {
      //******************** Print out EffectiveStress data ************
      temp = vtkData->GetTuple(i);
      std::cout << "EffectiveStress:" << i << " " << *temp << std::endl;
      //*************************************************************
     } 
  }

  if(cData->HasArray("EffectiveStrain")) {
    vtkData = cData->GetArray("EffectiveStrain");
    nelems = vtkData->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) {
      //******************** Print out EffectiveStrain data ************
      temp = vtkData->GetTuple(i);
      std::cout << "EffectiveStrain:" << i << " " << *temp << std::endl;
      //*************************************************************
     } 
  }

  if(cData->HasArray("Damage")) {
    vtkData = cData->GetArray("Damage");
    nelems = vtkData->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) {
      //******************** Print out Damage data ************
      temp = vtkData->GetTuple(i);
      std::cout << "Damage:" << i << " " << *temp << std::endl;
      //*************************************************************
     } 
  }
  
  if(cData->HasArray("Chi")) {
    vtkData = cData->GetArray("Chi");
    nelems = vtkData->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) {
      //******************** Print out Chi data ************
      temp = vtkData->GetTuple(i);
      std::cout << "Chi:" << i << " " << *temp << std::endl;
      //*************************************************************
     } 
  }

  if(cData->HasArray("CellProperty")) {
    vtkData = cData->GetArray("CellProperty");
    nelems = vtkData->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) {
      //******************** Print out CellProperty data ************
      temp = vtkData->GetTuple(i);
      std::cout << "CellProperty:" << i << " " << *temp << std::endl;
      //*************************************************************
     } 
  }

  if(cData->HasArray("F")) {
    vtkData = cData->GetArray("F");
    nelems = vtkData->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) {
      //******************** Print out F data ************
      temp = vtkData->GetTuple(i);
      std::cout << "----------F: " << i << " ---------- " << std::endl;
      for(int j = 0; j < 9; j++) {
        std::cout <<  "Component " << j << ": " << temp[j] << std::endl;
      }
      std::cout << std::endl;
      //*************************************************************
     } 
  }

  if(cData->HasArray("P")) {
    vtkData = cData->GetArray("P");
    nelems = vtkData->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) {
      //******************** Print out P data ************
      temp = vtkData->GetTuple(i);
      std::cout << "----------P: " << i << " ---------- " << std::endl;
      for(int j = 0; j < 9; j++) {
        std::cout <<  "Component " << j << ": " << temp[j] << std::endl;
      }
      std::cout << std::endl;
      //*************************************************************
     } 
  }

  if(cData->HasArray("W")) {
    vtkData = cData->GetArray("W");
    nelems = vtkData->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) {
      //******************** Print out W data ************
      temp = vtkData->GetTuple(i);
      std::cout << "W:" << i << " " << *temp << std::endl;
      //*************************************************************
     } 
  }


  //********************* Print out displacement data ***********
  //  vtkPointData *pData = grid->GetPointData(); --This is already initialized above
  if(pData->HasArray("Displacement")) {
    std::cout << "Displacement: " << std::endl;
    vtkDataArray *Displacement = pData->GetArray("Displacement");
    nelems = Displacement->GetNumberOfTuples();
    for(int i = 0; i < nelems; i++) {
      temp = Displacement->GetTuple(i);
      std::cout << "Tuple: " << i << std::endl << "Component1: " << temp[0] << std::endl << "Component2: " << temp[1] << std::endl << "Component3: " << temp[2] << std::endl << std::endl;
    }
  }
  //**************************************************************

  return EXIT_SUCCESS;
}

#else

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
  fprintf(stderr,"WARNING: calling VTK_IO function w/ no VTK support! (%s)",__func__);
  return NULL;
}

void PGFem3D_destroy_vtkUnstructuredGrid(void *grid)
{
  fprintf(stderr,"WARNING: calling VTK_IO function w/ no VTK support! (%s)",__func__);
}

void PGFem3D_write_vtkUnstructuredGrid(const char* filename,
				       const void *grid,
				       const int ascii)
{
  fprintf(stderr,"WARNING: calling VTK_IO function w/ no VTK support! (%s)",__func__);
}
				       
int read_VTK_file(char fn[], double *r)
{
  fprintf(stderr,"WARNING: calling VTK_IO function w/ no VTK support! (%s)",__func__);
  return EXIT_FAILURE;
}

#endif
