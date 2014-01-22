#include "Printing.h"

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

/* Titles */
void PrintTitleV1()
{
  PGFEM_printf ("\n");
  PGFEM_printf ("******************************\n");
  PGFEM_printf ("* Program PFEM3d ver. 1.0    *\n");
  PGFEM_printf ("* Karel Matous           (c) *\n");
  PGFEM_printf ("******************************\n");
  PGFEM_printf ("\n");
}

void PrintTitleV2()
{
  PGFEM_printf ("\n");
  PGFEM_printf ("****************************\n");
  PGFEM_printf ("* Program PFEM3d ver. 2.0  *\n");
  PGFEM_printf ("* Karel Matous         (c) *\n");
  PGFEM_printf ("****************************\n");
  PGFEM_printf ("\n");
  PGFEM_printf ("PFEM3d [input file] [output file] [OPTIONS]\n");
  PGFEM_printf ("\n");
  PGFEM_printf ("-cp  : Finite Strain Crystal Plasticity\n");
  PGFEM_printf ("-plc : Portevin-Le Chatelier effect\n");
  PGFEM_printf ("-fd   : Finite Strains elasticity, Venant-Kirchhoff material\n");
  PGFEM_printf ("-st # : Stabilized Finite Strains Analysis, P1/P1 elements; #stab. par.\n");
  PGFEM_printf ("-stc #: STB with Cohesive Damage Zones, Mooney-Rivlin material\n");
  PGFEM_printf ("-X   : Graphics Output\n");
  PGFEM_printf ("\n");
  PGFEM_printf ("* Others *\n");
  PGFEM_printf ("-sm  : Stress smoothing\n");
  PGFEM_printf ("-pr  : Periodic boundary conditions\n");
  PGFEM_printf ("\n");
}


/* File does not exist error messages */
void NoFileProc(char* filename, int proc)
{
  PGFEM_printf("Input file is not possible to open on processor [%d]\n",proc);
  PGFEM_printf("Input file is : %s\n",filename);
  PGFEM_printf("Check the input file and run program again\n");
}

void NoFile(char* filename)
{
  PGFEM_printf("Input file [%s] is not possible to open\n",filename);
  PGFEM_printf("Check the input file and run program again\n");
}

void NoFileOutProc(char* filename, int proc)
{
  PGFEM_printf("Output file is not possible to open on processor [%d]\n",proc);
  PGFEM_printf("Output file is : %s\n",filename);
  PGFEM_printf("Check the output file and run program again\n");
}

void NoFileOut(char* filename)
{
  PGFEM_printf("Output file [%s] is not possible to open\n",filename);
  PGFEM_printf("Check the output file and run program again\n");
}
