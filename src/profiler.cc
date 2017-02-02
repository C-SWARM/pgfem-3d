#include "profiler.h"

/* Macro to print both the timer name to a file and call TimerPrint
   function with proper arguments */
#define PRINT_TIMER(timer,file)   PGFEM_fprintf(file,"%% %s\n",#timer); \
                                  TimerPrint(&(timer),file)

inline void TimerInit(Timer *timer)
{
  timer->start = timer->end = timer->total = 0;
}

void TimerUpdate(Timer *timer)
{
  timer->total += timer->end - timer->start;
  timer->end = timer->start = 0;
}

void TimerPrint(Timer *timer, FILE *file)
{
  PGFEM_fprintf(file,"%.5f\n",timer->total);
}

void ProfilerInit(Profiler *times)
{
  TimerInit(&(times->Total));
  TimerInit(&(times->Setup));

  /* in.c */
  TimerInit(&(times->read_nodes));
  TimerInit(&(times->read_supports));
  TimerInit(&(times->read_elem));
  TimerInit(&(times->read_material));
  TimerInit(&(times->read_matgeom));
  TimerInit(&(times->read_nodal_load));
  TimerInit(&(times->read_elem_surface_load));

  /* incl.c */
  TimerInit(&(times->build_matgeom));
  TimerInit(&(times->build_node));
  TimerInit(&(times->build_elem));
  TimerInit(&(times->build_hommat));
  TimerInit(&(times->build_sig_el));
  TimerInit(&(times->build_eps_el));
  TimerInit(&(times->build_zatnode));
  TimerInit(&(times->build_zatelem));
  TimerInit(&(times->build_sig_il));
  TimerInit(&(times->build_eps_il));
  TimerInit(&(times->build_elem_inelas));
  TimerInit(&(times->build_pressure_nodes));
  TimerInit(&(times->build_crystal_plast));

  TimerInit(&(times->Build));
  TimerInit(&(times->METIS));
  TimerInit(&(times->Switch));
  TimerInit(&(times->Rebuild));
  TimerInit(&(times->NR));
  TimerInit(&(times->Assembly));
  TimerInit(&(times->HYPRE));
  TimerInit(&(times->Output));
  TimerInit(&(times->Stiffness));
  TimerInit(&(times->LocAssem));
}

void ProfilerPrint(Profiler *times, FILE *file)
{
  PRINT_TIMER(times->Total,file);
  PRINT_TIMER(times->Setup,file);

  /* in.c */
  PRINT_TIMER(times->read_nodes,file);
  PRINT_TIMER(times->read_supports,file);
  PRINT_TIMER(times->read_elem,file);
  PRINT_TIMER(times->read_material,file);
  PRINT_TIMER(times->read_matgeom,file);
  PRINT_TIMER(times->read_nodal_load,file);
  PRINT_TIMER(times->read_elem_surface_load,file);

  /* incl.c */
  PRINT_TIMER(times->build_matgeom,file);
  PRINT_TIMER(times->build_node,file);
  PRINT_TIMER(times->build_elem,file);
  PRINT_TIMER(times->build_hommat,file);
  PRINT_TIMER(times->build_sig_el,file);
  PRINT_TIMER(times->build_eps_el,file);
  PRINT_TIMER(times->build_zatnode,file);
  PRINT_TIMER(times->build_zatelem,file);
  PRINT_TIMER(times->build_sig_il,file);
  PRINT_TIMER(times->build_eps_il,file);
  PRINT_TIMER(times->build_elem_inelas,file);
  PRINT_TIMER(times->build_pressure_nodes,file);
  PRINT_TIMER(times->build_crystal_plast,file);
 
  PRINT_TIMER(times->Build,file);
  PRINT_TIMER(times->METIS,file);
  PRINT_TIMER(times->Switch,file);
  PRINT_TIMER(times->Rebuild,file);
  PRINT_TIMER(times->NR,file);
  PRINT_TIMER(times->Assembly,file);
  PRINT_TIMER(times->HYPRE,file);
  PRINT_TIMER(times->Output,file);
  PRINT_TIMER(times->Stiffness,file);
  PRINT_TIMER(times->LocAssem,file);
}
