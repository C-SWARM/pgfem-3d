#ifndef PROFILER_H
#define PROFILER_H

#ifndef PGFEM_IO_H
#include "PGFEM_io.h"
#endif

typedef struct TIMER{
  double start,end,total;
} Timer;

void TimerInit(Timer *timer);
void TimerUpdate(Timer *timer);
void TimerPrint(Timer *timer, FILE *file);

typedef struct PROFILER{
  Timer Total,
    Setup,
  /* in.c begin */
    read_nodes,
    read_supports,
    read_elem,
    read_material,
    read_matgeom,
    read_nodal_load,
    read_elem_surface_load,
  /* in.c end */
  /* incl.c begin */
    build_matgeom,
    build_node,
    build_elem,
    build_hommat,
    build_sig_el,
    build_eps_el,
    build_zatnode,
    build_zatelem,
    build_sig_il,
    build_eps_il,
    build_elem_inelas,
    build_pressure_nodes,
    build_crystal_plast,
  /* incl.c end */
    MPI,
    Build,
    METIS,
    Switch,
    Rebuild,
    NR,
    Assembly,
    HYPRE,
    Output,
    Stiffness,
    LocAssem;
} Profiler;

void ProfilerInit(Profiler *times);
void ProfilerPrint(Profiler *times, FILE *file);

#endif /* #ifndef _PROFILER_H_ */
