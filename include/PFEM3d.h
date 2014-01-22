/* HEADER */

#ifndef PFEM3d_h
#define PFEM3d_h

// Define all debug variables etc
#ifndef PFEM_DEBUG
#define PFEM_DEBUG       0
#endif

#ifndef PFEM_NEW_GRN
#define	PFEM_NEW_GRN     0
#endif

#ifndef PFEM_DEBUG_GRN
#define	PFEM_DEBUG_GRN   0
#endif

#ifndef PFEM_TRACE
#define	PFEM_TRACE       0
#endif

#ifndef PFEM_DEBUG_ALL
#define	PFEM_DEBUG_ALL   0
#endif

#ifndef PFEM_DEBUG_STIFFNESS
#define	PFEM_DEBUG_STIFFNESS 0
#endif

#ifndef PFEM_PRINT
#define	PFEM_PRINT       0
#endif

#ifndef PFEM_PRINT_MATRIX
#define	PFEM_PRINT_MATRIX       0
#endif

/* debugging print macro */
#if (PFEM_DEBUG)
#define DP(string) printf("[%d] %s, %s\n",myrank,__func__,(string))
#else
#define DP(string) {}
#endif

/* variables */
//extern double PII,VVolume;
/* extern int myrank,nproc; */
//extern long analysis,periodic,plc,coh,fil;

#endif /* HEADER */

//static const double PII = 3.141592653589793238462643;
