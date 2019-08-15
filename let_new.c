/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: let_new.instr (LET)
 * Date:       Tue Sep 18 00:58:23 2018
 * File:       let_new.c
 * Compile:    cc -o LET.exe let_new.c  -I@MCCODE_LIB@/libs/mcpl -L@MCCODE_LIB@/libs/mcpl -lmcpl
 * CFLAGS= -I@MCCODE_LIB@/libs/mcpl -L@MCCODE_LIB@/libs/mcpl -lmcpl
 */


#define MCCODE_STRING "McStas 2.4.1 - Jun. 26, 2017"
#define FLAVOR "mcstas"
#define FLAVOR_UPPER "MCSTAS"
#define MC_USE_DEFAULT_MAIN
#define MC_TRACE_ENABLED
#define MC_EMBEDDED_RUNTIME

#line 1 "mccode-r.h"
/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas 2.4.1
* Version: $Revision$
*
* Runtime system header for McStas/McXtrace.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas/McXtrace version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCCODE_R_H
#define MCCODE_R_H "$Revision$"

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <float.h>
#include <inttypes.h>

/* If the runtime is embedded in the simulation program, some definitions can
   be made static. */

#ifdef MC_EMBEDDED_RUNTIME
#define mcstatic static
#else
#define mcstatic
#endif

#ifdef __dest_os
#if (__dest_os == __mac_os)
#define MAC
#endif
#endif

#ifdef __FreeBSD__
#define NEED_STAT_H
#endif

#if defined(__APPLE__) && defined(__GNUC__)
#define NEED_STAT_H
#endif

#ifdef NEED_STAT_H
#include <sys/stat.h>
#endif

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !WIN32 */
#endif /* MC_PATHSEP_C */



/* the version string is replaced when building distribution with mkdist */
#ifndef MCCODE_STRING
#define MCCODE_STRING "McStas 2.4.1 - Jun. 26, 2017"
#endif

#ifndef MCCODE_DATE
#define MCCODE_DATE "Jun. 26, 2017"
#endif

#ifndef MCCODE_VERSION
#define MCCODE_VERSION "2.4.1"
#endif

#ifndef MCCODE_NAME
#define MCCODE_NAME "McStas"
#endif

#ifndef MCCODE_PARTICLE
#define MCCODE_PARTICLE "neutron"
#endif

#ifndef MCCODE_LIBENV
#define MCCODE_LIBENV "MCSTAS"
#endif

#ifndef FLAVOR_UPPER
#define FLAVOR_UPPER MCCODE_NAME
#endif

#ifdef MC_PORTABLE
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#ifdef MAC
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (USE_MPI == 0)
#undef USE_MPI
#endif

#ifdef USE_MPI  /* default is to disable signals with MPI, as MPICH uses them to communicate */
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (NOSIGNALS == 0)
#undef NOSIGNALS
#endif

/* Note: the enum instr_formal_types definition MUST be kept
   synchronized with the one in mccode.h and with the
   instr_formal_type_names array in cogen.c. */
enum instr_formal_types
  {
    instr_type_double, instr_type_int, instr_type_string
  };
struct mcinputtable_struct { /* defines instrument parameters */
  char *name; /* name of parameter */
  void *par;  /* pointer to instrument parameter (variable) */
  enum instr_formal_types type;
  char *val;  /* default value */
};

typedef double MCNUM;
typedef struct {MCNUM x, y, z;} Coords;
typedef MCNUM Rotation[3][3];

/* the following variables are defined in the McStas generated C code
   but should be defined externally in case of independent library usage */
#ifndef DANSE
extern struct mcinputtable_struct mcinputtable[]; /* list of instrument parameters */
extern int    mcnumipar;                          /* number of instrument parameters */
extern char   mcinstrument_name[], mcinstrument_source[]; /* instrument name and filename */
extern char  *mcinstrument_exe;                           /* executable path = argv[0] or NULL */
extern MCNUM  mccomp_storein[]; /* 11 coords * number of components in instrument */
extern MCNUM  mcAbsorbProp[];
extern MCNUM  mcScattered;      /* number of SCATTER calls in current component */
extern MCNUM  mcRestore;        /* Flag to indicate if neutron needs to be restored */
#ifndef MC_ANCIENT_COMPATIBILITY
extern int mctraceenabled, mcdefaultmain;
#endif
#endif


/* Useful macros ============================================================ */

/* MPI stuff */

#ifdef USE_MPI
#include "mpi.h"

#ifdef OMPI_MPI_H  /* openmpi does not use signals: we may install our sighandler */
#undef NOSIGNALS
#endif

/*
 * MPI_MASTER(i):
 * execution of i only on master node
 */
#define MPI_MASTER(statement) { \
  if(mpi_node_rank == mpi_node_root)\
  { statement; } \
}

#ifndef MPI_REDUCE_BLOCKSIZE
#define MPI_REDUCE_BLOCKSIZE 1000
#endif

int mc_MPI_Sum(double* buf, long count);
int mc_MPI_Send(void *sbuf, long count, MPI_Datatype dtype, int dest);
int mc_MPI_Recv(void *rbuf, long count, MPI_Datatype dtype, int source);

/* MPI_Finalize exits gracefully and should be preferred to MPI_Abort */
#define exit(code) do {                                   \
    MPI_Finalize();                                       \
    exit(code);                                           \
  } while(0)

#else /* !USE_MPI */
#define MPI_MASTER(instr) instr
#endif /* USE_MPI */

#ifdef USE_MPI
static int mpi_node_count;
#endif

#ifdef USE_THREADS  /* user want threads */
#error Threading (USE_THREADS) support has been removed for very poor efficiency. Use MPI/SSH grid instead.
#endif


void   mcset_ncount(unsigned long long count);    /* wrapper to get mcncount */
unsigned long long int mcget_ncount(void);            /* wrapper to set mcncount */
unsigned long long mcget_run_num(void);           /* wrapper to get mcrun_num=0:mcncount */


/* Following part is only embedded when not redundant with mccode.h ========= */

#ifndef MCCODE_H

#ifndef NOSIGNALS
#include <signal.h>
#define SIG_MESSAGE(msg) strcpy(mcsig_message, msg);
#else
#define SIG_MESSAGE(msg)
#endif /* !NOSIGNALS */

/* Useful macros and constants ============================================== */

#ifndef FLT_MAX
#define FLT_MAX         3.40282347E+38F /* max decimal value of a "float" */
#endif

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef SQR
#define SQR(x) ( (x) * (x) )
#endif
#ifndef SIGN
#define SIGN(x) (((x)>0.0)?(1):(-1))
#endif

#ifndef PI
# ifdef M_PI
#  define PI M_PI
# else
#  define PI 3.14159265358979323846
# endif
#endif

#define RAD2MIN  ((180*60)/PI)
#define MIN2RAD  (PI/(180*60))
#define DEG2RAD  (PI/180)
#define RAD2DEG  (180/PI)
#define FWHM2RMS 0.424660900144    /* Convert between full-width-half-max and */
#define RMS2FWHM 2.35482004503     /* root-mean-square (standard deviation) */
#define HBAR     1.05457168e-34    /* [Js] h bar Planck constant CODATA 2002 */
#define MNEUTRON 1.67492728e-27    /* [kg] mass of neutron CODATA 2002 */
#define GRAVITY  9.81              /* [m/s^2] gravitational acceleration */
#define NA       6.02214179e23     /* [#atoms/g .mole] Avogadro's number*/


/* wrapper to get absolute and relative position of comp */
/* mccomp_posa and mccomp_posr are defined in McStas generated C code */
#define POS_A_COMP_INDEX(index) \
    (mccomp_posa[index])
#define POS_R_COMP_INDEX(index) \
    (mccomp_posr[index])
/* number of SCATTER calls in current comp: mcScattered defined in generated C code */
#define SCATTERED mcScattered
/* Flag to indicate if neutron needs to be restored: mcRestore defined in generated C code */
#define RESTORE mcRestore


/* Retrieve component information from the kernel */
/* Name, position and orientation (both absolute and relative)  */
/* Any component: For "redundancy", see comment by KN */
#define tmp_name_comp(comp) #comp
#define NAME_COMP(comp) tmp_name_comp(comp)
#define tmp_pos_a_comp(comp) (mcposa ## comp)
#define POS_A_COMP(comp) tmp_pos_a_comp(comp)
#define tmp_pos_r_comp(comp) (mcposr ## comp)
#define POS_R_COMP(comp) tmp_pos_r_comp(comp)
#define tmp_rot_a_comp(comp) (mcrota ## comp)
#define ROT_A_COMP(comp) tmp_rot_a_comp(comp)
#define tmp_rot_r_comp(comp) (mcrotr ## comp)
#define ROT_R_COMP(comp) tmp_rot_r_comp(comp)

/* Current component name, index, position and orientation */
#define NAME_CURRENT_COMP  NAME_COMP(mccompcurname)
#define INDEX_CURRENT_COMP mccompcurindex
#define POS_A_CURRENT_COMP POS_A_COMP(mccompcurname)
#define POS_R_CURRENT_COMP POS_R_COMP(mccompcurname)
#define ROT_A_CURRENT_COMP ROT_A_COMP(mccompcurname)
#define ROT_R_CURRENT_COMP ROT_R_COMP(mccompcurname)

/* Note: The two-stage approach to MC_GETPAR is NOT redundant; without it,
* after #define C sample, MC_GETPAR(C,x) would refer to component C, not to
* component sample. Such are the joys of ANSI C.

* Anyway the usage of MCGETPAR requires that we use sometimes bare names...
*/
#define MC_GETPAR2(comp, par) (mcc ## comp ## _ ## par)
#define MC_GETPAR(comp, par) MC_GETPAR2(comp,par)

/* MCDISPLAY/trace and debugging message sent to stdout */
#ifdef MC_TRACE_ENABLED
#define DEBUG
#endif

#ifdef DEBUG
#define mcDEBUG_INSTR() if(!mcdotrace); else { printf("INSTRUMENT:\n"); printf("Instrument '%s' (%s)\n", mcinstrument_name, mcinstrument_source); }
#define mcDEBUG_COMPONENT(name,c,t) if(!mcdotrace); else {\
  printf("COMPONENT: \"%s\"\n" \
         "POS: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         name, c.x, c.y, c.z, t[0][0], t[0][1], t[0][2], \
         t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2]); \
  mcAccumulatedILength += coords_len(coords_sub(mcLastComp,c)); \
  printf("Component %30s AT (%g,%g,%g)    %g m from origin\n", name, c.x, c.y, c.z, mcAccumulatedILength); \
  mcLastComp=c;\
  }
#define mcDEBUG_INSTR_END() if(!mcdotrace); else printf("INSTRUMENT END:\n");
#define mcDEBUG_ENTER() if(!mcdotrace); else printf("ENTER:\n");
#define mcDEBUG_COMP(c) if(!mcdotrace); else printf("COMP: \"%s\"\n", c);
#define mcDEBUG_LEAVE() if(!mcdotrace); else printf("LEAVE:\n");
#define mcDEBUG_ABSORB() if(!mcdotrace); else printf("ABSORB:\n");
#else
#define mcDEBUG_INSTR()
#define mcDEBUG_COMPONENT(name,c,t)
#define mcDEBUG_INSTR_END()
#define mcDEBUG_ENTER()
#define mcDEBUG_COMP(c)
#define mcDEBUG_LEAVE()
#define mcDEBUG_ABSORB()
#endif

// mcDEBUG_STATE and mcDEBUG_SCATTER are defined by mcstas-r.h and mcxtrace-r.h



#ifdef TEST
#define test_printf printf
#else
#define test_printf while(0) printf
#endif

/* send MCDISPLAY message to stdout to show gemoetry */
void mcdis_magnify(char *what);
void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2);
void mcdis_dashed_linemcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n);
void mcdis_multiline(int count, ...);
void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height);
void mcdis_box(double x, double y, double z,
	       double width, double height, double length);
void mcdis_circle(char *plane, double x, double y, double z, double r);

/* selection of random number generator. default is MT */
#ifndef MC_RAND_ALG
#define MC_RAND_ALG 1
#endif

#if MC_RAND_ALG == 0
   /* Use system random() (not recommended). */
#  define MC_RAND_MAX RAND_MAX
#elif MC_RAND_ALG == 1
   /* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
#  define MC_RAND_MAX ((unsigned long)0xffffffff)
#  define random mt_random
#  define srandom mt_srandom
#elif MC_RAND_ALG == 2
   /* Algorithm used in McStas CVS-080208 and earlier (not recommended). */
#  define MC_RAND_MAX 0x7fffffff
#  define random mc_random
#  define srandom mc_srandom
#else
#  error "Bad value for random number generator choice."
#endif

typedef int mc_int32_t;
mc_int32_t mc_random(void);
void mc_srandom (unsigned int x);
unsigned long mt_random(void);
void mt_srandom (unsigned long x);

double rand01();
double randpm1();
double rand0max(double max);
double randminmax(double min, double max);

double randnorm(void);
double randtriangle(void);

#ifndef DANSE
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);
#endif

/* simple vector algebra ==================================================== */
#define vec_prod(x, y, z, x1, y1, z1, x2, y2, z2) \
	vec_prod_func(&x, &y, &z, x1, y1, z1, x2, y2, z2)
mcstatic inline void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1, double x2, double y2, double z2);

mcstatic inline double scalar_prod(
		double x1, double y1, double z1, double x2, double y2, double z2);

#define NORM(x,y,z) \
	norm_func(&x, &y, &z)
mcstatic inline void norm_func(double *x, double *y, double *z) {
	double temp = (*x * *x) + (*y * *y) + (*z * *z);
	if (temp != 0) {
		temp = sqrt(temp);
		*x /= temp;
		*y /= temp;
		*z /= temp;
	}
}
#define normal_vec(nx, ny, nz, x, y, z) \
    normal_vec_func(&(nx), &(ny), &(nz), x, y, z)
mcstatic inline void normal_vec_func(double *nx, double *ny, double *nz,
    double x, double y, double z);

/**
 * Rotate the vector vx,vy,vz psi radians around the vector ax,ay,az
 * and put the result in x,y,z.
 */
#define rotate(x, y, z, vx, vy, vz, phi, ax, ay, az) \
  do { \
    double mcrt_tmpx = (ax), mcrt_tmpy = (ay), mcrt_tmpz = (az); \
    double mcrt_vp, mcrt_vpx, mcrt_vpy, mcrt_vpz; \
    double mcrt_vnx, mcrt_vny, mcrt_vnz, mcrt_vn1x, mcrt_vn1y, mcrt_vn1z; \
    double mcrt_bx, mcrt_by, mcrt_bz; \
    double mcrt_cos, mcrt_sin; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vp = scalar_prod((vx), (vy), (vz), mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vpx = mcrt_vp*mcrt_tmpx; \
    mcrt_vpy = mcrt_vp*mcrt_tmpy; \
    mcrt_vpz = mcrt_vp*mcrt_tmpz; \
    mcrt_vnx = (vx) - mcrt_vpx; \
    mcrt_vny = (vy) - mcrt_vpy; \
    mcrt_vnz = (vz) - mcrt_vpz; \
    vec_prod(mcrt_bx, mcrt_by, mcrt_bz, \
             mcrt_tmpx, mcrt_tmpy, mcrt_tmpz, mcrt_vnx, mcrt_vny, mcrt_vnz); \
    mcrt_cos = cos((phi)); mcrt_sin = sin((phi)); \
    mcrt_vn1x = mcrt_vnx*mcrt_cos + mcrt_bx*mcrt_sin; \
    mcrt_vn1y = mcrt_vny*mcrt_cos + mcrt_by*mcrt_sin; \
    mcrt_vn1z = mcrt_vnz*mcrt_cos + mcrt_bz*mcrt_sin; \
    (x) = mcrt_vpx + mcrt_vn1x; \
    (y) = mcrt_vpy + mcrt_vn1y; \
    (z) = mcrt_vpz + mcrt_vn1z; \
  } while(0)

/**
 * Mirror (xyz) in the plane given by the point (rx,ry,rz) and normal (nx,ny,nz)
 *
 * TODO: This define is seemingly never used...
 */
#define mirror(x,y,z,rx,ry,rz,nx,ny,nz) \
  do { \
    double mcrt_tmpx= (nx), mcrt_tmpy = (ny), mcrt_tmpz = (nz); \
    double mcrt_tmpt; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_tmpt=scalar_prod((rx),(ry),(rz),mcrt_tmpx,mcrt_tmpy,mcrt_tmpz); \
    (x) = rx -2 * mcrt_tmpt*mcrt_rmpx; \
    (y) = ry -2 * mcrt_tmpt*mcrt_rmpy; \
    (z) = rz -2 * mcrt_tmpt*mcrt_rmpz; \
  } while (0)

Coords coords_set(MCNUM x, MCNUM y, MCNUM z);
Coords coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z);
Coords coords_add(Coords a, Coords b);
Coords coords_sub(Coords a, Coords b);
Coords coords_neg(Coords a);
Coords coords_scale(Coords b, double scale);
double coords_sp(Coords a, Coords b);
Coords coords_xp(Coords b, Coords c);
double coords_len(Coords a);
void   coords_print(Coords a);
mcstatic inline void coords_norm(Coords* c);

void rot_set_rotation(Rotation t, double phx, double phy, double phz);
int  rot_test_identity(Rotation t);
void rot_mul(Rotation t1, Rotation t2, Rotation t3);
void rot_copy(Rotation dest, Rotation src);
void rot_transpose(Rotation src, Rotation dst);
Coords rot_apply(Rotation t, Coords a);

void mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
    double *vx, double *vy, double *vz, double *sx, double *sy, double *sz);
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz);

double mcestimate_error(double N, double p1, double p2);
void mcreadparams(void);

/* this is now in mcstas-r.h and mcxtrace-r.h as the number of state parameters is no longer equal*/
/* void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);
*/
void mcgenstate(void);

/* trajectory/shape intersection routines */
int inside_rectangle(double, double, double, double);
int box_intersect(double *dt_in, double *dt_out, double x, double y, double z,
    double vx, double vy, double vz, double dx, double dy, double dz);
int cylinder_intersect(double *t0, double *t1, double x, double y, double z,
    double vx, double vy, double vz, double r, double h);
int sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r);
/* second order equation roots */
int solve_2nd_order(double *t1, double *t2,
    double A,  double B,  double C);

/* random vector generation to shape */
void randvec_target_circle(double *xo, double *yo, double *zo,
    double *solid_angle, double xi, double yi, double zi, double radius);
#define randvec_target_sphere randvec_target_circle
void randvec_target_rect_angular(double *xo, double *yo, double *zo,
    double *solid_angle,
               double xi, double yi, double zi, double height, double width, Rotation A);
#define randvec_target_rect(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9)  randvec_target_rect_real(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,0,0,0,1)
void randvec_target_rect_real(double *xo, double *yo, double *zo,
    double *solid_angle,
	       double xi, double yi, double zi, double height, double width, Rotation A,
			 double lx, double ly, double lz, int order);

/* this is the main() */
int mccode_main(int argc, char *argv[]);


#endif /* !MCCODE_H */

#ifndef MCCODE_R_IO_H
#define MCCODE_R_IO_H "$Revision$"

#if (USE_NEXUS == 0)
#undef USE_NEXUS
#endif

#ifndef CHAR_BUF_LENGTH
#define CHAR_BUF_LENGTH 1024
#endif

/* I/O section part ========================================================= */

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */


/* main DETECTOR structure which stores most information to write to data files */
struct mcdetector_struct {
  char   filename[CHAR_BUF_LENGTH];   /* file name of monitor */
  char   position[CHAR_BUF_LENGTH];   /* position of detector component */
  char   component[CHAR_BUF_LENGTH];  /* component instance name */
  char   instrument[CHAR_BUF_LENGTH]; /* instrument name */
  char   type[CHAR_BUF_LENGTH];       /* data type, e.g. 0d, 1d, 2d, 3d */
  char   user[CHAR_BUF_LENGTH];       /* user name, e.g. HOME */
  char   date[CHAR_BUF_LENGTH];       /* date of simulation end/write time */
  char   title[CHAR_BUF_LENGTH];      /* title of detector */
  char   xlabel[CHAR_BUF_LENGTH];     /* X axis label */
  char   ylabel[CHAR_BUF_LENGTH];     /* Y axis label */
  char   zlabel[CHAR_BUF_LENGTH];     /* Z axis label */
  char   xvar[CHAR_BUF_LENGTH];       /* X variable name */
  char   yvar[CHAR_BUF_LENGTH];       /* Y variable name */
  char   zvar[CHAR_BUF_LENGTH];       /* Z variable name */
  char   ncount[CHAR_BUF_LENGTH];     /* number of events initially generated */
  char   limits[CHAR_BUF_LENGTH];     /* X Y Z limits, e.g. [xmin xmax ymin ymax zmin zmax] */
  char   variables[CHAR_BUF_LENGTH];  /* variables written into data block */
  char   statistics[CHAR_BUF_LENGTH]; /* center, mean and half width along axis */
  char   signal[CHAR_BUF_LENGTH];     /* min max and mean of signal (data block) */
  char   values[CHAR_BUF_LENGTH];     /* integrated values e.g. [I I_err N] */
  double xmin,xmax;                   /* min max of axes */
  double ymin,ymax;
  double zmin,zmax;
  double intensity;                   /* integrated values for data block */
  double error;
  double events;
  double min;                         /* statistics for data block */
  double max;
  double mean;
  double centerX;                     /* statistics for axes */
  double halfwidthX;
  double centerY;
  double halfwidthY;
  int    rank;                        /* dimensionaly of monitor, e.g. 0 1 2 3 */
  char   istransposed;                /* flag to transpose matrix for some formats */

  long   m,n,p;                       /* dimensions of data block and along axes */
  long   date_l;                      /* same as date, but in sec since 1970 */

  double *p0, *p1, *p2;               /* pointers to saved data, NULL when freed */
  char   format[CHAR_BUF_LENGTH];    /* format for file generation */
};

typedef struct mcdetector_struct MCDETECTOR;

static   char *mcdirname             = NULL;      /* name of output directory */
static   char *mcsiminfo_name        = "mccode";  /* default output sim file name */
char    *mcformat                    = NULL;      /* NULL (default) or a specific format */

/* file I/O definitions and function prototypes */

#ifndef MC_EMBEDDED_RUNTIME /* the mcstatic variables (from mccode-r.c) */
extern FILE * mcsiminfo_file;     /* handle to the output siminfo file */
extern int    mcgravitation;      /* flag to enable gravitation */
extern int    mcdotrace;          /* flag to print MCDISPLAY messages */
#else
mcstatic FILE *mcsiminfo_file        = NULL;
#endif

/* I/O function prototypes ================================================== */

/* output functions */
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2, char *c, Coords pos);
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
                  char *xvar, double x1, double x2, long n,
                  double *p0, double *p1, double *p2, char *f, char *c, Coords pos);
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2, long m,
                  long n, double *p0, double *p1, double *p2, char *f,
                  char *c, Coords pos);
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa);

/* wrappers to output functions, that automatically set NAME and POSITION */
#define DETECTOR_OUT(p0,p1,p2) mcdetector_out_0D(NAME_CURRENT_COMP,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_0D(t,p0,p1,p2) mcdetector_out_0D(t,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f) \
     mcdetector_out_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f) \
     mcdetector_out_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)

#ifdef USE_NEXUS
#include "napi.h"
NXhandle nxhandle;
#endif

#endif /* ndef MCCODE_R_IO_H */

#endif /* MCCODE_R_H */
/* End of file "mccode-r.h". */

#line 691 "let_new.c"

#line 1 "mcstas-r.h"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system header for McStas.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#define MCSTAS_R_H "$Revision$"

/* Following part is only embedded when not redundent with mcstas.h ========= */

#ifndef MCCODE_H

#define AA2MS    629.622368        /* Convert k[1/AA] to v[m/s] */
#define MS2AA    1.58825361e-3     /* Convert v[m/s] to k[1/AA] */
#define K2V      AA2MS
#define V2K      MS2AA
#define Q2V      AA2MS
#define V2Q      MS2AA
#define SE2V     437.393377        /* Convert sqrt(E)[meV] to v[m/s] */
#define VS2E     5.22703725e-6     /* Convert (v[m/s])**2 to E[meV] */

#define SCATTER do {mcDEBUG_SCATTER(mcnlx, mcnly, mcnlz, mcnlvx, mcnlvy, mcnlvz, \
        mcnlt,mcnlsx,mcnlsy,mcnlsz, mcnlp); mcScattered++;} while(0)
#define ABSORB do {mcDEBUG_STATE(mcnlx, mcnly, mcnlz, mcnlvx, mcnlvy, mcnlvz, \
        mcnlt,mcnlsx,mcnlsy,mcnlsz, mcnlp); mcDEBUG_ABSORB(); MAGNET_OFF; goto mcabsorb;} while(0)

#define STORE_NEUTRON(index, x, y, z, vx, vy, vz, t, sx, sy, sz, p) \
  mcstore_neutron(mccomp_storein,index, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
#define RESTORE_NEUTRON(index, x, y, z, vx, vy, vz, t, sx, sy, sz, p) \
  mcrestore_neutron(mccomp_storein,index, &x, &y, &z, &vx, &vy, &vz, &t, &sx, &sy, &sz, &p);

#define MAGNET_ON \
  do { \
    mcMagnet = 1; \
  } while(0)

#define MAGNET_OFF \
  do { \
    mcMagnet = 0; \
  } while(0)

#define ALLOW_BACKPROP \
  do { \
    mcallowbackprop = 1; \
  } while(0)

#define DISALLOW_BACKPROP \
  do { \
    mcallowbackprop = 0; \
  } while(0)

#define PROP_MAGNET(dt) \
  do { \
  }while (0)
    /* change coordinates from local system to magnet system */
/*    Rotation rotLM, rotTemp; \
      Coords   posLM = coords_sub(POS_A_CURRENT_COMP, mcMagnetPos); \
      rot_transpose(ROT_A_CURRENT_COMP, rotTemp); \
      rot_mul(rotTemp, mcMagnetRot, rotLM); \
      mcMagnetPrecession(mcnlx, mcnly, mcnlz, mcnlt, mcnlvx, mcnlvy, mcnlvz, \
               &mcnlsx, &mcnlsy, &mcnlsz, dt, posLM, rotLM); \
      } while(0)
*/

#define mcPROP_DT(dt) \
  do { \
    if (mcMagnet && dt > 0) PROP_MAGNET(dt);\
    mcnlx += mcnlvx*(dt); \
    mcnly += mcnlvy*(dt); \
    mcnlz += mcnlvz*(dt); \
    mcnlt += (dt); \
    if (isnan(p) || isinf(p)) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
  } while(0)

/* ADD: E. Farhi, Aug 6th, 2001 PROP_GRAV_DT propagation with acceleration */
#define PROP_GRAV_DT(dt, Ax, Ay, Az) \
  do { \
    if(dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
    if (mcMagnet) printf("Spin precession gravity\n"); \
    mcnlx  += mcnlvx*(dt) + (Ax)*(dt)*(dt)/2; \
    mcnly  += mcnlvy*(dt) + (Ay)*(dt)*(dt)/2; \
    mcnlz  += mcnlvz*(dt) + (Az)*(dt)*(dt)/2; \
    mcnlvx += (Ax)*(dt); \
    mcnlvy += (Ay)*(dt); \
    mcnlvz += (Az)*(dt); \
    mcnlt  += (dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_DT(dt) \
  do { \
    if(dt < 0) { RESTORE=1; goto mcabsorbComp; }; \
    if (mcgravitation) { Coords mcLocG; double mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    PROP_GRAV_DT(dt, mc_gx, mc_gy, mc_gz); } \
    else mcPROP_DT(dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_Z0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gz/2, -mcnlvz, -mcnlz); \
    if (mc_ret && mc_dt>=0) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); mcnlz=0;}\
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_Z0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_Z0 \
  do { \
    double mc_dt; \
    if(mcnlvz == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnlz/mcnlvz; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnlz = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_X0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gx/2, -mcnlvx, -mcnlx); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_X0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_X0 \
  do { \
    double mc_dt; \
    if(mcnlvx == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnlx/mcnlvx; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnlx = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_Y0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gy/2, -mcnlvy, -mcnly); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_Y0; \
    DISALLOW_BACKPROP;\
  } while(0)


#define mcPROP_Y0 \
  do { \
    double mc_dt; \
    if(mcnlvy == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnly/mcnlvy; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnly = 0; \
    DISALLOW_BACKPROP; \
  } while(0)

/*moved from mccode-r.h*/
void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);

#ifdef DEBUG

#define mcDEBUG_STATE(x,y,z,vx,vy,vz,t,sx,sy,sz,p) if(!mcdotrace); else \
  printf("STATE: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);
#define mcDEBUG_SCATTER(x,y,z,vx,vy,vz,t,sx,sy,sz,p) if(!mcdotrace); else \
  printf("SCATTER: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);

#else

#define mcDEBUG_STATE(x,y,z,vx,vy,vz,t,sx,sy,sz,p)
#define mcDEBUG_SCATTER(x,y,z,vx,vy,vz,t,sx,sy,sz,p)

#endif

#endif /* !MCCODE_H */

#endif /* MCSTAS_R_H */
/* End of file "mcstas-r.h". */

#line 924 "let_new.c"

#line 1 "mccode-r.c"
/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y/McXtrace X.Y
* Version: $Revision$
*
* Runtime system for McStas and McXtrace.
* Embedded within instrument in runtime mode.
* Contains SECTIONS:
*   MPI handling (sum, send, recv)
*   format definitions
*   I/O
*   mcdisplay support
*   random numbers
*   coordinates handling
*   vectors math (solve 2nd order, normals, randvec...)
*   parameter handling
*   signal and main handlers
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/


/** Include header files to avoid implicit declarations (not allowed on LLVM) */
#include <ctype.h>
#include <sys/types.h>

// UNIX specific headers (non-Windows)
#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#include <sys/stat.h>
#endif


#ifndef DANSE
#ifdef MC_ANCIENT_COMPATIBILITY
int mctraceenabled = 0;
int mcdefaultmain  = 0;
#endif
/* else defined directly in the McCode generated C code */

static   long mcseed                 = 0; /* seed for random generator */
static   long mcstartdate            = 0; /* start simulation time */
static   int  mcdisable_output_files = 0; /* --no-output-files */
mcstatic int  mcgravitation          = 0; /* use gravitation flag, for PROP macros */
int      mcMagnet                    = 0; /* magnet stack flag */
mcstatic int  mcdotrace              = 0; /* flag for --trace and messages for DISPLAY */
int      mcallowbackprop             = 0;         /* flag to enable negative/backprop */

/* Number of particle histories to simulate. */
#ifdef NEUTRONICS
mcstatic unsigned long long int mcncount             = 1;
mcstatic unsigned long long int mcrun_num            = 0;
#else
mcstatic unsigned long long int mcncount             = 1000000;
mcstatic unsigned long long int mcrun_num            = 0;
#endif /* NEUTRONICS */

#else
#include "mcstas-globals.h"
#endif /* !DANSE */

/* SECTION: MPI handling ==================================================== */

#ifdef USE_MPI
/* MPI rank */
static int mpi_node_rank;
static int mpi_node_root = 0;


/*******************************************************************************
* mc_MPI_Reduce: Gathers arrays from MPI nodes using Reduce function.
*******************************************************************************/
int mc_MPI_Sum(double *sbuf, long count)
{
  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to reduce */
  else {
    /* we must cut the buffer into blocks not exceeding the MPI max buffer size of 32000 */
    long   offset=0;
    double *rbuf=NULL;
    int    length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */
    int    i=0;
    rbuf = calloc(count, sizeof(double));
    if (!rbuf)
      exit(-fprintf(stderr, "Error: Out of memory %li (mc_MPI_Sum)\n", count*sizeof(double)));
    while (offset < count) {
      if (!length || offset+length > count-1) length=count-offset;
      else length=MPI_REDUCE_BLOCKSIZE;
      if (MPI_Allreduce((double*)(sbuf+offset), (double*)(rbuf+offset),
              length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
        return MPI_ERR_COUNT;
      offset += length;
    }

    for (i=0; i<count; i++) sbuf[i] = rbuf[i];
    free(rbuf);
  }
  return MPI_SUCCESS;
} /* mc_MPI_Sum */

/*******************************************************************************
* mc_MPI_Send: Send array to MPI node by blocks to avoid buffer limit
*******************************************************************************/
int mc_MPI_Send(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int dest)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to send */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Send((void*)(sbuf+offset*dsize), length, dtype, dest, tag++, MPI_COMM_WORLD) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Send */

/*******************************************************************************
* mc_MPI_Recv: Receives arrays from MPI nodes by blocks to avoid buffer limit
*             the buffer must have been allocated previously.
*******************************************************************************/
int mc_MPI_Recv(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int source)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to recv */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Recv((void*)(sbuf+offset*dsize), length, dtype, source, tag++,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Recv */

#endif /* USE_MPI */

/* SECTION: parameters handling ============================================= */

/* Instrument input parameter type handling. */
/*******************************************************************************
* mcparm_double: extract double value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_double(char *s, void *vptr)
{
  char *p;
  double *v = (double *)vptr;

  if (!s) { *v = 0; return(1); }
  *v = strtod(s, &p);
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_double: display parameter type double
*******************************************************************************/
static char *
mcparminfo_double(char *parmname)
{
  return "double";
}

/*******************************************************************************
* mcparmerror_double: display error message when failed extract double
*******************************************************************************/
static void
mcparmerror_double(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for floating point parameter %s (mcparmerror_double)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_double: convert double to string
*******************************************************************************/
static void
mcparmprinter_double(char *f, void *vptr)
{
  double *v = (double *)vptr;
  sprintf(f, "%g", *v);
}

/*******************************************************************************
* mcparm_int: extract int value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_int(char *s, void *vptr)
{
  char *p;
  int *v = (int *)vptr;
  long x;

  if (!s) { *v = 0; return(1); }
  *v = 0;
  x = strtol(s, &p, 10);
  if(x < INT_MIN || x > INT_MAX)
    return 0;                        /* Under/overflow */
  *v = x;
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_int: display parameter type int
*******************************************************************************/
static char *
mcparminfo_int(char *parmname)
{
  return "int";
}

/*******************************************************************************
* mcparmerror_int: display error message when failed extract int
*******************************************************************************/
static void
mcparmerror_int(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for integer parameter %s (mcparmerror_int)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_int: convert int to string
*******************************************************************************/
static void
mcparmprinter_int(char *f, void *vptr)
{
  int *v = (int *)vptr;
  sprintf(f, "%d", *v);
}

/*******************************************************************************
* mcparm_string: extract char* value from 's' into 'vptr' (copy)
*******************************************************************************/
static int
mcparm_string(char *s, void *vptr)
{
  char **v = (char **)vptr;
  if (!s) { *v = NULL; return(1); }
  *v = (char *)malloc(strlen(s) + 1);
  if(*v == NULL)
  {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcparm_string).\n", (long)strlen(s) + 1));
  }
  strcpy(*v, s);
  return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_string: display parameter type string
*******************************************************************************/
static char *
mcparminfo_string(char *parmname)
{
  return "string";
}

/*******************************************************************************
* mcparmerror_string: display error message when failed extract string
*******************************************************************************/
static void
mcparmerror_string(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for string parameter %s (mcparmerror_string)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_string: convert string to string (including esc chars)
*******************************************************************************/
static void
mcparmprinter_string(char *f, void *vptr)
{
  char **v = (char **)vptr;
  char *p;

  if (!*v) { *f='\0'; return; }
  strcpy(f, "");
  for(p = *v; *p != '\0'; p++)
  {
    switch(*p)
    {
      case '\n':
        strcat(f, "\\n");
        break;
      case '\r':
        strcat(f, "\\r");
        break;
      case '"':
        strcat(f, "\\\"");
        break;
      case '\\':
        strcat(f, "\\\\");
        break;
      default:
        strncat(f, p, 1);
    }
  }
  /* strcat(f, "\""); */
} /* mcparmprinter_string */

/* now we may define the parameter structure, using previous functions */
static struct
  {
    int (*getparm)(char *, void *);
    char * (*parminfo)(char *);
    void (*error)(char *, char *);
    void (*printer)(char *, void *);
} mcinputtypes[] = {
  {
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }, {
    mcparm_int, mcparminfo_int, mcparmerror_int,
    mcparmprinter_int
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
  }
};

/*******************************************************************************
* mcestimate_error: compute sigma from N,p,p2 in Gaussian large numbers approx
*******************************************************************************/
double mcestimate_error(double N, double p1, double p2)
{
  double pmean, n1;
  if(N <= 1)
    return p1;
  pmean = p1 / N;
  n1 = N - 1;
  /* Note: underflow may cause p2 to become zero; the fabs() below guards
     against this. */
  return sqrt((N/n1)*fabs(p2 - pmean*pmean));
}

double (*mcestimate_error_p)
  (double V2, double psum, double p2sum)=mcestimate_error;

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */

#ifndef MCCODE_R_IO_C
#define MCCODE_R_IO_C "$Revision$"

/* SECTION: file i/o handling ================================================ */

#ifndef HAVE_STRCASESTR
// from msysgit: https://code.google.com/p/msysgit/source/browse/compat/strcasestr.c
char *strcasestr(const char *haystack, const char *needle)
{
  int nlen = strlen(needle);
  int hlen = strlen(haystack) - nlen + 1;
  int i;

  for (i = 0; i < hlen; i++) {
    int j;
    for (j = 0; j < nlen; j++) {
            unsigned char c1 = haystack[i+j];
            unsigned char c2 = needle[j];
            if (toupper(c1) != toupper(c2))
                    goto next;
    }
    return (char *) haystack + i;
  next:
    ;
  }
  return NULL;
}


#endif
#ifndef HAVE_STRCASECMP
int strcasecmp( const char *s1, const char *s2 )
{
  int c1, c2;
  do {
    c1 = tolower( (unsigned char) *s1++ );
    c2 = tolower( (unsigned char) *s2++ );
  } while (c1 == c2 && c1 != 0);
  return c2 > c1 ? -1 : c1 > c2;
}
#endif

/*******************************************************************************
* mcfull_file: allocates a full file name=mcdirname+file. Catenate extension if missing.
*******************************************************************************/
char *mcfull_file(char *name, char *ext)
{
  int   dirlen=0;
  char *mem   =NULL;

  dirlen = mcdirname ? strlen(mcdirname) : 0;
  mem = (char*)malloc(dirlen + strlen(name) + CHAR_BUF_LENGTH);
  if(!mem) {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcfull_file)\n", (long)(dirlen + strlen(name) + 256)));
  }
  strcpy(mem, "");

  /* prepend directory name to path if name does not contain a path */
  if (dirlen > 0 && !strchr(name, MC_PATHSEP_C)) {
    strcat(mem, mcdirname);
    strcat(mem, MC_PATHSEP_S);
  } /* dirlen */

  strcat(mem, name);
  if (!strchr(name, '.') && ext && strlen(ext))
  { /* add extension if not in file name already */
    strcat(mem, ".");
    strcat(mem, ext);
  }
  return(mem);
} /* mcfull_file */

/*******************************************************************************
* mcnew_file: opens a new file within mcdirname if non NULL
*             the file is opened in "a" (append, create if does not exist)
*             the extension 'ext' is added if the file name does not include one.
*             the last argument is set to 0 if file did not exist, else to 1.
*******************************************************************************/
FILE *mcnew_file(char *name, char *ext, int *exists)
{
  char *mem;
  FILE *file=NULL;

  if (!name || strlen(name) == 0 || mcdisable_output_files) return(NULL);
  
  mem  = mcfull_file(name, ext); /* create mcdirname/name.ext */
  
  /* check for existence */
  file = fopen(mem, "r"); /* for reading -> fails if does not exist */
  if (file) {
    fclose(file);
    *exists=1;
  } else
    *exists=0;
  
  /* open the file for writing/appending */
#ifdef USE_NEXUS
  if (mcformat && strcasestr(mcformat, "NeXus")) {
    /* NXhandle nxhandle is defined in the .h with USE_NEXUS */
    NXaccess mode = (*exists ? NXACC_CREATE5 | NXACC_RDWR : NXACC_CREATE5);
      
    if (NXopen(mem, mode, &nxhandle) != NX_OK)
      file = NULL;
    else
      file = (FILE*)&nxhandle; /* to make it non NULL */
  } else
#endif
    file = fopen(mem, "a+"); 
    
  if(!file)
    fprintf(stderr, "Warning: could not open output file '%s' for %s (mcnew_file)\n", 
      mem, *exists ? "append" : "create");
  free(mem);

  return file;
} /* mcnew_file */

/*******************************************************************************
* mcdetector_statistics: compute detector statistics, error bars, [x I I_err N] 1D
* RETURN:            updated detector structure
* Used by: mcdetector_import
*******************************************************************************/
MCDETECTOR mcdetector_statistics(
  MCDETECTOR detector)
{

  if (!detector.p1 || !detector.m || !detector.filename)
    return(detector);
  
  /* compute statistics and update MCDETECTOR structure ===================== */
  double sum_z  = 0, min_z  = 0, max_z  = 0;
  double fmon_x =0,  smon_x = 0, fmon_y =0, smon_y=0, mean_z=0;
  double Nsum=0, P2sum=0;

  double sum_xz = 0, sum_yz = 0, sum_x = 0, sum_y = 0, sum_x2z = 0, sum_y2z = 0;
  int    i,j;
  char   hasnan=0, hasinf=0;
  char   israw = ((char*)strcasestr(detector.format,"raw") != NULL);
  double *this_p1=NULL; /* new 1D McCode array [x I E N]. Freed after writing data */

  /* if McCode/PGPLOT and rank==1 we create a new m*4 data block=[x I E N] */
  if (detector.rank == 1 && strcasestr(detector.format,"McCode")) {
    this_p1 = (double *)calloc(detector.m*detector.n*detector.p*4, sizeof(double));
    if (!this_p1)
      exit(-fprintf(stderr, "Error: Out of memory creating %li 1D " MCCODE_STRING " data set for file '%s' (mcdetector_import)\n",
        detector.m*detector.n*detector.p*4*sizeof(double*), detector.filename));
  }

  max_z = min_z = detector.p1[0];
  
  /* compute sum and moments (not for lists) */
  if (!strcasestr(detector.format,"list") && detector.m)
  for(j = 0; j < detector.n*detector.p; j++)
  {
    for(i = 0; i < detector.m; i++)
    {
      double x,y,z;
      double N, E;
      long   index= !detector.istransposed ? i*detector.n*detector.p + j : i+j*detector.m;
      char   hasnaninf=0;

      if (detector.m) 
        x = detector.xmin + (i + 0.5)/detector.m*(detector.xmax - detector.xmin); 
      else x = 0;
      if (detector.n && detector.p) 
        y = detector.ymin + (j + 0.5)/detector.n/detector.p*(detector.ymax - detector.ymin); 
      else y = 0;
      z = detector.p1[index];
      N = detector.p0 ? detector.p0[index] : 1;
      E = detector.p2 ? detector.p2[index] : 0;
      if (detector.p2 && !israw) 
        detector.p2[index] = (*mcestimate_error_p)(detector.p0[index],detector.p1[index],detector.p2[index]); /* set sigma */
      
      if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
        /* fill-in 1D McCode array [x I E N] */
        this_p1[index*4]   = x;
        this_p1[index*4+1] = z;
        this_p1[index*4+2] = detector.p2 ? detector.p2[index] : 0;
        this_p1[index*4+3] = N;
      }
      
      if (isnan(z) || isnan(E) || isnan(N)) hasnaninf=hasnan=1;
      if (isinf(z) || isinf(E) || isinf(N)) hasnaninf=hasinf=1;

      /* compute stats integrals */
      if (!hasnaninf) {
        sum_xz += x*z;
        sum_yz += y*z;
        sum_x  += x;
        sum_y  += y;
        sum_z  += z;
        sum_x2z += x*x*z;
        sum_y2z += y*y*z;
        if (z > max_z) max_z = z;
        if (z < min_z) min_z = z;

        Nsum += N;
        P2sum += E;
      }

    }
  } /* for j */

  /* compute 1st and 2nd moments. For lists, sum_z=0 so this is skipped. */
  if (sum_z && detector.n*detector.m*detector.p)
  {
    fmon_x = sum_xz/sum_z;
    fmon_y = sum_yz/sum_z;
    smon_x = sum_x2z/sum_z-fmon_x*fmon_x; smon_x = smon_x > 0 ? sqrt(smon_x) : 0;
    smon_y = sum_y2z/sum_z-fmon_y*fmon_y; smon_y = smon_y > 0 ? sqrt(smon_y) : 0;
    mean_z = sum_z/detector.n/detector.m/detector.p;
  }
  /* store statistics into detector */
  detector.intensity = sum_z;
  detector.error     = Nsum ? (*mcestimate_error_p)(Nsum, sum_z, P2sum) : 0;
  detector.events    = Nsum;
  detector.min       = min_z;
  detector.max       = max_z;
  detector.mean      = mean_z;
  detector.centerX   = fmon_x;
  detector.halfwidthX= smon_x;
  detector.centerY   = fmon_y;
  detector.halfwidthY= smon_y;

  /* if McCode/PGPLOT and rank==1 replace p1 with new m*4 1D McCode and clear others */
  if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
    
    detector.p1 = this_p1;
    detector.n  = detector.m; detector.m  = 4;
    detector.p0 = detector.p2 = NULL;
    detector.istransposed = 1;
  }

  if (detector.n*detector.m*detector.p > 1)
    snprintf(detector.signal, CHAR_BUF_LENGTH, 
      "Min=%g; Max=%g; Mean=%g;", detector.min, detector.max, detector.mean);
  else
    strcpy(detector.signal, "None");
  snprintf(detector.values, CHAR_BUF_LENGTH,
    "%g %g %g", detector.intensity, detector.error, detector.events);

  switch (detector.rank) {
    case 1:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g;",
      detector.centerX, detector.halfwidthX); break;
    case 2:
    case 3:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g; Y0=%g; dY=%g;",
      detector.centerX, detector.halfwidthX, detector.centerY, detector.halfwidthY);
      break;
    default: strcpy(detector.statistics, "None");
  }
  
  if (hasnan)
    printf("WARNING: Nan detected in component/file %s %s\n", 
      detector.component, strlen(detector.filename) ? detector.filename : "");
  if (hasinf)
    printf("WARNING: Inf detected in component/file %s %s\n", 
      detector.component, strlen(detector.filename) ? detector.filename : "");
  
  return(detector);
  
} /* mcdetector_statistics */

/*******************************************************************************
* mcdetector_import: build detector structure, merge non-lists from MPI
*                    compute basic stat, write "Detector:" line
* RETURN:            detector structure. Invalid data if detector.p1 == NULL
*                    Invalid detector sets m=0 and filename=""
*                    Simulation data  sets m=0 and filename=mcsiminfo_name
* This function is equivalent to the old 'mcdetector_out', returning a structure
*******************************************************************************/
MCDETECTOR mcdetector_import(
  char *format,
  char *component, char *title,
  long m, long n,  long p,
  char *xlabel, char *ylabel, char *zlabel,
  char *xvar, char *yvar, char *zvar,
  double x1, double x2, double y1, double y2, double z1, double z2,
  char *filename,
  double *p0, double *p1, double *p2,
  Coords position)
{
  time_t t;       /* for detector.date */
  long   date_l;  /* date as a long number */
  char   istransposed=0;
  char   c[CHAR_BUF_LENGTH]; /* temp var for signal label */

  MCDETECTOR detector;

  /* build MCDETECTOR structure ============================================= */
  /* make sure we do not have NULL for char fields */

  /* these also apply to simfile */
  strncpy (detector.filename,  filename ? filename : "",        CHAR_BUF_LENGTH);
  strncpy (detector.format,    format   ? format   : "McCode" , CHAR_BUF_LENGTH);
  /* add extension if missing */
  if (strlen(detector.filename) && !strchr(detector.filename, '.'))
  { /* add extension if not in file name already */
    strcat(detector.filename, ".dat");
  }
  strncpy (detector.component, component ? component : MCCODE_STRING " component", CHAR_BUF_LENGTH);

  snprintf(detector.instrument, CHAR_BUF_LENGTH, "%s (%s)", mcinstrument_name, mcinstrument_source);
  snprintf(detector.user, CHAR_BUF_LENGTH,      "%s on %s",
        getenv("USER") ? getenv("USER") : MCCODE_NAME,
        getenv("HOST") ? getenv("HOST") : "localhost");
  time(&t);         /* get current write time */
  date_l = (long)t; /* same but as a long */
  snprintf(detector.date, CHAR_BUF_LENGTH, "%s", ctime(&t));
  if (strlen(detector.date))   detector.date[strlen(detector.date)-1] = '\0'; /* remove last \n in date */
  detector.date_l = date_l;

  if (!mcget_run_num() || mcget_run_num() >= mcget_ncount())
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%llu", mcget_ncount()
#ifdef USE_MPI
*mpi_node_count
#endif
  );
  else
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%g/%g", (double)mcget_run_num(), (double)mcget_ncount());

  detector.p0         = p0;
  detector.p1         = p1;
  detector.p2         = p2;

  /* handle transposition (not for NeXus) */
  if (!strcasestr(detector.format, "NeXus")) {
    if (m<0 || n<0 || p<0)             istransposed = !istransposed;
    if (strcasestr(detector.format, "transpose")) istransposed = !istransposed;
    if (istransposed) { /* do the swap once for all */
      long i=m; m=n; n=i;
    }
  }

  m=abs(m); n=abs(n); p=abs(p); /* make sure dimensions are positive */
  detector.istransposed = istransposed;

  /* determine detector rank (dimensionality) */
  if (!m || !n || !p || !p1) detector.rank = 4; /* invalid: exit with m=0 filename="" */
  else if (m*n*p == 1)       detector.rank = 0; /* 0D */
  else if (n == 1 || m == 1) detector.rank = 1; /* 1D */
  else if (p == 1)           detector.rank = 2; /* 2D */
  else                       detector.rank = 3; /* 3D */

  /* from rank, set type */
  switch (detector.rank) {
    case 0:  strcpy(detector.type,  "array_0d"); m=n=p=1; break;
    case 1:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_1d(%ld)", m*n*p); m *= n*p; n=p=1; break;
    case 2:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_2d(%ld, %ld)", m, n*p); n *= p; p=1; break;
    case 3:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_3d(%ld, %ld, %ld)", m, n, p); break;
    default: m=0; strcpy(detector.type, ""); strcpy(detector.filename, "");/* invalid */
  }

  detector.m    = m;
  detector.n    = n;
  detector.p    = p;

  /* these only apply to detector files ===================================== */

  snprintf(detector.position, CHAR_BUF_LENGTH, "%g %g %g", position.x, position.y, position.z);
  /* may also store actual detector orientation in the future */

  strncpy(detector.title,      title && strlen(title) ? title : component,       CHAR_BUF_LENGTH);
  strncpy(detector.xlabel,     xlabel && strlen(xlabel) ? xlabel : "X", CHAR_BUF_LENGTH); /* axis labels */
  strncpy(detector.ylabel,     ylabel && strlen(ylabel) ? ylabel : "Y", CHAR_BUF_LENGTH);
  strncpy(detector.zlabel,     zlabel && strlen(zlabel) ? zlabel : "Z", CHAR_BUF_LENGTH);
  strncpy(detector.xvar,       xvar && strlen(xvar) ? xvar :       "x", CHAR_BUF_LENGTH); /* axis variables */
  strncpy(detector.yvar,       yvar && strlen(yvar) ? yvar :       detector.xvar, CHAR_BUF_LENGTH);
  strncpy(detector.zvar,       zvar && strlen(zvar) ? zvar :       detector.yvar, CHAR_BUF_LENGTH);

  /* set "variables" as e.g. "I I_err N" */
  strcpy(c, "I ");
  if (strlen(detector.zvar))      strncpy(c, detector.zvar,32);
  else if (strlen(detector.yvar)) strncpy(c, detector.yvar,32);
  else if (strlen(detector.xvar)) strncpy(c, detector.xvar,32);

  if (detector.rank == 1)
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s %s_err N", detector.xvar, c, c);
  else
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s_err N", c, c);

  /* limits */
  detector.xmin = x1;
  detector.xmax = x2;
  detector.ymin = y1;
  detector.ymax = y2;
  detector.zmin = z1;
  detector.zmax = z2;
  if (abs(detector.rank) == 1)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g", x1, x2);
  else if (detector.rank == 2)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g", x1, x2, y1, y2);
  else
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g %g %g", x1, x2, y1, y2, z1, z2);

  /* if MPI and nodes_nb > 1: reduce data sets when using MPI =============== */
#ifdef USE_MPI
  if (!strcasestr(detector.format,"list") && mpi_node_count > 1 && m) {
    /* we save additive data: reduce everything into mpi_node_root */
    if (p0) mc_MPI_Sum(p0, m*n*p);
    if (p1) mc_MPI_Sum(p1, m*n*p);
    if (p2) mc_MPI_Sum(p2, m*n*p);
    if (!p0) {  /* additive signal must be then divided by the number of nodes */
      int i;
      for (i=0; i<m*n*p; i++) {
        p1[i] /= mpi_node_count;
        if (p2) p2[i] /= mpi_node_count;
      }
    }
  }
#endif /* USE_MPI */

  /* compute statistics, Nsum, intensity, Error bars */
  detector = mcdetector_statistics(detector);

#ifdef USE_MPI
  /* slaves are done */
  if(mpi_node_rank != mpi_node_root) {
    return detector;
  }
#endif

  /* output "Detector:" line ================================================ */
  /* when this is a detector written by a component (not the SAVE from instrument),
     not an event lists */
  if (!m) return(detector);
  if (!strcasestr(detector.format,"list")) {
    if (!strcmp(detector.component, mcinstrument_name)) {
      if (strlen(detector.filename))  /* we name it from its filename, or from its title */
        strncpy(c, detector.filename, CHAR_BUF_LENGTH);
      else
        snprintf(c, CHAR_BUF_LENGTH, "%s", mcinstrument_name);
    } else
      strncpy(c, detector.component, CHAR_BUF_LENGTH);  /* usual detectors written by components */

    printf("Detector: %s_I=%g %s_ERR=%g %s_N=%g",
           c, detector.intensity,
           c, detector.error,
           c, detector.events);
    printf(" \"%s\"\n", strlen(detector.filename) ? detector.filename : detector.component);
  }
  

  return(detector);
} /* mcdetector_import */

/* end MCDETECTOR import section ============================================ */

















/* ========================================================================== */

/*                               ASCII output                                 */
/*     The SIM file is YAML based, the data files have '#' headers            */

/* ========================================================================== */


/*******************************************************************************
* mcinfo_out: output instrument tags/info (only in SIM)
* Used in: mcsiminfo_init (ascii), mcinfo(stdout)
*******************************************************************************/
static void mcinfo_out(char *pre, FILE *f)
{
  char Parameters[CHAR_BUF_LENGTH] = "";
  int  i;

  if (!f || mcdisable_output_files) return;

  /* create parameter string ================================================ */
  for(i = 0; i < mcnumipar; i++)
  {
    char ThisParam[CHAR_BUF_LENGTH];
    if (strlen(mcinputtable[i].name) > CHAR_BUF_LENGTH) break;
    snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
            (*mcinputtypes[mcinputtable[i].type].parminfo)
                (mcinputtable[i].name));
    strcat(Parameters, ThisParam);
    if (strlen(Parameters) >= CHAR_BUF_LENGTH-64) break;
  }

  /* output data ============================================================ */
  if (f != stdout)
    fprintf(f, "%sFile: %s%c%s\n",    pre, mcdirname, MC_PATHSEP_C, mcsiminfo_name);
  else
    fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);

  fprintf(f, "%sSource: %s\n",   pre, mcinstrument_source);
  fprintf(f, "%sParameters: %s\n",    pre, Parameters);
  
  fprintf(f, "%sTrace_enabled: %s\n", pre, mctraceenabled ? "yes" : "no");
  fprintf(f, "%sDefault_main: %s\n",  pre, mcdefaultmain ?  "yes" : "no");
  fprintf(f, "%sEmbedded_runtime: %s\n", pre, 
#ifdef MC_EMBEDDED_RUNTIME
         "yes"
#else
         "no"
#endif
         );

  fflush(f);
} /* mcinfo_out */

/*******************************************************************************
* mcruninfo_out: output simulation tags/info (both in SIM and data files)
* Used in: mcsiminfo_init (ascii case), mcdetector_out_xD_ascii
*******************************************************************************/
static void mcruninfo_out(char *pre, FILE *f)
{
  int i;
  char Parameters[CHAR_BUF_LENGTH];

  if (!f || mcdisable_output_files) return;

  fprintf(f, "%sFormat: %s%s\n",      pre, 
    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME,
    mcformat && strcasestr(mcformat,"McCode") ? " with text headers" : "");
  fprintf(f, "%sURL: %s\n",         pre, "http://www.mccode.org");
  fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);
  fprintf(f, "%sInstrument: %s\n", pre, mcinstrument_source);
  fprintf(f, "%sNcount: %llu\n",        pre, mcget_ncount());
  fprintf(f, "%sTrace: %s\n",       pre, mcdotrace ? "yes" : "no");
  fprintf(f, "%sGravitation: %s\n", pre, mcgravitation ? "yes" : "no");
  snprintf(Parameters, CHAR_BUF_LENGTH, "%ld", mcseed);
  fprintf(f, "%sSeed: %s\n",        pre, Parameters);
  fprintf(f, "%sDirectory: %s\n",        pre, mcdirname ? mcdirname : ".");
#ifdef USE_MPI
  if (mpi_node_count > 1)
    fprintf(f, "%sNodes: %i\n",        pre, mpi_node_count);
#endif

  /* output parameter string ================================================ */
  for(i = 0; i < mcnumipar; i++) {
    if (mcget_run_num() || (mcinputtable[i].val && strlen(mcinputtable[i].val))) {
      if (mcinputtable[i].par == NULL)
        strncpy(Parameters, (mcinputtable[i].val ? mcinputtable[i].val : ""), CHAR_BUF_LENGTH);
      else
        (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);

      fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
    }
  }
  fflush(f);
} /* mcruninfo_out */

/*******************************************************************************
* mcsiminfo_out:    wrapper to fprintf(mcsiminfo_file)
*******************************************************************************/
void mcsiminfo_out(char *format, ...)
{
  va_list ap;

  if(mcsiminfo_file && !mcdisable_output_files)
  {
    va_start(ap, format);
    vfprintf(mcsiminfo_file, format, ap);
    va_end(ap);
  }
} /* mcsiminfo_out */


/*******************************************************************************
* mcdatainfo_out: output detector header
*   mcdatainfo_out(prefix, file_handle, detector) writes info to data file
*******************************************************************************/
static void
mcdatainfo_out(char *pre, FILE *f, MCDETECTOR detector)
{
  if (!f || !detector.m || mcdisable_output_files) return;
  
  /* output data ============================================================ */
  fprintf(f, "%sDate: %s (%li)\n",       pre, detector.date, detector.date_l);
  fprintf(f, "%stype: %s\n",       pre, detector.type);
  fprintf(f, "%sSource: %s\n",     pre, detector.instrument);
  fprintf(f, "%scomponent: %s\n",  pre, detector.component);
  fprintf(f, "%sposition: %s\n",   pre, detector.position);

  fprintf(f, "%stitle: %s\n",      pre, detector.title);
  fprintf(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
             "%sNcount: %s\n" : 
             "%sratio: %s\n",  pre, detector.ncount);

  if (strlen(detector.filename)) {
    fprintf(f, "%sfilename: %s\n", pre, detector.filename);
  }

  fprintf(f, "%sstatistics: %s\n", pre, detector.statistics);
  fprintf(f, "%ssignal: %s\n",     pre, detector.signal);
  fprintf(f, "%svalues: %s\n",     pre, detector.values);

  if (detector.rank >= 1)
  {
    fprintf(f, "%sxvar: %s\n",     pre, detector.xvar);
    fprintf(f, "%syvar: %s\n",     pre, detector.yvar);
    fprintf(f, "%sxlabel: %s\n",   pre, detector.xlabel);
    fprintf(f, "%sylabel: %s\n",   pre, detector.ylabel);
    if (detector.rank > 1) {
      fprintf(f, "%szvar: %s\n",   pre, detector.zvar);
      fprintf(f, "%szlabel: %s\n", pre, detector.zlabel);
    }
  }

  fprintf(f, 
    abs(detector.rank)==1 ?
             "%sxlimits: %s\n" : 
             "%sxylimits: %s\n", pre, detector.limits);
  fprintf(f, "%svariables: %s\n", pre, 
    strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
    
  fflush(f);

} /* mcdatainfo_out */

/* mcdetector_out_array_ascii: output a single array to a file
 *   m: columns
 *   n: rows
 *   p: array
 *   f: file handle (already opened)
 */
static void mcdetector_out_array_ascii(long m, long n, double *p, FILE *f, char istransposed)
{
  if(f)
  {
    int i,j;
    for(j = 0; j < n; j++)
    {
      for(i = 0; i < m; i++)
      {
          fprintf(f, "%.10g ", p[!istransposed ? i*n + j : j*m+i]);
      }
      fprintf(f,"\n");
    }
  }
} /* mcdetector_out_array_ascii */

/*******************************************************************************
* mcdetector_out_0D_ascii: called by mcdetector_out_0D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_0D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;
  
  /* Write data set information to simulation description file. */
  MPI_MASTER(
    mcsiminfo_out("\nbegin data\n"); // detector.component
    mcdatainfo_out("  ", mcsiminfo_file, detector);
    mcsiminfo_out("end data\n");
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.component, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* write I I_err N */
      fprintf(outfile, "%g %g %g\n", 
        detector.intensity, detector.error, detector.events);
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
} /* mcdetector_out_0D_ascii */

/*******************************************************************************
* mcdetector_out_1D_ascii: called by mcdetector_out_1D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_1D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  MPI_MASTER(
    /* Write data set information to simulation description file. */
    mcsiminfo_out("\nbegin data\n"); // detector.filename
    mcdatainfo_out("  ", mcsiminfo_file, detector);
    mcsiminfo_out("end data\n");
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* output the 1D array columns */
      mcdetector_out_array_ascii(detector.m, detector.n, detector.p1, outfile, detector.istransposed);
      
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
  
}  /* mcdetector_out_1D_ascii */

/*******************************************************************************
* mcdetector_out_2D_ascii: called by mcdetector_out_2D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_2D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;
  
  MPI_MASTER(
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write header only if file has just been created (not appending) */
      if (!exists) {
        /* Write data set information to simulation description file. */
        mcsiminfo_out("\nbegin data\n"); // detector.filename
        mcdatainfo_out("  ", mcsiminfo_file, detector);
        mcsiminfo_out("end data\n");
      
        mcruninfo_out( "# ", outfile);
        mcdatainfo_out("# ", outfile,   detector);
        fprintf(outfile, "# Data [%s/%s] %s:\n", detector.component, detector.filename, detector.zvar);
      }
      mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1, 
        outfile, detector.istransposed);
      if (detector.p2) {
        fprintf(outfile, "# Errors [%s/%s] %s_err:\n", detector.component, detector.filename, detector.zvar);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p2, 
          outfile, detector.istransposed);
      }
      if (detector.p0) {
        fprintf(outfile, "# Events [%s/%s] N:\n", detector.component, detector.filename);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p0, 
          outfile, detector.istransposed);
      }
      fclose(outfile);
      
      if (!exists) {
        if (strcasestr(detector.format, "list"))
          printf("Events:   \"%s\"\n",  
            strlen(detector.filename) ? detector.filename : detector.component);
      }
    } /* if outfile */
  ); /* MPI_MASTER */
#ifdef USE_MPI
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    int node_i=0;
    /* loop along MPI nodes to write sequentially */
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      /* MPI: slaves wait for the master to write its block, then append theirs */
      MPI_Barrier(MPI_COMM_WORLD);
      if (node_i != mpi_node_root && node_i == mpi_node_rank) {
        if(strlen(detector.filename) && !mcdisable_output_files)	/* Don't write if filename is NULL */
          outfile = mcnew_file(detector.filename, "dat", &exists);
        if (!exists)
          fprintf(stderr, "Warning: [MPI node %i] file '%s' does not exist yet, "
                          "MASTER should have opened it before.\n",
            mpi_node_rank, detector.filename);
        if(outfile) {
          mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1, 
            outfile, detector.istransposed);
          fclose(outfile);
        }
      }
    }
  } /* if strcasestr list */
#endif
  return(detector);
} /* mcdetector_out_2D_ascii */

/*******************************************************************************
* strcpy_valid: makes a valid string for variable names.
*   copy 'original' into 'valid', replacing invalid characters by '_'
*   char arrays must be pre-allocated
*******************************************************************************/
static char *strcpy_valid(char *valid, char *original)
{
  long i;
  int  n=32; /* max length of valid names */

  if (original == NULL || !strlen(original)) return(NULL);

  if (n > strlen(original)) n = strlen(original);
  else original += strlen(original)-n;
  strncpy(valid, original, n);

  for (i=0; i < n; i++)
  {
    if ( (valid[i] > 122)
      || (valid[i] < 32)
      || (strchr("!\"#$%&'()*+,-.:;<=>?@[\\]^`/ \n\r\t", valid[i]) != NULL) )
    {
      if (i) valid[i] = '_'; else valid[i] = 'm';
    }
  }
  valid[i] = '\0';

  return(valid);
} /* strcpy_valid */

/* end ascii output section ================================================= */







#ifdef USE_NEXUS

/* ========================================================================== */

/*                               NeXus output                                 */

/* ========================================================================== */

#define nxprintf(...)    nxstr('d', __VA_ARGS__)
#define nxprintattr(...) nxstr('a', __VA_ARGS__)

/*******************************************************************************
* nxstr: output a tag=value data set (char) in NeXus/current group
*   when 'format' is larger that 1024 chars it is used as value for the 'tag'
*   else the value is assembled with format and following arguments.
*   type='d' -> data set
*        'a' -> attribute for current data set
*******************************************************************************/
static int nxstr(char type, NXhandle *f, char *tag, char *format, ...)
{
  va_list ap;
  char value[CHAR_BUF_LENGTH];
  int  i;
  int  ret=NX_OK;
  
  if (!tag || !format || !strlen(tag) || !strlen(format)) return(NX_OK);
  
  /* assemble the value string */
  if (strlen(format) < CHAR_BUF_LENGTH) {
    va_start(ap, format);
    ret = vsnprintf(value, CHAR_BUF_LENGTH, format, ap);
    va_end(ap);
  
    i = strlen(value);
  } else {
    i = strlen(format);
  }

  if (type == 'd') {
    /* open/put/close data set */
    if (NXmakedata (f, tag, NX_CHAR, 1, &i) != NX_OK) return(NX_ERROR);
    NXopendata (f, tag);
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputdata  (f, value);
    else
      ret = NXputdata  (f, format);
    NXclosedata(f);
  } else {
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputattr  (f, tag, value, strlen(value), NX_CHAR);
    else
      ret = NXputattr  (f, tag, format, strlen(format), NX_CHAR);
  }
  
  return(ret);
  
} /* nxstr */

/*******************************************************************************
* mcinfo_readfile: read a full file into a string buffer which is allocated
*   Think to free the buffer after use.
* Used in: mcinfo_out_nexus (nexus)
*******************************************************************************/
char *mcinfo_readfile(char *filename)
{
  FILE *f = fopen(filename, "r");
  if (!f) return(NULL);
  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  rewind(f);
  char *string = malloc(fsize + 1);
  if (string) {
    int n = fread(string, fsize, 1, f);
    fclose(f);

    string[fsize] = 0;
  }
  return(string);
}

/*******************************************************************************
* mcinfo_out: output instrument/simulation groups in NeXus file
* Used in: mcsiminfo_init (nexus)
*******************************************************************************/
static void mcinfo_out_nexus(NXhandle f)
{
  FILE  *fid;     /* for intrument source code/C/IDF */
  char  *buffer=NULL;
  time_t t     =time(NULL); /* for date */
  char   entry0[CHAR_BUF_LENGTH];
  int    count=0;
  char   name[CHAR_BUF_LENGTH];
  char   class[CHAR_BUF_LENGTH];
  
  if (!f || mcdisable_output_files) return;
  
  /* write NeXus NXroot attributes */
  /* automatically added: file_name, HDF5_Version, file_time, NeXus_version */ 
  nxprintattr(f, "creator",   "%s generated with " MCCODE_STRING, mcinstrument_name);
  
  /* count the number of existing NXentry and create the next one */
  NXgetgroupinfo(f, &count, name, class);
  sprintf(entry0, "entry%i", count+1);

  /* create the main NXentry (mandatory in NeXus) */
  if (NXmakegroup(f, entry0, "NXentry") == NX_OK) 
  if (NXopengroup(f, entry0, "NXentry") == NX_OK) {
    
    nxprintf(nxhandle, "program_name", MCCODE_STRING);
    nxprintf(f, "start_time", ctime(&t));
    nxprintf(f, "title", "%s%s%s simulation generated by instrument %s", 
      mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name,
      mcinstrument_name);
    nxprintattr(f, "program_name", MCCODE_STRING);
    nxprintattr(f, "instrument",   mcinstrument_name);
    nxprintattr(f, "simulation",   "%s%s%s",
        mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);

    /* write NeXus instrument group */
    if (NXmakegroup(f, "instrument", "NXinstrument") == NX_OK)
    if (NXopengroup(f, "instrument", "NXinstrument") == NX_OK) {
      int   i;
      char *string=NULL;

      /* write NeXus parameters(types) data =================================== */
      string = (char*)malloc(CHAR_BUF_LENGTH);
      if (string) {
        strcpy(string, "");
        for(i = 0; i < mcnumipar; i++)
        {
          char ThisParam[CHAR_BUF_LENGTH];
          snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
                  (*mcinputtypes[mcinputtable[i].type].parminfo)
                      (mcinputtable[i].name));
          if (strlen(string) + strlen(ThisParam) < CHAR_BUF_LENGTH)
            strcat(string, ThisParam);
        }
        nxprintattr(f, "Parameters",    string);
        free(string);
      }
        
      nxprintattr(f, "name",          mcinstrument_name);
      nxprintf   (f, "name",          mcinstrument_name);
      nxprintattr(f, "Source",        mcinstrument_source);
      
      nxprintattr(f, "Trace_enabled", mctraceenabled ? "yes" : "no");
      nxprintattr(f, "Default_main",  mcdefaultmain ?  "yes" : "no");
      nxprintattr(f, "Embedded_runtime",  
  #ifdef MC_EMBEDDED_RUNTIME
           "yes"
  #else
           "no"
  #endif
           );
           
      /* add instrument source code when available */
      buffer = mcinfo_readfile(mcinstrument_source);
      if (buffer && strlen(buffer)) {
        long length=strlen(buffer);
        nxprintf (f, "description", buffer);
        NXopendata(f,"description");
        nxprintattr(f, "file_name", mcinstrument_source);
        nxprintattr(f, "file_size", "%li", length);
        nxprintattr(f, "MCCODE_STRING", MCCODE_STRING);
        NXclosedata(f);
        nxprintf (f,"instrument_source", "%s " MCCODE_NAME " " MCCODE_PARTICLE " Monte Carlo simulation", mcinstrument_name);
        free(buffer);
      } else
        nxprintf (f, "description", "File %s not found (instrument description %s is missing)", 
          mcinstrument_source, mcinstrument_name);
      
      /* add Mantid/IDF.xml when available */
      char *IDFfile=NULL;
      IDFfile = (char*)malloc(CHAR_BUF_LENGTH);
      sprintf(IDFfile,"%s%s",mcinstrument_source,".xml");
      buffer = mcinfo_readfile(IDFfile);
      if (buffer && strlen(buffer)) {
        NXmakegroup (nxhandle, "instrument_xml", "NXnote");
        NXopengroup (nxhandle, "instrument_xml", "NXnote");
        nxprintf(f, "data", buffer);
        nxprintf(f, "description", "IDF.xml file found with instrument %s", mcinstrument_source);
        nxprintf(f, "type", "text/xml");
        NXclosegroup(f); /* instrument_xml */
        free(buffer);
      }
      free(IDFfile);
      NXclosegroup(f); /* instrument */
    } /* NXinstrument */

    /* write NeXus simulation group */
    if (NXmakegroup(f, "simulation", "NXnote") == NX_OK)
    if (NXopengroup(f, "simulation", "NXnote") == NX_OK) {

      nxprintattr(f, "name",   "%s%s%s",
        mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);
      
      nxprintf   (f, "name",      "%s",     mcsiminfo_name);
      nxprintattr(f, "Format",    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME);
      nxprintattr(f, "URL",       "http://www.mccode.org");
      nxprintattr(f, "program",   MCCODE_STRING);
      nxprintattr(f, "Instrument",mcinstrument_source);
      nxprintattr(f, "Trace",     mcdotrace ?     "yes" : "no");
      nxprintattr(f, "Gravitation",mcgravitation ? "yes" : "no");
      nxprintattr(f, "Seed",      "%li", mcseed);
      nxprintattr(f, "Directory", mcdirname);
    #ifdef USE_MPI
      if (mpi_node_count > 1)
        nxprintf(f, "Nodes", "%i",        mpi_node_count);
    #endif
    
      /* output parameter string ================================================ */
      if (NXmakegroup(f, "Param", "NXparameters") == NX_OK)
      if (NXopengroup(f, "Param", "NXparameters") == NX_OK) {
        int i;
        char string[CHAR_BUF_LENGTH];
        for(i = 0; i < mcnumipar; i++) {
          if (mcget_run_num() || (mcinputtable[i].val && strlen(mcinputtable[i].val))) {
            if (mcinputtable[i].par == NULL)
              strncpy(string, (mcinputtable[i].val ? mcinputtable[i].val : ""), CHAR_BUF_LENGTH);
            else
              (*mcinputtypes[mcinputtable[i].type].printer)(string, mcinputtable[i].par);

            nxprintf(f,  mcinputtable[i].name, "%s", string);
            nxprintattr(f, mcinputtable[i].name, string);
          }
        }
        NXclosegroup(f); /* Param */
      } /* NXparameters */
      
      NXclosegroup(f); /* simulation */
    } /* NXsimulation */
    
    /* create a group to hold all monitors */
    NXmakegroup(f, "data", "NXdetector");

    /* leave the NXentry opened (closed at exit) */
  } /* NXentry */
} /* mcinfo_out_nexus */

/*******************************************************************************
* mcdatainfo_out_nexus: output detector header
*   mcdatainfo_out_nexus(detector) create group and write info to NeXus data file
*   open data:NXdetector then filename:NXdata and write headers/attributes
*   requires: NXentry to be opened
*******************************************************************************/
static void
mcdatainfo_out_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[32];
  if (!f || !detector.m || mcdisable_output_files) return;
  
  strcpy_valid(data_name, 
    detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (mcsiminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* create and open the data group */
    /* this may fail when appending to list -> ignore/skip */
    NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
    
    if (NXmakegroup(f, data_name, "NXdata") == NX_OK)
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {
    
      /* output metadata (as attributes) ======================================== */
      nxprintattr(f, "Date",       detector.date);
      nxprintattr(f, "type",       detector.type);
      nxprintattr(f, "Source",     detector.instrument);
      nxprintattr(f, "component",  detector.component);
      nxprintattr(f, "position",   detector.position);

      nxprintattr(f, "title",      detector.title);
      nxprintattr(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
                 "Ncount" : 
                 "ratio",  detector.ncount);

      if (strlen(detector.filename)) {
        nxprintattr(f, "filename", detector.filename);
      }

      nxprintattr(f, "statistics", detector.statistics);
      nxprintattr(f, "signal",     detector.signal);
      nxprintattr(f, "values",     detector.values);

      if (detector.rank >= 1)
      {
        nxprintattr(f, "xvar",     detector.xvar);
        nxprintattr(f, "yvar",     detector.yvar);
        nxprintattr(f, "xlabel",   detector.xlabel);
        nxprintattr(f, "ylabel",   detector.ylabel);
        if (detector.rank > 1) {
          nxprintattr(f, "zvar",   detector.zvar);
          nxprintattr(f, "zlabel", detector.zlabel);
        }
      }

      nxprintattr(f, abs(detector.rank)==1 ?
                 "xlimits" : 
                 "xylimits", detector.limits);
      nxprintattr(f, "variables", 
        strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
      nxprintf(f, "distance", detector.position);
      nxprintf(f, "acquisition_mode",
        strcasestr(detector.format, "list") ? "event" : "summed");
        
      NXclosegroup(f);
    } /* NXdata (filename) */
    NXMEnableErrorReporting();  /* re-enable NeXus error messages */
    NXclosegroup(f);
  } /* NXdetector (data) */
  
} /* mcdatainfo_out_nexus */

/*******************************************************************************
* mcdetector_out_axis_nexus: write detector axis into current NXdata
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_axis_nexus(NXhandle f, char *label, char *var, int rank, long length, double min, double max)
{
  if (!f || length <= 1 || mcdisable_output_files || max == min) return(NX_OK);
  else {
    double axis[length];
    char valid[32];
    int dim=(int)length;
    int i;
    int nprimary=1;
    /* create an axis from [min:max] */
    for(i = 0; i < length; i++)
      axis[i] = min+(max-min)*(i+0.5)/length;
    /* create the data set */
    strcpy_valid(valid, label);
    NXcompmakedata(f, valid, NX_FLOAT64, 1, &dim, NX_COMP_LZW, &dim);
    /* open it */
    if (NXopendata(f, valid) != NX_OK) {
      fprintf(stderr, "Warning: could not open axis rank %i '%s' (NeXus)\n",
        rank, valid);
      return(NX_ERROR);
    }
    /* put the axis and its attributes */
    NXputdata  (f, axis);
    nxprintattr(f, "long_name",  label);
    nxprintattr(f, "short_name", var);
    NXputattr  (f, "axis",       &rank,     1, NX_INT32);
    nxprintattr(f, "units",      var);
    NXputattr  (f, "primary",    &nprimary, 1, NX_INT32);
    NXclosedata(f);
    
    return(NX_OK);
  }
} /* mcdetector_out_axis_nexus */

/*******************************************************************************
* mcdetector_out_array_nexus: write detector array into current NXdata (1D,2D)
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_array_nexus(NXhandle f, char *part, double *data, MCDETECTOR detector)
{
  
  int dims[3]={detector.m,detector.n,detector.p};  /* number of elements to write */
  int signal=1;
  int exists=0;
  int current_dims[3]={0,0,0};
  int ret=NX_OK;
  
  if (!f || !data || !detector.m || mcdisable_output_files) return(NX_OK);
  
  /* when this is a list, we set 1st dimension to NX_UNLIMITED for creation */
  if (strcasestr(detector.format, "list")) dims[0] = NX_UNLIMITED;
  
  /* create the data set in NXdata group */
  NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
  /* NXcompmakedata fails with NX_UNLIMITED */
  if (strcasestr(detector.format, "list"))
    ret = NXmakedata(    f, part, NX_FLOAT64, detector.rank, dims);
  else
    ret = NXcompmakedata(f, part, NX_FLOAT64, detector.rank, dims, NX_COMP_LZW, dims);
  if (ret != NX_OK) {
    /* failed: data set already exists */
    int datatype=0;
    int rank=0;
    exists=1;
    /* inquire current size of data set (nb of events stored) */
    NXopendata(f, part);
    NXgetinfo(f, &rank, current_dims, &datatype);
    NXclosedata(f);
  }
  NXMEnableErrorReporting();  /* re-enable NeXus error messages */
  dims[0] = detector.m; /* restore actual dimension from data writing */
  
  /* open the data set */
  if (NXopendata(f, part) == NX_ERROR) {
    fprintf(stderr, "Warning: could not open DataSet %s '%s' (NeXus)\n",
      part, detector.title);
    return(NX_ERROR);
  }
  if (strcasestr(detector.format, "list")) {
    current_dims[1] = current_dims[2] = 0; /* set starting location for writing slab */
    NXputslab(f, data, current_dims, dims);
    if (!exists)
      printf("Events:   \"%s\"\n",  
        strlen(detector.filename) ? detector.filename : detector.component);
  } else {
    NXputdata (f, data);
  }
  
  if (strstr(part,"data") || strstr(part, "events")) {
    NXputattr(f, "signal", &signal, 1, NX_INT32);
    nxprintattr(f, "short_name", detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);
  }
  nxprintattr(f, "long_name", "%s '%s'", part, detector.title);
  NXclosedata(f);
  
  return(NX_OK);
} /* mcdetector_out_array_nexus */

/*******************************************************************************
* mcdetector_out_data_nexus: write detector axes+data into current NXdata
*   The data:NXdetector is opened, then filename:NXdata
*   requires: NXentry to be opened
*******************************************************************************/
int mcdetector_out_data_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[32];
  
  if (!f || !detector.m || mcdisable_output_files) return(NX_OK);
  
  strcpy_valid(data_name, 
    detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (mcsiminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* the NXdata group has been created in mcdatainfo_out_nexus */
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {
  
      /* write axes, for histogram data sets, not for lists */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_axis_nexus(f, detector.xlabel, detector.xvar, 
          1, detector.m, detector.xmin, detector.xmax);
          
        mcdetector_out_axis_nexus(f, detector.ylabel, detector.yvar, 
          2, detector.n, detector.ymin, detector.ymax);
          
        mcdetector_out_axis_nexus(f, detector.zlabel, detector.zvar, 
          3, detector.p, detector.zmin, detector.zmax);

      } /* !list */
      
      /* write the actual data (appended if already exists) */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_array_nexus(f, "data", detector.p1, detector);
        mcdetector_out_array_nexus(f, "errors", detector.p2, detector);
        mcdetector_out_array_nexus(f, "ncount", detector.p0, detector);
      } else
        mcdetector_out_array_nexus(  f, "events", detector.p1, detector);
      
      NXclosegroup(f);
    } /* NXdata */
    NXclosegroup(f);
  } /* NXdetector */
  
  return(NX_OK);
} /* mcdetector_out_array_nexus */

#ifdef USE_MPI
/*******************************************************************************
* mcdetector_out_list_slaves: slaves send their list data to master which writes
*   requires: NXentry to be opened
* WARNING: this method has a flaw: it requires all nodes to flush the lists
*   the same number of times. In case one node is just below the buffer size
*   when finishing (e.g. monitor_nd), it may not trigger save but others may. 
*   Then the number of recv/send is not constant along nodes, and simulation stalls.  
*******************************************************************************/
MCDETECTOR mcdetector_out_list_slaves(MCDETECTOR detector)
{
  int     node_i=0;
  MPI_MASTER(
	     printf("\n** MPI master gathering slave node list data ** \n");
  );
  
  if (mpi_node_rank != mpi_node_root) {
    /* MPI slave: slaves send their data to master: 2 MPI_Send calls */
    /* m, n, p must be sent first, since all slaves do not have the same number of events */
    int mnp[3]={detector.m,detector.n,detector.p};

    if (mc_MPI_Send(mnp, 3, MPI_INT, mpi_node_root)!= MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send mnp list error (mcdetector_out_list_slaves)\n", mpi_node_rank);
    if (!detector.p1
     || mc_MPI_Send(detector.p1, mnp[0]*mnp[1]*mnp[2], MPI_DOUBLE, mpi_node_root) != MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send p1 list error: mnp=%i (mcdetector_out_list_slaves)\n", mpi_node_rank, abs(mnp[0]*mnp[1]*mnp[2]));
    /* slaves are done: sent mnp and p1 */
    return (detector);
  } /* end slaves */

  /* MPI master: receive data from slaves sequentially: 2 MPI_Recv calls */

  if (mpi_node_rank == mpi_node_root) {
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      double *this_p1=NULL;                               /* buffer to hold the list from slaves */
      int     mnp[3]={0,0,0};  /* size of this buffer */
      if (node_i != mpi_node_root) { /* get data from slaves */
	if (mc_MPI_Recv(mnp, 3, MPI_INT, node_i) != MPI_SUCCESS)
	  fprintf(stderr, "Warning: master from proc %i: "
		  "MPI_Recv mnp list error (mcdetector_write_data)\n", node_i);
	if (mnp[0]*mnp[1]*mnp[2]) {
	  this_p1 = (double *)calloc(mnp[0]*mnp[1]*mnp[2], sizeof(double));
	  if (!this_p1 || mc_MPI_Recv(this_p1, abs(mnp[0]*mnp[1]*mnp[2]), MPI_DOUBLE, node_i)!= MPI_SUCCESS)
	    fprintf(stderr, "Warning: master from proc %i: "
		    "MPI_Recv p1 list error: mnp=%i (mcdetector_write_data)\n", node_i, mnp[0]*mnp[1]*mnp[2]);
	  else {
	    printf(". MPI master writing data for slave node %i\n",node_i);
	    detector.p1 = this_p1;
	    detector.m  = mnp[0]; detector.n  = mnp[1]; detector.p  = mnp[2];
	    
	    mcdetector_out_data_nexus(nxhandle, detector);
	  }
	}
      } /* if not master */
    } /* for */
  MPI_MASTER(
	     printf("\n** Done ** \n");
  );   
  }
}
#endif

MCDETECTOR mcdetector_out_0D_nexus(MCDETECTOR detector)
{
  /* Write data set information to NeXus file. */
  MPI_MASTER(
    mcdatainfo_out_nexus(nxhandle, detector);
  );
  
  return(detector);
} /* mcdetector_out_0D_ascii */

MCDETECTOR mcdetector_out_1D_nexus(MCDETECTOR detector)
{
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  return(detector);
} /* mcdetector_out_1D_ascii */

MCDETECTOR mcdetector_out_2D_nexus(MCDETECTOR detector)
{
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  
#ifdef USE_MPI // and USE_NEXUS
  /* NeXus: slave nodes have master write their lists */
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    mcdetector_out_list_slaves(detector);
  }
#endif /* USE_MPI */

  return(detector);
} /* mcdetector_out_2D_nexus */

#endif /* USE_NEXUS*/








/* ========================================================================== */

/*                            Main input functions                            */
/*            DETECTOR_OUT_xD function calls -> ascii or NeXus                */

/* ========================================================================== */

/*******************************************************************************
* mcsiminfo_init:   open SIM and write header
*******************************************************************************/
FILE *mcsiminfo_init(FILE *f)
{
  int exists=0;
  int index;
  
  /* check format */      
  if (!mcformat || !strlen(mcformat) 
   || !strcasecmp(mcformat, "MCSTAS") || !strcasecmp(mcformat, "MCXTRACE") 
   || !strcasecmp(mcformat, "PGPLOT") || !strcasecmp(mcformat, "GNUPLOT") || !strcasecmp(mcformat, "MCCODE")
   || !strcasecmp(mcformat, "MATLAB")) {
    mcformat="McCode";
#ifdef USE_NEXUS
  } else if (strcasestr(mcformat, "NeXus")) {
    /* Do nothing */
#endif
  } else {
    fprintf(stderr,
	    "Warning: You have requested the output format %s which is unsupported by this binary. Resetting to standard %s format.\n",mcformat ,"McCode");
    mcformat="McCode";
  }
  
  /* open the SIM file if not defined yet */
  if (mcsiminfo_file || mcdisable_output_files) 
    return (mcsiminfo_file);
    
#ifdef USE_NEXUS
  /* only master writes NeXus header: calls NXopen(nxhandle) */
  if (mcformat && strcasestr(mcformat, "NeXus")) {
	  MPI_MASTER(
	  mcsiminfo_file = mcnew_file(mcsiminfo_name, "h5", &exists);
    if(!mcsiminfo_file)
      fprintf(stderr,
	      "Warning: could not open simulation description file '%s'\n",
	      mcsiminfo_name);
	  else
	    mcinfo_out_nexus(nxhandle);
	  );
    return(mcsiminfo_file); /* points to nxhandle */
  }
#endif
  
  /* write main description file (only MASTER) */
  MPI_MASTER(

  mcsiminfo_file = mcnew_file(mcsiminfo_name, "sim", &exists);
  if(!mcsiminfo_file)
    fprintf(stderr,
	    "Warning: could not open simulation description file '%s'\n",
	    mcsiminfo_name);
  else
  {
    /* write SIM header */
    time_t t=time(NULL);
    mcsiminfo_out("%s simulation description file for %s.\n", 
      MCCODE_NAME, mcinstrument_name);
    mcsiminfo_out("Date:    %s", ctime(&t)); /* includes \n */
    mcsiminfo_out("Program: %s\n\n", MCCODE_STRING);
    
    mcsiminfo_out("begin instrument: %s\n", mcinstrument_name);
    mcinfo_out(   "  ", mcsiminfo_file);
    mcsiminfo_out("end instrument\n");

    mcsiminfo_out("\nbegin simulation: %s\n", mcdirname);
    mcruninfo_out("  ", mcsiminfo_file);
    mcsiminfo_out("end simulation\n");

  }
  return (mcsiminfo_file);
  
  ); /* MPI_MASTER */
  
} /* mcsiminfo_init */

/*******************************************************************************
*   mcsiminfo_close:  close SIM
*******************************************************************************/
void mcsiminfo_close()
{
  MPI_MASTER(
  if(mcsiminfo_file && !mcdisable_output_files) {
#ifdef USE_NEXUS
    if (mcformat && strcasestr(mcformat, "NeXus")) {
      time_t t=time(NULL);
      nxprintf(nxhandle, "end_time", ctime(&t));
      nxprintf(nxhandle, "duration", "%li", (long)t-mcstartdate);
      NXclosegroup(nxhandle); /* NXentry */
      NXclose(&nxhandle);
    } else
#endif
      fclose(mcsiminfo_file);
    );
    mcsiminfo_file = NULL;
  }
} /* mcsiminfo_close */

/*******************************************************************************
* mcdetector_out_0D: wrapper for 0D (single value).
*   Output single detector/monitor data (p0, p1, p2).
*   Title is t, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2,
                         char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI reduce) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : MCCODE_STRING " data"),
    1, 1, 1,
    "I", "", "",
    "I", "", "",
    0, 0, 0, 0, 0, 0, "",
    &p0, &p1, &p2, posa); /* write Detector: line */

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_0D_nexus(detector));
  else
#endif
    return(mcdetector_out_0D_ascii(detector));
    
} /* mcdetector_out_0D */



/*******************************************************************************
* mcdetector_out_1D: wrapper for 1D.
*   Output 1d detector data (p0, p1, p2) for n bins linearly
*   distributed across the range x1..x2 (x1 is lower limit of first
*   bin, x2 is upper limit of last bin). Title is t, axis labels are xl
*   and yl. File name is f, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
        char *xvar, double x1, double x2,
        long n,
        double *p0, double *p1, double *p2, char *f,
        char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : MCCODE_STRING " 1D data"),
    n, 1, 1,
    xl, yl, (n > 1 ? "Signal per bin" : " Signal"),
    xvar, "(I,I_err)", "I",
    x1, x2, 0, 0, 0, 0, f,
    p0, p1, p2, posa); /* write Detector: line */
  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_1D_nexus(detector));
  else
#endif
    return(mcdetector_out_1D_ascii(detector));
  
} /* mcdetector_out_1D */

/*******************************************************************************
* mcdetector_out_2D: wrapper for 2D.
*   special case for list: master creates file first, then slaves append their blocks without header
*******************************************************************************/
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2,
                  long m, long n,
                  double *p0, double *p1, double *p2, char *f,
                  char *c, Coords posa)
{
  char xvar[CHAR_BUF_LENGTH];
  char yvar[CHAR_BUF_LENGTH];
  
  /* create short axes labels */
  if (xl && strlen(xl)) { strncpy(xvar, xl, CHAR_BUF_LENGTH); xvar[2]='\0'; }
  else strcpy(xvar, "x");
  if (yl && strlen(yl)) { strncpy(yvar, yl, CHAR_BUF_LENGTH); yvar[2]='\0'; }
  else strcpy(yvar, "y");

  MCDETECTOR detector;

  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  if (abs(m) == 1) {/* n>1 on Y, m==1 on X: 1D, no X axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      n, 1, 1,
      yl, "", "Signal per bin",
      yvar, "(I,Ierr)", "I",
      y1, y2, x1, x2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  } else if (abs(n)==1) {/* m>1 on X, n==1 on Y: 1D, no Y axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      m, 1, 1,
      xl, "", "Signal per bin",
      xvar, "(I,Ierr)", "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }else {
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 2D data"),
      m, n, 1,
      xl, yl, "Signal per bin",
      xvar, yvar, "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }

  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_2D_nexus(detector));
  else
#endif
    return(mcdetector_out_2D_ascii(detector));
  
} /* mcdetector_out_2D */

/*******************************************************************************
* mcdetector_out_list: wrapper for list output (calls out_2D with mcformat+"list").
*   m=number of events, n=size of each event
*******************************************************************************/
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa)
{
  char       format_new[CHAR_BUF_LENGTH];
  char      *format_org;
  MCDETECTOR detector;
  
  format_org = mcformat;
  strcpy(format_new, mcformat);
  strcat(format_new, " list");
  mcformat = format_new;

  detector = mcdetector_out_2D(t, xl, yl,
                  1,abs(m),1,abs(n),
                  m,n,
                  NULL, p1, NULL, f,
                  c, posa);
  
  mcformat = format_org;
  return(detector);
}

/*******************************************************************************
 * mcuse_dir: set data/sim storage directory and create it,
 * or exit with error if exists
 ******************************************************************************/
static void
mcuse_dir(char *dir)
{
  if (!dir || !strlen(dir)) return;
#ifdef MC_PORTABLE
  fprintf(stderr, "Error: "
          "Directory output cannot be used with portable simulation (mcuse_dir)\n");
  exit(1);
#else  /* !MC_PORTABLE */
  /* handle file://directory URL type */
  if (strncmp(dir, "file://", strlen("file://")))
    mcdirname = dir;
  else
    mcdirname = dir+strlen("file://");
  
  
  
  MPI_MASTER(
    if(mkdir(mcdirname, 0777)) {
#ifndef DANSE
      fprintf(stderr, "Error: unable to create directory '%s' (mcuse_dir)\n", dir);
      fprintf(stderr, "(Maybe the directory already exists?)\n");
#endif
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    exit(-1);
    }
  ); /* MPI_MASTER */
  
  /* remove trailing PATHSEP (if any) */
  while (strlen(mcdirname) && mcdirname[strlen(mcdirname) - 1] == MC_PATHSEP_C)
    mcdirname[strlen(mcdirname) - 1]='\0';
#endif /* !MC_PORTABLE */
} /* mcuse_dir */

/*******************************************************************************
* mcinfo: display instrument simulation info to stdout and exit
*******************************************************************************/
static void
mcinfo(void)
{
  fprintf(stdout, "begin instrument: %s\n", mcinstrument_name);
  mcinfo_out("  ", stdout);
  fprintf(stdout, "end instrument\n");
  fprintf(stdout, "begin simulation: %s\n", mcdirname ? mcdirname : ".");
  mcruninfo_out("  ", stdout);
  fprintf(stdout, "end simulation\n");
  exit(0); /* includes MPI_Finalize in MPI mode */
} /* mcinfo */

#endif /* ndef MCCODE_R_IO_C */

/* end of the I/O section =================================================== */







/*******************************************************************************
* mcset_ncount: set total number of rays to generate
*******************************************************************************/
void mcset_ncount(unsigned long long int count)
{
  mcncount = count;
}

/* mcget_ncount: get total number of rays to generate */
unsigned long long int mcget_ncount(void)
{
  return mcncount;
}

/* mcget_run_num: get curent number of rays in TRACE */
unsigned long long int mcget_run_num(void)
{
  return mcrun_num;
}

/* mcsetn_arg: get ncount from a string argument */
static void
mcsetn_arg(char *arg)
{
  mcset_ncount((long long int) strtod(arg, NULL));
}

/* mcsetseed: set the random generator seed from a string argument */
static void
mcsetseed(char *arg)
{
  mcseed = atol(arg);
  if(mcseed) {
    srandom(mcseed);
  } else {
    fprintf(stderr, "Error: seed must not be zero (mcsetseed)\n");
    exit(1);
  }
}

/* Following part is only embedded when not redundent with mccode-r.h ========= */

#ifndef MCCODE_H

/* SECTION: MCDISPLAY support. =============================================== */

/*******************************************************************************
* Just output MCDISPLAY keywords to be caught by an external plotter client.
*******************************************************************************/

void mcdis_magnify(char *what){
  printf("MCDISPLAY: magnify('%s')\n", what);
}

void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2){
  printf("MCDISPLAY: multiline(2,%g,%g,%g,%g,%g,%g)\n",
         x1,y1,z1,x2,y2,z2);
}

void mcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n){
  int i;
  const double dx = (x2-x1)/(2*n+1);
  const double dy = (y2-y1)/(2*n+1);
  const double dz = (z2-z1)/(2*n+1);

  for(i = 0; i < n+1; i++)
    mcdis_line(x1 + 2*i*dx,     y1 + 2*i*dy,     z1 + 2*i*dz,
	       x1 + (2*i+1)*dx, y1 + (2*i+1)*dy, z1 + (2*i+1)*dz);
}

void mcdis_multiline(int count, ...){
  va_list ap;
  double x,y,z;

  printf("MCDISPLAY: multiline(%d", count);
  va_start(ap, count);
  while(count--)
    {
    x = va_arg(ap, double);
    y = va_arg(ap, double);
    z = va_arg(ap, double);
    printf(",%g,%g,%g", x, y, z);
    }
  va_end(ap);
  printf(")\n");
}

void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height){
  /* draws a rectangle in the plane           */
  /* x is ALWAYS width and y is ALWAYS height */
  if (strcmp("xy", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y - height/2, z,
		    x + width/2, y - height/2, z,
		    x + width/2, y + height/2, z,
		    x - width/2, y + height/2, z,
		    x - width/2, y - height/2, z);
  } else if (strcmp("xz", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y, z - height/2,
		    x + width/2, y, z - height/2,
		    x + width/2, y, z + height/2,
		    x - width/2, y, z + height/2,
		    x - width/2, y, z - height/2);
  } else if (strcmp("yz", plane)==0) {
    mcdis_multiline(5,
		    x, y - height/2, z - width/2,
		    x, y - height/2, z + width/2,
		    x, y + height/2, z + width/2,
		    x, y + height/2, z - width/2,
		    x, y - height/2, z - width/2);
  } else {

    fprintf(stderr, "Error: Definition of plane %s unknown\n", plane);
    exit(1);
  }
}

/*  draws a box with center at (x, y, z) and
    width (deltax), height (deltay), length (deltaz) */
void mcdis_box(double x, double y, double z,
	       double width, double height, double length){

  mcdis_rectangle("xy", x, y, z-length/2, width, height);
  mcdis_rectangle("xy", x, y, z+length/2, width, height);
  mcdis_line(x-width/2, y-height/2, z-length/2,
	     x-width/2, y-height/2, z+length/2);
  mcdis_line(x-width/2, y+height/2, z-length/2,
	     x-width/2, y+height/2, z+length/2);
  mcdis_line(x+width/2, y-height/2, z-length/2,
	     x+width/2, y-height/2, z+length/2);
  mcdis_line(x+width/2, y+height/2, z-length/2,
	     x+width/2, y+height/2, z+length/2);
}

void mcdis_circle(char *plane, double x, double y, double z, double r){
  printf("MCDISPLAY: circle('%s',%g,%g,%g,%g)\n", plane, x, y, z, r);
}

/* SECTION: coordinates handling ============================================ */

/*******************************************************************************
* Since we use a lot of geometric calculations using Cartesian coordinates,
* we collect some useful routines here. However, it is also permissible to
* work directly on the underlying struct coords whenever that is most
* convenient (that is, the type Coords is not abstract).
*
* Coordinates are also used to store rotation angles around x/y/z axis.
*
* Since coordinates are used much like a basic type (such as double), the
* structure itself is passed and returned, rather than a pointer.
*
* At compile-time, the values of the coordinates may be unknown (for example
* a motor position). Hence coordinates are general expressions and not simple
* numbers. For this we used the type Coords_exp which has three CExp
* fields. For runtime (or calculations possible at compile time), we use
* Coords which contains three double fields.
*******************************************************************************/

/* coords_set: Assign coordinates. */
Coords
coords_set(MCNUM x, MCNUM y, MCNUM z)
{
  Coords a;

  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

/* coords_get: get coordinates. Required when 'x','y','z' are #defined as ray pars */
Coords
coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z)
{
  *x = a.x;
  *y = a.y;
  *z = a.z;
  return a;
}

/* coords_add: Add two coordinates. */
Coords
coords_add(Coords a, Coords b)
{
  Coords c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_sub: Subtract two coordinates. */
Coords
coords_sub(Coords a, Coords b)
{
  Coords c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_neg: Negate coordinates. */
Coords
coords_neg(Coords a)
{
  Coords b;

  b.x = -a.x;
  b.y = -a.y;
  b.z = -a.z;
  return b;
}

/* coords_scale: Scale a vector. */
Coords coords_scale(Coords b, double scale) {
  Coords a;

  a.x = b.x*scale;
  a.y = b.y*scale;
  a.z = b.z*scale;
  return a;
}

/* coords_sp: Scalar product: a . b */
double coords_sp(Coords a, Coords b) {
  double value;

  value = a.x*b.x + a.y*b.y + a.z*b.z;
  return value;
}

/* coords_xp: Cross product: a = b x c. */
Coords coords_xp(Coords b, Coords c) {
  Coords a;

  a.x = b.y*c.z - c.y*b.z;
  a.y = b.z*c.x - c.z*b.x;
  a.z = b.x*c.y - c.x*b.y;
  return a;
}

/* coords_len: Gives length of coords set. */
double coords_len(Coords a) {
  return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

/* coords_mirror: Mirror a in plane (through the origin) defined by normal n*/
Coords coords_mirror(Coords a, Coords n) {
  double t = scalar_prod(n.x, n.y, n.z, n.x, n.y, n.z);
  Coords b;
  if (t!=1) {
    t = sqrt(t);
    n.x /= t;
    n.y /= t;
    n.z /= t;
  }
  t=scalar_prod(a.x, a.y, a.z, n.x, n.y, n.z);
  b.x = a.x-2*t*n.x;
  b.y = a.y-2*t*n.y;
  b.z = a.z-2*t*n.z;
  return b;
}

/* coords_print: Print out vector values. */
void coords_print(Coords a) {

  fprintf(stdout, "(%f, %f, %f)\n", a.x, a.y, a.z);
  return;
}

mcstatic inline void coords_norm(Coords* c) {
	double temp = coords_sp(*c,*c);

	// Skip if we will end dividing by zero
	if (temp == 0) return;

	temp = sqrt(temp);

	c->x /= temp;
	c->y /= temp;
	c->z /= temp;
}

/*******************************************************************************
* The Rotation type implements a rotation transformation of a coordinate
* system in the form of a double[3][3] matrix.
*
* Contrary to the Coords type in coords.c, rotations are passed by
* reference. Functions that yield new rotations do so by writing to an
* explicit result parameter; rotations are not returned from functions. The
* reason for this is that arrays cannot by returned from functions (though
* structures can; thus an alternative would have been to wrap the
* double[3][3] array up in a struct). Such are the ways of C programming.
*
* A rotation represents the tranformation of the coordinates of a vector when
* changing between coordinate systems that are rotated with respect to each
* other. For example, suppose that coordinate system Q is rotated 45 degrees
* around the Z axis with respect to coordinate system P. Let T be the
* rotation transformation representing a 45 degree rotation around Z. Then to
* get the coordinates of a vector r in system Q, apply T to the coordinates
* of r in P. If r=(1,0,0) in P, it will be (sqrt(1/2),-sqrt(1/2),0) in
* Q. Thus we should be careful when interpreting the sign of rotation angles:
* they represent the rotation of the coordinate systems, not of the
* coordinates (which has opposite sign).
*******************************************************************************/

/*******************************************************************************
* rot_set_rotation: Get transformation for rotation first phx around x axis,
* then phy around y, then phz around z.
*******************************************************************************/
void
rot_set_rotation(Rotation t, double phx, double phy, double phz)
{
  if ((phx == 0) && (phy == 0) && (phz == 0)) {
    t[0][0] = 1.0;
    t[0][1] = 0.0;
    t[0][2] = 0.0;
    t[1][0] = 0.0;
    t[1][1] = 1.0;
    t[1][2] = 0.0;
    t[2][0] = 0.0;
    t[2][1] = 0.0;
    t[2][2] = 1.0;
  } else {
    double cx = cos(phx);
    double sx = sin(phx);
    double cy = cos(phy);
    double sy = sin(phy);
    double cz = cos(phz);
    double sz = sin(phz);

    t[0][0] = cy*cz;
    t[0][1] = sx*sy*cz + cx*sz;
    t[0][2] = sx*sz - cx*sy*cz;
    t[1][0] = -cy*sz;
    t[1][1] = cx*cz - sx*sy*sz;
    t[1][2] = sx*cz + cx*sy*sz;
    t[2][0] = sy;
    t[2][1] = -sx*cy;
    t[2][2] = cx*cy;
  }
}

/*******************************************************************************
* rot_test_identity: Test if rotation is identity
*******************************************************************************/
int
rot_test_identity(Rotation t)
{
  return (t[0][0] + t[1][1] + t[2][2] == 3);
}

/*******************************************************************************
* rot_mul: Matrix multiplication of transformations (this corresponds to
* combining transformations). After rot_mul(T1, T2, T3), doing T3 is
* equal to doing first T2, then T1.
* Note that T3 must not alias (use the same array as) T1 or T2.
*******************************************************************************/
void
rot_mul(Rotation t1, Rotation t2, Rotation t3)
{
  if (rot_test_identity(t1)) {
    rot_copy(t3, t2);
  } else if (rot_test_identity(t2)) {
    rot_copy(t3, t1);
  } else {
    int i,j;
    for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
	t3[i][j] = t1[i][0]*t2[0][j] + t1[i][1]*t2[1][j] + t1[i][2]*t2[2][j];
  }
}

/*******************************************************************************
* rot_copy: Copy a rotation transformation (arrays cannot be assigned in C).
*******************************************************************************/
void
rot_copy(Rotation dest, Rotation src)
{
  int i,j;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      dest[i][j] = src[i][j];
}

/*******************************************************************************
* rot_transpose: Matrix transposition, which is inversion for Rotation matrices
*******************************************************************************/
void
rot_transpose(Rotation src, Rotation dst)
{
  dst[0][0] = src[0][0];
  dst[0][1] = src[1][0];
  dst[0][2] = src[2][0];
  dst[1][0] = src[0][1];
  dst[1][1] = src[1][1];
  dst[1][2] = src[2][1];
  dst[2][0] = src[0][2];
  dst[2][1] = src[1][2];
  dst[2][2] = src[2][2];
}

/*******************************************************************************
* rot_apply: returns t*a
*******************************************************************************/
Coords
rot_apply(Rotation t, Coords a)
{
  Coords b;
  if (rot_test_identity(t)) {
    return a;
  } else {
    b.x = t[0][0]*a.x + t[0][1]*a.y + t[0][2]*a.z;
    b.y = t[1][0]*a.x + t[1][1]*a.y + t[1][2]*a.z;
    b.z = t[2][0]*a.x + t[2][1]*a.y + t[2][2]*a.z;
    return b;
  }
}

/**
 * Pretty-printing of rotation matrices.
 */
void rot_print(Rotation rot) {
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[0][0], rot[0][1], rot[0][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[1][0], rot[1][1], rot[1][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n\n",
			rot[2][0], rot[2][1], rot[2][2]);
}

/**
 * Vector product: used by vec_prod (mccode-r.h). Use coords_xp for Coords.
 */
mcstatic inline void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
    *x = (y1)*(z2) - (y2)*(z1);
    *y = (z1)*(x2) - (z2)*(x1);
    *z = (x1)*(y2) - (x2)*(y1);
}

/**
 * Scalar product: use coords_sp for Coords.
 */
mcstatic inline double scalar_prod(
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
	return ((x1 * x2) + (y1 * y2) + (z1 * z2));
}

/*******************************************************************************
* mccoordschange: applies rotation to (x y z) and (vx vy vz) and Spin (sx,sy,sz)
*******************************************************************************/
void
mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *x;
  b.y = *y;
  b.z = *z;
  c = rot_apply(t, b);
  b = coords_add(c, a);
  *x = b.x;
  *y = b.y;
  *z = b.z;

  if ( (vz && vy  && vx) && (*vz != 0.0 || *vx != 0.0 || *vy != 0.0) ) mccoordschange_polarisation(t, vx, vy, vz);

  if ( (sz && sy  && sx) && (*sz != 0.0 || *sx != 0.0 || *sy != 0.0) ) mccoordschange_polarisation(t, sx, sy, sz);

}

/*******************************************************************************
* mccoordschange_polarisation: applies rotation to vector (sx sy sz)
*******************************************************************************/
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *sx;
  b.y = *sy;
  b.z = *sz;
  c = rot_apply(t, b);
  *sx = c.x;
  *sy = c.y;
  *sz = c.z;
}

/* SECTION: vector math  ==================================================== */

/* normal_vec_func: Compute normal vector to (x,y,z). */
mcstatic inline void normal_vec_func(double *nx, double *ny, double *nz,
                double x, double y, double z)
{
  double ax = fabs(x);
  double ay = fabs(y);
  double az = fabs(z);
  double l;
  if(x == 0 && y == 0 && z == 0)
  {
    *nx = 0;
    *ny = 0;
    *nz = 0;
    return;
  }
  if(ax < ay)
  {
    if(ax < az)
    {                           /* Use X axis */
      l = sqrt(z*z + y*y);
      *nx = 0;
      *ny = z/l;
      *nz = -y/l;
      return;
    }
  }
  else
  {
    if(ay < az)
    {                           /* Use Y axis */
      l = sqrt(z*z + x*x);
      *nx = z/l;
      *ny = 0;
      *nz = -x/l;
      return;
    }
  }
  /* Use Z axis */
  l = sqrt(y*y + x*x);
  *nx = y/l;
  *ny = -x/l;
  *nz = 0;
} /* normal_vec */

/*******************************************************************************
 * solve_2nd_order: second order equation solve: A*t^2 + B*t + C = 0
 * solve_2nd_order(&t1, NULL, A,B,C)
 *   returns 0 if no solution was found, or set 't1' to the smallest positive
 *   solution.
 * solve_2nd_order(&t1, &t2, A,B,C)
 *   same as with &t2=NULL, but also returns the second solution.
 * EXAMPLE usage for intersection of a trajectory with a plane in gravitation
 * field (gx,gy,gz):
 * The neutron starts at point r=(x,y,z) with velocityv=(vx vy vz). The plane
 * has a normal vector n=(nx,ny,nz) and contains the point W=(wx,wy,wz).
 * The problem consists in solving the 2nd order equation:
 *      1/2.n.g.t^2 + n.v.t + n.(r-W) = 0
 * so that A = 0.5 n.g; B = n.v; C = n.(r-W);
 * Without acceleration, t=-n.(r-W)/n.v
 ******************************************************************************/
int solve_2nd_order(double *t1, double *t2,
                  double A,  double B,  double C)
{
  int ret=0;

  if (!t1) return 0;
  *t1 = 0;
  if (t2) *t2=0;

  if (fabs(A) < 1E-10) /* approximate to linear equation: A ~ 0 */
  {
    if (B) {  *t1 = -C/B; ret=1; if (t2) *t2=*t1; }
    /* else no intersection: A=B=0 ret=0 */
  }
  else
  {
    double D;
    D = B*B - 4*A*C;
    if (D >= 0) /* Delta > 0: two solutions */
    {
      double sD, dt1, dt2;
      sD = sqrt(D);
      dt1 = (-B + sD)/2/A;
      dt2 = (-B - sD)/2/A;
      /* we identify very small values with zero */
      if (fabs(dt1) < 1e-10) dt1=0.0;
      if (fabs(dt2) < 1e-10) dt2=0.0;

      /* now we choose the smallest positive solution */
      if      (dt1<=0.0 && dt2>0.0) ret=2; /* dt2 positive */
      else if (dt2<=0.0 && dt1>0.0) ret=1; /* dt1 positive */
      else if (dt1> 0.0 && dt2>0.0)
      {  if (dt1 < dt2) ret=1; else ret=2; } /* all positive: min(dt1,dt2) */
      /* else two solutions are negative. ret=-1 */
      if (ret==1) { *t1 = dt1;  if (t2) *t2=dt2; }
      else        { *t1 = dt2;  if (t2) *t2=dt1; }
      ret=2;  /* found 2 solutions and t1 is the positive one */
    } /* else Delta <0: no intersection. ret=0 */
  }
  return(ret);
} /* solve_2nd_order */

/*******************************************************************************
 * randvec_target_circle: Choose random direction towards target at (x,y,z)
 * with given radius.
 * If radius is zero, choose random direction in full 4PI, no target.
 ******************************************************************************/
void
randvec_target_circle(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double radius)
{
  double l2, phi, theta, nx, ny, nz, xt, yt, zt, xu, yu, zu;

  if(radius == 0.0)
  {
    /* No target, choose uniformly a direction in full 4PI solid angle. */
    theta = acos (1 - rand0max(2));
    phi = rand0max(2 * PI);
    if(solid_angle)
      *solid_angle = 4*PI;
    nx = 1;
    ny = 0;
    nz = 0;
    yi = sqrt(xi*xi+yi*yi+zi*zi);
    zi = 0;
    xi = 0;
  }
  else
  {
    double costheta0;
    l2 = xi*xi + yi*yi + zi*zi; /* sqr Distance to target. */
    costheta0 = sqrt(l2/(radius*radius+l2));
    if (radius < 0) costheta0 *= -1;
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
        *solid_angle = 2*PI*(1 - costheta0);
    }

    /* Now choose point uniformly on circle surface within angle theta0 */
    theta = acos (1 - rand0max(1 - costheta0)); /* radius on circle */
    phi = rand0max(2 * PI); /* rotation on circle at given radius */
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around a
       perpendicular axis u=i x n and then angle phi around i. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, xu, yu, zu);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xi, yi, zi);
} /* randvec_target_circle */

/*******************************************************************************
 * randvec_target_rect_angular: Choose random direction towards target at
 * (xi,yi,zi) with given ANGULAR dimension height x width. height=phi_x=[0,PI],
 * width=phi_y=[0,2*PI] (radians)
 * If height or width is zero, choose random direction in full 4PI, no target.
 *******************************************************************************/
void
randvec_target_rect_angular(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double width, double height, Rotation A)
{
  double theta, phi, nx, ny, nz, xt, yt, zt, xu, yu, zu;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
      *solid_angle = 2*fabs(width*sin(height/2));
    }

    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Now choose point uniformly on the unit sphere segment with angle theta/phi */
    phi   = width*randpm1()/2.0;
    theta = asin(randpm1()*sin(height/2.0));
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around
       n, and then phi around u. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, nx, ny, nz);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xu,  yu,  zu);

  /* Go back to local coordinate system */
  tmp = coords_set(*xo, *yo, *zo);
  tmp = rot_apply(A, tmp);
  coords_get(tmp, &*xo, &*yo, &*zo);

} /* randvec_target_rect_angular */

/*******************************************************************************
 * randvec_target_rect_real: Choose random direction towards target at (xi,yi,zi)
 * with given dimension height x width (in meters !).
 *
 * Local emission coordinate is taken into account and corrected for 'order' times.
 * (See remarks posted to mcstas-users by George Apostolopoulus <gapost@ipta.demokritos.gr>)
 *
 * If height or width is zero, choose random direction in full 4PI, no target.
 *
 * Traditionally, this routine had the name randvec_target_rect - this is now a
 * a define (see mcstas-r.h) pointing here. If you use the old rouine, you are NOT
 * taking the local emmission coordinate into account.
*******************************************************************************/

void
randvec_target_rect_real(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi,
               double width, double height, Rotation A,
               double lx, double ly, double lz, int order)
{
  double dx, dy, dist, dist_p, nx, ny, nz, mx, my, mz, n_norm, m_norm;
  double cos_theta;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {

    /* Now choose point uniformly on rectangle within width x height */
    dx = width*randpm1()/2.0;
    dy = height*randpm1()/2.0;

    /* Determine distance to target plane*/
    dist = sqrt(xi*xi + yi*yi + zi*zi);
    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Determine vector normal to trajectory axis (z) and gravity [0 1 0] */
    vec_prod(nx, ny, nz, xi, yi, zi, 0, 1, 0);

    /* This now defines the x-axis, normalize: */
    n_norm=sqrt(nx*nx + ny*ny + nz*nz);
    nx = nx/n_norm;
    ny = ny/n_norm;
    nz = nz/n_norm;

    /* Now, determine our y-axis (vertical in many cases...) */
    vec_prod(mx, my, mz, xi, yi, zi, nx, ny, nz);
    m_norm=sqrt(mx*mx + my*my + mz*mz);
    mx = mx/m_norm;
    my = my/m_norm;
    mz = mz/m_norm;

    /* Our output, random vector can now be defined by linear combination: */

    *xo = xi + dx * nx + dy * mx;
    *yo = yi + dx * ny + dy * my;
    *zo = zi + dx * nz + dy * mz;

    /* Go back to local coordinate system */
    tmp = coords_set(*xo, *yo, *zo);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &*xo, &*yo, &*zo);

    /* Go back to local coordinate system */
    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    if (solid_angle) {
      /* Calculate vector from local point to remote random point */
      lx = *xo - lx;
      ly = *yo - ly;
      lz = *zo - lz;
      dist_p = sqrt(lx*lx + ly*ly + lz*lz);

      /* Adjust the 'solid angle' */
      /* 1/r^2 to the chosen point times cos(\theta) between the normal */
      /* vector of the target rectangle and direction vector of the chosen point. */
      cos_theta = (xi * lx + yi * ly + zi * lz) / (dist * dist_p);
      *solid_angle = width * height / (dist_p * dist_p);
      int counter;
      for (counter = 0; counter < order; counter++) {
	*solid_angle = *solid_angle * cos_theta;
      }
    }
  }
} /* randvec_target_rect_real */

/* SECTION: random numbers ================================================== */

/*
 * Copyright (c) 1983 Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the University of California, Berkeley.  The name of the
 * University may not be used to endorse or promote products derived
 * from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 * This is derived from the Berkeley source:
 *        @(#)random.c        5.5 (Berkeley) 7/6/88
 * It was reworked for the GNU C Library by Roland McGrath.
 * Rewritten to use reentrant functions by Ulrich Drepper, 1995.
 */

/*******************************************************************************
* Modified for McStas from glibc 2.0.7pre1 stdlib/random.c and
* stdlib/random_r.c.
*
* This way random() is more than four times faster compared to calling
* standard glibc random() on ix86 Linux, probably due to multithread support,
* ELF shared library overhead, etc. It also makes McStas generated
* simulations more portable (more likely to behave identically across
* platforms, important for parrallel computations).
*******************************************************************************/


#define        TYPE_3                3
#define        BREAK_3                128
#define        DEG_3                31
#define        SEP_3                3

static mc_int32_t randtbl[DEG_3 + 1] =
  {
    TYPE_3,

    -1726662223, 379960547, 1735697613, 1040273694, 1313901226,
    1627687941, -179304937, -2073333483, 1780058412, -1989503057,
    -615974602, 344556628, 939512070, -1249116260, 1507946756,
    -812545463, 154635395, 1388815473, -1926676823, 525320961,
    -1009028674, 968117788, -123449607, 1284210865, 435012392,
    -2017506339, -911064859, -370259173, 1132637927, 1398500161,
    -205601318,
  };

static mc_int32_t *fptr = &randtbl[SEP_3 + 1];
static mc_int32_t *rptr = &randtbl[1];
static mc_int32_t *state = &randtbl[1];
#define rand_deg DEG_3
#define rand_sep SEP_3
static mc_int32_t *end_ptr = &randtbl[sizeof (randtbl) / sizeof (randtbl[0])];

mc_int32_t
mc_random (void)
{
  mc_int32_t result;

  *fptr += *rptr;
  /* Chucking least random bit.  */
  result = (*fptr >> 1) & 0x7fffffff;
  ++fptr;
  if (fptr >= end_ptr)
  {
    fptr = state;
    ++rptr;
  }
  else
  {
    ++rptr;
    if (rptr >= end_ptr)
      rptr = state;
  }
  return result;
}

void
mc_srandom (unsigned int x)
{
  /* We must make sure the seed is not 0.  Take arbitrarily 1 in this case.  */
  state[0] = x ? x : 1;
  {
    long int i;
    for (i = 1; i < rand_deg; ++i)
    {
      /* This does:
         state[i] = (16807 * state[i - 1]) % 2147483647;
         but avoids overflowing 31 bits.  */
      long int hi = state[i - 1] / 127773;
      long int lo = state[i - 1] % 127773;
      long int test = 16807 * lo - 2836 * hi;
      state[i] = test + (test < 0 ? 2147483647 : 0);
    }
    fptr = &state[rand_sep];
    rptr = &state[0];
    for (i = 0; i < 10 * rand_deg; ++i)
      random ();
  }
}

/* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
/* See http://www.math.keio.ac.jp/~matumoto/emt.html for original source. */


/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using mt_srandom(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

#include <stdio.h>

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void mt_srandom(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;
    mt_srandom(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long mt_random(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if mt_srandom() has not been called, */
            mt_srandom(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK

/* End of "Mersenne Twister". */

/* End of McCode random number routine. */

/* randnorm: generate a random number from normal law */
double
randnorm(void)
{
  static double v1, v2, s;
  static int phase = 0;
  double X, u1, u2;

  if(phase == 0)
  {
    do
    {
      u1 = rand01();
      u2 = rand01();
      v1 = 2*u1 - 1;
      v2 = 2*u2 - 1;
      s = v1*v1 + v2*v2;
    } while(s >= 1 || s == 0);

    X = v1*sqrt(-2*log(s)/s);
  }
  else
  {
    X = v2*sqrt(-2*log(s)/s);
  }

  phase = 1 - phase;
  return X;
}

/**
 * Generate a random number from -1 to 1 with triangle distribution
 */
double randtriangle(void) {
	double randnum = rand01();
	if (randnum>0.5) return(1-sqrt(2*(randnum-0.5)));
	else return(sqrt(2*randnum)-1);
}

/**
 * Random number between 0.0 and 1.0 (including?)
 */
double rand01() {
	double randnum;
	randnum = (double) random();
	randnum /= (double) MC_RAND_MAX + 1;
	return randnum;
}

/**
 * Return a random number between 1 and -1
 */
double randpm1() {
	double randnum;
	randnum = (double) random();
	randnum /= ((double) MC_RAND_MAX + 1) / 2;
	randnum -= 1;
	return randnum;
}

/**
 * Return a random number between 0 and max.
 */
double rand0max(double max) {
	double randnum;
	randnum = (double) random();
	randnum /= ((double) MC_RAND_MAX + 1) / max;
	return randnum;
}

/**
 * Return a random number between min and max.
 */
double randminmax(double min, double max) {
	return rand0max(max - min) + max;
}

/* SECTION: main and signal handlers ======================================== */

/*******************************************************************************
* mchelp: displays instrument executable help with possible options
*******************************************************************************/
static void
mchelp(char *pgmname)
{
  int i;

  fprintf(stderr, "%s (%s) instrument simulation, generated with " MCCODE_STRING " (" MCCODE_DATE ")\n", mcinstrument_name, mcinstrument_source);
  fprintf(stderr, "Usage: %s [options] [parm=value ...]\n", pgmname);
  fprintf(stderr,
"Options are:\n"
"  -s SEED   --seed=SEED      Set random seed (must be != 0)\n"
"  -n COUNT  --ncount=COUNT   Set number of " MCCODE_PARTICLE "s to simulate.\n"
"  -d DIR    --dir=DIR        Put all data files in directory DIR.\n"
"  -t        --trace          Enable trace of " MCCODE_PARTICLE "s through instrument.\n"
"  -g        --gravitation    Enable gravitation for all trajectories.\n"
"  --no-output-files          Do not write any data files.\n"
"  -h        --help           Show this help message.\n"
"  -i        --info           Detailed instrument information.\n"
"  --format=FORMAT            Output data files using FORMAT="
   FLAVOR_UPPER
#ifdef USE_NEXUS
   " NEXUS"
#endif
"\n\n"
);
#ifdef USE_MPI
  fprintf(stderr,
  "This instrument has been compiled with MPI support.\n  Use 'mpirun %s [options] [parm=value ...]'.\n", pgmname);
#endif
  if(mcnumipar > 0)
  {
    fprintf(stderr, "Instrument parameters are:\n");
    for(i = 0; i < mcnumipar; i++)
      if (mcinputtable[i].val && strlen(mcinputtable[i].val))
        fprintf(stderr, "  %-16s(%s) [default='%s']\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name),
        mcinputtable[i].val);
      else
        fprintf(stderr, "  %-16s(%s)\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name));
  }

#ifndef NOSIGNALS
  fprintf(stderr, "Known signals are: "
#ifdef SIGUSR1
  "USR1 (status) "
#endif
#ifdef SIGUSR2
  "USR2 (save) "
#endif
#ifdef SIGBREAK
  "BREAK (save) "
#endif
#ifdef SIGTERM
  "TERM (save and exit)"
#endif
  "\n");
#endif /* !NOSIGNALS */
} /* mchelp */


/* mcshowhelp: show help and exit with 0 */
static void
mcshowhelp(char *pgmname)
{
  mchelp(pgmname);
  exit(0);
}

/* mcusage: display usage when error in input arguments and exit with 1 */
static void
mcusage(char *pgmname)
{
  fprintf(stderr, "Error: incorrect command line arguments\n");
  mchelp(pgmname);
  exit(1);
}

/* mcenabletrace: enable trace/mcdisplay or error if requires recompile */
static void
mcenabletrace(void)
{
 if(mctraceenabled)
  mcdotrace = 1;
 else
 {
   fprintf(stderr,
           "Error: trace not enabled (mcenabletrace)\n"
           "Please re-run the " MCCODE_NAME " compiler "
                   "with the --trace option, or rerun the\n"
           "C compiler with the MC_TRACE_ENABLED macro defined.\n");
   exit(1);
 }
}

/*******************************************************************************
* mcreadparams: request parameters from the prompt (or use default)
*******************************************************************************/
void
mcreadparams(void)
{
  int i,j,status;
  static char buf[CHAR_BUF_LENGTH];
  char *p;
  int len;

  MPI_MASTER(printf("Instrument parameters for %s (%s)\n",
                    mcinstrument_name, mcinstrument_source));

  for(i = 0; mcinputtable[i].name != 0; i++)
  {
    do
    {
      MPI_MASTER(
                 if (mcinputtable[i].val && strlen(mcinputtable[i].val))
                   printf("Set value of instrument parameter %s (%s) [default='%s']:\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name), mcinputtable[i].val);
                 else
                   printf("Set value of instrument parameter %s (%s):\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name));
                 fflush(stdout);
                 );
#ifdef USE_MPI
      if(mpi_node_rank == mpi_node_root)
        {
          p = fgets(buf, CHAR_BUF_LENGTH, stdin);
          if(p == NULL)
            {
              fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
              exit(1);
            }
        }
      else
        p = buf;
      MPI_Bcast(buf, CHAR_BUF_LENGTH, MPI_CHAR, mpi_node_root, MPI_COMM_WORLD);
#else /* !USE_MPI */
      p = fgets(buf, CHAR_BUF_LENGTH, stdin);
      if(p == NULL)
        {
          fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
          exit(1);
        }
#endif /* USE_MPI */
      len = strlen(buf);
      if (!len || (len == 1 && (buf[0] == '\n' || buf[0] == '\r')))
      {
        if (mcinputtable[i].val && strlen(mcinputtable[i].val)) {
          strncpy(buf, mcinputtable[i].val, CHAR_BUF_LENGTH);  /* use default value */
          len = strlen(buf);
        }
      }
      for(j = 0; j < 2; j++)
      {
        if(len > 0 && (buf[len - 1] == '\n' || buf[len - 1] == '\r'))
        {
          len--;
          buf[len] = '\0';
        }
      }

      status = (*mcinputtypes[mcinputtable[i].type].getparm)
                   (buf, mcinputtable[i].par);
      if(!status)
      {
        (*mcinputtypes[mcinputtable[i].type].error)(mcinputtable[i].name, buf);
        if (!mcinputtable[i].val || strlen(mcinputtable[i].val)) {
          fprintf(stderr, "       Change %s default value in instrument definition.\n", mcinputtable[i].name);
          exit(1);
        }
      }
    } while(!status);
  }
} /* mcreadparams */

/*******************************************************************************
* mcparseoptions: parse command line arguments (options, parameters)
*******************************************************************************/
void
mcparseoptions(int argc, char *argv[])
{
  int i, j;
  char *p;
  int paramset = 0, *paramsetarray;
  char *usedir=NULL;

  /* Add one to mcnumipar to avoid allocating zero size memory block. */
  paramsetarray = (int*)malloc((mcnumipar + 1)*sizeof(*paramsetarray));
  if(paramsetarray == NULL)
  {
    fprintf(stderr, "Error: insufficient memory (mcparseoptions)\n");
    exit(1);
  }
  for(j = 0; j < mcnumipar; j++)
    {
      paramsetarray[j] = 0;
      if (mcinputtable[j].val != NULL && strlen(mcinputtable[j].val))
      {
        int  status;
        char buf[CHAR_BUF_LENGTH];
        strncpy(buf, mcinputtable[j].val, CHAR_BUF_LENGTH);
        status = (*mcinputtypes[mcinputtable[j].type].getparm)
                   (buf, mcinputtable[j].par);
        if(!status) fprintf(stderr, "Invalid '%s' default value %s in instrument definition (mcparseoptions)\n", mcinputtable[j].name, buf);
        else paramsetarray[j] = 1;
      } else {
        (*mcinputtypes[mcinputtable[j].type].getparm)
          (NULL, mcinputtable[j].par);
        paramsetarray[j] = 0;
      }
    }
  for(i = 1; i < argc; i++)
  {
    if(!strcmp("-s", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("-s", argv[i], 2))
      mcsetseed(&argv[i][2]);
    else if(!strcmp("--seed", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("--seed=", argv[i], 7))
      mcsetseed(&argv[i][7]);
    else if(!strcmp("-n", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("-n", argv[i], 2))
      mcsetn_arg(&argv[i][2]);
    else if(!strcmp("--ncount", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("--ncount=", argv[i], 9))
      mcsetn_arg(&argv[i][9]);
    else if(!strcmp("-d", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];  /* will create directory after parsing all arguments (end of this function) */
    else if(!strncmp("-d", argv[i], 2))
      usedir=&argv[i][2];
    else if(!strcmp("--dir", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];
    else if(!strncmp("--dir=", argv[i], 6))
      usedir=&argv[i][6];
    else if(!strcmp("-h", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("--help", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("-i", argv[i])) {
      mcformat=FLAVOR_UPPER;
      mcinfo();
    }
    else if(!strcmp("--info", argv[i]))
      mcinfo();
    else if(!strcmp("-t", argv[i]))
      mcenabletrace();
    else if(!strcmp("--trace", argv[i]))
      mcenabletrace();
    else if(!strcmp("--gravitation", argv[i]))
      mcgravitation = 1;
    else if(!strcmp("-g", argv[i]))
      mcgravitation = 1;
    else if(!strncmp("--format=", argv[i], 9)) {
      mcformat=&argv[i][9];
    }
    else if(!strcmp("--format", argv[i]) && (i + 1) < argc) {
      mcformat=argv[++i];
    }
    else if(!strcmp("--no-output-files", argv[i]))
      mcdisable_output_files = 1;
    else if(argv[i][0] != '-' && (p = strchr(argv[i], '=')) != NULL)
    {
      *p++ = '\0';

      for(j = 0; j < mcnumipar; j++)
        if(!strcmp(mcinputtable[j].name, argv[i]))
        {
          int status;
          status = (*mcinputtypes[mcinputtable[j].type].getparm)(p,
                        mcinputtable[j].par);
          if(!status || !strlen(p))
          {
            (*mcinputtypes[mcinputtable[j].type].error)
              (mcinputtable[j].name, p);
            exit(1);
          }
          paramsetarray[j] = 1;
          paramset = 1;
          break;
        }
      if(j == mcnumipar)
      {                                /* Unrecognized parameter name */
        fprintf(stderr, "Error: unrecognized parameter %s (mcparseoptions)\n", argv[i]);
        exit(1);
      }
    }
    else if(argv[i][0] == '-') {
      fprintf(stderr, "Error: unrecognized option argument %s (mcparseoptions). Ignored.\n", argv[i++]);
    }
    else {
      fprintf(stderr, "Error: unrecognized argument %s (mcparseoptions). Aborting.\n", argv[i]);
      mcusage(argv[0]);
    }
  }
  if(!paramset)
    mcreadparams();                /* Prompt for parameters if not specified. */
  else
  {
    for(j = 0; j < mcnumipar; j++)
      if(!paramsetarray[j])
      {
        fprintf(stderr, "Error: Instrument parameter %s left unset (mcparseoptions)\n",
                mcinputtable[j].name);
        exit(1);
      }
  }
  free(paramsetarray);
#ifdef USE_MPI
  if (mcdotrace) mpi_node_count=1; /* disable threading when in trace mode */
#endif
  if (usedir && strlen(usedir)) mcuse_dir(usedir);
} /* mcparseoptions */

#ifndef NOSIGNALS
mcstatic char  mcsig_message[256];


/*******************************************************************************
* sighandler: signal handler that makes simulation stop, and save results
*******************************************************************************/
void sighandler(int sig)
{
  /* MOD: E. Farhi, Sep 20th 2001: give more info */
  time_t t1, t0;
#define SIG_SAVE 0
#define SIG_TERM 1
#define SIG_STAT 2
#define SIG_ABRT 3

  printf("\n# " MCCODE_STRING ": [pid %i] Signal %i detected", getpid(), sig);
#ifdef USE_MPI
  printf(" [proc %i]", mpi_node_rank);
#endif
#if defined(SIGUSR1) && defined(SIGUSR2) && defined(SIGKILL)
  if (!strcmp(mcsig_message, "sighandler") && (sig != SIGUSR1) && (sig != SIGUSR2))
  {
    printf("\n# Fatal : unrecoverable loop ! Suicide (naughty boy).\n");
    kill(0, SIGKILL); /* kill myself if error occurs within sighandler: loops */
  }
#endif
  switch (sig) {
#ifdef SIGINT
    case SIGINT : printf(" SIGINT (interrupt from terminal, Ctrl-C)"); sig = SIG_TERM; break;
#endif
#ifdef SIGILL
    case SIGILL  : printf(" SIGILL (Illegal instruction)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGFPE
    case SIGFPE  : printf(" SIGFPE (Math Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGSEGV
    case SIGSEGV : printf(" SIGSEGV (Mem Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGTERM
    case SIGTERM : printf(" SIGTERM (Termination)"); sig = SIG_TERM; break;
#endif
#ifdef SIGABRT
    case SIGABRT : printf(" SIGABRT (Abort)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGQUIT
    case SIGQUIT : printf(" SIGQUIT (Quit from terminal)"); sig = SIG_TERM; break;
#endif
#ifdef SIGTRAP
    case SIGTRAP : printf(" SIGTRAP (Trace trap)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGPIPE
    case SIGPIPE : printf(" SIGPIPE (Broken pipe)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGUSR1
    case SIGUSR1 : printf(" SIGUSR1 (Display info)"); sig = SIG_STAT; break;
#endif
#ifdef SIGUSR2
    case SIGUSR2 : printf(" SIGUSR2 (Save simulation)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGHUP
    case SIGHUP  : printf(" SIGHUP (Hangup/update)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGBUS
    case SIGBUS  : printf(" SIGBUS (Bus error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGURG
    case SIGURG  : printf(" SIGURG (Urgent socket condition)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGBREAK
    case SIGBREAK: printf(" SIGBREAK (Break signal, Ctrl-Break)"); sig = SIG_SAVE; break;
#endif
    default : printf(" (look at signal list for signification)"); sig = SIG_ABRT; break;
  }
  printf("\n");
  printf("# Simulation: %s (%s) \n", mcinstrument_name, mcinstrument_source);
  printf("# Breakpoint: %s ", mcsig_message);
  if (strstr(mcsig_message, "Save") && (sig == SIG_SAVE))
    sig = SIG_STAT;
  SIG_MESSAGE("sighandler");
  if (mcget_ncount() == 0)
    printf("(0 %%)\n" );
  else
  {
    printf("%.2f %% (%10.1f/%10.1f)\n", 100.0*mcget_run_num()/mcget_ncount(), 1.0*mcget_run_num(), 1.0*mcget_ncount());
  }
  t0 = (time_t)mcstartdate;
  t1 = time(NULL);
  printf("# Date:      %s", ctime(&t1));
  printf("# Started:   %s", ctime(&t0));

  if (sig == SIG_STAT)
  {
    printf("# " MCCODE_STRING ": Resuming simulation (continue)\n");
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_SAVE)
  {
    printf("# " MCCODE_STRING ": Saving data and resume simulation (continue)\n");
    mcsave(NULL);
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_TERM)
  {
    printf("# " MCCODE_STRING ": Finishing simulation (save results and exit)\n");
    mcfinally();
    exit(0);
  }
  else
  {
    fflush(stdout);
    perror("# Last I/O Error");
    printf("# " MCCODE_STRING ": Simulation stop (abort).\n");
// This portion of the signal handling only works on UNIX
#if defined(__unix__) || defined(__APPLE__)
    signal(sig, SIG_DFL); /* force to use default sighandler now */
    kill(getpid(), sig);  /* and trigger it with the current signal */
#endif
    exit(-1);
  }
#undef SIG_SAVE
#undef SIG_TERM
#undef SIG_STAT
#undef SIG_ABRT

} /* sighandler */
#endif /* !NOSIGNALS */

/*******************************************************************************
* mccode_main: McCode main() function.
*******************************************************************************/
int mccode_main(int argc, char *argv[])
{
/*  double run_num = 0; */
  time_t  t;
#ifdef USE_MPI
  char mpi_node_name[MPI_MAX_PROCESSOR_NAME];
  int  mpi_node_name_len;
#endif /* USE_MPI */

#ifdef MAC
  argc = ccommand(&argv);
#endif

#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_node_count); /* get number of nodes */
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_node_rank);
  MPI_Comm_set_name(MPI_COMM_WORLD, mcinstrument_name);
  MPI_Get_processor_name(mpi_node_name, &mpi_node_name_len);
#endif /* USE_MPI */

t = time(NULL);
mcseed = (long)t+(long)getpid();

#ifdef USE_MPI
/* *** print number of nodes *********************************************** */
  if (mpi_node_count > 1) {
    MPI_MASTER(
    printf("Simulation '%s' (%s): running on %i nodes (master is '%s', MPI version %i.%i).\n",
      mcinstrument_name, mcinstrument_source, mpi_node_count, mpi_node_name, MPI_VERSION, MPI_SUBVERSION);
    );
  }
#endif /* USE_MPI */
  
  mcstartdate = (long)t;  /* set start date before parsing options and creating sim file */

/* *** parse options ******************************************************* */
  SIG_MESSAGE("main (Start)");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;
  mcinstrument_exe = argv[0]; /* store the executable path */
  /* read simulation parameters and options */
  mcparseoptions(argc, argv); /* sets output dir and format */
  
#ifdef USE_MPI
  if (mpi_node_count > 1) {
    /* share the same seed, then adapt random seed for each node */
    MPI_Bcast(&mcseed, 1, MPI_LONG, 0, MPI_COMM_WORLD); /* root sends its seed to slaves */
    mcseed += mpi_node_rank; /* make sure we use different seeds per node */
  }
#endif
  srandom(mcseed);

/* *** install sig handler, but only once !! after parameters parsing ******* */
#ifndef NOSIGNALS
#ifdef SIGQUIT
  if (signal( SIGQUIT ,sighandler) == SIG_IGN)
    signal( SIGQUIT,SIG_IGN);   /* quit (ASCII FS) */
#endif
#ifdef SIGABRT
  if (signal( SIGABRT ,sighandler) == SIG_IGN)
    signal( SIGABRT,SIG_IGN);   /* used by abort, replace SIGIOT in the future */
#endif
#ifdef SIGTERM
  if (signal( SIGTERM ,sighandler) == SIG_IGN)
    signal( SIGTERM,SIG_IGN);   /* software termination signal from kill */
#endif
#ifdef SIGUSR1
  if (signal( SIGUSR1 ,sighandler) == SIG_IGN)
    signal( SIGUSR1,SIG_IGN);   /* display simulation status */
#endif
#ifdef SIGUSR2
  if (signal( SIGUSR2 ,sighandler) == SIG_IGN)
    signal( SIGUSR2,SIG_IGN);
#endif
#ifdef SIGHUP
  if (signal( SIGHUP ,sighandler) == SIG_IGN)
    signal( SIGHUP,SIG_IGN);
#endif
#ifdef SIGILL
  if (signal( SIGILL ,sighandler) == SIG_IGN)
    signal( SIGILL,SIG_IGN);    /* illegal instruction (not reset when caught) */
#endif
#ifdef SIGFPE
  if (signal( SIGFPE ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* floating point exception */
#endif
#ifdef SIGBUS
  if (signal( SIGBUS ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* bus error */
#endif
#ifdef SIGSEGV
  if (signal( SIGSEGV ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);   /* segmentation violation */
#endif
#endif /* !NOSIGNALS */
  mcsiminfo_init(NULL); /* open SIM */
  SIG_MESSAGE("main (Init)");
  mcinit();
#ifndef NOSIGNALS
#ifdef SIGINT
  if (signal( SIGINT ,sighandler) == SIG_IGN)
    signal( SIGINT,SIG_IGN);    /* interrupt (rubout) only after INIT */
#endif
#endif /* !NOSIGNALS */

/* ================ main particle generation/propagation loop ================ */
#if defined (USE_MPI)
  /* sliced Ncount on each MPI node */
  mcncount = mpi_node_count > 1 ?
    floor(mcncount / mpi_node_count) :
    mcncount; /* number of rays per node */
#endif

/* main particle event loop */
while(mcrun_num < mcncount || mcrun_num < mcget_ncount())
  {
#ifndef NEUTRONICS
    mcgenstate();
#endif
    /* old init: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
    mcraytrace();
    mcrun_num++;
  }

#ifdef USE_MPI
 /* merge run_num from MPI nodes */
  if (mpi_node_count > 1) {
  double mcrun_num_double = (double)mcrun_num;
  mc_MPI_Sum(&mcrun_num_double, 1);
  mcrun_num = (unsigned long long)mcrun_num_double;
  }
#endif

/* save/finally executed by master node/thread */
  mcfinally();

#ifdef USE_MPI
  MPI_Finalize();
#endif /* USE_MPI */

  return 0;
} /* mccode_main */

#ifdef NEUTRONICS
/*Main neutronics function steers the McStas calls, initializes parameters etc */
/* Only called in case NEUTRONICS = TRUE */
void neutronics_main_(float *inx, float *iny, float *inz, float *invx, float *invy, float *invz, float *intime, float *insx, float *insy, float *insz, float *inw, float *outx, float *outy, float *outz, float *outvx, float *outvy, float *outvz, float *outtime, float *outsx, float *outsy, float *outsz, float *outwgt)
{

  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  /* External code governs iteration - McStas is iterated once per call to neutronics_main. I.e. below counter must be initiancated for each call to neutronics_main*/
  mcrun_num=0;

  time_t t;
  t = (time_t)mcstartdate;
  mcstartdate = t;  /* set start date before parsing options and creating sim file */
  mcinit();

  /* *** parse options *** */
  SIG_MESSAGE("main (Start)");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;

  /* Set neutron state based on input from neutronics code */
  mcsetstate(*inx,*iny,*inz,*invx,*invy,*invz,*intime,*insx,*insy,*insz,*inw);

  /* main neutron event loop - runs only one iteration */

  //mcstas_raytrace(&mcncount); /* prior to McStas 1.12 */

  mcallowbackprop = 1; //avoid absorbtion from negative dt
  int argc=1;
  char *argv[0];
  int dummy = mccode_main(argc, argv);

  *outx =  mcnx;
  *outy =  mcny;
  *outz =  mcnz;
  *outvx =  mcnvx;
  *outvy =  mcnvy;
  *outvz =  mcnvz;
  *outtime =  mcnt;
  *outsx =  mcnsx;
  *outsy =  mcnsy;
  *outsz =  mcnsz;
  *outwgt =  mcnp;

  return;
} /* neutronics_main */

#endif /*NEUTRONICS*/

#endif /* !MCCODE_H */
/* End of file "mccode-r.c". */
/* End of file "mccode-r.c". */

#line 4856 "let_new.c"

#line 1 "mcstas-r.c"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system for McStas.
* Embedded within instrument in runtime mode.
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#include "mcstas-r.h"
#endif
#ifdef DANSE
#include "mcstas-globals.h"
#endif

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/

/*the magnet stack*/
#ifdef MC_POL_COMPAT
void (*mcMagnetPrecession) (double, double, double, double, double, double,
    double, double*, double*, double*, double, Coords, Rotation)=NULL;
Coords   mcMagnetPos;
Rotation mcMagnetRot;
double*  mcMagnetData                = NULL;
/* mcMagneticField(x, y, z, t, Bx, By, Bz) */
int (*mcMagneticField) (double, double, double, double,
    double*, double*, double*, void *) = NULL;
#endif

#ifndef MCSTAS_H

/*******************************************************************************
* mcstore_neutron: stores neutron coodinates into global array (per component)
*******************************************************************************/
void
mcstore_neutron(MCNUM *s, int index, double x, double y, double z,
               double vx, double vy, double vz, double t,
               double sx, double sy, double sz, double p)
{
    double *dptr = &s[11*index];
    *dptr++  = x;
    *dptr++  = y ;
    *dptr++  = z ;
    *dptr++  = vx;
    *dptr++  = vy;
    *dptr++  = vz;
    *dptr++  = t ;
    *dptr++  = sx;
    *dptr++  = sy;
    *dptr++  = sz;
    *dptr    = p ;
} /* mcstore_neutron */

/*******************************************************************************
* mcrestore_neutron: restores neutron coodinates from global array
*******************************************************************************/
void
mcrestore_neutron(MCNUM *s, int index, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *t,
               double *sx, double *sy, double *sz, double *p)
{
    double *dptr = &s[11*index];
    *x  =  *dptr++;
    *y  =  *dptr++;
    *z  =  *dptr++;
    *vx =  *dptr++;
    *vy =  *dptr++;
    *vz =  *dptr++;
    *t  =  *dptr++;
    *sx =  *dptr++;
    *sy =  *dptr++;
    *sz =  *dptr++;
    *p  =  *dptr;
} /* mcrestore_neutron */

/*******************************************************************************
* mcsetstate: transfer parameters into global McStas variables 
*******************************************************************************/
void
mcsetstate(double x, double y, double z, double vx, double vy, double vz,
           double t, double sx, double sy, double sz, double p)
{
  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  mcnx = x;
  mcny = y;
  mcnz = z;
  mcnvx = vx;
  mcnvy = vy;
  mcnvz = vz;
  mcnt = t;
  mcnsx = sx;
  mcnsy = sy;
  mcnsz = sz;
  mcnp = p;
} /* mcsetstate */

/*******************************************************************************
* mcgenstate: set default neutron parameters 
*******************************************************************************/
void
mcgenstate(void)
{
  mcsetstate(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
  /* old initialisation: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
}

/* intersection routines ==================================================== */

/*******************************************************************************
* inside_rectangle: Check if (x,y) is inside rectangle (xwidth, yheight) 
* return 0 if outside and 1 if inside 
*******************************************************************************/
int inside_rectangle(double x, double y, double xwidth, double yheight)
{
  if (x>-xwidth/2 && x<xwidth/2 && y>-yheight/2 && y<yheight/2)
    return 1;
  else
    return 0;
}

/*******************************************************************************
 * box_intersect: compute time intersection with a box
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times dt_in and dt_out
 * This function written by Stine Nyborg, 1999. 
 *******************************************************************************/
int box_intersect(double *dt_in, double *dt_out,
                  double x, double y, double z,
                  double vx, double vy, double vz,
                  double dx, double dy, double dz)
{
  double x_in, y_in, z_in, tt, t[6], a, b;
  int i, count, s;

      /* Calculate intersection time for each of the six box surface planes
       *  If the box surface plane is not hit, the result is zero.*/

  if(vx != 0)
   {
    tt = -(dx/2 + x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[0] = tt;
    else
      t[0] = 0;

    tt = (dx/2 - x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[1] = tt;
    else
      t[1] = 0;
   }
  else
    t[0] = t[1] = 0;

  if(vy != 0)
   {
    tt = -(dy/2 + y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[2] = tt;
    else
      t[2] = 0;

    tt = (dy/2 - y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[3] = tt;
    else
      t[3] = 0;
   }
  else
    t[2] = t[3] = 0;

  if(vz != 0)
   {
    tt = -(dz/2 + z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[4] = tt;
    else
      t[4] = 0;

    tt = (dz/2 - z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[5] = tt;
    else
      t[5] = 0;
   }
  else
    t[4] = t[5] = 0;

  /* The intersection is evaluated and *dt_in and *dt_out are assigned */

  a = b = s = 0;
  count = 0;

  for( i = 0; i < 6; i = i + 1 )
    if( t[i] == 0 )
      s = s+1;
    else if( count == 0 )
    {
      a = t[i];
      count = 1;
    }
    else
    {
      b = t[i];
      count = 2;
    }

  if ( a == 0 && b == 0 )
    return 0;
  else if( a < b )
  {
    *dt_in = a;
    *dt_out = b;
    return 1;
  }
  else
  {
    *dt_in = b;
    *dt_out = a;
    return 1;
  }

} /* box_intersect */

/*******************************************************************************
 * cylinder_intersect: compute intersection with a cylinder
 * returns 0 when no intersection is found
 *      or 2/4/8/16 bits depending on intersection,
 *     and resulting times t0 and t1
 * Written by: EM,NB,ABA 4.2.98 
  *******************************************************************************/
int
cylinder_intersect(double *t0, double *t1, double x, double y, double z,
                   double vx, double vy, double vz, double r, double h)
{
  double D, t_in, t_out, y_in, y_out;
  int ret=1;

  D = (2*vx*x + 2*vz*z)*(2*vx*x + 2*vz*z)
    - 4*(vx*vx + vz*vz)*(x*x + z*z - r*r);

  if (D>=0)
  {
    if (vz*vz + vx*vx) {
      t_in  = (-(2*vz*z + 2*vx*x) - sqrt(D))/(2*(vz*vz + vx*vx));
      t_out = (-(2*vz*z + 2*vx*x) + sqrt(D))/(2*(vz*vz + vx*vx));
    } else if (vy) { /* trajectory parallel to cylinder axis */
      t_in = (-h/2-y)/vy;
      t_out = (h/2-y)/vy;
      if (t_in>t_out){
        double tmp=t_in;
        t_in=t_out;t_out=tmp;
      }
    } else return 0;
    y_in = vy*t_in + y;
    y_out =vy*t_out + y;

    if ( (y_in > h/2 && y_out > h/2) || (y_in < -h/2 && y_out < -h/2) )
      return 0;
    else
    {
      if (y_in > h/2)
        { t_in = ((h/2)-y)/vy; ret += 2; }
      else if (y_in < -h/2)
        { t_in = ((-h/2)-y)/vy; ret += 4; }
      if (y_out > h/2)
        { t_out = ((h/2)-y)/vy; ret += 8; }
      else if (y_out < -h/2)
        { t_out = ((-h/2)-y)/vy; ret += 16; }
    }
    *t0 = t_in;
    *t1 = t_out;
    return ret;
  }
  else
  {
    *t0 = *t1 = 0;
    return 0;
  }
} /* cylinder_intersect */


/*******************************************************************************
 * sphere_intersect: Calculate intersection between a line and a sphere.
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times t0 and t1 
 *******************************************************************************/
int
sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r)
{
  double A, B, C, D, v;

  v = sqrt(vx*vx + vy*vy + vz*vz);
  A = v*v;
  B = 2*(x*vx + y*vy + z*vz);
  C = x*x + y*y + z*z - r*r;
  D = B*B - 4*A*C;
  if(D < 0)
    return 0;
  D = sqrt(D);
  *t0 = (-B - D) / (2*A);
  *t1 = (-B + D) / (2*A);
  return 1;
} /* sphere_intersect */

/*******************************************************************************
 * plane_intersect: Calculate intersection between a plane and a line.
 * returns 0 when no intersection is found (i.e. line is parallel to the plane)
 * returns 1 or -1 when intersection time is positive and negative respectively
 *******************************************************************************/
int
plane_intersect(double *t, double x, double y, double z,
                 double vx, double vy, double vz, double nx, double ny, double nz, double wx, double wy, double wz)
{
  double s;
  if (fabs(s=scalar_prod(nx,ny,nz,vx,vy,vz))<FLT_EPSILON) return 0;
  *t = - scalar_prod(nx,ny,nz,x-wx,y-wy,z-wz)/s;
  if (*t<0) return -1;
  else return 1;
} /* plane_intersect */

#endif /* !MCSTAS_H */
/* End of file "mcstas-r.c". */

#line 5216 "let_new.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\"
int mcdefaultmain = 1;
char mcinstrument_name[] = "LET";
char mcinstrument_source[] = "let_new.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'Guide_channeled'. */
#line 76 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
/*****************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/ref-lib.h
*
* %Identification
* Written by: Peter Christiansen
* Date: August, 2006
* Origin: RISOE
* Release: McStas 1.10
* Version: $Revision$
*
* Commonly used reflection functions are declared in this file which
* are used by some guide and mirror components.
*
* Depends on read_table-lib
*
* Usage: within SHARE
* %include "ref-lib"
*
****************************************************************************/
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.h
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions.
*
* This library may be used directly as an external library. It has no dependency
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#define READ_TABLE_LIB_H "$Revision$"

#define READ_TABLE_STEPTOL  0.04 /* tolerancy for constant step approx */

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#ifdef MAC
#define MC_PATHSEP_C ':'
#define MC_PATHSEP_S ":"
#else  /* !MAC */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MC_PATHSEP_C */

#ifndef MCSTAS
#ifdef WIN32
#define MCSTAS "C:\\mcstas\\lib"
#else  /* !WIN32 */
#ifdef MAC
#define MCSTAS ":mcstas:lib" /* ToDo: What to put here? */
#else  /* !MAC */
#define MCSTAS "/usr/local/lib/mcstas"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MCSTAS */

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

  typedef struct struct_table
  {
    char    filename[1024];
    long    filesize;
    char   *header;  /* text header, e.g. comments */
    double *data;    /* vector { x[0], y[0], ... x[n-1], y[n-1]... } */
    double  min_x;   /* min value of first column */
    double  max_x;   /* max value of first column */
    double  step_x;  /* minimal step value of first column */
    long    rows;    /* number of rows in matrix block */
    long    columns; /* number of columns in matrix block */

    long    begin;   /* start fseek index of block */
    long    end;     /* stop  fseek index of block */
    long    block_number;  /* block index. 0 is catenation of all */
    long    array_length;  /* number of elements in the t_Table array */
    char    monotonic;     /* true when 1st column/vector data is monotonic */
    char    constantstep;  /* true when 1st column/vector data has constant step */
    char    method[32];    /* interpolation method: nearest, linear */
  } t_Table;

typedef struct t_Read_table_file_item {
    int ref_count;
    t_Table *table_ref;
} t_Read_table_file_item;

typedef enum enum_Read_table_file_actions {STORE,FIND,GC}  t_Read_table_file_actions;

/* read_table-lib function prototypes */
/* ========================================================================= */

/* 'public' functions */
long     Table_Read              (t_Table *Table, char *File, long block_number);
long     Table_Read_Offset       (t_Table *Table, char *File, long block_number,
                                  long *offset, long max_lines);
long     Table_Read_Offset_Binary(t_Table *Table, char *File, char *Type,
                                  long *Offset, long Rows, long Columns);
long     Table_Rebin(t_Table *Table); /* rebin table with regular 1st column and interpolate all columns 2:end */
long     Table_Info (t_Table Table);
double   Table_Index(t_Table Table,   long i, long j); /* get indexed value */
double   Table_Value(t_Table Table, double X, long j); /* search X in 1st column and return interpolated value in j-column */
t_Table *Table_Read_Array(char *File, long *blocks);
void     Table_Free_Array(t_Table *Table);
long     Table_Info_Array(t_Table *Table);
int      Table_SetElement(t_Table *Table, long i, long j, double value);
long     Table_Init(t_Table *Table, long rows, long columns); /* create a Table */
double   Table_Value2d(t_Table Table, double X, double Y);    /* same as Table_Index with non-integer indices and 2d interpolation */
MCDETECTOR Table_Write(t_Table Table, char*file, char*xl, char*yl, 
           double x1, double x2, double y1, double y2); /* write Table to disk */
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier);
t_Table *Table_File_List_find(char *name, int block, int offset);
int Table_File_List_gc(t_Table *tab);
void *Table_File_List_store(t_Table *tab);

#define Table_ParseHeader(header, ...) \
  Table_ParseHeader_backend(header,__VA_ARGS__,NULL);

char **Table_ParseHeader_backend(char *header, ...);

/* private functions */
void Table_Free(t_Table *Table);
long Table_Read_Handle(t_Table *Table, FILE *fid, long block_number, long max_lines, char *name);
static void Table_Stat(t_Table *Table);
double Table_Interp1d(double x, double x1, double y1, double x2, double y2);
double Table_Interp1d_nearest(double x, double x1, double y1, double x2, double y2);
double Table_Interp2d(double x, double y, double x1, double y1, double x2, double y2,
double z11, double z12, double z21, double z22);

#endif

/* end of read_table-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.c
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas CVS_090504
* Version: $Revision: 5052 $
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions. Embedded within instrument in runtime mode.
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#endif


/*******************************************************************************
 * void *Table_File_List_Handler(action, item, item_modifier)
 *   ACTION: handle file entries in the read_table-lib file list. If a file is read - it is supposed to be
 *   stored in a list such that we can avoid reading the same file many times.
 *   input  action: FIND, STORE, GC. check if file exists in the list, store an item in the list, or check if it can be garbage collected.
 *   input item: depends on the action.
 *    FIND)  item is a filename, and item_modifier is the block number
 *    STORE) item is the Table to store - item_modifier is ignored
 *    GC)    item is the Table to check. If it has a ref_count >1 then this is simply decremented.
 *   return  depends on the action
 *    FIND)  return a reference to a table+ref_count item if found - NULL otherwise. I.e. NULL means the file has not been read before and must be read again.
 *    STORE) return NULL always
 *    GC)    return NULL if no garbage collection is needed, return an adress to the t_Table which should be garbage collected. 0x1 is returned if
 *           the item is not found in the list
*******************************************************************************/
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier){

    /* logic here is Read_Table should include a call to FIND. If found the return value shoud just be used as
     * if the table had been read. If not found then read the table and STORE.
     * Table_Free should include a call to GC. If this returns non-NULL then we shoudl proceed with freeing the memory
     * associated with the table item - otherwise do nothing since there are more references that may need it.*/ 

    static t_Read_table_file_item read_table_file_list[1024];  
    static int read_table_file_count=0;

    t_Read_table_file_item *tr;
    switch(action){
        case FIND:
            /*interpret data item as a filename, if it is found return a pointer to the table and increment refcount.
             * if not found return the item itself*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                int i=*((int*) item_modifier);
                int j=*( ((int*) item_modifier)+1);
                if ( !strcmp(tr->table_ref->filename,(char *) item) &&
                        tr->table_ref->block_number==i && tr->table_ref->begin==j ){
                    tr->ref_count++;
                    return (void *) tr;
                }
                tr++;
            }
            return NULL;
        case STORE:
            /*find an available slot and store references to table there*/
            tr=&(read_table_file_list[read_table_file_count++]);
            tr->table_ref=(t_Table *)calloc(1,sizeof(t_Table));
            /*copy the contents of the table handle*/
            *(tr->table_ref)= *((t_Table *) item);
            tr->ref_count++;
            return NULL;
        case GC:
            /* Should this item be garbage collected (freed) - if so scratch the entry and return the address of the item - 
             * else decrement ref_count and return NULL.
             * A non-NULL return expects the item to actually be freed afterwards.*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                if ( tr->table_ref->data ==((t_Table *)item)->data && 
                        tr->table_ref->block_number == ((t_Table *)item)->block_number){
                    /*matching item found*/
                    if (tr->ref_count>1){
                        /*the item is found - no garbage collection needed*/
                        tr->ref_count--;
                        return NULL;
                    }else{
                        /* The item is found - move remaining list items up one slot,
                         * and return the table for garbage collection by caller*/
                        while (tr->table_ref!=NULL){
                            *tr=*(tr+1);
                            tr++;
                        }
                        read_table_file_count--;
                        return (t_Table *) item;
                    }
                }
                tr++;
            }
            return (void *)0x1 ;/*item not found*/ 
    } 

}

/* Access functions to the handler*/

/********************************************
 * t_Table *Table_File_List_find(char *name, int block, int offset)
 * input name: filename to search for in the file list
 * input block: data block in the file as each file may contain more than 1 data block.
 * return a ref. to a table if it is found (you may use this pointer and skip reading the file), NULL otherwise (i.e. go ahead and read the file)
*********************************************/
t_Table *Table_File_List_find(char *name, int block, int offset){
    int vars[2]={block,offset};
    t_Read_table_file_item *item = Table_File_List_Handler(FIND,name, vars);
    if (item == NULL){
        return NULL;
    }else{
        return item->table_ref;
    }
}
/********************************************
 * int Table_File_List_gc(t_Table *tab)
 * input tab: the table to check for references.
 * return 0: no garbage collection needed
 *        1: Table's data and header (at least) should be freed.
*********************************************/
int Table_File_List_gc(t_Table *tab){
    void *rval=Table_File_List_Handler(GC,tab,0);
    if (rval==NULL) return 0;
    else return 1;
}


/*****************************************************************************
 * void *Table_File_List_store(t_Table *tab)
 * input tab: pointer to table to store.
 * return None. 
*******************************************************************************/
void *Table_File_List_store(t_Table *tab){
    Table_File_List_Handler(STORE,tab,0);
}


/*******************************************************************************
* FILE *Open_File(char *name, char *Mode, char *path)
*   ACTION: search for a file and open it. Optionally return the opened path.
*   input   name:  file name from which table should be extracted
*           mode: "r", "w", "a" or any valid fopen mode
*           path:  NULL or a pointer to at least 1024 allocated chars
*   return  initialized file handle or NULL in case of error
*******************************************************************************/

  FILE *Open_File(char *File, const char *Mode, char *Path)
  {
    char path[1024];
    FILE *hfile = NULL;
    
    if (!File || File[0]=='\0')                     return(NULL);
    if (!strcmp(File,"NULL") || !strcmp(File,"0"))  return(NULL);
    
    /* search in current or full path */
    strncpy(path, File, 1024);
    hfile = fopen(path, Mode);
    if(!hfile)
    {
      char dir[1024];

      if (!hfile && mcinstrument_source && strlen(mcinstrument_source)) /* search in instrument source location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_source, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_source;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_source, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile && mcinstrument_exe && strlen(mcinstrument_exe)) /* search in PWD instrument executable location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_exe, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_exe;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_exe, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile) /* search in HOME or . */
      {
        strcpy(dir, getenv("HOME") ? getenv("HOME") : ".");
        snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MCSTAS/data */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "data", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MVCSTAS/contrib */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "contrib", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if(!hfile)
      {
        fprintf(stderr, "Error: Could not open input file '%s' (Open_File)\n", File);
        return (NULL);
      }
    }
    if (Path) strncpy(Path, path, 1024);
    return(hfile);
  } /* end Open_File */

/*******************************************************************************
* long Read_Table(t_Table *Table, char *name, int block_number)
*   ACTION: read a single Table from a text file
*   input   Table: pointer to a t_Table structure
*           name:  file name from which table should be extracted
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* File is opened, read and closed
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebinned with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read(t_Table *Table, char *File, long block_number)
  { /* reads all or a single data block from 'file' and returns a Table structure  */
    return(Table_Read_Offset(Table, File, block_number, NULL, 0));
  } /* end Table_Read */

/*******************************************************************************
* long Table_Read_Offset(t_Table *Table, char *name, int block_number, long *offset
*                        long max_rows)
*   ACTION: read a single Table from a text file, starting at offset
*     Same as Table_Read(..) except:
*   input   offset:    pointer to an offset (*offset should be 0 at start)
*           max_rows: max number of data rows to read from file (0 means all)
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset(t_Table *Table, char *File,
                         long block_number, long *offset,
                         long max_rows)
  { /* reads all/a data block in 'file' and returns a Table structure  */
    FILE *hfile;
    long  nelements=0;
    long  begin=0;
    long  filesize=0;
    char  name[1024];
    char  path[1024];
    struct stat stfile;

    /*Need to be able to store the pointer*/
    if (!Table) return(-1);
    
    //if (offset && *offset) snprintf(name, 1024, "%s@%li", File, *offset);
    //else                   
    strncpy(name, File, 1024);
    if(offset && *offset){
        begin=*offset;
    }
    /* Check if the table has already been read from file.
     * If so just reuse the table, if not (this is flagged by returning NULL
     * set up a new table and read the data into it */
    t_Table *tab_p= Table_File_List_find(name,block_number,begin);
    if ( tab_p!=NULL ){
        /*table was found in the Table_File_List*/
        printf("Reusing input file '%s' (Table_Read_Offset)\n", name);
        *Table=*tab_p;
        return Table->rows*Table->columns;
    }

    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
      printf("Opening input file '%s' (Table_Read_Offset)\n", path);
      );
    }
    
    /* read file state */
    stat(path,&stfile); filesize = stfile.st_size;
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    
    Table_Init(Table, 0, 0);

    /* read file content and set the Table */
    nelements = Table_Read_Handle(Table, hfile, block_number, max_rows, name);
    Table->begin = begin;
    Table->end   = ftell(hfile);
    Table->filesize = (filesize>0 ? filesize : 0);
    Table_Stat(Table);
    
    Table_File_List_store(Table);

    if (offset) *offset=Table->end;
    fclose(hfile);
    return(nelements);

  } /* end Table_Read_Offset */

/*******************************************************************************
* long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
*                               long *offset, long rows, long columns)
*   ACTION: read a single Table from a binary file, starting at offset
*     Same as Table_Read_Offset(..) except that it handles binary files.
*   input   type: may be "float"/NULL or "double"
*           offset: pointer to an offset (*offset should be 0 at start)
*           rows   : number of rows (0 means read all)
*           columns: number of columns
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
                                long *offset, long rows, long columns)
  { /* reads all/a data block in binary 'file' and returns a Table structure  */
    long    nelements, sizeofelement;
    long    filesize;
    FILE   *hfile;
    char    path[1024];
    struct stat stfile;
    double *data;
    long    i;
    long    begin;

    if (!Table) return(-1);

    Table_Init(Table, 0, 0);
    
    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
      printf("Opening input file '%s' (Table_Read, Binary)\n", path);
      );
    }
    
    /* read file state */
    stat(File,&stfile);
    filesize = stfile.st_size;
    Table->filesize=filesize;
    
    /* read file content */
    if (type && !strcmp(type,"double")) sizeofelement = sizeof(double);
    else  sizeofelement = sizeof(float);
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    if (rows && filesize > sizeofelement*columns*rows)
      nelements = columns*rows;
    else nelements = (long)(filesize/sizeofelement);
    if (!nelements || filesize <= *offset) return(0);
    data    = (double*)malloc(nelements*sizeofelement);
    if (!data) {
      fprintf(stderr,"Error: allocating %ld elements for %s file '%s'. Too big (Table_Read_Offset_Binary).\n", nelements, type, File);
      exit(-1);
    }
    nelements = fread(data, sizeofelement, nelements, hfile);

    if (!data || !nelements)
    {
      fprintf(stderr,"Error: reading %ld elements from %s file '%s' (Table_Read_Offset_Binary)\n", nelements, type, File);
      exit(-1);
    }
    Table->begin   = begin;
    Table->end     = ftell(hfile);
    if (offset) *offset=Table->end;
    fclose(hfile);
    data = (double*)realloc(data, (double)nelements*sizeofelement);
    /* copy file data into Table */
    if (type && !strcmp(type,"double")) Table->data = data;
    else {
      float  *s;
      double *dataf;
      s     = (float*)data;
      dataf = (double*)malloc(sizeof(double)*nelements);
      for (i=0; i<nelements; i++)
        dataf[i]=s[i];
      free(data);
      Table->data = dataf;
    }
    strncpy(Table->filename, File, 1024);
    Table->rows    = nelements/columns;
    Table->columns = columns;
    Table->array_length = 1;
    Table->block_number = 1;

    Table_Stat(Table);

    return(nelements);
  } /* end Table_Read_Offset_Binary */

/*******************************************************************************
* long Table_Read_Handle(t_Table *Table, FILE *fid, int block_number, long max_rows, char *name)
*   ACTION: read a single Table from a text file handle (private)
*   input   Table:pointer to a t_Table structure
*           fid:  pointer to FILE handle
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*           max_rows: if non 0, only reads that number of lines
*   return  initialized single Table t_Table structure containing data, header, ...
*           modified Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebined with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read_Handle(t_Table *Table, FILE *hfile,
                         long block_number, long max_rows, char *name)
  { /* reads all/a data block from 'file' handle and returns a Table structure  */
    double *Data;
    char *Header              = NULL;
    long  malloc_size         = CHAR_BUF_LENGTH;
    long  malloc_size_h       = 4096;
    long  Rows = 0,   Columns = 0;
    long  count_in_array      = 0;
    long  count_in_header     = 0;
    long  block_Current_index = 0;
    char  flag_End_row_loop   = 0;

    if (!Table) return(-1);
    Table_Init(Table, 0, 0);
    if (name && name[0]!='\0') strncpy(Table->filename, name, 1024);

    if(!hfile) {
       fprintf(stderr, "Error: File handle is NULL (Table_Read_Handle).\n");
       return (-1);
    }
    Header = (char*)  calloc(malloc_size_h, sizeof(char));
    Data   = (double*)calloc(malloc_size,   sizeof(double));
    if ((Header == NULL) || (Data == NULL)) {
       fprintf(stderr, "Error: Could not allocate Table and Header (Table_Read_Handle).\n");
       return (-1);
    }

    int flag_In_array = 0;
    do { /* while (!flag_End_row_loop) */
      char  line[1024*CHAR_BUF_LENGTH];
      long  back_pos=0;   /* ftell start of line */

      back_pos = ftell(hfile);
      if (fgets(line, 1024*CHAR_BUF_LENGTH, hfile) != NULL) { /* analyse line */
        /* first skip blank and tabulation characters */
        int i = strspn(line, " \t");

        /* handle comments: stored in header */
        if (NULL != strchr("#%;/", line[i]))
        { /* line is a comment */
          count_in_header += strlen(line);
          if (count_in_header >= malloc_size_h) {
            /* if succeed and in array : add (and realloc if necessary) */
            malloc_size_h = count_in_header+4096;
            Header        = (char*)realloc(Header, malloc_size_h*sizeof(char));
          }
          strncat(Header, line, 4096);
          flag_In_array=0;
          /* exit line and file if passed desired block */
          if (block_number > 0 && block_number == block_Current_index) {
            flag_End_row_loop = 1;
          }

          /* Continue with next line */
          continue;
        }

        /* get the number of columns splitting line with strtok */
        char  *lexeme;
        char  flag_End_Line = 0;
        long  block_Num_Columns = 0;
        const char seps[] = " ,;\t\n\r";

        lexeme = strtok(line, seps);
        while (!flag_End_Line) {
          if ((lexeme != NULL) && (lexeme[0] != '\0')) {
            /* reading line: the token is not empty */
            double X;
            int    count=1;
            /* test if we have 'NaN','Inf' */
            if (!strncasecmp(lexeme,"NaN",3))
              X = 0;
            else if (!strncasecmp(lexeme,"Inf",3) || !strncasecmp(lexeme,"+Inf",4))
              X = FLT_MAX;
            else if (!strncasecmp(lexeme,"-Inf",4))
              X = -FLT_MAX;
            else
              count = sscanf(lexeme,"%lg",&X);
            if (count == 1) {
              /* reading line: the token is a number in the line */
              if (!flag_In_array) {
                /* reading num: not already in a block: starts a new data block */
                block_Current_index++;
                flag_In_array    = 1;
                block_Num_Columns= 0;
                if (block_number > 0) {
                  /* initialise a new data block */
                  Rows = 0;
                  count_in_array = 0;
                } /* else append */
              }
              /* reading num: all blocks or selected block */
              if (flag_In_array && (block_number == 0 ||
                  block_number == block_Current_index)) {
                /* starting block: already the desired number of rows ? */
                if (block_Num_Columns == 0 &&
                    max_rows > 0 && Rows >= max_rows) {
                  flag_End_Line      = 1;
                  flag_End_row_loop  = 1;
                  flag_In_array      = 0;
                  /* reposition to begining of line (ignore line) */
                  fseek(hfile, back_pos, SEEK_SET);
                } else { /* store into data array */
                  if (count_in_array >= malloc_size) {
                    /* realloc data buffer if necessary */
                    malloc_size = count_in_array+CHAR_BUF_LENGTH;
                    Data = (double*) realloc(Data, malloc_size*sizeof(double));
                    if (Data == NULL) {
                      fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Handle).\n",
                              malloc_size*sizeof(double));
                      return (-1);
                    }
                  }
                  if (0 == block_Num_Columns) Rows++;
                  Data[count_in_array] = X;
                  count_in_array++;
                  block_Num_Columns++;
                }
              } /* reading num: end if flag_In_array */
            } /* end reading num: end if sscanf lexeme -> numerical */
            else {
              /* reading line: the token is not numerical in that line. end block */
              if (block_Current_index == block_number) {
                flag_End_Line = 1;
                flag_End_row_loop = 1;
              } else {
                flag_In_array = 0;
                flag_End_Line = 1;
              }
            }
          }
          else {
            /* no more tokens in line */
            flag_End_Line = 1;
            if (block_Num_Columns > 0) Columns = block_Num_Columns;
          }

          // parse next token
          lexeme = strtok(NULL, seps);

        } /* while (!flag_End_Line) */
      } /* end: if fgets */
      else flag_End_row_loop = 1; /* else fgets : end of file */

    } while (!flag_End_row_loop); /* end while flag_End_row_loop */

    Table->block_number = block_number;
    Table->array_length = 1;

    // shrink header to actual size (plus terminating 0-byte)
    if (count_in_header) {
      Header = (char*)realloc(Header, count_in_header*sizeof(char) + 1);
    }
    Table->header = Header;

    if (count_in_array*Rows*Columns == 0)
    {
      Table->rows         = 0;
      Table->columns      = 0;
      free(Data);
      return (0);
    }
    if (Rows * Columns != count_in_array)
    {
      fprintf(stderr, "Warning: Read_Table :%s %s Data has %li values that should be %li x %li\n",
        (Table->filename ? Table->filename : ""),
        (!block_number ? " catenated" : ""),
        count_in_array, Rows, Columns);
      Columns = count_in_array; Rows = 1;
    }
    Data     = (double*)realloc(Data, count_in_array*sizeof(double));
    Table->data         = Data;
    Table->rows         = Rows;
    Table->columns      = Columns;

    return (count_in_array);

  } /* end Table_Read_Handle */

/*******************************************************************************
* long Table_Rebin(t_Table *Table)
*   ACTION: rebin a single Table, sorting 1st column in ascending order
*   input   Table: single table containing data.
*                  The data block is reallocated in this process
*   return  updated Table with increasing, evenly spaced first column (index 0)
*           number of data elements (-1: error, 0:empty data)
*******************************************************************************/
  long Table_Rebin(t_Table *Table)
  {
    double new_step=0;
    long   i;
    /* performs linear interpolation on X axis (0-th column) */

    if (!Table) return(-1);
    if (!Table->data 
    || Table->rows*Table->columns == 0 || !Table->step_x)
      return(0);
    Table_Stat(Table); /* recompute statitstics and minimal step */
    new_step = Table->step_x; /* minimal step in 1st column */

    if (!(Table->constantstep)) /* not already evenly spaced */
    {
      long Length_Table;
      double *New_Table;

      Length_Table = ceil(fabs(Table->max_x - Table->min_x)/new_step)+1;
      New_Table    = (double*)malloc(Length_Table*Table->columns*sizeof(double));

      for (i=0; i < Length_Table; i++)
      {
        long   j;
        double X;
        X = Table->min_x + i*new_step;
        New_Table[i*Table->columns] = X;
        for (j=1; j < Table->columns; j++)
          New_Table[i*Table->columns+j]
                = Table_Value(*Table, X, j);
      } /* end for i */

      Table->rows = Length_Table;
      Table->step_x = new_step;
      Table->max_x = Table->min_x + (Length_Table-1)*new_step; 
      /*max might not be the same anymore
       * Use Length_Table -1 since the first and laset rows are the limits of the defined interval.*/
      free(Table->data);
      Table->data = New_Table;
      Table->constantstep=1;
    } /* end else (!constantstep) */
    return (Table->rows*Table->columns);
  } /* end Table_Rebin */

/*******************************************************************************
* double Table_Index(t_Table Table, long i, long j)
*   ACTION: read an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*   return  Value = data[i][j]
* Returns Value from the i-th row, j-th column of Table
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

double Table_Index(t_Table Table, long i, long j)
{
  long AbsIndex;

  if (Table.rows == 1 || Table.columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table.columns*Table.rows - 1);
    i = 0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table.rows - 1);
    j = MIN(MAX(0, j), Table.columns - 1);
  }

  /* handle vectors specifically */
  AbsIndex = i*(Table.columns)+j;

  if (Table.data != NULL)
    return (Table.data[AbsIndex]);
  else
    return 0;
} /* end Table_Index */

/*******************************************************************************
* void Table_SetElement(t_Table *Table, long i, long j, double value)
*   ACTION: set an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*           value = data[i][j]
* Returns 0 in case of error
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/
int Table_SetElement(t_Table *Table, long i, long j,
                     double value)
{
  long AbsIndex;

  if (Table->rows == 1 || Table->columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table->columns*Table->rows - 1); i=0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table->rows - 1);
    j = MIN(MAX(0, j), Table->columns - 1);
  }

  AbsIndex = i*(Table->columns)+j;
  if (Table->data != NULL) {
    Table->data[AbsIndex] = value;
    return 1;
  }

  return 0;
} /* end Table_SetElement */

/*******************************************************************************
* double Table_Value(t_Table Table, double X, long j)
*   ACTION: read column [j] of a single Table at row which 1st column is X
*   input   Table: table containing data.
*           X : data value in the first column (index 0)
*           j : index of column from which is extracted the Value (0:Columns-1)
*   return  Value = data[index for X][j] with linear interpolation
* Returns Value from the j-th column of Table corresponding to the
* X value for the 1st column (index 0)
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
double Table_Value(t_Table Table, double X, long j)
{
  long   Index = -1;
  double X1=0, Y1=0, X2=0, Y2=0;
  double ret=0;

  if (X > Table.max_x) return Table_Index(Table,Table.rows-1  ,j);
  if (X < Table.min_x) return Table_Index(Table,0  ,j);

  // Use constant-time lookup when possible
  if(Table.constantstep) {
    Index = (long)floor(
              (X - Table.min_x) / (Table.max_x - Table.min_x) * (Table.rows-1));
    X1 = Table_Index(Table,Index  ,0);
    X2 = Table_Index(Table,Index+1,0);
  }
  // Use binary search on large, monotonic tables
  else if(Table.monotonic && Table.rows > 100) {
    long left = Table.min_x;
    long right = Table.max_x;

    while (!((X1 <= X) && (X < X2)) && (right - left > 1)) {
      Index = (left + right) / 2;

      X1 = Table_Index(Table, Index-1, 0);
      X2 = Table_Index(Table, Index,   0);

      if (X < X1) {
        right = Index;
      } else {
        left  = Index;
      }
    }
  }

  // Fall back to linear search, if no-one else has set X1, X2 correctly
  if (!((X1 <= X) && (X < X2))) {
    /* look for index surrounding X in the table -> Index */
    for (Index=1; Index < Table.rows-1; Index++) {
        X1 = Table_Index(Table, Index-1,0);
        X2 = Table_Index(Table, Index  ,0);
        if ((X1 <= X) && (X < X2)) break;
      } /* end for Index */
  }

  Y1 = Table_Index(Table,Index-1,j);
  Y2 = Table_Index(Table,Index  ,j);

  if (!strcmp(Table.method,"linear")) {
    ret = Table_Interp1d(X, X1,Y1, X2,Y2);
  }
  else if (!strcmp(Table.method,"nearest")) {
    ret = Table_Interp1d_nearest(X, X1,Y1, X2,Y2);
  }

  return ret;
} /* end Table_Value */

/*******************************************************************************
* double Table_Value2d(t_Table Table, double X, double Y)
*   ACTION: read element [X,Y] of a matrix Table
*   input   Table: table containing data.
*           X : row index, may be non integer
*           Y : column index, may be non integer
*   return  Value = data[index X][index Y] with bi-linear interpolation
* Returns Value for the indices [X,Y]
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
  double Table_Value2d(t_Table Table, double X, double Y)
  {
    long   x1,x2,y1,y2;
    double z11,z12,z21,z22;
    double ret=0;

    x1 = (long)floor(X);
    y1 = (long)floor(Y);

    if (x1 > Table.rows-1 || x1 < 0) {
      x2 = x1;
    } else {
      x2 = x1 + 1;
    }

    if (y1 > Table.columns-1 || y1 < 0) {
      y2 = y1;
    } else {
      y2 = y1 + 1;
    }

    z11 = Table_Index(Table, x1, y1);

    if (y2 != y1) z12=Table_Index(Table, x1, y2); else z12 = z11;
    if (x2 != x1) z21=Table_Index(Table, x2, y1); else z21 = z11;
    if (y2 != y1) z22=Table_Index(Table, x2, y2); else z22 = z21;

    if (!strcmp(Table.method,"linear"))
      ret = Table_Interp2d(X,Y, x1,y1,x2,y2, z11,z12,z21,z22);
    else {
      if (fabs(X-x1) < fabs(X-x2)) {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z11; else ret = z12;
      } else {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z21; else ret = z22;
      }
    }
    return ret;
  } /* end Table_Value2d */


/*******************************************************************************
* void Table_Free(t_Table *Table)
*   ACTION: free a single Table
*   return: empty Table
*******************************************************************************/
  void Table_Free(t_Table *Table)
  {
    if( !Table_File_List_gc(Table) ){
       return;
    } 
    if (!Table) return;
    if (Table->data   != NULL) free(Table->data);
    if (Table->header != NULL) free(Table->header);
    Table->data   = NULL;
    Table->header = NULL;
  } /* end Table_Free */

/******************************************************************************
* void Table_Info(t_Table Table)
*    ACTION: print informations about a single Table
*******************************************************************************/
  long Table_Info(t_Table Table)
  {
    char buffer[256];
    long ret=0;

    if (!Table.block_number) strcpy(buffer, "catenated");
    else sprintf(buffer, "block %li", Table.block_number);
    printf("Table from file '%s' (%s)",
      Table.filename ? Table.filename : "", buffer);
    if ((Table.data != NULL) && (Table.rows*Table.columns))
    {
      printf(" is %li x %li ", Table.rows, Table.columns);
      if (Table.rows*Table.columns > 1)
        printf("(x=%g:%g)", Table.min_x, Table.max_x);
      else printf("(x=%g) ", Table.min_x);
      ret = Table.rows*Table.columns;
      if (Table.monotonic)    printf(", monotonic");
      if (Table.constantstep) printf(", constant step");
      printf(". interpolation: %s\n", Table.method);
    }
    else printf(" is empty.\n");

    if (Table.header && strlen(Table.header)) {
      char *header;
      int  i;
      header = malloc(80);
      if (!header) return(ret);
      for (i=0; i<80; header[i++]=0);
      strncpy(header, Table.header, 75);
      if (strlen(Table.header) > 75) {
        strcat( header, " ...");
      }
      for (i=0; i<strlen(header); i++)
        if (header[i] == '\n' || header[i] == '\r') header[i] = ';';
      printf("  '%s'\n", header);
      free(header);
    }

    return(ret);
  } /* end Table_Info */

/******************************************************************************
* long Table_Init(t_Table *Table, m, n)
*   ACTION: initialise a Table to empty m by n table
*   return: empty Table
******************************************************************************/
long Table_Init(t_Table *Table, long rows, long columns)
{
  double *data=NULL;
  long   i;

  if (!Table) return(0);

  Table->header  = NULL;
  Table->filename[0]= '\0';
  Table->filesize= 0;
  Table->min_x   = 0;
  Table->max_x   = 0;
  Table->step_x  = 0;
  Table->block_number = 0;
  Table->array_length = 0;
  Table->monotonic    = 0;
  Table->constantstep = 0;
  Table->begin   = 0;
  Table->end     = 0;
  strcpy(Table->method,"linear");

  if (rows*columns >= 1) {
    data    = (double*)malloc(rows*columns*sizeof(double));
    if (data) for (i=0; i < rows*columns; data[i++]=0);
    else {
      fprintf(stderr,"Error: allocating %ld double elements."
                     "Too big (Table_Init).\n", rows*columns);
      rows = columns = 0;
    }
  }
  Table->rows    = (rows >= 1 ? rows : 0);
  Table->columns = (columns >= 1 ? columns : 0);
  Table->data    = data;
  return(Table->rows*Table->columns);
} /* end Table_Init */

/******************************************************************************
* long Table_Write(t_Table Table, char *file, x1,x2, y1,y2)
*   ACTION: write a Table to disk (ascii).
*     when x1=x2=0 or y1=y2=0, the table default limits are used.
*   return: 0=all is fine, non-0: error
*******************************************************************************/
MCDETECTOR Table_Write(t_Table Table, char *file, char *xl, char *yl, 
  double x1, double x2, double y1, double y2)
{
  long    i =0;
  MCDETECTOR detector;

  if ((Table.data == NULL) && (Table.rows*Table.columns)) {
    detector.m = 0;
    return(detector); /* Table is empty - nothing to do */
  }
  if (!x1 && !x2) {
    x1 = Table.min_x;
    x2 = Table.max_x;
  }
  if (!y1 && !y2) {
    y1 = 1;
    y2 = Table.columns;
  }

  /* transfer content of the Table into a 2D detector */
  Coords coords = { 0, 0, 0};

  if (Table.rows == 1 || Table.columns == 1) {
    detector = mcdetector_out_1D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      "x", x1, x2,
                      Table.rows * Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  } else {
    detector = mcdetector_out_2D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      x1, x2, y1, y2,
                      Table.rows, Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  }
  return(detector);
}

/******************************************************************************
* void Table_Stat(t_Table *Table)
*   ACTION: computes min/max/mean step of 1st column for a single table (private)
*   return: updated Table
*******************************************************************************/
  static void Table_Stat(t_Table *Table)
  {
    long   i;
    double max_x, min_x;
    double row=1;
    char   monotonic=1;
    char   constantstep=1;
    double step=0;
    long n;

    if (!Table) return;
    if (!Table->rows || !Table->columns) return;
    if (Table->rows == 1) row=0; // single row
    max_x = -FLT_MAX;
    min_x =  FLT_MAX;
    n     = (row ? Table->rows : Table->columns);
    /* get min and max of first column/vector */
    for (i=0; i < n; i++)
    {
      double X;
      X = (row ? Table_Index(*Table,i  ,0)
                               : Table_Index(*Table,0, i));
      if (X < min_x) min_x = X;
      if (X > max_x) max_x = X;
    } /* for */
    
    /* test for monotonicity and constant step if the table is an XY or single vector */
    if (n > 1) {
      /* mean step */
      step = (max_x - min_x)/(n-1);
      /* now test if table is monotonic on first column, and get minimal step size */
      for (i=0; i < n-1; i++) {
        double X, diff;;
        X    = (row ? Table_Index(*Table,i  ,0)
                    : Table_Index(*Table,0,  i));
        diff = (row ? Table_Index(*Table,i+1,0)
                    : Table_Index(*Table,0,  i+1)) - X;
        if (diff && fabs(diff) < fabs(step)) step = diff;
        /* change sign ? */
        if ((max_x - min_x)*diff < 0 && monotonic)
          monotonic = 0;
      } /* end for */
      
      /* now test if steps are constant within READ_TABLE_STEPTOL */
      if(!step){
        /*means there's a disconitnuity -> not constantstep*/
        constantstep=0;
      }else if (monotonic) {
        for (i=0; i < n-1; i++) {
          double X, diff;
          X    = (row ? Table_Index(*Table,i  ,0)
              : Table_Index(*Table,0,  i));
          diff = (row ? Table_Index(*Table,i+1,0)
              : Table_Index(*Table,0,  i+1)) - X;
          if ( fabs(step)*(1+READ_TABLE_STEPTOL) < fabs(diff) ||
                fabs(diff) < fabs(step)*(1-READ_TABLE_STEPTOL) )
          { constantstep = 0; break; }
        }
      }

    }
    Table->step_x= step;
    Table->max_x = max_x;
    Table->min_x = min_x;
    Table->monotonic = monotonic;
    Table->constantstep = constantstep;
  } /* end Table_Stat */

/******************************************************************************
* t_Table *Table_Read_Array(char *File, long *blocks)
*   ACTION: read as many data blocks as available, iteratively from file
*   return: initialized t_Table array, last element is an empty Table.
*           the number of extracted blocks in non NULL pointer *blocks
*******************************************************************************/
  t_Table *Table_Read_Array(char *File, long *blocks)
  {
    t_Table *Table_Array=NULL;
    long offset=0;
    long block_number=0;
    long allocated=256;
    long nelements=1;

    /* fisrt allocate an initial empty t_Table array */
    Table_Array = (t_Table *)malloc(allocated*sizeof(t_Table));
    if (!Table_Array) {
      fprintf(stderr, "Error: Can not allocate memory %li (Table_Read_Array).\n",
         allocated*sizeof(t_Table));
      *blocks = 0;
      return (NULL);
    }

    while (nelements > 0)
    {
      t_Table Table;

      /* if ok, set t_Table block number else exit loop */
      block_number++;
      Table.block_number = block_number;
      
      /* access file at offset and get following block. Block number is from the set offset
       * hence the hardcoded 1 - i.e. the next block counted from offset.*/
      nelements = Table_Read_Offset(&Table, File, 1, &offset,0);
      /* if t_Table array is not long enough, expand and realocate */
      if (block_number >= allocated-1) {
        allocated += 256;
        Table_Array = (t_Table *)realloc(Table_Array,
           allocated*sizeof(t_Table));
        if (!Table_Array) {
          fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Array).\n",
              allocated*sizeof(t_Table));
          *blocks = 0;
          return (NULL);
        }
      }
      /* store it into t_Table array */
      //snprintf(Table.filename, 1024, "%s#%li", File, block_number-1);
      Table_Array[block_number-1] = Table;
      /* continues until we find an empty block */
    }
    /* send back number of extracted blocks */
    if (blocks) *blocks = block_number-1;

    /* now store total number of elements in Table array */
    for (offset=0; offset < block_number;
      Table_Array[offset++].array_length = block_number-1);

    return(Table_Array);
  } /* end Table_Read_Array */
/*******************************************************************************
* void Table_Free_Array(t_Table *Table)
*   ACTION: free a Table array
*******************************************************************************/
  void Table_Free_Array(t_Table *Table)
  {
    long index=0;
    if (!Table) return;
    while (Table[index].data || Table[index].header){
            Table_Free(&Table[index]);
            index++;
    }
    free(Table);
  } /* end Table_Free_Array */

/******************************************************************************
* long Table_Info_Array(t_Table *Table)
*    ACTION: print informations about a Table array
*    return: number of elements in the Table array
*******************************************************************************/
  long Table_Info_Array(t_Table *Table)
  {
    long index=0;

    if (!Table) return(-1);
    while (index < Table[index].array_length
       && (Table[index].data || Table[index].header)
       && (Table[index].rows*Table[index].columns) ) {
      Table_Info(Table[index]);
      index++;
    }
    printf("This Table array contains %li elements\n", index);
    return(index);
  } /* end Table_Info_Array */

/******************************************************************************
* char **Table_ParseHeader(char *header, symbol1, symbol2, ..., NULL)
*    ACTION: search for char* symbols in header and return their value or NULL
*            the search is not case sensitive.
*            Last argument MUST be NULL
*    return: array of char* with line following each symbol, or NULL if not found
*******************************************************************************/
#ifndef MyNL_ARGMAX
#define MyNL_ARGMAX 50
#endif

char **Table_ParseHeader_backend(char *header, ...){
  va_list ap;
  char exit_flag=0;
  int counter   =0;
  char **ret    =NULL;
  if (!header || header[0]=='\0') return(NULL);

  ret = (char**)calloc(MyNL_ARGMAX, sizeof(char*));
  if (!ret) {
    printf("Table_ParseHeader: Cannot allocate %i values array for Parser (Table_ParseHeader).\n",
      MyNL_ARGMAX);
    return(NULL);
  }
  for (counter=0; counter < MyNL_ARGMAX; ret[counter++] = NULL);
  counter=0;

  va_start(ap, header);
  while(!exit_flag && counter < MyNL_ARGMAX-1)
  {
    char *arg_char=NULL;
    char *pos     =NULL;
    /* get variable argument value as a char */
    arg_char = va_arg(ap, char *);
    if (!arg_char || arg_char[0]=='\0'){
      exit_flag = 1; break;
    }
    /* search for the symbol in the header */
    pos = (char*)strcasestr(header, arg_char);
    if (pos) {
      char *eol_pos;
      eol_pos = strchr(pos+strlen(arg_char), '\n');
      if (!eol_pos)
        eol_pos = strchr(pos+strlen(arg_char), '\r');
      if (!eol_pos)
        eol_pos = pos+strlen(pos)-1;
      ret[counter] = (char*)malloc(eol_pos - pos);
      if (!ret[counter]) {
        printf("Table_ParseHeader: Cannot allocate value[%i] array for Parser searching for %s (Table_ParseHeader).\n",
          counter, arg_char);
        exit_flag = 1; break;
      }
      strncpy(ret[counter], pos+strlen(arg_char), eol_pos - pos - strlen(arg_char));
      ret[counter][eol_pos - pos - strlen(arg_char)]='\0';
    }
    counter++;
  }
  va_end(ap);
  return(ret);
} /* Table_ParseHeader */

/******************************************************************************
* double Table_Interp1d(x, x1, y1, x2, y2)
*    ACTION: interpolates linearly at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d(double x,
  double x1, double y1,
  double x2, double y2)
{
  double slope;
  if (x2 == x1) return (y1+y2)/2;
  if (y1 == y2) return  y1;
  slope = (y2 - y1)/(x2 - x1);
  return y1+slope*(x - x1);
} /* Table_Interp1d */

/******************************************************************************
* double Table_Interp1d_nearest(x, x1, y1, x2, y2)
*    ACTION: table lookup with nearest method at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d_nearest(double x,
  double x1, double y1,
  double x2, double y2)
{
  if (fabs(x-x1) < fabs(x-x2)) return (y1);
  else return(y2);
} /* Table_Interp1d_nearest */

/******************************************************************************
* double Table_Interp2d(x,y, x1,y1, x2,y2, z11,z12,z21,z22)
*    ACTION: interpolates bi-linearly at (x,y) between z1=f(x1,y1) and z2=f(x2,y2)
*    return: z=f(x,y) value
*    x,y |   x1   x2
*    ----------------
*     y1 |   z11  z21
*     y2 |   z12  z22
*******************************************************************************/
double Table_Interp2d(double x, double y,
  double x1, double y1,
  double x2, double y2,
  double z11, double z12, double z21, double z22)
{
  double ratio_x, ratio_y;
  if (x2 == x1) return Table_Interp1d(y, y1,z11, y2,z12);
  if (y1 == y2) return Table_Interp1d(x, x1,z11, x2,z21);

  ratio_y = (y - y1)/(y2 - y1);
  ratio_x = (x - x1)/(x2 - x1);
  return (1-ratio_x)*(1-ratio_y)*z11 + ratio_x*(1-ratio_y)*z21
    + ratio_x*ratio_y*z22         + (1-ratio_x)*ratio_y*z12;
} /* Table_Interp2d */

/* end of read_table-lib.c */


#ifndef REF_LIB_H
#define REF_LIB_H "$Revision$"

void StdReflecFunc(double, double*, double*);
void TableReflecFunc(double, t_Table*, double*);

#endif

/* end of ref-lib.h */
/****************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/ref-lib.c
*
* %Identification
* Written by: Peter Christiansen
* Date: August, 2006
* Origin: RISOE
* Release: McStas 1.10
* Version: $Revision$
*
* Commonly used reflection functions are declared in this file which
* are used by some guide and mirror components.
*
* Variable names have prefix 'mc_ref_' for 'McStas Reflection' 
* to avoid conflicts
*
* Usage: within SHARE
* %include "ref-lib"
*
****************************************************************************/

#ifndef REF_LIB_H
#include "ref-lib.h"
#endif

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#include "read_table-lib.c"
#endif

/****************************************************************************
* void StdReflecFunc(double q, double *par, double *r)
* 
* The McStas standard analytic parametrization of the reflectivity.
* The parameters are:
* R0:      [1]    Low-angle reflectivity
* Qc:      [AA-1] Critical scattering vector
* alpha:   [AA]   Slope of reflectivity
* m:       [1]    m-value of material. Zero means completely absorbing.
* W:       [AA-1] Width of supermirror cut-off
*****************************************************************************/
void StdReflecFunc(double mc_pol_q, double *mc_pol_par, double *mc_pol_r) {
    double R0    = mc_pol_par[0];
    double Qc    = mc_pol_par[1];
    double alpha = mc_pol_par[2];
    double m     = mc_pol_par[3];
    double W     = mc_pol_par[4];
    double beta  = 0;
    mc_pol_q     = fabs(mc_pol_q);
    double arg;
        
    /* Simpler parametrization from Henrik Jacobsen uses these values that depend on m only.
       double m_value=m*0.9853+0.1978;
       double W=-0.0002*m_value+0.0022;
       double alpha=0.2304*m_value+5.0944;
       double beta=-7.6251*m_value+68.1137; 
       If W and alpha are set to 0, use Henrik's approach for estimating these parameters
       and apply the formulation:
       arg = R0*0.5*(1-tanh(arg))*(1-alpha*(q-Qc)+beta*(q-Qc)*(q-Qc));
    */  
    if (W==0 && alpha==0) {
      m=m*0.9853+0.1978;
      W=-0.0002*m+0.0022;
      alpha=0.2304*m+5.0944;
      beta=-7.6251*m+68.1137;
      if (m<=3) {
	alpha=m;
	beta=0;
      }
    }
    
    arg = W > 0 ? (mc_pol_q - m*Qc)/W : 11;

    if (arg > 10 || m <= 0 || Qc <=0 || R0 <= 0) {
      *mc_pol_r = 0;
      return;
    }
    
    if (m < 1) { Qc *= m; m=1; }
    
    if(mc_pol_q <= Qc) {      
      *mc_pol_r = R0;
      return;
    }
    
    
    *mc_pol_r = R0*0.5*(1 - tanh(arg))*(1 - alpha*(mc_pol_q - Qc) + beta*(mc_pol_q - Qc)*(mc_pol_q - Qc));
    
    return;
  }

/****************************************************************************
* void TableReflecFunc(double q, t_Table *par, double *r) {
* 
* Looks up the reflectivity in a table using the routines in read_table-lib.
*****************************************************************************/
void TableReflecFunc(double mc_pol_q, t_Table *mc_pol_par, double *mc_pol_r) {
    
  *mc_pol_r = Table_Value(*mc_pol_par, mc_pol_q, 1);
  if(*mc_pol_r>1)
    *mc_pol_r = 1;
  return;
}

/* end of ref-lib.c */

#line 6812 "let_new.c"

/* Shared user declarations for all components 'MCPL_output_horace'. */
#line 63 "MCPL_output_horace.comp"
#include <mcpl.h>
#ifndef MCPL_FILE_EXISTS
#define MCPL_FILE_EXISTS
int mcpl_file_exist (char *filename)
  {
    struct stat   buffer;
    return (stat (filename, &buffer) == 0);
  }
#endif
#line 6825 "let_new.c"

/* Instrument parameters. */
MCNUM mcipEi;
MCNUM mcipfreq;
char* mcipchopper;
char* mcipoutput_filename;
char* mcipsample;

#define mcNUMIPAR 5
int mcnumipar = 5;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "Ei", &mcipEi, instr_type_double, "1.89", 
  "freq", &mcipfreq, instr_type_double, "100", 
  "chopper", &mcipchopper, instr_type_string, "high resolution", 
  "output_filename", &mcipoutput_filename, instr_type_string, "mcstas.mcpl", 
  "sample", &mcipsample, instr_type_string, "plate", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  LET
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaLET coords_set(0,0,0)
#define Ei mcipEi
#define freq mcipfreq
#define chopper mcipchopper
#define output_filename mcipoutput_filename
#define sample mcipsample
#line 5 "let_new.instr"
double v_foc;
double phase_time, E_min, E_max;
double PR, Res5, E_foc;
double pha_offset = 0;
double jitter = 0;
double slot_size;
int plate;
#line 6863 "let_new.c"
#undef sample
#undef output_filename
#undef chopper
#undef freq
#undef Ei
#undef mcposaLET
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*23];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[23];
Coords mccomp_posr[23];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[23];
MCNUM  mcPCounter[23];
MCNUM  mcP2Counter[23];
#define mcNUMCOMP 22 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[23];
/* Flag true when previous component acted on the neutron (SCATTER) */
MCNUM mcScattered=0;
/* Flag true when neutron should be restored (RESTORE) */
MCNUM mcRestore=0;
/* Declarations of component definition and setting parameters. */

/* Setting parameters for component 'Source' [2]. */
char mccSource_Face[16384];
MCNUM mccSource_E0;
MCNUM mccSource_E1;
MCNUM mccSource_tally;
MCNUM mccSource_modPosition;
MCNUM mccSource_modXsize;
MCNUM mccSource_modZsize;
MCNUM mccSource_xw;
MCNUM mccSource_yh;
MCNUM mccSource_dist;

/* Setting parameters for component 'shutter' [3]. */
MCNUM mccshutter_w1;
MCNUM mccshutter_h1;
MCNUM mccshutter_w2;
MCNUM mccshutter_h2;
MCNUM mccshutter_l;
MCNUM mccshutter_R0;
MCNUM mccshutter_Qc;
MCNUM mccshutter_alpha;
MCNUM mccshutter_m;
MCNUM mccshutter_nslit;
MCNUM mccshutter_d;
MCNUM mccshutter_Qcx;
MCNUM mccshutter_Qcy;
MCNUM mccshutter_alphax;
MCNUM mccshutter_alphay;
MCNUM mccshutter_W;
MCNUM mccshutter_mx;
MCNUM mccshutter_my;
MCNUM mccshutter_nu;
MCNUM mccshutter_phase;

/* Setting parameters for component 'insert' [4]. */
MCNUM mccinsert_w1;
MCNUM mccinsert_h1;
MCNUM mccinsert_w2;
MCNUM mccinsert_h2;
MCNUM mccinsert_l;
MCNUM mccinsert_R0;
MCNUM mccinsert_Qc;
MCNUM mccinsert_alpha;
MCNUM mccinsert_m;
MCNUM mccinsert_nslit;
MCNUM mccinsert_d;
MCNUM mccinsert_Qcx;
MCNUM mccinsert_Qcy;
MCNUM mccinsert_alphax;
MCNUM mccinsert_alphay;
MCNUM mccinsert_W;
MCNUM mccinsert_mx;
MCNUM mccinsert_my;
MCNUM mccinsert_nu;
MCNUM mccinsert_phase;

/* Setting parameters for component 'afterinsert' [5]. */
MCNUM mccafterinsert_w1;
MCNUM mccafterinsert_h1;
MCNUM mccafterinsert_w2;
MCNUM mccafterinsert_h2;
MCNUM mccafterinsert_l;
MCNUM mccafterinsert_R0;
MCNUM mccafterinsert_Qc;
MCNUM mccafterinsert_alpha;
MCNUM mccafterinsert_m;
MCNUM mccafterinsert_nslit;
MCNUM mccafterinsert_d;
MCNUM mccafterinsert_Qcx;
MCNUM mccafterinsert_Qcy;
MCNUM mccafterinsert_alphax;
MCNUM mccafterinsert_alphay;
MCNUM mccafterinsert_W;
MCNUM mccafterinsert_mx;
MCNUM mccafterinsert_my;
MCNUM mccafterinsert_nu;
MCNUM mccafterinsert_phase;

/* Setting parameters for component 'Res1' [6]. */
MCNUM mccRes1_theta_0;
MCNUM mccRes1_radius;
MCNUM mccRes1_yheight;
MCNUM mccRes1_nu;
MCNUM mccRes1_nslit;
MCNUM mccRes1_jitter;
MCNUM mccRes1_delay;
MCNUM mccRes1_isfirst;
MCNUM mccRes1_n_pulse;
MCNUM mccRes1_abs_out;
MCNUM mccRes1_phase;
MCNUM mccRes1_xwidth;
MCNUM mccRes1_verbose;

/* Setting parameters for component 'Res1_counter' [7]. */
MCNUM mccRes1_counter_theta_0;
MCNUM mccRes1_counter_radius;
MCNUM mccRes1_counter_yheight;
MCNUM mccRes1_counter_nu;
MCNUM mccRes1_counter_nslit;
MCNUM mccRes1_counter_jitter;
MCNUM mccRes1_counter_delay;
MCNUM mccRes1_counter_isfirst;
MCNUM mccRes1_counter_n_pulse;
MCNUM mccRes1_counter_abs_out;
MCNUM mccRes1_counter_phase;
MCNUM mccRes1_counter_xwidth;
MCNUM mccRes1_counter_verbose;

/* Setting parameters for component 'betweenchop1and2' [8]. */
MCNUM mccbetweenchop1and2_w1;
MCNUM mccbetweenchop1and2_h1;
MCNUM mccbetweenchop1and2_w2;
MCNUM mccbetweenchop1and2_h2;
MCNUM mccbetweenchop1and2_l;
MCNUM mccbetweenchop1and2_R0;
MCNUM mccbetweenchop1and2_Qc;
MCNUM mccbetweenchop1and2_alpha;
MCNUM mccbetweenchop1and2_m;
MCNUM mccbetweenchop1and2_nslit;
MCNUM mccbetweenchop1and2_d;
MCNUM mccbetweenchop1and2_Qcx;
MCNUM mccbetweenchop1and2_Qcy;
MCNUM mccbetweenchop1and2_alphax;
MCNUM mccbetweenchop1and2_alphay;
MCNUM mccbetweenchop1and2_W;
MCNUM mccbetweenchop1and2_mx;
MCNUM mccbetweenchop1and2_my;
MCNUM mccbetweenchop1and2_nu;
MCNUM mccbetweenchop1and2_phase;

/* Setting parameters for component 'betweenchop2and3' [9]. */
MCNUM mccbetweenchop2and3_w1;
MCNUM mccbetweenchop2and3_h1;
MCNUM mccbetweenchop2and3_w2;
MCNUM mccbetweenchop2and3_h2;
MCNUM mccbetweenchop2and3_l;
MCNUM mccbetweenchop2and3_R0;
MCNUM mccbetweenchop2and3_Qc;
MCNUM mccbetweenchop2and3_alpha;
MCNUM mccbetweenchop2and3_m;
MCNUM mccbetweenchop2and3_nslit;
MCNUM mccbetweenchop2and3_d;
MCNUM mccbetweenchop2and3_Qcx;
MCNUM mccbetweenchop2and3_Qcy;
MCNUM mccbetweenchop2and3_alphax;
MCNUM mccbetweenchop2and3_alphay;
MCNUM mccbetweenchop2and3_W;
MCNUM mccbetweenchop2and3_mx;
MCNUM mccbetweenchop2and3_my;
MCNUM mccbetweenchop2and3_nu;
MCNUM mccbetweenchop2and3_phase;

/* Setting parameters for component 'PR' [10]. */
MCNUM mccPR_theta_0;
MCNUM mccPR_radius;
MCNUM mccPR_yheight;
MCNUM mccPR_nu;
MCNUM mccPR_nslit;
MCNUM mccPR_jitter;
MCNUM mccPR_delay;
MCNUM mccPR_isfirst;
MCNUM mccPR_n_pulse;
MCNUM mccPR_abs_out;
MCNUM mccPR_phase;
MCNUM mccPR_xwidth;
MCNUM mccPR_verbose;

/* Setting parameters for component 'betweenchop3and4' [11]. */
MCNUM mccbetweenchop3and4_w1;
MCNUM mccbetweenchop3and4_h1;
MCNUM mccbetweenchop3and4_w2;
MCNUM mccbetweenchop3and4_h2;
MCNUM mccbetweenchop3and4_l;
MCNUM mccbetweenchop3and4_R0;
MCNUM mccbetweenchop3and4_Qc;
MCNUM mccbetweenchop3and4_alpha;
MCNUM mccbetweenchop3and4_m;
MCNUM mccbetweenchop3and4_nslit;
MCNUM mccbetweenchop3and4_d;
MCNUM mccbetweenchop3and4_Qcx;
MCNUM mccbetweenchop3and4_Qcy;
MCNUM mccbetweenchop3and4_alphax;
MCNUM mccbetweenchop3and4_alphay;
MCNUM mccbetweenchop3and4_W;
MCNUM mccbetweenchop3and4_mx;
MCNUM mccbetweenchop3and4_my;
MCNUM mccbetweenchop3and4_nu;
MCNUM mccbetweenchop3and4_phase;

/* Setting parameters for component 'CR' [12]. */
MCNUM mccCR_theta_0;
MCNUM mccCR_radius;
MCNUM mccCR_yheight;
MCNUM mccCR_nu;
MCNUM mccCR_nslit;
MCNUM mccCR_jitter;
MCNUM mccCR_delay;
MCNUM mccCR_isfirst;
MCNUM mccCR_n_pulse;
MCNUM mccCR_abs_out;
MCNUM mccCR_phase;
MCNUM mccCR_xwidth;
MCNUM mccCR_verbose;

/* Setting parameters for component 'betweenchop4andatt' [13]. */
MCNUM mccbetweenchop4andatt_w1;
MCNUM mccbetweenchop4andatt_h1;
MCNUM mccbetweenchop4andatt_w2;
MCNUM mccbetweenchop4andatt_h2;
MCNUM mccbetweenchop4andatt_l;
MCNUM mccbetweenchop4andatt_R0;
MCNUM mccbetweenchop4andatt_Qc;
MCNUM mccbetweenchop4andatt_alpha;
MCNUM mccbetweenchop4andatt_m;
MCNUM mccbetweenchop4andatt_nslit;
MCNUM mccbetweenchop4andatt_d;
MCNUM mccbetweenchop4andatt_Qcx;
MCNUM mccbetweenchop4andatt_Qcy;
MCNUM mccbetweenchop4andatt_alphax;
MCNUM mccbetweenchop4andatt_alphay;
MCNUM mccbetweenchop4andatt_W;
MCNUM mccbetweenchop4andatt_mx;
MCNUM mccbetweenchop4andatt_my;
MCNUM mccbetweenchop4andatt_nu;
MCNUM mccbetweenchop4andatt_phase;

/* Setting parameters for component 'guide_move' [14]. */
MCNUM mccguide_move_w1;
MCNUM mccguide_move_h1;
MCNUM mccguide_move_w2;
MCNUM mccguide_move_h2;
MCNUM mccguide_move_l;
MCNUM mccguide_move_R0;
MCNUM mccguide_move_Qc;
MCNUM mccguide_move_alpha;
MCNUM mccguide_move_m;
MCNUM mccguide_move_nslit;
MCNUM mccguide_move_d;
MCNUM mccguide_move_Qcx;
MCNUM mccguide_move_Qcy;
MCNUM mccguide_move_alphax;
MCNUM mccguide_move_alphay;
MCNUM mccguide_move_W;
MCNUM mccguide_move_mx;
MCNUM mccguide_move_my;
MCNUM mccguide_move_nu;
MCNUM mccguide_move_phase;

/* Setting parameters for component 'funnel' [15]. */
MCNUM mccfunnel_w1;
MCNUM mccfunnel_h1;
MCNUM mccfunnel_w2;
MCNUM mccfunnel_h2;
MCNUM mccfunnel_l;
MCNUM mccfunnel_R0;
MCNUM mccfunnel_Qc;
MCNUM mccfunnel_alpha;
MCNUM mccfunnel_m;
MCNUM mccfunnel_nslit;
MCNUM mccfunnel_d;
MCNUM mccfunnel_Qcx;
MCNUM mccfunnel_Qcy;
MCNUM mccfunnel_alphax;
MCNUM mccfunnel_alphay;
MCNUM mccfunnel_W;
MCNUM mccfunnel_mx;
MCNUM mccfunnel_my;
MCNUM mccfunnel_nu;
MCNUM mccfunnel_phase;

/* Setting parameters for component 'C5' [16]. */
MCNUM mccC5_theta_0;
MCNUM mccC5_radius;
MCNUM mccC5_yheight;
MCNUM mccC5_nu;
MCNUM mccC5_nslit;
MCNUM mccC5_jitter;
MCNUM mccC5_delay;
MCNUM mccC5_isfirst;
MCNUM mccC5_n_pulse;
MCNUM mccC5_abs_out;
MCNUM mccC5_phase;
MCNUM mccC5_xwidth;
MCNUM mccC5_verbose;

/* Setting parameters for component 'C5_counter' [17]. */
MCNUM mccC5_counter_theta_0;
MCNUM mccC5_counter_radius;
MCNUM mccC5_counter_yheight;
MCNUM mccC5_counter_nu;
MCNUM mccC5_counter_nslit;
MCNUM mccC5_counter_jitter;
MCNUM mccC5_counter_delay;
MCNUM mccC5_counter_isfirst;
MCNUM mccC5_counter_n_pulse;
MCNUM mccC5_counter_abs_out;
MCNUM mccC5_counter_phase;
MCNUM mccC5_counter_xwidth;
MCNUM mccC5_counter_verbose;

/* Setting parameters for component 'endguide' [18]. */
MCNUM mccendguide_w1;
MCNUM mccendguide_h1;
MCNUM mccendguide_w2;
MCNUM mccendguide_h2;
MCNUM mccendguide_l;
MCNUM mccendguide_R0;
MCNUM mccendguide_Qc;
MCNUM mccendguide_alpha;
MCNUM mccendguide_m;
MCNUM mccendguide_nslit;
MCNUM mccendguide_d;
MCNUM mccendguide_Qcx;
MCNUM mccendguide_Qcy;
MCNUM mccendguide_alphax;
MCNUM mccendguide_alphay;
MCNUM mccendguide_W;
MCNUM mccendguide_mx;
MCNUM mccendguide_my;
MCNUM mccendguide_nu;
MCNUM mccendguide_phase;

/* Setting parameters for component 'snot' [19]. */
MCNUM mccsnot_w1;
MCNUM mccsnot_h1;
MCNUM mccsnot_w2;
MCNUM mccsnot_h2;
MCNUM mccsnot_l;
MCNUM mccsnot_R0;
MCNUM mccsnot_Qc;
MCNUM mccsnot_alpha;
MCNUM mccsnot_m;
MCNUM mccsnot_nslit;
MCNUM mccsnot_d;
MCNUM mccsnot_Qcx;
MCNUM mccsnot_Qcy;
MCNUM mccsnot_alphax;
MCNUM mccsnot_alphay;
MCNUM mccsnot_W;
MCNUM mccsnot_mx;
MCNUM mccsnot_my;
MCNUM mccsnot_nu;
MCNUM mccsnot_phase;

/* Setting parameters for component 'detslit' [20]. */
MCNUM mccdetslit_xmin;
MCNUM mccdetslit_xmax;
MCNUM mccdetslit_ymin;
MCNUM mccdetslit_ymax;
MCNUM mccdetslit_radius;
MCNUM mccdetslit_xwidth;
MCNUM mccdetslit_yheight;

/* Definition parameters for component 'mcplout' [21]. */
#define mccmcplout_polarisationuse 0
#define mccmcplout_doubleprec 0
#define mccmcplout_verbose 0
#define mccmcplout_radius 0.021
#define mccmcplout_thickness 0.005
#define mccmcplout_yheight 0.05
#define mccmcplout_restore_neutron 1
#define mccmcplout_userflag 0
#define mccmcplout_plate plate
/* Setting parameters for component 'mcplout' [21]. */
char mccmcplout_filename[16384];
char mccmcplout_userflagcomment[16384];
MCNUM mccmcplout_merge_mpi;
MCNUM mccmcplout_keep_mpi_unmerged;

/* User component declarations. */

/* User declarations for component 'Origin' [1]. */
#define mccompcurname  Origin
#define mccompcurtype  Arm
#define mccompcurindex 1
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Source' [2]. */
#define mccompcurname  Source
#define mccompcurtype  ViewModISISver1
#define mccompcurindex 2
#define Face mccSource_Face
#define E0 mccSource_E0
#define E1 mccSource_E1
#define tally mccSource_tally
#define modPosition mccSource_modPosition
#define modXsize mccSource_modXsize
#define modZsize mccSource_modZsize
#define xw mccSource_xw
#define yh mccSource_yh
#define dist mccSource_dist
#line 90 "ViewModISISver1.comp"

#include <ctype.h>

typedef struct
{
  int nEnergy;        ///< Number of energy bins
  int nTime;          ///< number of time bins

  double XAxis;
  double ZAxis;

  double rdumMid;     ///< tally time Window mid point
  double timeOffset;  ///< Time separation
  double* TimeBin;    ///< Time bins
  double* EnergyBin;  ///< Energy bins

  double** Flux;       ///< Flux per bin (integrated)
  double* EInt;        ///< Integrated Energy point
  double Total;        ///< Integrated Total

} Source;

  /* New functions */

  int cmdnumberD(char *,double*);
  int cmdnumberI(char *,int*,const int);
  double polInterp(double*,double*,int,double);
  int findRDUM(char*);
  int findTimeOffset(char*,double*);
  FILE* openFile(char*);
  FILE* openFileTest(char*);
  void processHeader(FILE*);
  void readHtable(FILE*,const double,const double);
  int timeStart(char*);
  int timeEnd(char*);
  int energyBin(char*,double,double,double*,double*);
  int notComment(char*);
  void orderEnergy(double*,double*);
  void cutToNumber(char*);
  double convertEnergy(double);
  double EtoLambda(double);
  double calcRDum(double*);

  double strArea();

  /* global variables */

  double p_in;            /* Polorization term (from McSTAS) */
  int Tnpts;              /* Number of points in parameteriation */

  double scaleSize;        /* correction for the actual area of the moderator viewed */

  double angleArea;       /* Area seen by the window  */

  double Nsim;		/* Total number of neutrons to be simulated */

  int Ncount;          /* Number of neutron simulate so far*/
  Source TS;

  /* runtime variables*/
  double fullAngle;
  double rtE0,rtE1;       /* runtime Energy minima and maxima so we can use angstroms as negative input */


double**
matrix(const int m,const int n)
 /*!
   Determine a double matrix
 */
{
  int i;
  double* pv;
  double** pd;

  if (m<1) return 0;
  if (n<1) return 0;
  pv = (double*) malloc(m*n*sizeof(double));
  pd = (double**) malloc(m*sizeof(double*));
  if (!pd)
    {
      fprintf(stderr,"No room for matrix!\n");
      exit(1);
    }
  for (i=0;i<m;i++)
    pd[i]=pv + (i*n);
  return pd;
}


double
polInterp(double* X,double* Y,int Psize,double Aim)
  /*!
    returns the interpolated polynomial between Epnts
    and the integration
    \param X :: X coordinates
    \param Y :: Y coordinates
    \param Psize :: number of valid point in array to use
    \param Aim :: Aim point to intepolate result (X coordinate)
    \returns Energy point
  */
{
  double out,errOut;         /* out put variables */
  double C[Psize],D[Psize];
  double testDiff,diff;

  double w,den,ho,hp;           /* intermediate variables */
  int i,m,ns;


  ns=0;
  diff=fabs(Aim-X[0]);
  C[0]=Y[0];
  D[0]=Y[0];
  for(i=1;i<Psize;i++)
    {
      testDiff=fabs(Aim-X[i]);
      if (diff>testDiff)
	{
	  ns=i;
	  diff=testDiff;
	}
      C[i]=Y[i];
      D[i]=Y[i];
    }

  out=Y[ns];
  ns--;              /* Now can be -1 !!!! */

  for(m=1;m<Psize;m++)
    {
      for(i=0;i<Psize-m;i++)
	{
	  ho=X[i]-Aim;
	  hp=X[i+m]-Aim;
	  w=C[i+1]-D[i];
	  /*	  den=ho-hp;  -- test !=0.0 */
	  den=w/(ho-hp);
	  D[i]=hp*den;
	  C[i]=ho*den;
	}

      errOut= (2*(ns+1)<(Psize-m)) ? C[ns+1] : D[ns--];
      out+=errOut;
    }
  return out;
}

int
binSearch(int Npts,double* AR,double V)
  /*!
    Object is to find the point in
    array AR, closest to the value V
    Checked for ordered array returns lower of backeting objects
  */
{
  int klo,khi,k;
  if (Npts<=0)
    return 0;
  if (V>AR[Npts-1])
    return Npts;

  if(AR[0]>0.0)AR[0]=0.0;

  if (V<AR[0])
    {
      // if(AR[0]>0.0)AR[0]=0.0;
      fprintf(stderr,"here");
      return 0;
    }
  klo=0;
  khi= Npts-1;
  while (khi-klo >1)
    {
      k=(khi+klo) >> 1;    // quick division by 2
      if (AR[k]>V)
	khi=k;
      else
	klo=k;
    }
  return khi;
}

int
cmdnumberD(char *mc,double* num)
 /*!
   \returns 1 on success 0 on failure
 */
{
  int i,j;
  char* ss;
  char **endptr;
  double nmb;
  int len;

  len=strlen(mc);
  j=0;

  for(i=0;i<len && mc[i] &&
	(mc[i]=='\t' || mc[i]==' '  || mc[i]==',');i++);
  if(i==len || !mc[i]) return 0;
  ss=malloc(sizeof(char)*(len+1));

  for(;i<len && mc[i]!='\n' && mc[i]
	&& mc[i]!='\t' && mc[i]!=' ' && mc[i]!=',';i++)
    {
      ss[j]=mc[i];
      j++;
    }
  if (!j)
    {
      free(ss);
      return 0;         //This should be impossible
    }
  ss[j]=0;
  endptr=malloc(sizeof(char*));
  nmb = strtod(ss,endptr);
  if (*endptr != ss+j)
    {
      free(endptr);
      free(ss);
      return 0;
    }
  *num = (double) nmb;
  for(j=0;j<i && mc[j];j++)
    mc[j]=' ';
  free(endptr);
  free(ss);
  return 1;
}

int
notComment(char* Line)
 /*!
   \returns 0 on a comment, 1 on a non-comment
 */
{
  int len,i;

  len=strlen(Line);
  for(i=0;i<len && isspace(Line[i]);i++);

  if (!Line[i] || Line[i]=='c' || Line[i]=='C' ||
      Line[i]=='!' || Line[i]=='#')
    return 0;
  return 1;
}

void
orderEnergy(double* A,double *B)
{
  double tmp;
  if (*A>*B)
    {
      tmp=*A;
      *A=*B;
      *B=tmp;
    }
  return;
}

int
timeStart(char* Line)
 /*!
   Search for a word time at the start of
   the line.
   \param Line :: Line to search
   \returns 1 on success 0 on failure
 */
{
  int len,i;

  len=strlen(Line);
  for(i=0;i<len && isspace(Line[i]);i++);
  if (len-i<4) return 0;
  return (strncmp(Line+i,"time",4)) ? 0 : 1;
}

void
cutToNumber(char* Line)
{
  int len,i;

  len=strlen(Line);
  for(i=0;i<len && !isdigit(Line[i]) && Line[i]!='-';i++)
    Line[i]=' ';
  return;
}

int
findTimeOffset(char* Line,double* TO)
  /*!
    Determine the time offset from the file
    \return 1 on success
   */
{
  int len,i;
  double D;

  len=strlen(Line);
  for(i=0;i<len && Line[i]!='T';i++);
  if (len-i<11) return 0;

  if (strncmp(Line+i,"TimeOffset",11) &&
      cmdnumberD(Line+i+11,&D))
    {
      *TO=D/100.0;
      return 1;
    }
  return 0;
}

int
findRDUM(char* Line)
 /*!
   Search for a word rdum at the start of
   the line.
   \param Line :: Line to search
   \returns 1 on success 0 on failure
 */
{
  int len,i;

  len=strlen(Line);
  for(i=0;i<len && Line[i]!='r';i++);
  if (len-i<4) return 0;
  return (strncmp(Line+i,"rdum",4)) ? 0 : 1;
}

double 
convertEnergy(double E)
  /*!
    Convert the energy from eV [not change] or 
    from negative -ve which is angstrom
  */
{
  return (E>0.0) ? E : 81.793936/(E*E);
}

double 
EtoLambda(double E)
  /*!
    Convert the energy from eV [not change] o
    to lambda [A]
  */
{
  return sqrt(81.793936/E);
}
  

int
timeEnd(char* Line)
 /*!
   Search for a word time at the start of
   the line.
   \param Line :: Line to search
   \returns 1 on success 0 on failure
 */
{
  int len,i;

  len=strlen(Line);
  for(i=0;i<len && isspace(Line[i]);i++);
  if (len-i<5) return 0;
  return (strncmp(Line+i,"total",5)) ? 0 : 1;
}

int
energyBin(char* Line,double Einit,double Eend,double* Ea,double* Eb)
     /*!
       Search for a word "energy bin:" at the start of
       the line. Then separte off the energy bin values
       \param Line :: Line to search
       \param Ea :: first energy bin [meV]
       \param Eb :: second energy bin [meV]
       \returns 1 on success 0 on failure
     */
{
  int len,i;
  double A,B;

  len=strlen(Line);
  for(i=0;i<len && isspace(Line[i]);i++);
  if (len-i<11) return 0;


  if (strncmp(Line+i,"energy bin:",11))
    return 0;

  i+=11;
  if (!cmdnumberD(Line+i,&A))
    return 0;
  // remove 'to'
  for(;i<len-1 && Line[i]!='o';i++);
  i++;
  if (!cmdnumberD(Line+i,&B))
    return 0;
  A*=1e9;
  B*=1e9;
  *Ea=A;
  *Eb=B;
  if (*Eb>Einit && *Ea<Eend)
    return 1;
  return 0;
}

double
calcFraction(double EI,double EE,double Ea,double Eb)
 /*!
   Calculate the fraction of the bin between Ea -> Eb
   that is encompassed by EI->EE
 */
{
  double frac;
  double dRange;

  if (EI>Eb)
    return 0.0;
  if (EE<Ea)
    return 0.0;

  dRange=Eb-Ea;
  frac=(EI>Ea) ? (Eb-EI)/dRange : 1.0;


  frac-=(EE<Eb) ? (Eb-EE)/dRange : 0.0;

  //  if(frac != 1.0)
  //  fprintf(stderr,"frac %g, Ea %g,Eb %g, EI %g, EE %g\n",frac,Ea,Eb,EI,EE);

  return frac;
}

double
calcRDum(double* RPts)
  /*!
    Caluclate the mean distance from the window centre point
   */
{
  double sum,zMid;
  int i;
  extern Source TS;

  sum=0.0;
  for(i=0;i<4;i++)
    {
      fprintf(stderr,"RDUM %g %g %g\n",RPts[i*3],RPts[i*3+1],RPts[i*3+2]);
      sum+=sqrt(RPts[i*3]*RPts[i*3]+
		RPts[i*3+1]*RPts[i*3+1]);
    }
  sum/=400.0;   // Convert to metres from cm
  return sum;
}

void
processHeader(FILE* TFile)
  /*!
    Determine the window and the time offset
    \param TFile :: file [rewound on exit]
  */
{
  char ss[255];          /* BIG space for line */
  int rdumCnt;
  double Ea,Eb;
  double RPts[12];
  int timeOffsetFlag;
  int i,j;
  double D;
  

  extern Source TS;
  
  rdumCnt=0;
  timeOffsetFlag=0;

  while(fgets(ss,255,TFile))
    {
      if (!timeOffsetFlag)
	timeOffsetFlag=findTimeOffset(ss,&TS.timeOffset);

      if (!rdumCnt)
	rdumCnt=findRDUM(ss);
      if (rdumCnt && rdumCnt<5)
	{
	  cutToNumber(ss);
	  for(i=0;i<3 && cmdnumberD(ss,&D);i++)
	    RPts[(rdumCnt-1)*3+i]=D;
	  rdumCnt++;
	}
      // EXIT CONDITION:
      if (rdumCnt*timeOffsetFlag==5)
	{
	  for(j=0;j<3;j++)
	    {
	      TS.ZAxis+=(RPts[3+j]-RPts[j])*(RPts[3+j]-RPts[j]);
	      TS.XAxis+=(RPts[6+j]-RPts[3+j])*(RPts[6+j]-RPts[3+j]);
	    }	


	  TS.ZAxis=sqrt(TS.ZAxis)/100.0;   // convert to metres from cm
	  TS.XAxis=sqrt(TS.XAxis)/100.0;
              if (!modPosition)
                 {
	          TS.ZAxis=modZsize;
	          TS.XAxis=modXsize;
                  }
	  fprintf(stderr,"Time off sec == %g \n",TS.timeOffset);
	  fprintf(stderr,"mod size z == %g \n",TS.ZAxis);
	  // TS.rdumMid=calcRDum(RPts);
	  TS.rdumMid=TS.timeOffset; // Goran
	  return;
	}
    }
  return;
}



void
readHtable(FILE* TFile,const double Einit,const double Eend)
/*!
  Process a general h.o file to create an integrated
  table of results from Einit -> Eend
  \param Einit :: inital Energy
  \param Eend  :: final energy
*/
{
  char ss[255];          /* BIG space for line */
  double Ea,Eb;
  double T,D;
  double Efrac;          // Fraction of an Energy Bin
  int Ftime;             // time Flag
  int eIndex;             // energy Index
  int tIndex;             // time Index
  double Tsum;           // Running integration
  double Efraction;      // Amount to use for an energy/time bin

  extern Source TS;

  int i;
  /*!
    Status Flag::
    Ftime=1 :: [time ] Reading Time : Data : Err [Exit on Total]

    Double Read File to determine how many bins and
    memery size
  */
  Ea=0.0;
  Eb=0.0;
  fprintf(stderr,"Energy == %g %g\n",Einit,Eend);
  eIndex= -1;
  Ftime=0;
  tIndex=0;
  TS.nTime=0;
  TS.nEnergy=0;
  
  processHeader(TFile);

  // Read file and get time bins:
  while(fgets(ss,255,TFile) && Eend>Ea)
    {
      if (notComment(ss))
	{
          if (!Ftime)
	    {
	      // find :: energy bin: <Number> to <Number> 
	      if (energyBin(ss,Einit,Eend,&Ea,&Eb))  
		{
		  if (eIndex==0)
		    TS.nTime=tIndex;
		  eIndex++;
		}
	      else if (timeStart(ss))    // find :: <spc> time 
		{
		  Ftime=1;
		  tIndex=0;
		}
	    }
	  else  // In the time data section
	    {
	      if (timeEnd(ss))     // Found "total"
		Ftime=0;
	      else
		{
		  // Need to read the line in the case of first run
		  if (TS.nTime==0)
		    {
		      if (cmdnumberD(ss,&T) &&
			  cmdnumberD(ss,&D))
			tIndex++;
		    }
		}
	    }
	}
    }
  // 
  // Plus 2 since we have a 0 counter and we have missed the last line.
  // 
  TS.nEnergy=eIndex+2;
  if (!TS.nTime && tIndex)  // SMALL TEST if reading first energy bin
    TS.nTime=tIndex;
  // printf("tIndex %d %d %d %d \n",tIndex,eIndex,TS.nEnergy,TS.nTime);

  /* SECOND TIME THROUGH:: */
  rewind(TFile);

  TS.Flux=matrix(TS.nEnergy,TS.nTime);
  TS.EInt=(double*) malloc(TS.nEnergy*sizeof(double));     // Runing integral of energy 
  TS.TimeBin=(double*) malloc(TS.nTime*sizeof(double));
  TS.EnergyBin=(double*) malloc(TS.nEnergy*sizeof(double));
  fprintf(stderr,"NUMBER %d %d \n",TS.nEnergy,TS.nTime);

  Tsum=0.0;
  Ea=0.0;
  Eb=0.0;
  eIndex=-1;
  Ftime=0;
  tIndex=0;
  TS.EInt[0]=0.0;

  // Read file and get time bins
  while(fgets(ss,255,TFile) && Eend>Ea)
    {
      if (notComment(ss))
	{
          if (!Ftime)
	    {
	      if (energyBin(ss,Einit,Eend,&Ea,&Eb))
		{
		  eIndex++;
		  TS.EnergyBin[eIndex]=(Einit>Ea) ? Einit : Ea;
		  Efraction=calcFraction(Einit,Eend,Ea,Eb);
		  Ftime++;
		}
	    }
	  else if (Ftime==1)
	    {
	      if (timeStart(ss))
		{
		  Ftime=2;
		  tIndex=0;
		}
	    }
	  else           // In the time section
	    {
	      if (timeEnd(ss))     // Found "total"
		{
		  Ftime=0;
		  TS.EInt[eIndex+1]=Tsum;

		  fprintf(stderr,"Energy %g %g \n",EtoLambda(Ea),(TS.EInt[eIndex+1]-TS.EInt[eIndex])*3.744905847e14);
		}
	      else
		{
		  // Need to read the line in the case of first run
		  if (cmdnumberD(ss,&T) &&
		      cmdnumberD(ss,&D))
		    {
		      TS.TimeBin[tIndex]=T/1e8;     // convert Time into second (from shakes)
		      Tsum+=D*Efraction;
		      TS.Flux[eIndex][tIndex]=Tsum;
		      tIndex++;
		    }
		}
	    }
	}
    }
  TS.EnergyBin[eIndex+1]=Eend;
  TS.Total=Tsum;

  // printf("tIndex %d %d %d \n",tIndex,eIndex,TS.nTime);
  // printf("Tsum %g \n",Tsum);
  // fprintf(stderr,"ebin1 ebinN %g %g\n",TS.EnergyBin[0],TS.EnergyBin[TS.nEnergy-1]);

  return;
}

void
getPoint(double* TV,double* EV,double* lim1, double* lim2)
 /*!
   Calculate the Time and Energy
   by sampling the file.
   Uses TS table to find the point
   \param TV ::
   \param EV ::
   \param lim1 ::
   \param lim2 ::
 */
{
  int i;

  extern Source TS;
  double R0,R1,R,Rend;
  int Epnt;       ///< Points to the next higher index of the neutron integral
  int Tpnt;
  int iStart,iEnd;
  double TRange,Tspread;
  double Espread,Estart;
  double *EX;

  // So that lowPoly+highPoly==maxPoly
  const int maxPoly=6;
  const int highPoly=maxPoly/2;
  const int lowPoly=maxPoly-highPoly;

  // static int testVar=0;

  R0=rand01();
  Rend=R=TS.Total*R0;
  // This gives Eint[Epnt-1] > R > Eint[Epnt]
  Epnt=binSearch(TS.nEnergy-1,TS.EInt,R);
  Tpnt=binSearch(TS.nTime-1,TS.Flux[Epnt-1],R);

  if(Epnt<0 || R<TS.Flux[Epnt-1][Tpnt-1] || R >TS.Flux[Epnt-1][Tpnt] )
    {
      fprintf(stderr,"outside bin limits Tpnt/Epnt problem  %12.6e %12.6e %12.6e \n",
      TS.Flux[Epnt-1][Tpnt-1],R,TS.Flux[Epnt-1][Tpnt]);
    }

  if(Epnt == 0)
    {
      Estart=0.0;
      Espread=TS.EInt[0];
      *EV=TS.EnergyBin[1];
    }
  else
    {
      Estart=TS.EInt[Epnt-1];
      Espread=TS.EInt[Epnt]-TS.EInt[Epnt-1];
      *EV=(Epnt>TS.nEnergy-1) ? TS.EnergyBin[Epnt+1] : TS.EnergyBin[Epnt];
    }

  if (Tpnt==0 || Epnt==0 || Tpnt==TS.nTime)
    {
      fprintf(stderr,"BIG ERROR WITH Tpnt: %d and Epnt: %d\n",Tpnt,Epnt);
      exit(1);
    }
  *TV=TS.TimeBin[Tpnt-1];
  TRange=TS.TimeBin[Tpnt]-TS.TimeBin[Tpnt-1];
  Tspread=TS.Flux[Epnt-1][Tpnt]-TS.Flux[Epnt-1][Tpnt-1];
  R-=TS.Flux[Epnt-1][Tpnt-1];

  R/=Tspread;
  *TV+=TRange*R;


  R1=TS.EInt[Epnt-1]+Espread*rand01();
  iStart=Epnt>lowPoly ? Epnt-lowPoly : 0;                  // max(Epnt-halfPoly,0)
  iEnd=TS.nEnergy>Epnt+highPoly ? Epnt+highPoly : TS.nEnergy-1;  // min(nEnergy-1,Epnt+highPoly

  *EV=polInterp(TS.EInt+iStart,TS.EnergyBin+iStart,1+iEnd-iStart,R1);

  if(*TV < TS.TimeBin[Tpnt-1] || *TV > TS.TimeBin[Tpnt])
    {
      fprintf(stderr,"%d Tpnt %d Tval %g Epnt %d \n",TS.nTime,Tpnt,*TV,Epnt);
      fprintf(stderr,"TBoundary == %12.6e,%g , %12.6e \n\n",TS.TimeBin[Tpnt-1],*TV,TS.TimeBin[Tpnt]);
    }
  

  if(*EV < *lim1 || *EV > *lim2)
    {
      fprintf(stderr,"outside boundaries\n Epnt= %d, Tpnt= %d binlo %g|%g| binhi %g \n",
	      Epnt,Tpnt,TS.EnergyBin[Epnt-1],*EV,TS.EnergyBin[Epnt]);

      fprintf(stderr,"TS == %g %g :: %d %d \n",TS.EInt[Epnt-1],TS.EInt[Epnt],iStart,iEnd);
      fprintf(stderr,"Points (%g) == ",R1);
      for(i=0;i<iEnd-iStart;i++)
	fprintf(stderr,"Time start %g %g",*(TS.EInt+i+iStart),*(TS.EnergyBin+iStart+i));
      fprintf(stderr,"\n");
      exit(1);
    }
  return;
}

int
cmdnumberI(char *mc,int* num,const int len)
  /*!
    \param mc == character string to use
    \param num :: Place to put output
    \param len == length of the character string to process
    returns 1 on success and 0 on failure
  */
{
  int i,j;
  char* ss;
  char **endptr;
  double nmb;

  if (len<1)
    return 0;
  j=0;
  
  for(i=0;i<len && mc[i] &&
	(mc[i]=='\t' || mc[i]==' '  || mc[i]==',');i++);
  if(i==len || !mc[i]) return 0;
  ss=malloc(sizeof(char)*(len+1));
  /*  char *ss=new char[len+1]; */
  for(;i<len && mc[i]!='\n' && mc[i]
	&& mc[i]!='\t' && mc[i]!=' ' && mc[i]!=',';i++)
    {
      ss[j]=mc[i];
      j++;
    }
  if (!j)
    {
      free(ss);
      return 0;         //This should be impossible
    }
  ss[j]=0;
  endptr=malloc(sizeof(char*));
  nmb = strtod(ss,endptr);
  if (*endptr != ss+j)
    {
      free(endptr);
      free(ss);
      return 0;
    }
  *num = (double) nmb;
  for(j=0;j<i && mc[j];j++)
    mc[j]=' ';
  free(endptr);
  free(ss);
  return 1;
}


FILE* 
openFile(char* FileName)
/*
  Tries to open the file:
  (i) In current working directory 
  (ii) In MC_Path directory + ISIS_tables extension
*/
{
  FILE* efile=0;
  int fLen;
  char extFileName[1024];
  char testFileName[2048];
  char mct[2048];
  
  fLen=strlen(FileName);

  if (fLen>512) 
    {
      fprintf(stderr,"Filename excessivley long [EXIT]:\n %s\n",FileName);
      exit(1);
    }
  
  
  strcpy(extFileName,FileName);
  strcpy(extFileName+fLen,".mcstas");
  

  /* Is the file located in working dir? */
  efile=fopen(FileName,"r");
  if (efile) return efile; 
  
  efile=fopen(extFileName,"r");
  if (efile) return efile; 
  
  /* Is the file in a local 'tablefiles' folder? */
  sprintf(testFileName,"%c%s%s", MC_PATHSEP_C,"ISIS_tables",FileName);
  efile=fopen(testFileName,"r");
  if (efile) return efile; 


  /* Is MCTABLES set, files located there? */
  if (getenv("MCTABLES")) 
    {
      strcpy(mct, getenv("MCTABLES"));
      sprintf(testFileName, "%s%c%s", mct, MC_PATHSEP_C, FileName);
      efile=fopen(testFileName,"r");
      if (efile) return efile;
    }

  fprintf(stderr,"Searching -- %s, MCTABLES directory\n",testFileName);
  fprintf(stderr,"\nPlease check your McStas installation and/or MCTABLES environment variable!\n");
  exit(1);
  return efile;
}

double strArea()
{
  /*
    Returns the mean Str view of the viewport
    This integrates over each point on the window xw to yh
    View port is symmetric so use only 1/4 of the view
    for the calcuation.
    Control Values TS.XAxis TS.ZAxis xw yh
  */
  
  double A;
  double Vx,Vy;        // view temp points
  double Mx,My;        // moderator x,y
  double D2;           // Distance ^2
  double projArea;     // viewport projection to moderator
  int i,j,aa,bb;       // loop variables
  int NStep;

  NStep=50;
  D2=dist*dist;
  A=0.0;
  fprintf(stderr,"TS axis == %g %g\n",TS.XAxis,TS.ZAxis);
  fprintf(stderr,"AW axis == %g %g\n",xw,yh);
  for(i=0;i<NStep;i++)              // Mod X
    {
      Mx=(i+0.5)*TS.XAxis/(NStep*2.0);
      for(j=0;j<NStep;j++)         // Mod Y
	{
	  My=(j+0.5)*TS.ZAxis/(NStep*2.0);
	  // Position on moderator == (Mx,My)
	  for(aa= -NStep;aa<=NStep;aa++)  //view port
	    for(bb= -NStep;bb<=NStep;bb++)
	      {
		Vx=aa*xw/(2.0*NStep+1.0);
		Vy=bb*yh/(2.0*NStep+1.0);
		A+=1.0/((Mx-Vx)*(Mx-Vx)+(My-Vy)*(My-Vy)+D2);
	      }
	}
    }
  //change to Mx*My 
  A*= (TS.XAxis*TS.ZAxis)/(NStep*NStep*(2.0*NStep+1.0)*(2.0*NStep+1.0));
  // Correct for the area of the viewport. (tables are per cm^2)
  A*=xw*yh*10000;

  //   testing only - Goran              
    // if (!modPosition)
        // {
            // projArea=xw*yh/tally/tally*(tally+dist)*(tally+dist);
	    // A*=TS.XAxis*TS.ZAxis/(TS.XAxis*TS.ZAxis-projArea);
         // }  

  fprintf(stderr,"Viewport == %g %g Moderator size == (%g * %g) m^2 \n",xw,yh,TS.XAxis,TS.ZAxis);
  fprintf(stderr,"Dist == %g (metres) \n",dist);
  fprintf(stderr,"Viewport Solid angle == %g str\n",A/(xw*yh*10000.0));
  fprintf(stderr,"Solid angle used == %g str\n",A);
  return A;
}

#line 8225 "let_new.c"
#undef dist
#undef yh
#undef xw
#undef modZsize
#undef modXsize
#undef modPosition
#undef tally
#undef E1
#undef E0
#undef Face
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'shutter' [3]. */
#define mccompcurname  shutter
#define mccompcurtype  Guide_channeled
#define mccompcurindex 3
#define w1c mccshutter_w1c
#define w2c mccshutter_w2c
#define ww mccshutter_ww
#define hh mccshutter_hh
#define whalf mccshutter_whalf
#define hhalf mccshutter_hhalf
#define lwhalf mccshutter_lwhalf
#define lhhalf mccshutter_lhhalf
#define w1 mccshutter_w1
#define h1 mccshutter_h1
#define w2 mccshutter_w2
#define h2 mccshutter_h2
#define l mccshutter_l
#define R0 mccshutter_R0
#define Qc mccshutter_Qc
#define alpha mccshutter_alpha
#define m mccshutter_m
#define nslit mccshutter_nslit
#define d mccshutter_d
#define Qcx mccshutter_Qcx
#define Qcy mccshutter_Qcy
#define alphax mccshutter_alphax
#define alphay mccshutter_alphay
#define W mccshutter_W
#define mx mccshutter_mx
#define my mccshutter_my
#define nu mccshutter_nu
#define phase mccshutter_phase
#line 80 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 8278 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'insert' [4]. */
#define mccompcurname  insert
#define mccompcurtype  Guide_channeled
#define mccompcurindex 4
#define w1c mccinsert_w1c
#define w2c mccinsert_w2c
#define ww mccinsert_ww
#define hh mccinsert_hh
#define whalf mccinsert_whalf
#define hhalf mccinsert_hhalf
#define lwhalf mccinsert_lwhalf
#define lhhalf mccinsert_lhhalf
#define w1 mccinsert_w1
#define h1 mccinsert_h1
#define w2 mccinsert_w2
#define h2 mccinsert_h2
#define l mccinsert_l
#define R0 mccinsert_R0
#define Qc mccinsert_Qc
#define alpha mccinsert_alpha
#define m mccinsert_m
#define nslit mccinsert_nslit
#define d mccinsert_d
#define Qcx mccinsert_Qcx
#define Qcy mccinsert_Qcy
#define alphax mccinsert_alphax
#define alphay mccinsert_alphay
#define W mccinsert_W
#define mx mccinsert_mx
#define my mccinsert_my
#define nu mccinsert_nu
#define phase mccinsert_phase
#line 80 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 8349 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'afterinsert' [5]. */
#define mccompcurname  afterinsert
#define mccompcurtype  Guide_channeled
#define mccompcurindex 5
#define w1c mccafterinsert_w1c
#define w2c mccafterinsert_w2c
#define ww mccafterinsert_ww
#define hh mccafterinsert_hh
#define whalf mccafterinsert_whalf
#define hhalf mccafterinsert_hhalf
#define lwhalf mccafterinsert_lwhalf
#define lhhalf mccafterinsert_lhhalf
#define w1 mccafterinsert_w1
#define h1 mccafterinsert_h1
#define w2 mccafterinsert_w2
#define h2 mccafterinsert_h2
#define l mccafterinsert_l
#define R0 mccafterinsert_R0
#define Qc mccafterinsert_Qc
#define alpha mccafterinsert_alpha
#define m mccafterinsert_m
#define nslit mccafterinsert_nslit
#define d mccafterinsert_d
#define Qcx mccafterinsert_Qcx
#define Qcy mccafterinsert_Qcy
#define alphax mccafterinsert_alphax
#define alphay mccafterinsert_alphay
#define W mccafterinsert_W
#define mx mccafterinsert_mx
#define my mccafterinsert_my
#define nu mccafterinsert_nu
#define phase mccafterinsert_phase
#line 80 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 8420 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Res1' [6]. */
#define mccompcurname  Res1
#define mccompcurtype  DiskChopper
#define mccompcurindex 6
#define Tg mccRes1_Tg
#define To mccRes1_To
#define delta_y mccRes1_delta_y
#define height mccRes1_height
#define omega mccRes1_omega
#define theta_0 mccRes1_theta_0
#define radius mccRes1_radius
#define yheight mccRes1_yheight
#define nu mccRes1_nu
#define nslit mccRes1_nslit
#define jitter mccRes1_jitter
#define delay mccRes1_delay
#define isfirst mccRes1_isfirst
#define n_pulse mccRes1_n_pulse
#define abs_out mccRes1_abs_out
#define phase mccRes1_phase
#define xwidth mccRes1_xwidth
#define verbose mccRes1_verbose
#line 63 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 8477 "let_new.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'Res1_counter' [7]. */
#define mccompcurname  Res1_counter
#define mccompcurtype  DiskChopper
#define mccompcurindex 7
#define Tg mccRes1_counter_Tg
#define To mccRes1_counter_To
#define delta_y mccRes1_counter_delta_y
#define height mccRes1_counter_height
#define omega mccRes1_counter_omega
#define theta_0 mccRes1_counter_theta_0
#define radius mccRes1_counter_radius
#define yheight mccRes1_counter_yheight
#define nu mccRes1_counter_nu
#define nslit mccRes1_counter_nslit
#define jitter mccRes1_counter_jitter
#define delay mccRes1_counter_delay
#define isfirst mccRes1_counter_isfirst
#define n_pulse mccRes1_counter_n_pulse
#define abs_out mccRes1_counter_abs_out
#define phase mccRes1_counter_phase
#define xwidth mccRes1_counter_xwidth
#define verbose mccRes1_counter_verbose
#line 63 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 8524 "let_new.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'betweenchop1and2' [8]. */
#define mccompcurname  betweenchop1and2
#define mccompcurtype  Guide_channeled
#define mccompcurindex 8
#define w1c mccbetweenchop1and2_w1c
#define w2c mccbetweenchop1and2_w2c
#define ww mccbetweenchop1and2_ww
#define hh mccbetweenchop1and2_hh
#define whalf mccbetweenchop1and2_whalf
#define hhalf mccbetweenchop1and2_hhalf
#define lwhalf mccbetweenchop1and2_lwhalf
#define lhhalf mccbetweenchop1and2_lhhalf
#define w1 mccbetweenchop1and2_w1
#define h1 mccbetweenchop1and2_h1
#define w2 mccbetweenchop1and2_w2
#define h2 mccbetweenchop1and2_h2
#define l mccbetweenchop1and2_l
#define R0 mccbetweenchop1and2_R0
#define Qc mccbetweenchop1and2_Qc
#define alpha mccbetweenchop1and2_alpha
#define m mccbetweenchop1and2_m
#define nslit mccbetweenchop1and2_nslit
#define d mccbetweenchop1and2_d
#define Qcx mccbetweenchop1and2_Qcx
#define Qcy mccbetweenchop1and2_Qcy
#define alphax mccbetweenchop1and2_alphax
#define alphay mccbetweenchop1and2_alphay
#define W mccbetweenchop1and2_W
#define mx mccbetweenchop1and2_mx
#define my mccbetweenchop1and2_my
#define nu mccbetweenchop1and2_nu
#define phase mccbetweenchop1and2_phase
#line 80 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 8585 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'betweenchop2and3' [9]. */
#define mccompcurname  betweenchop2and3
#define mccompcurtype  Guide_channeled
#define mccompcurindex 9
#define w1c mccbetweenchop2and3_w1c
#define w2c mccbetweenchop2and3_w2c
#define ww mccbetweenchop2and3_ww
#define hh mccbetweenchop2and3_hh
#define whalf mccbetweenchop2and3_whalf
#define hhalf mccbetweenchop2and3_hhalf
#define lwhalf mccbetweenchop2and3_lwhalf
#define lhhalf mccbetweenchop2and3_lhhalf
#define w1 mccbetweenchop2and3_w1
#define h1 mccbetweenchop2and3_h1
#define w2 mccbetweenchop2and3_w2
#define h2 mccbetweenchop2and3_h2
#define l mccbetweenchop2and3_l
#define R0 mccbetweenchop2and3_R0
#define Qc mccbetweenchop2and3_Qc
#define alpha mccbetweenchop2and3_alpha
#define m mccbetweenchop2and3_m
#define nslit mccbetweenchop2and3_nslit
#define d mccbetweenchop2and3_d
#define Qcx mccbetweenchop2and3_Qcx
#define Qcy mccbetweenchop2and3_Qcy
#define alphax mccbetweenchop2and3_alphax
#define alphay mccbetweenchop2and3_alphay
#define W mccbetweenchop2and3_W
#define mx mccbetweenchop2and3_mx
#define my mccbetweenchop2and3_my
#define nu mccbetweenchop2and3_nu
#define phase mccbetweenchop2and3_phase
#line 80 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 8656 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'PR' [10]. */
#define mccompcurname  PR
#define mccompcurtype  DiskChopper
#define mccompcurindex 10
#define Tg mccPR_Tg
#define To mccPR_To
#define delta_y mccPR_delta_y
#define height mccPR_height
#define omega mccPR_omega
#define theta_0 mccPR_theta_0
#define radius mccPR_radius
#define yheight mccPR_yheight
#define nu mccPR_nu
#define nslit mccPR_nslit
#define jitter mccPR_jitter
#define delay mccPR_delay
#define isfirst mccPR_isfirst
#define n_pulse mccPR_n_pulse
#define abs_out mccPR_abs_out
#define phase mccPR_phase
#define xwidth mccPR_xwidth
#define verbose mccPR_verbose
#line 63 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 8713 "let_new.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'betweenchop3and4' [11]. */
#define mccompcurname  betweenchop3and4
#define mccompcurtype  Guide_channeled
#define mccompcurindex 11
#define w1c mccbetweenchop3and4_w1c
#define w2c mccbetweenchop3and4_w2c
#define ww mccbetweenchop3and4_ww
#define hh mccbetweenchop3and4_hh
#define whalf mccbetweenchop3and4_whalf
#define hhalf mccbetweenchop3and4_hhalf
#define lwhalf mccbetweenchop3and4_lwhalf
#define lhhalf mccbetweenchop3and4_lhhalf
#define w1 mccbetweenchop3and4_w1
#define h1 mccbetweenchop3and4_h1
#define w2 mccbetweenchop3and4_w2
#define h2 mccbetweenchop3and4_h2
#define l mccbetweenchop3and4_l
#define R0 mccbetweenchop3and4_R0
#define Qc mccbetweenchop3and4_Qc
#define alpha mccbetweenchop3and4_alpha
#define m mccbetweenchop3and4_m
#define nslit mccbetweenchop3and4_nslit
#define d mccbetweenchop3and4_d
#define Qcx mccbetweenchop3and4_Qcx
#define Qcy mccbetweenchop3and4_Qcy
#define alphax mccbetweenchop3and4_alphax
#define alphay mccbetweenchop3and4_alphay
#define W mccbetweenchop3and4_W
#define mx mccbetweenchop3and4_mx
#define my mccbetweenchop3and4_my
#define nu mccbetweenchop3and4_nu
#define phase mccbetweenchop3and4_phase
#line 80 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 8774 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'CR' [12]. */
#define mccompcurname  CR
#define mccompcurtype  DiskChopper
#define mccompcurindex 12
#define Tg mccCR_Tg
#define To mccCR_To
#define delta_y mccCR_delta_y
#define height mccCR_height
#define omega mccCR_omega
#define theta_0 mccCR_theta_0
#define radius mccCR_radius
#define yheight mccCR_yheight
#define nu mccCR_nu
#define nslit mccCR_nslit
#define jitter mccCR_jitter
#define delay mccCR_delay
#define isfirst mccCR_isfirst
#define n_pulse mccCR_n_pulse
#define abs_out mccCR_abs_out
#define phase mccCR_phase
#define xwidth mccCR_xwidth
#define verbose mccCR_verbose
#line 63 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 8831 "let_new.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'betweenchop4andatt' [13]. */
#define mccompcurname  betweenchop4andatt
#define mccompcurtype  Guide_channeled
#define mccompcurindex 13
#define w1c mccbetweenchop4andatt_w1c
#define w2c mccbetweenchop4andatt_w2c
#define ww mccbetweenchop4andatt_ww
#define hh mccbetweenchop4andatt_hh
#define whalf mccbetweenchop4andatt_whalf
#define hhalf mccbetweenchop4andatt_hhalf
#define lwhalf mccbetweenchop4andatt_lwhalf
#define lhhalf mccbetweenchop4andatt_lhhalf
#define w1 mccbetweenchop4andatt_w1
#define h1 mccbetweenchop4andatt_h1
#define w2 mccbetweenchop4andatt_w2
#define h2 mccbetweenchop4andatt_h2
#define l mccbetweenchop4andatt_l
#define R0 mccbetweenchop4andatt_R0
#define Qc mccbetweenchop4andatt_Qc
#define alpha mccbetweenchop4andatt_alpha
#define m mccbetweenchop4andatt_m
#define nslit mccbetweenchop4andatt_nslit
#define d mccbetweenchop4andatt_d
#define Qcx mccbetweenchop4andatt_Qcx
#define Qcy mccbetweenchop4andatt_Qcy
#define alphax mccbetweenchop4andatt_alphax
#define alphay mccbetweenchop4andatt_alphay
#define W mccbetweenchop4andatt_W
#define mx mccbetweenchop4andatt_mx
#define my mccbetweenchop4andatt_my
#define nu mccbetweenchop4andatt_nu
#define phase mccbetweenchop4andatt_phase
#line 80 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 8892 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'guide_move' [14]. */
#define mccompcurname  guide_move
#define mccompcurtype  Guide_channeled
#define mccompcurindex 14
#define w1c mccguide_move_w1c
#define w2c mccguide_move_w2c
#define ww mccguide_move_ww
#define hh mccguide_move_hh
#define whalf mccguide_move_whalf
#define hhalf mccguide_move_hhalf
#define lwhalf mccguide_move_lwhalf
#define lhhalf mccguide_move_lhhalf
#define w1 mccguide_move_w1
#define h1 mccguide_move_h1
#define w2 mccguide_move_w2
#define h2 mccguide_move_h2
#define l mccguide_move_l
#define R0 mccguide_move_R0
#define Qc mccguide_move_Qc
#define alpha mccguide_move_alpha
#define m mccguide_move_m
#define nslit mccguide_move_nslit
#define d mccguide_move_d
#define Qcx mccguide_move_Qcx
#define Qcy mccguide_move_Qcy
#define alphax mccguide_move_alphax
#define alphay mccguide_move_alphay
#define W mccguide_move_W
#define mx mccguide_move_mx
#define my mccguide_move_my
#define nu mccguide_move_nu
#define phase mccguide_move_phase
#line 80 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 8963 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'funnel' [15]. */
#define mccompcurname  funnel
#define mccompcurtype  Guide_channeled
#define mccompcurindex 15
#define w1c mccfunnel_w1c
#define w2c mccfunnel_w2c
#define ww mccfunnel_ww
#define hh mccfunnel_hh
#define whalf mccfunnel_whalf
#define hhalf mccfunnel_hhalf
#define lwhalf mccfunnel_lwhalf
#define lhhalf mccfunnel_lhhalf
#define w1 mccfunnel_w1
#define h1 mccfunnel_h1
#define w2 mccfunnel_w2
#define h2 mccfunnel_h2
#define l mccfunnel_l
#define R0 mccfunnel_R0
#define Qc mccfunnel_Qc
#define alpha mccfunnel_alpha
#define m mccfunnel_m
#define nslit mccfunnel_nslit
#define d mccfunnel_d
#define Qcx mccfunnel_Qcx
#define Qcy mccfunnel_Qcy
#define alphax mccfunnel_alphax
#define alphay mccfunnel_alphay
#define W mccfunnel_W
#define mx mccfunnel_mx
#define my mccfunnel_my
#define nu mccfunnel_nu
#define phase mccfunnel_phase
#line 80 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 9034 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'C5' [16]. */
#define mccompcurname  C5
#define mccompcurtype  DiskChopper
#define mccompcurindex 16
#define Tg mccC5_Tg
#define To mccC5_To
#define delta_y mccC5_delta_y
#define height mccC5_height
#define omega mccC5_omega
#define theta_0 mccC5_theta_0
#define radius mccC5_radius
#define yheight mccC5_yheight
#define nu mccC5_nu
#define nslit mccC5_nslit
#define jitter mccC5_jitter
#define delay mccC5_delay
#define isfirst mccC5_isfirst
#define n_pulse mccC5_n_pulse
#define abs_out mccC5_abs_out
#define phase mccC5_phase
#define xwidth mccC5_xwidth
#define verbose mccC5_verbose
#line 63 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 9091 "let_new.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'C5_counter' [17]. */
#define mccompcurname  C5_counter
#define mccompcurtype  DiskChopper
#define mccompcurindex 17
#define Tg mccC5_counter_Tg
#define To mccC5_counter_To
#define delta_y mccC5_counter_delta_y
#define height mccC5_counter_height
#define omega mccC5_counter_omega
#define theta_0 mccC5_counter_theta_0
#define radius mccC5_counter_radius
#define yheight mccC5_counter_yheight
#define nu mccC5_counter_nu
#define nslit mccC5_counter_nslit
#define jitter mccC5_counter_jitter
#define delay mccC5_counter_delay
#define isfirst mccC5_counter_isfirst
#define n_pulse mccC5_counter_n_pulse
#define abs_out mccC5_counter_abs_out
#define phase mccC5_counter_phase
#define xwidth mccC5_counter_xwidth
#define verbose mccC5_counter_verbose
#line 63 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
double Tg,To,delta_y,height,omega;
#line 9138 "let_new.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'endguide' [18]. */
#define mccompcurname  endguide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 18
#define w1c mccendguide_w1c
#define w2c mccendguide_w2c
#define ww mccendguide_ww
#define hh mccendguide_hh
#define whalf mccendguide_whalf
#define hhalf mccendguide_hhalf
#define lwhalf mccendguide_lwhalf
#define lhhalf mccendguide_lhhalf
#define w1 mccendguide_w1
#define h1 mccendguide_h1
#define w2 mccendguide_w2
#define h2 mccendguide_h2
#define l mccendguide_l
#define R0 mccendguide_R0
#define Qc mccendguide_Qc
#define alpha mccendguide_alpha
#define m mccendguide_m
#define nslit mccendguide_nslit
#define d mccendguide_d
#define Qcx mccendguide_Qcx
#define Qcy mccendguide_Qcy
#define alphax mccendguide_alphax
#define alphay mccendguide_alphay
#define W mccendguide_W
#define mx mccendguide_mx
#define my mccendguide_my
#define nu mccendguide_nu
#define phase mccendguide_phase
#line 80 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 9199 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'snot' [19]. */
#define mccompcurname  snot
#define mccompcurtype  Guide_channeled
#define mccompcurindex 19
#define w1c mccsnot_w1c
#define w2c mccsnot_w2c
#define ww mccsnot_ww
#define hh mccsnot_hh
#define whalf mccsnot_whalf
#define hhalf mccsnot_hhalf
#define lwhalf mccsnot_lwhalf
#define lhhalf mccsnot_lhhalf
#define w1 mccsnot_w1
#define h1 mccsnot_h1
#define w2 mccsnot_w2
#define h2 mccsnot_h2
#define l mccsnot_l
#define R0 mccsnot_R0
#define Qc mccsnot_Qc
#define alpha mccsnot_alpha
#define m mccsnot_m
#define nslit mccsnot_nslit
#define d mccsnot_d
#define Qcx mccsnot_Qcx
#define Qcy mccsnot_Qcy
#define alphax mccsnot_alphax
#define alphay mccsnot_alphay
#define W mccsnot_W
#define mx mccsnot_mx
#define my mccsnot_my
#define nu mccsnot_nu
#define phase mccsnot_phase
#line 80 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
double w1c;
double w2c;
double ww, hh;
double whalf, hhalf;
double lwhalf, lhhalf;
#line 9270 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'detslit' [20]. */
#define mccompcurname  detslit
#define mccompcurtype  Slit
#define mccompcurindex 20
#define xmin mccdetslit_xmin
#define xmax mccdetslit_xmax
#define ymin mccdetslit_ymin
#define ymax mccdetslit_ymax
#define radius mccdetslit_radius
#define xwidth mccdetslit_xwidth
#define yheight mccdetslit_yheight
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'mcplout' [21]. */
#define mccompcurname  mcplout
#define mccompcurtype  MCPL_output_horace
#define mccompcurindex 21
#define polarisationuse mccmcplout_polarisationuse
#define doubleprec mccmcplout_doubleprec
#define verbose mccmcplout_verbose
#define radius mccmcplout_radius
#define thickness mccmcplout_thickness
#define yheight mccmcplout_yheight
#define restore_neutron mccmcplout_restore_neutron
#define userflag mccmcplout_userflag
#define plate mccmcplout_plate
#define outputfile mccmcplout_outputfile
#define particle mccmcplout_particle
#define Particle mccmcplout_Particle
#define userflagenabled mccmcplout_userflagenabled
#define filename mccmcplout_filename
#define userflagcomment mccmcplout_userflagcomment
#define merge_mpi mccmcplout_merge_mpi
#define keep_mpi_unmerged mccmcplout_keep_mpi_unmerged
#line 76 "MCPL_output_horace.comp"
    mcpl_outfile_t outputfile;
    mcpl_particle_t *particle,Particle;
    int userflagenabled;
#line 9350 "let_new.c"
#undef keep_mpi_unmerged
#undef merge_mpi
#undef userflagcomment
#undef filename
#undef userflagenabled
#undef Particle
#undef particle
#undef outputfile
#undef plate
#undef userflag
#undef restore_neutron
#undef yheight
#undef thickness
#undef radius
#undef verbose
#undef doubleprec
#undef polarisationuse
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

Coords mcposaOrigin, mcposrOrigin;
Rotation mcrotaOrigin, mcrotrOrigin;
Coords mcposaSource, mcposrSource;
Rotation mcrotaSource, mcrotrSource;
Coords mcposashutter, mcposrshutter;
Rotation mcrotashutter, mcrotrshutter;
Coords mcposainsert, mcposrinsert;
Rotation mcrotainsert, mcrotrinsert;
Coords mcposaafterinsert, mcposrafterinsert;
Rotation mcrotaafterinsert, mcrotrafterinsert;
Coords mcposaRes1, mcposrRes1;
Rotation mcrotaRes1, mcrotrRes1;
Coords mcposaRes1_counter, mcposrRes1_counter;
Rotation mcrotaRes1_counter, mcrotrRes1_counter;
Coords mcposabetweenchop1and2, mcposrbetweenchop1and2;
Rotation mcrotabetweenchop1and2, mcrotrbetweenchop1and2;
Coords mcposabetweenchop2and3, mcposrbetweenchop2and3;
Rotation mcrotabetweenchop2and3, mcrotrbetweenchop2and3;
Coords mcposaPR, mcposrPR;
Rotation mcrotaPR, mcrotrPR;
Coords mcposabetweenchop3and4, mcposrbetweenchop3and4;
Rotation mcrotabetweenchop3and4, mcrotrbetweenchop3and4;
Coords mcposaCR, mcposrCR;
Rotation mcrotaCR, mcrotrCR;
Coords mcposabetweenchop4andatt, mcposrbetweenchop4andatt;
Rotation mcrotabetweenchop4andatt, mcrotrbetweenchop4andatt;
Coords mcposaguide_move, mcposrguide_move;
Rotation mcrotaguide_move, mcrotrguide_move;
Coords mcposafunnel, mcposrfunnel;
Rotation mcrotafunnel, mcrotrfunnel;
Coords mcposaC5, mcposrC5;
Rotation mcrotaC5, mcrotrC5;
Coords mcposaC5_counter, mcposrC5_counter;
Rotation mcrotaC5_counter, mcrotrC5_counter;
Coords mcposaendguide, mcposrendguide;
Rotation mcrotaendguide, mcrotrendguide;
Coords mcposasnot, mcposrsnot;
Rotation mcrotasnot, mcrotrsnot;
Coords mcposadetslit, mcposrdetslit;
Rotation mcrotadetslit, mcrotrdetslit;
Coords mcposamcplout, mcposrmcplout;
Rotation mcrotamcplout, mcrotrmcplout;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  LET
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposaLET coords_set(0,0,0)
#define Ei mcipEi
#define freq mcipfreq
#define chopper mcipchopper
#define output_filename mcipoutput_filename
#define sample mcipsample
#line 15 "let_new.instr"
{
E_foc = Ei;
v_foc = SE2V*sqrt(E_foc);
fprintf(stderr, "\nv_foc: %f\n\n",v_foc);
phase_time = 23.5 / v_foc + pha_offset;
E_min = 2.28e-3 * 23.5 / (phase_time+1.e-3); E_min *= E_min;
E_max = 2.28e-3 * 23.5 / (phase_time-1.e-3); E_max *= E_max;
PR = freq / 2.;
Res5 = freq;
plate = strcmp(sample, "plate")==0 ? 1 : 0;
for (int ii=0; chopper[ii]; ii++)
  chopper[ii] = tolower(chopper[ii]);
if(strstr(chopper, "res"))
  slot_size = 0.015;
else if(strstr(chopper, "flux"))
  slot_size = 0.031;
else
  slot_size = 0.020;
fprintf(stderr, "Using mode %s with slot size = %f mm\n", chopper, slot_size);
}
#line 9450 "let_new.c"
#undef sample
#undef output_filename
#undef chopper
#undef freq
#undef Ei
#undef mcposaLET
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname
  /* Computation of coordinate transformations. */
  {
    Coords mctc1, mctc2, mcLastComp;
    Rotation mctr1;
    double mcAccumulatedILength = 0;
    /* Initialize "last" component origin as (0,0,0) */
    mcLastComp = coords_set(0,0,0);

    mcDEBUG_INSTR()
  /* Component initializations. */
    /* Component Origin. */
  /* Setting parameters for component Origin. */
  SIG_MESSAGE("Origin (Init:SetPar)");

  SIG_MESSAGE("Origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaOrigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9479 "let_new.c"
  rot_copy(mcrotrOrigin, mcrotaOrigin);
  mcposaOrigin = coords_set(
#line 39 "let_new.instr"
    0,
#line 39 "let_new.instr"
    0,
#line 39 "let_new.instr"
    0);
#line 9488 "let_new.c"
  mctc1 = coords_neg(mcposaOrigin);
  mcposrOrigin = rot_apply(mcrotaOrigin, mctc1);
  mcDEBUG_COMPONENT("Origin", mcposaOrigin, mcrotaOrigin)
  mccomp_posa[1] = mcposaOrigin;
  mccomp_posr[1] = mcposrOrigin;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component Source. */
  /* Setting parameters for component Source. */
  SIG_MESSAGE("Source (Init:SetPar)");
#line 42 "let_new.instr"
  if("Let_Base.mcstas") strncpy(mccSource_Face, "Let_Base.mcstas" ? "Let_Base.mcstas" : "", 16384); else mccSource_Face[0]='\0';
#line 42 "let_new.instr"
  mccSource_E0 = E_min;
#line 42 "let_new.instr"
  mccSource_E1 = E_max;
#line 43 "let_new.instr"
  mccSource_tally = 8.4;
#line 43 "let_new.instr"
  mccSource_modPosition = 0;
#line 43 "let_new.instr"
  mccSource_modXsize = 0.12;
#line 43 "let_new.instr"
  mccSource_modZsize = 0.12;
#line 42 "let_new.instr"
  mccSource_xw = 0.04;
#line 42 "let_new.instr"
  mccSource_yh = 0.09;
#line 42 "let_new.instr"
  mccSource_dist = 1.68;
#line 9519 "let_new.c"

  SIG_MESSAGE("Source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9526 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaSource);
  rot_transpose(mcrotaOrigin, mctr1);
  rot_mul(mcrotaSource, mctr1, mcrotrSource);
  mctc1 = coords_set(
#line 44 "let_new.instr"
    0,
#line 44 "let_new.instr"
    0,
#line 44 "let_new.instr"
    0);
#line 9537 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaSource = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaOrigin, mcposaSource);
  mcposrSource = rot_apply(mcrotaSource, mctc1);
  mcDEBUG_COMPONENT("Source", mcposaSource, mcrotaSource)
  mccomp_posa[2] = mcposaSource;
  mccomp_posr[2] = mcposrSource;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component shutter. */
  /* Setting parameters for component shutter. */
  SIG_MESSAGE("shutter (Init:SetPar)");
#line 47 "let_new.instr"
  mccshutter_w1 = 0.04;
#line 47 "let_new.instr"
  mccshutter_h1 = 0.09;
#line 47 "let_new.instr"
  mccshutter_w2 = 0.04;
#line 47 "let_new.instr"
  mccshutter_h2 = 0.09;
#line 47 "let_new.instr"
  mccshutter_l = 1.980;
#line 71 "let_new.instr"
  mccshutter_R0 = 0.995;
#line 71 "let_new.instr"
  mccshutter_Qc = 0;
#line 71 "let_new.instr"
  mccshutter_alpha = 0;
#line 71 "let_new.instr"
  mccshutter_m = 0;
#line 71 "let_new.instr"
  mccshutter_nslit = 1;
#line 71 "let_new.instr"
  mccshutter_d = 0.0005;
#line 72 "let_new.instr"
  mccshutter_Qcx = 0.0218;
#line 72 "let_new.instr"
  mccshutter_Qcy = 0.0218;
#line 47 "let_new.instr"
  mccshutter_alphax = 4.38;
#line 47 "let_new.instr"
  mccshutter_alphay = 4.38;
#line 48 "let_new.instr"
  mccshutter_W = 3e-3;
#line 48 "let_new.instr"
  mccshutter_mx = 2;
#line 48 "let_new.instr"
  mccshutter_my = 3;
#line 72 "let_new.instr"
  mccshutter_nu = 0;
#line 72 "let_new.instr"
  mccshutter_phase = 0;
#line 9591 "let_new.c"

  SIG_MESSAGE("shutter (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9598 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotashutter);
  rot_transpose(mcrotaSource, mctr1);
  rot_mul(mcrotashutter, mctr1, mcrotrshutter);
  mctc1 = coords_set(
#line 49 "let_new.instr"
    0,
#line 49 "let_new.instr"
    0,
#line 49 "let_new.instr"
    1.680);
#line 9609 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposashutter = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaSource, mcposashutter);
  mcposrshutter = rot_apply(mcrotashutter, mctc1);
  mcDEBUG_COMPONENT("shutter", mcposashutter, mcrotashutter)
  mccomp_posa[3] = mcposashutter;
  mccomp_posr[3] = mcposrshutter;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component insert. */
  /* Setting parameters for component insert. */
  SIG_MESSAGE("insert (Init:SetPar)");
#line 52 "let_new.instr"
  mccinsert_w1 = 0.04;
#line 52 "let_new.instr"
  mccinsert_h1 = 0.09;
#line 52 "let_new.instr"
  mccinsert_w2 = 0.04;
#line 52 "let_new.instr"
  mccinsert_h2 = 0.09;
#line 52 "let_new.instr"
  mccinsert_l = 2.500;
#line 71 "let_new.instr"
  mccinsert_R0 = 0.995;
#line 71 "let_new.instr"
  mccinsert_Qc = 0;
#line 71 "let_new.instr"
  mccinsert_alpha = 0;
#line 71 "let_new.instr"
  mccinsert_m = 0;
#line 71 "let_new.instr"
  mccinsert_nslit = 1;
#line 71 "let_new.instr"
  mccinsert_d = 0.0005;
#line 72 "let_new.instr"
  mccinsert_Qcx = 0.0218;
#line 72 "let_new.instr"
  mccinsert_Qcy = 0.0218;
#line 52 "let_new.instr"
  mccinsert_alphax = 4.38;
#line 52 "let_new.instr"
  mccinsert_alphay = 4.38;
#line 53 "let_new.instr"
  mccinsert_W = 3e-3;
#line 53 "let_new.instr"
  mccinsert_mx = 2;
#line 53 "let_new.instr"
  mccinsert_my = 3;
#line 72 "let_new.instr"
  mccinsert_nu = 0;
#line 72 "let_new.instr"
  mccinsert_phase = 0;
#line 9663 "let_new.c"

  SIG_MESSAGE("insert (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9670 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotainsert);
  rot_transpose(mcrotashutter, mctr1);
  rot_mul(mcrotainsert, mctr1, mcrotrinsert);
  mctc1 = coords_set(
#line 54 "let_new.instr"
    0,
#line 54 "let_new.instr"
    0,
#line 54 "let_new.instr"
    3.740);
#line 9681 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposainsert = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposashutter, mcposainsert);
  mcposrinsert = rot_apply(mcrotainsert, mctc1);
  mcDEBUG_COMPONENT("insert", mcposainsert, mcrotainsert)
  mccomp_posa[4] = mcposainsert;
  mccomp_posr[4] = mcposrinsert;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component afterinsert. */
  /* Setting parameters for component afterinsert. */
  SIG_MESSAGE("afterinsert (Init:SetPar)");
#line 57 "let_new.instr"
  mccafterinsert_w1 = 0.04;
#line 57 "let_new.instr"
  mccafterinsert_h1 = 0.09;
#line 57 "let_new.instr"
  mccafterinsert_w2 = 0.04;
#line 57 "let_new.instr"
  mccafterinsert_h2 = 0.09;
#line 57 "let_new.instr"
  mccafterinsert_l = 1.514;
#line 71 "let_new.instr"
  mccafterinsert_R0 = 0.995;
#line 71 "let_new.instr"
  mccafterinsert_Qc = 0;
#line 71 "let_new.instr"
  mccafterinsert_alpha = 0;
#line 71 "let_new.instr"
  mccafterinsert_m = 0;
#line 71 "let_new.instr"
  mccafterinsert_nslit = 1;
#line 71 "let_new.instr"
  mccafterinsert_d = 0.0005;
#line 72 "let_new.instr"
  mccafterinsert_Qcx = 0.0218;
#line 72 "let_new.instr"
  mccafterinsert_Qcy = 0.0218;
#line 57 "let_new.instr"
  mccafterinsert_alphax = 4.38;
#line 57 "let_new.instr"
  mccafterinsert_alphay = 4.38;
#line 58 "let_new.instr"
  mccafterinsert_W = 3e-3;
#line 58 "let_new.instr"
  mccafterinsert_mx = 2;
#line 58 "let_new.instr"
  mccafterinsert_my = 3;
#line 72 "let_new.instr"
  mccafterinsert_nu = 0;
#line 72 "let_new.instr"
  mccafterinsert_phase = 0;
#line 9735 "let_new.c"

  SIG_MESSAGE("afterinsert (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9742 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaafterinsert);
  rot_transpose(mcrotainsert, mctr1);
  rot_mul(mcrotaafterinsert, mctr1, mcrotrafterinsert);
  mctc1 = coords_set(
#line 59 "let_new.instr"
    0,
#line 59 "let_new.instr"
    0,
#line 59 "let_new.instr"
    6.300);
#line 9753 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaafterinsert = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposainsert, mcposaafterinsert);
  mcposrafterinsert = rot_apply(mcrotaafterinsert, mctc1);
  mcDEBUG_COMPONENT("afterinsert", mcposaafterinsert, mcrotaafterinsert)
  mccomp_posa[5] = mcposaafterinsert;
  mccomp_posr[5] = mcposrafterinsert;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component Res1. */
  /* Setting parameters for component Res1. */
  SIG_MESSAGE("Res1 (Init:SetPar)");
#line 56 "let_new.instr"
  mccRes1_theta_0 = 0;
#line 62 "let_new.instr"
  mccRes1_radius = 0.279;
#line 62 "let_new.instr"
  mccRes1_yheight = 0.09;
#line 62 "let_new.instr"
  mccRes1_nu = PR / 2.0;
#line 62 "let_new.instr"
  mccRes1_nslit = 6;
#line 62 "let_new.instr"
  mccRes1_jitter = jitter;
#line 62 "let_new.instr"
  mccRes1_delay = 7.83 / v_foc + pha_offset;
#line 56 "let_new.instr"
  mccRes1_isfirst = 0;
#line 56 "let_new.instr"
  mccRes1_n_pulse = 1;
#line 56 "let_new.instr"
  mccRes1_abs_out = 1;
#line 56 "let_new.instr"
  mccRes1_phase = 0;
#line 62 "let_new.instr"
  mccRes1_xwidth = 0.04;
#line 56 "let_new.instr"
  mccRes1_verbose = 0;
#line 9793 "let_new.c"

  SIG_MESSAGE("Res1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9800 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaRes1);
  rot_transpose(mcrotaafterinsert, mctr1);
  rot_mul(mcrotaRes1, mctr1, mcrotrRes1);
  mctc1 = coords_set(
#line 63 "let_new.instr"
    0,
#line 63 "let_new.instr"
    0,
#line 63 "let_new.instr"
    7.83);
#line 9811 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaRes1 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaafterinsert, mcposaRes1);
  mcposrRes1 = rot_apply(mcrotaRes1, mctc1);
  mcDEBUG_COMPONENT("Res1", mcposaRes1, mcrotaRes1)
  mccomp_posa[6] = mcposaRes1;
  mccomp_posr[6] = mcposrRes1;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component Res1_counter. */
  /* Setting parameters for component Res1_counter. */
  SIG_MESSAGE("Res1_counter (Init:SetPar)");
#line 56 "let_new.instr"
  mccRes1_counter_theta_0 = 0;
#line 66 "let_new.instr"
  mccRes1_counter_radius = 0.279;
#line 66 "let_new.instr"
  mccRes1_counter_yheight = 0.09;
#line 66 "let_new.instr"
  mccRes1_counter_nu = - PR / 2.0;
#line 66 "let_new.instr"
  mccRes1_counter_nslit = 6;
#line 66 "let_new.instr"
  mccRes1_counter_jitter = jitter;
#line 66 "let_new.instr"
  mccRes1_counter_delay = 7.83 / v_foc + pha_offset;
#line 56 "let_new.instr"
  mccRes1_counter_isfirst = 0;
#line 56 "let_new.instr"
  mccRes1_counter_n_pulse = 1;
#line 56 "let_new.instr"
  mccRes1_counter_abs_out = 1;
#line 56 "let_new.instr"
  mccRes1_counter_phase = 0;
#line 66 "let_new.instr"
  mccRes1_counter_xwidth = 0.04;
#line 56 "let_new.instr"
  mccRes1_counter_verbose = 0;
#line 9851 "let_new.c"

  SIG_MESSAGE("Res1_counter (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9858 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaRes1_counter);
  rot_transpose(mcrotaRes1, mctr1);
  rot_mul(mcrotaRes1_counter, mctr1, mcrotrRes1_counter);
  mctc1 = coords_set(
#line 67 "let_new.instr"
    0,
#line 67 "let_new.instr"
    0,
#line 67 "let_new.instr"
    7.83 + 2e-6);
#line 9869 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaRes1_counter = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaRes1, mcposaRes1_counter);
  mcposrRes1_counter = rot_apply(mcrotaRes1_counter, mctc1);
  mcDEBUG_COMPONENT("Res1_counter", mcposaRes1_counter, mcrotaRes1_counter)
  mccomp_posa[7] = mcposaRes1_counter;
  mccomp_posr[7] = mcposrRes1_counter;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component betweenchop1and2. */
  /* Setting parameters for component betweenchop1and2. */
  SIG_MESSAGE("betweenchop1and2 (Init:SetPar)");
#line 70 "let_new.instr"
  mccbetweenchop1and2_w1 = 0.04;
#line 70 "let_new.instr"
  mccbetweenchop1and2_h1 = 0.09;
#line 70 "let_new.instr"
  mccbetweenchop1and2_w2 = 0.04;
#line 70 "let_new.instr"
  mccbetweenchop1and2_h2 = 0.09;
#line 70 "let_new.instr"
  mccbetweenchop1and2_l = 0.312;
#line 71 "let_new.instr"
  mccbetweenchop1and2_R0 = 0.995;
#line 71 "let_new.instr"
  mccbetweenchop1and2_Qc = 0;
#line 71 "let_new.instr"
  mccbetweenchop1and2_alpha = 0;
#line 71 "let_new.instr"
  mccbetweenchop1and2_m = 0;
#line 71 "let_new.instr"
  mccbetweenchop1and2_nslit = 1;
#line 71 "let_new.instr"
  mccbetweenchop1and2_d = 0.0005;
#line 72 "let_new.instr"
  mccbetweenchop1and2_Qcx = 0.0218;
#line 72 "let_new.instr"
  mccbetweenchop1and2_Qcy = 0.0218;
#line 70 "let_new.instr"
  mccbetweenchop1and2_alphax = 4.38;
#line 70 "let_new.instr"
  mccbetweenchop1and2_alphay = 4.38;
#line 71 "let_new.instr"
  mccbetweenchop1and2_W = 3e-3;
#line 71 "let_new.instr"
  mccbetweenchop1and2_mx = 2;
#line 71 "let_new.instr"
  mccbetweenchop1and2_my = 3;
#line 72 "let_new.instr"
  mccbetweenchop1and2_nu = 0;
#line 72 "let_new.instr"
  mccbetweenchop1and2_phase = 0;
#line 9923 "let_new.c"

  SIG_MESSAGE("betweenchop1and2 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 9930 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotabetweenchop1and2);
  rot_transpose(mcrotaRes1_counter, mctr1);
  rot_mul(mcrotabetweenchop1and2, mctr1, mcrotrbetweenchop1and2);
  mctc1 = coords_set(
#line 72 "let_new.instr"
    0,
#line 72 "let_new.instr"
    0,
#line 72 "let_new.instr"
    7.852);
#line 9941 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabetweenchop1and2 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaRes1_counter, mcposabetweenchop1and2);
  mcposrbetweenchop1and2 = rot_apply(mcrotabetweenchop1and2, mctc1);
  mcDEBUG_COMPONENT("betweenchop1and2", mcposabetweenchop1and2, mcrotabetweenchop1and2)
  mccomp_posa[8] = mcposabetweenchop1and2;
  mccomp_posr[8] = mcposrbetweenchop1and2;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component betweenchop2and3. */
  /* Setting parameters for component betweenchop2and3. */
  SIG_MESSAGE("betweenchop2and3 (Init:SetPar)");
#line 79 "let_new.instr"
  mccbetweenchop2and3_w1 = 0.04;
#line 79 "let_new.instr"
  mccbetweenchop2and3_h1 = 0.09;
#line 79 "let_new.instr"
  mccbetweenchop2and3_w2 = 0.04;
#line 79 "let_new.instr"
  mccbetweenchop2and3_h2 = 0.09;
#line 79 "let_new.instr"
  mccbetweenchop2and3_l = 3.499;
#line 71 "let_new.instr"
  mccbetweenchop2and3_R0 = 0.995;
#line 71 "let_new.instr"
  mccbetweenchop2and3_Qc = 0;
#line 71 "let_new.instr"
  mccbetweenchop2and3_alpha = 0;
#line 71 "let_new.instr"
  mccbetweenchop2and3_m = 0;
#line 71 "let_new.instr"
  mccbetweenchop2and3_nslit = 1;
#line 71 "let_new.instr"
  mccbetweenchop2and3_d = 0.0005;
#line 72 "let_new.instr"
  mccbetweenchop2and3_Qcx = 0.0218;
#line 72 "let_new.instr"
  mccbetweenchop2and3_Qcy = 0.0218;
#line 79 "let_new.instr"
  mccbetweenchop2and3_alphax = 4.38;
#line 79 "let_new.instr"
  mccbetweenchop2and3_alphay = 4.38;
#line 80 "let_new.instr"
  mccbetweenchop2and3_W = 3e-3;
#line 80 "let_new.instr"
  mccbetweenchop2and3_mx = 2;
#line 80 "let_new.instr"
  mccbetweenchop2and3_my = 3;
#line 72 "let_new.instr"
  mccbetweenchop2and3_nu = 0;
#line 72 "let_new.instr"
  mccbetweenchop2and3_phase = 0;
#line 9995 "let_new.c"

  SIG_MESSAGE("betweenchop2and3 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10002 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotabetweenchop2and3);
  rot_transpose(mcrotabetweenchop1and2, mctr1);
  rot_mul(mcrotabetweenchop2and3, mctr1, mcrotrbetweenchop2and3);
  mctc1 = coords_set(
#line 81 "let_new.instr"
    0,
#line 81 "let_new.instr"
    0,
#line 81 "let_new.instr"
    8.236);
#line 10013 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabetweenchop2and3 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposabetweenchop1and2, mcposabetweenchop2and3);
  mcposrbetweenchop2and3 = rot_apply(mcrotabetweenchop2and3, mctc1);
  mcDEBUG_COMPONENT("betweenchop2and3", mcposabetweenchop2and3, mcrotabetweenchop2and3)
  mccomp_posa[9] = mcposabetweenchop2and3;
  mccomp_posr[9] = mcposrbetweenchop2and3;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component PR. */
  /* Setting parameters for component PR. */
  SIG_MESSAGE("PR (Init:SetPar)");
#line 56 "let_new.instr"
  mccPR_theta_0 = 0;
#line 84 "let_new.instr"
  mccPR_radius = 0.29;
#line 84 "let_new.instr"
  mccPR_yheight = 0.1;
#line 84 "let_new.instr"
  mccPR_nu = PR;
#line 84 "let_new.instr"
  mccPR_nslit = 2;
#line 84 "let_new.instr"
  mccPR_jitter = jitter;
#line 84 "let_new.instr"
  mccPR_delay = 11.75 / v_foc + pha_offset;
#line 56 "let_new.instr"
  mccPR_isfirst = 0;
#line 56 "let_new.instr"
  mccPR_n_pulse = 1;
#line 56 "let_new.instr"
  mccPR_abs_out = 1;
#line 56 "let_new.instr"
  mccPR_phase = 0;
#line 84 "let_new.instr"
  mccPR_xwidth = 0.058;
#line 56 "let_new.instr"
  mccPR_verbose = 0;
#line 10053 "let_new.c"

  SIG_MESSAGE("PR (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10060 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaPR);
  rot_transpose(mcrotabetweenchop2and3, mctr1);
  rot_mul(mcrotaPR, mctr1, mcrotrPR);
  mctc1 = coords_set(
#line 85 "let_new.instr"
    0,
#line 85 "let_new.instr"
    0,
#line 85 "let_new.instr"
    11.75);
#line 10071 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaPR = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposabetweenchop2and3, mcposaPR);
  mcposrPR = rot_apply(mcrotaPR, mctc1);
  mcDEBUG_COMPONENT("PR", mcposaPR, mcrotaPR)
  mccomp_posa[10] = mcposaPR;
  mccomp_posr[10] = mcposrPR;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
    /* Component betweenchop3and4. */
  /* Setting parameters for component betweenchop3and4. */
  SIG_MESSAGE("betweenchop3and4 (Init:SetPar)");
#line 89 "let_new.instr"
  mccbetweenchop3and4_w1 = 0.04;
#line 89 "let_new.instr"
  mccbetweenchop3and4_h1 = 0.09;
#line 89 "let_new.instr"
  mccbetweenchop3and4_w2 = 0.04;
#line 89 "let_new.instr"
  mccbetweenchop3and4_h2 = 0.09;
#line 89 "let_new.instr"
  mccbetweenchop3and4_l = 3.886;
#line 71 "let_new.instr"
  mccbetweenchop3and4_R0 = 0.995;
#line 71 "let_new.instr"
  mccbetweenchop3and4_Qc = 0;
#line 71 "let_new.instr"
  mccbetweenchop3and4_alpha = 0;
#line 71 "let_new.instr"
  mccbetweenchop3and4_m = 0;
#line 71 "let_new.instr"
  mccbetweenchop3and4_nslit = 1;
#line 71 "let_new.instr"
  mccbetweenchop3and4_d = 0.0005;
#line 72 "let_new.instr"
  mccbetweenchop3and4_Qcx = 0.0218;
#line 72 "let_new.instr"
  mccbetweenchop3and4_Qcy = 0.0218;
#line 89 "let_new.instr"
  mccbetweenchop3and4_alphax = 4.38;
#line 89 "let_new.instr"
  mccbetweenchop3and4_alphay = 4.38;
#line 90 "let_new.instr"
  mccbetweenchop3and4_W = 3e-3;
#line 90 "let_new.instr"
  mccbetweenchop3and4_mx = 2;
#line 90 "let_new.instr"
  mccbetweenchop3and4_my = 3;
#line 72 "let_new.instr"
  mccbetweenchop3and4_nu = 0;
#line 72 "let_new.instr"
  mccbetweenchop3and4_phase = 0;
#line 10125 "let_new.c"

  SIG_MESSAGE("betweenchop3and4 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10132 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotabetweenchop3and4);
  rot_transpose(mcrotaPR, mctr1);
  rot_mul(mcrotabetweenchop3and4, mctr1, mcrotrbetweenchop3and4);
  mctc1 = coords_set(
#line 91 "let_new.instr"
    0,
#line 91 "let_new.instr"
    0,
#line 91 "let_new.instr"
    11.765);
#line 10143 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabetweenchop3and4 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaPR, mcposabetweenchop3and4);
  mcposrbetweenchop3and4 = rot_apply(mcrotabetweenchop3and4, mctc1);
  mcDEBUG_COMPONENT("betweenchop3and4", mcposabetweenchop3and4, mcrotabetweenchop3and4)
  mccomp_posa[11] = mcposabetweenchop3and4;
  mccomp_posr[11] = mcposrbetweenchop3and4;
  mcNCounter[11]  = mcPCounter[11] = mcP2Counter[11] = 0;
  mcAbsorbProp[11]= 0;
    /* Component CR. */
  /* Setting parameters for component CR. */
  SIG_MESSAGE("CR (Init:SetPar)");
#line 56 "let_new.instr"
  mccCR_theta_0 = 0;
#line 94 "let_new.instr"
  mccCR_radius = 0.29;
#line 94 "let_new.instr"
  mccCR_yheight = 0.1;
#line 94 "let_new.instr"
  mccCR_nu = Res5 / 2.0;
#line 94 "let_new.instr"
  mccCR_nslit = 6;
#line 94 "let_new.instr"
  mccCR_jitter = jitter;
#line 94 "let_new.instr"
  mccCR_delay = 15.66 / v_foc + pha_offset;
#line 56 "let_new.instr"
  mccCR_isfirst = 0;
#line 56 "let_new.instr"
  mccCR_n_pulse = 1;
#line 56 "let_new.instr"
  mccCR_abs_out = 1;
#line 56 "let_new.instr"
  mccCR_phase = 0;
#line 94 "let_new.instr"
  mccCR_xwidth = 0.054;
#line 56 "let_new.instr"
  mccCR_verbose = 0;
#line 10183 "let_new.c"

  SIG_MESSAGE("CR (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10190 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaCR);
  rot_transpose(mcrotabetweenchop3and4, mctr1);
  rot_mul(mcrotaCR, mctr1, mcrotrCR);
  mctc1 = coords_set(
#line 95 "let_new.instr"
    0,
#line 95 "let_new.instr"
    0,
#line 95 "let_new.instr"
    15.66);
#line 10201 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaCR = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposabetweenchop3and4, mcposaCR);
  mcposrCR = rot_apply(mcrotaCR, mctc1);
  mcDEBUG_COMPONENT("CR", mcposaCR, mcrotaCR)
  mccomp_posa[12] = mcposaCR;
  mccomp_posr[12] = mcposrCR;
  mcNCounter[12]  = mcPCounter[12] = mcP2Counter[12] = 0;
  mcAbsorbProp[12]= 0;
    /* Component betweenchop4andatt. */
  /* Setting parameters for component betweenchop4andatt. */
  SIG_MESSAGE("betweenchop4andatt (Init:SetPar)");
#line 99 "let_new.instr"
  mccbetweenchop4andatt_w1 = 0.04;
#line 99 "let_new.instr"
  mccbetweenchop4andatt_h1 = 0.09;
#line 99 "let_new.instr"
  mccbetweenchop4andatt_w2 = 0.04;
#line 99 "let_new.instr"
  mccbetweenchop4andatt_h2 = 0.0639;
#line 99 "let_new.instr"
  mccbetweenchop4andatt_l = 5.807;
#line 71 "let_new.instr"
  mccbetweenchop4andatt_R0 = 0.995;
#line 71 "let_new.instr"
  mccbetweenchop4andatt_Qc = 0;
#line 71 "let_new.instr"
  mccbetweenchop4andatt_alpha = 0;
#line 71 "let_new.instr"
  mccbetweenchop4andatt_m = 0;
#line 71 "let_new.instr"
  mccbetweenchop4andatt_nslit = 1;
#line 71 "let_new.instr"
  mccbetweenchop4andatt_d = 0.0005;
#line 72 "let_new.instr"
  mccbetweenchop4andatt_Qcx = 0.0218;
#line 72 "let_new.instr"
  mccbetweenchop4andatt_Qcy = 0.0218;
#line 99 "let_new.instr"
  mccbetweenchop4andatt_alphax = 4.38;
#line 99 "let_new.instr"
  mccbetweenchop4andatt_alphay = 4.38;
#line 100 "let_new.instr"
  mccbetweenchop4andatt_W = 3e-3;
#line 100 "let_new.instr"
  mccbetweenchop4andatt_mx = 2;
#line 100 "let_new.instr"
  mccbetweenchop4andatt_my = 3;
#line 72 "let_new.instr"
  mccbetweenchop4andatt_nu = 0;
#line 72 "let_new.instr"
  mccbetweenchop4andatt_phase = 0;
#line 10255 "let_new.c"

  SIG_MESSAGE("betweenchop4andatt (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10262 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotabetweenchop4andatt);
  rot_transpose(mcrotaCR, mctr1);
  rot_mul(mcrotabetweenchop4andatt, mctr1, mcrotrbetweenchop4andatt);
  mctc1 = coords_set(
#line 101 "let_new.instr"
    0,
#line 101 "let_new.instr"
    0,
#line 101 "let_new.instr"
    15.681);
#line 10273 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposabetweenchop4andatt = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaCR, mcposabetweenchop4andatt);
  mcposrbetweenchop4andatt = rot_apply(mcrotabetweenchop4andatt, mctc1);
  mcDEBUG_COMPONENT("betweenchop4andatt", mcposabetweenchop4andatt, mcrotabetweenchop4andatt)
  mccomp_posa[13] = mcposabetweenchop4andatt;
  mccomp_posr[13] = mcposrbetweenchop4andatt;
  mcNCounter[13]  = mcPCounter[13] = mcP2Counter[13] = 0;
  mcAbsorbProp[13]= 0;
    /* Component guide_move. */
  /* Setting parameters for component guide_move. */
  SIG_MESSAGE("guide_move (Init:SetPar)");
#line 104 "let_new.instr"
  mccguide_move_w1 = 0.04;
#line 104 "let_new.instr"
  mccguide_move_h1 = 0.0639;
#line 104 "let_new.instr"
  mccguide_move_w2 = 0.0311;
#line 104 "let_new.instr"
  mccguide_move_h2 = 0.05718;
#line 104 "let_new.instr"
  mccguide_move_l = 0.8823;
#line 71 "let_new.instr"
  mccguide_move_R0 = 0.995;
#line 71 "let_new.instr"
  mccguide_move_Qc = 0;
#line 71 "let_new.instr"
  mccguide_move_alpha = 0;
#line 71 "let_new.instr"
  mccguide_move_m = 0;
#line 71 "let_new.instr"
  mccguide_move_nslit = 1;
#line 71 "let_new.instr"
  mccguide_move_d = 0.0005;
#line 72 "let_new.instr"
  mccguide_move_Qcx = 0.0218;
#line 72 "let_new.instr"
  mccguide_move_Qcy = 0.0218;
#line 104 "let_new.instr"
  mccguide_move_alphax = 4.38;
#line 104 "let_new.instr"
  mccguide_move_alphay = 4.38;
#line 105 "let_new.instr"
  mccguide_move_W = 3e-3;
#line 105 "let_new.instr"
  mccguide_move_mx = 4;
#line 105 "let_new.instr"
  mccguide_move_my = 4;
#line 72 "let_new.instr"
  mccguide_move_nu = 0;
#line 72 "let_new.instr"
  mccguide_move_phase = 0;
#line 10327 "let_new.c"

  SIG_MESSAGE("guide_move (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10334 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaguide_move);
  rot_transpose(mcrotabetweenchop4andatt, mctr1);
  rot_mul(mcrotaguide_move, mctr1, mcrotrguide_move);
  mctc1 = coords_set(
#line 106 "let_new.instr"
    0,
#line 106 "let_new.instr"
    0,
#line 106 "let_new.instr"
    21.489);
#line 10345 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaguide_move = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposabetweenchop4andatt, mcposaguide_move);
  mcposrguide_move = rot_apply(mcrotaguide_move, mctc1);
  mcDEBUG_COMPONENT("guide_move", mcposaguide_move, mcrotaguide_move)
  mccomp_posa[14] = mcposaguide_move;
  mccomp_posr[14] = mcposrguide_move;
  mcNCounter[14]  = mcPCounter[14] = mcP2Counter[14] = 0;
  mcAbsorbProp[14]= 0;
    /* Component funnel. */
  /* Setting parameters for component funnel. */
  SIG_MESSAGE("funnel (Init:SetPar)");
#line 109 "let_new.instr"
  mccfunnel_w1 = 0.031;
#line 109 "let_new.instr"
  mccfunnel_h1 = 0.05711;
#line 109 "let_new.instr"
  mccfunnel_w2 = 0.02;
#line 109 "let_new.instr"
  mccfunnel_h2 = 0.04868;
#line 109 "let_new.instr"
  mccfunnel_l = 1.117;
#line 71 "let_new.instr"
  mccfunnel_R0 = 0.995;
#line 71 "let_new.instr"
  mccfunnel_Qc = 0;
#line 71 "let_new.instr"
  mccfunnel_alpha = 0;
#line 71 "let_new.instr"
  mccfunnel_m = 0;
#line 71 "let_new.instr"
  mccfunnel_nslit = 1;
#line 71 "let_new.instr"
  mccfunnel_d = 0.0005;
#line 72 "let_new.instr"
  mccfunnel_Qcx = 0.0218;
#line 72 "let_new.instr"
  mccfunnel_Qcy = 0.0218;
#line 109 "let_new.instr"
  mccfunnel_alphax = 4.38;
#line 109 "let_new.instr"
  mccfunnel_alphay = 4.38;
#line 110 "let_new.instr"
  mccfunnel_W = 3e-3;
#line 110 "let_new.instr"
  mccfunnel_mx = 4;
#line 110 "let_new.instr"
  mccfunnel_my = 4;
#line 72 "let_new.instr"
  mccfunnel_nu = 0;
#line 72 "let_new.instr"
  mccfunnel_phase = 0;
#line 10399 "let_new.c"

  SIG_MESSAGE("funnel (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10406 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotafunnel);
  rot_transpose(mcrotaguide_move, mctr1);
  rot_mul(mcrotafunnel, mctr1, mcrotrfunnel);
  mctc1 = coords_set(
#line 111 "let_new.instr"
    0,
#line 111 "let_new.instr"
    0,
#line 111 "let_new.instr"
    22.373);
#line 10417 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposafunnel = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaguide_move, mcposafunnel);
  mcposrfunnel = rot_apply(mcrotafunnel, mctc1);
  mcDEBUG_COMPONENT("funnel", mcposafunnel, mcrotafunnel)
  mccomp_posa[15] = mcposafunnel;
  mccomp_posr[15] = mcposrfunnel;
  mcNCounter[15]  = mcPCounter[15] = mcP2Counter[15] = 0;
  mcAbsorbProp[15]= 0;
    /* Component C5. */
  /* Setting parameters for component C5. */
  SIG_MESSAGE("C5 (Init:SetPar)");
#line 56 "let_new.instr"
  mccC5_theta_0 = 0;
#line 114 "let_new.instr"
  mccC5_radius = 0.279;
#line 114 "let_new.instr"
  mccC5_yheight = 0.0487;
#line 114 "let_new.instr"
  mccC5_nu = Res5;
#line 114 "let_new.instr"
  mccC5_nslit = 2;
#line 114 "let_new.instr"
  mccC5_jitter = jitter;
#line 114 "let_new.instr"
  mccC5_delay = 23.5 / v_foc + pha_offset;
#line 56 "let_new.instr"
  mccC5_isfirst = 0;
#line 56 "let_new.instr"
  mccC5_n_pulse = 1;
#line 56 "let_new.instr"
  mccC5_abs_out = 1;
#line 56 "let_new.instr"
  mccC5_phase = 0;
#line 114 "let_new.instr"
  mccC5_xwidth = slot_size;
#line 56 "let_new.instr"
  mccC5_verbose = 0;
#line 10457 "let_new.c"

  SIG_MESSAGE("C5 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10464 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaC5);
  rot_transpose(mcrotafunnel, mctr1);
  rot_mul(mcrotaC5, mctr1, mcrotrC5);
  mctc1 = coords_set(
#line 115 "let_new.instr"
    0,
#line 115 "let_new.instr"
    0,
#line 115 "let_new.instr"
    23.5);
#line 10475 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaC5 = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposafunnel, mcposaC5);
  mcposrC5 = rot_apply(mcrotaC5, mctc1);
  mcDEBUG_COMPONENT("C5", mcposaC5, mcrotaC5)
  mccomp_posa[16] = mcposaC5;
  mccomp_posr[16] = mcposrC5;
  mcNCounter[16]  = mcPCounter[16] = mcP2Counter[16] = 0;
  mcAbsorbProp[16]= 0;
    /* Component C5_counter. */
  /* Setting parameters for component C5_counter. */
  SIG_MESSAGE("C5_counter (Init:SetPar)");
#line 56 "let_new.instr"
  mccC5_counter_theta_0 = 0;
#line 118 "let_new.instr"
  mccC5_counter_radius = 0.279;
#line 118 "let_new.instr"
  mccC5_counter_yheight = 0.0487;
#line 118 "let_new.instr"
  mccC5_counter_nu = - Res5;
#line 118 "let_new.instr"
  mccC5_counter_nslit = 2;
#line 118 "let_new.instr"
  mccC5_counter_jitter = jitter;
#line 118 "let_new.instr"
  mccC5_counter_delay = 23.51 / v_foc + pha_offset;
#line 56 "let_new.instr"
  mccC5_counter_isfirst = 0;
#line 56 "let_new.instr"
  mccC5_counter_n_pulse = 1;
#line 56 "let_new.instr"
  mccC5_counter_abs_out = 1;
#line 56 "let_new.instr"
  mccC5_counter_phase = 0;
#line 118 "let_new.instr"
  mccC5_counter_xwidth = slot_size;
#line 56 "let_new.instr"
  mccC5_counter_verbose = 0;
#line 10515 "let_new.c"

  SIG_MESSAGE("C5_counter (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10522 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaC5_counter);
  rot_transpose(mcrotaC5, mctr1);
  rot_mul(mcrotaC5_counter, mctr1, mcrotrC5_counter);
  mctc1 = coords_set(
#line 119 "let_new.instr"
    0,
#line 119 "let_new.instr"
    0,
#line 119 "let_new.instr"
    23.51);
#line 10533 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaC5_counter = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaC5, mcposaC5_counter);
  mcposrC5_counter = rot_apply(mcrotaC5_counter, mctc1);
  mcDEBUG_COMPONENT("C5_counter", mcposaC5_counter, mcrotaC5_counter)
  mccomp_posa[17] = mcposaC5_counter;
  mccomp_posr[17] = mcposrC5_counter;
  mcNCounter[17]  = mcPCounter[17] = mcP2Counter[17] = 0;
  mcAbsorbProp[17]= 0;
    /* Component endguide. */
  /* Setting parameters for component endguide. */
  SIG_MESSAGE("endguide (Init:SetPar)");
#line 122 "let_new.instr"
  mccendguide_w1 = 0.02;
#line 122 "let_new.instr"
  mccendguide_h1 = 0.0484;
#line 122 "let_new.instr"
  mccendguide_w2 = 0.020;
#line 122 "let_new.instr"
  mccendguide_h2 = 0.04;
#line 122 "let_new.instr"
  mccendguide_l = 1.1;
#line 71 "let_new.instr"
  mccendguide_R0 = 0.995;
#line 71 "let_new.instr"
  mccendguide_Qc = 0;
#line 71 "let_new.instr"
  mccendguide_alpha = 0;
#line 71 "let_new.instr"
  mccendguide_m = 0;
#line 71 "let_new.instr"
  mccendguide_nslit = 1;
#line 71 "let_new.instr"
  mccendguide_d = 0.0005;
#line 72 "let_new.instr"
  mccendguide_Qcx = 0.0218;
#line 72 "let_new.instr"
  mccendguide_Qcy = 0.0218;
#line 122 "let_new.instr"
  mccendguide_alphax = 4.38;
#line 122 "let_new.instr"
  mccendguide_alphay = 4.38;
#line 123 "let_new.instr"
  mccendguide_W = 3e-3;
#line 123 "let_new.instr"
  mccendguide_mx = 4;
#line 123 "let_new.instr"
  mccendguide_my = 4;
#line 72 "let_new.instr"
  mccendguide_nu = 0;
#line 72 "let_new.instr"
  mccendguide_phase = 0;
#line 10587 "let_new.c"

  SIG_MESSAGE("endguide (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10594 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotaendguide);
  rot_transpose(mcrotaC5_counter, mctr1);
  rot_mul(mcrotaendguide, mctr1, mcrotrendguide);
  mctc1 = coords_set(
#line 124 "let_new.instr"
    0,
#line 124 "let_new.instr"
    0,
#line 124 "let_new.instr"
    23.52);
#line 10605 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaendguide = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaC5_counter, mcposaendguide);
  mcposrendguide = rot_apply(mcrotaendguide, mctc1);
  mcDEBUG_COMPONENT("endguide", mcposaendguide, mcrotaendguide)
  mccomp_posa[18] = mcposaendguide;
  mccomp_posr[18] = mcposrendguide;
  mcNCounter[18]  = mcPCounter[18] = mcP2Counter[18] = 0;
  mcAbsorbProp[18]= 0;
    /* Component snot. */
  /* Setting parameters for component snot. */
  SIG_MESSAGE("snot (Init:SetPar)");
#line 127 "let_new.instr"
  mccsnot_w1 = 0.02;
#line 127 "let_new.instr"
  mccsnot_h1 = 0.04;
#line 127 "let_new.instr"
  mccsnot_w2 = 0.02;
#line 127 "let_new.instr"
  mccsnot_h2 = 0.04;
#line 127 "let_new.instr"
  mccsnot_l = 0.23;
#line 71 "let_new.instr"
  mccsnot_R0 = 0.995;
#line 71 "let_new.instr"
  mccsnot_Qc = 0;
#line 71 "let_new.instr"
  mccsnot_alpha = 0;
#line 71 "let_new.instr"
  mccsnot_m = 0;
#line 71 "let_new.instr"
  mccsnot_nslit = 1;
#line 71 "let_new.instr"
  mccsnot_d = 0.0005;
#line 72 "let_new.instr"
  mccsnot_Qcx = 0.0218;
#line 72 "let_new.instr"
  mccsnot_Qcy = 0.0218;
#line 127 "let_new.instr"
  mccsnot_alphax = 4.38;
#line 127 "let_new.instr"
  mccsnot_alphay = 4.38;
#line 128 "let_new.instr"
  mccsnot_W = 3e-3;
#line 128 "let_new.instr"
  mccsnot_mx = 2;
#line 128 "let_new.instr"
  mccsnot_my = 4;
#line 72 "let_new.instr"
  mccsnot_nu = 0;
#line 72 "let_new.instr"
  mccsnot_phase = 0;
#line 10659 "let_new.c"

  SIG_MESSAGE("snot (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10666 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotasnot);
  rot_transpose(mcrotaendguide, mctr1);
  rot_mul(mcrotasnot, mctr1, mcrotrsnot);
  mctc1 = coords_set(
#line 129 "let_new.instr"
    0,
#line 129 "let_new.instr"
    0,
#line 129 "let_new.instr"
    24.622);
#line 10677 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasnot = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposaendguide, mcposasnot);
  mcposrsnot = rot_apply(mcrotasnot, mctc1);
  mcDEBUG_COMPONENT("snot", mcposasnot, mcrotasnot)
  mccomp_posa[19] = mcposasnot;
  mccomp_posr[19] = mcposrsnot;
  mcNCounter[19]  = mcPCounter[19] = mcP2Counter[19] = 0;
  mcAbsorbProp[19]= 0;
    /* Component detslit. */
  /* Setting parameters for component detslit. */
  SIG_MESSAGE("detslit (Init:SetPar)");
#line 131 "let_new.instr"
  mccdetslit_xmin = -0.025;
#line 131 "let_new.instr"
  mccdetslit_xmax = 0.025;
#line 131 "let_new.instr"
  mccdetslit_ymin = -0.025;
#line 131 "let_new.instr"
  mccdetslit_ymax = 0.025;
#line 43 "let_new.instr"
  mccdetslit_radius = 0;
#line 43 "let_new.instr"
  mccdetslit_xwidth = 0;
#line 43 "let_new.instr"
  mccdetslit_yheight = 0;
#line 10705 "let_new.c"

  SIG_MESSAGE("detslit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10712 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotadetslit);
  rot_transpose(mcrotasnot, mctr1);
  rot_mul(mcrotadetslit, mctr1, mcrotrdetslit);
  mctc1 = coords_set(
#line 132 "let_new.instr"
    0,
#line 132 "let_new.instr"
    0,
#line 132 "let_new.instr"
    25.0);
#line 10723 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposadetslit = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposasnot, mcposadetslit);
  mcposrdetslit = rot_apply(mcrotadetslit, mctc1);
  mcDEBUG_COMPONENT("detslit", mcposadetslit, mcrotadetslit)
  mccomp_posa[20] = mcposadetslit;
  mccomp_posr[20] = mcposrdetslit;
  mcNCounter[20]  = mcPCounter[20] = mcP2Counter[20] = 0;
  mcAbsorbProp[20]= 0;
    /* Component mcplout. */
  /* Setting parameters for component mcplout. */
  SIG_MESSAGE("mcplout (Init:SetPar)");
#line 134 "let_new.instr"
  if("mcstas.mcpl") strncpy(mccmcplout_filename, "mcstas.mcpl" ? "mcstas.mcpl" : "", 16384); else mccmcplout_filename[0]='\0';
#line 57 "let_new.instr"
  if("") strncpy(mccmcplout_userflagcomment, "" ? "" : "", 16384); else mccmcplout_userflagcomment[0]='\0';
#line 57 "let_new.instr"
  mccmcplout_merge_mpi = 1;
#line 57 "let_new.instr"
  mccmcplout_keep_mpi_unmerged = 0;
#line 10745 "let_new.c"

  SIG_MESSAGE("mcplout (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10752 "let_new.c"
  rot_mul(mctr1, mcrotaOrigin, mcrotamcplout);
  rot_transpose(mcrotadetslit, mctr1);
  rot_mul(mcrotamcplout, mctr1, mcrotrmcplout);
  mctc1 = coords_set(
#line 135 "let_new.instr"
    0,
#line 135 "let_new.instr"
    0,
#line 135 "let_new.instr"
    25.0);
#line 10763 "let_new.c"
  rot_transpose(mcrotaOrigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposamcplout = coords_add(mcposaOrigin, mctc2);
  mctc1 = coords_sub(mcposadetslit, mcposamcplout);
  mcposrmcplout = rot_apply(mcrotamcplout, mctc1);
  mcDEBUG_COMPONENT("mcplout", mcposamcplout, mcrotamcplout)
  mccomp_posa[21] = mcposamcplout;
  mccomp_posr[21] = mcposrmcplout;
  mcNCounter[21]  = mcPCounter[21] = mcP2Counter[21] = 0;
  mcAbsorbProp[21]= 0;
  /* Component initializations. */
  /* Initializations for component Origin. */
  SIG_MESSAGE("Origin (Init)");

  /* Initializations for component Source. */
  SIG_MESSAGE("Source (Init)");
#define mccompcurname  Source
#define mccompcurtype  ViewModISISver1
#define mccompcurindex 2
#define Face mccSource_Face
#define E0 mccSource_E0
#define E1 mccSource_E1
#define tally mccSource_tally
#define modPosition mccSource_modPosition
#define modXsize mccSource_modXsize
#define modZsize mccSource_modZsize
#define xw mccSource_xw
#define yh mccSource_yh
#define dist mccSource_dist
#line 1031 "ViewModISISver1.comp"
{
  /* READ IN THE ENERGY FILE */
  FILE* Tfile;
  

  Nsim=mcget_ncount();  // Number of points requested.

  Tfile=openFile(Face);	              // Get open file
  rtE0=convertEnergy(E0);
  rtE1=convertEnergy(E1);
  orderEnergy(&rtE0,&rtE1);

  readHtable(Tfile,rtE0,rtE1);
  fclose(Tfile);

  /**********************************************************************/

  Tnpts=0;
  Ncount=0;

  fprintf(stderr,"Face == %s \n",Face);
  fprintf(stderr,"Number of Energy Points == %d\n",TS.nEnergy);
  if (dist<0)
    {
      dist=TS.rdumMid;
      fprintf(stderr,"Setting distance to moderator surface == %g\n",
	      dist);
    }
  /* Do solid angle correction */
  angleArea= strArea();

  fprintf(stderr,"Totals:: %g %d %d \n",TS.Total,TS.nEnergy,TS.nTime);

}
#line 10828 "let_new.c"
#undef dist
#undef yh
#undef xw
#undef modZsize
#undef modXsize
#undef modPosition
#undef tally
#undef E1
#undef E0
#undef Face
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component shutter. */
  SIG_MESSAGE("shutter (Init)");
#define mccompcurname  shutter
#define mccompcurtype  Guide_channeled
#define mccompcurindex 3
#define w1c mccshutter_w1c
#define w2c mccshutter_w2c
#define ww mccshutter_ww
#define hh mccshutter_hh
#define whalf mccshutter_whalf
#define hhalf mccshutter_hhalf
#define lwhalf mccshutter_lwhalf
#define lhhalf mccshutter_lhhalf
#define w1 mccshutter_w1
#define h1 mccshutter_h1
#define w2 mccshutter_w2
#define h2 mccshutter_h2
#define l mccshutter_l
#define R0 mccshutter_R0
#define Qc mccshutter_Qc
#define alpha mccshutter_alpha
#define m mccshutter_m
#define nslit mccshutter_nslit
#define d mccshutter_d
#define Qcx mccshutter_Qcx
#define Qcy mccshutter_Qcy
#define alphax mccshutter_alphax
#define alphay mccshutter_alphay
#define W mccshutter_W
#define mx mccshutter_mx
#define my mccshutter_my
#define nu mccshutter_nu
#define phase mccshutter_phase
#line 88 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
if (!w2) w2=w1;
  if (!h2) h2=h1;
  if (nslit <= 0 || W <=0)
  { fprintf(stderr,"Guide_channeled: %s: nslit and W must be positive\n", NAME_CURRENT_COMP);
    exit(-1); }
  w1c = (w1 + d)/(double)nslit;
  w2c = (w2 + d)/(double)nslit;
  ww = .5*(w2c - w1c);
  hh = .5*(h2 - h1);
  whalf = .5*(w1c - d);
  hhalf = .5*h1;
  lwhalf = l*whalf;
  lhhalf = l*hhalf;

  if (m)     { mx=my=m; }
  if (Qc)    { Qcx=Qcy=Qc; }
  if (alpha) { alphax=alphay=alpha; }

  if ((nslit > 1) && (w1 != w2))
  {
    fprintf(stderr,"WARNING: Guide_channeled: %s:"
    "This component does not work with multichannel focusing guide\n"
    "Use Guide_gravity for that.\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_channeled: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (mcgravitation) fprintf(stderr,"WARNING: Guide_channeled: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  if (nu != 0 || phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_channeled: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_channeled: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, nu, phase);
    }
}
#line 10917 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component insert. */
  SIG_MESSAGE("insert (Init)");
#define mccompcurname  insert
#define mccompcurtype  Guide_channeled
#define mccompcurindex 4
#define w1c mccinsert_w1c
#define w2c mccinsert_w2c
#define ww mccinsert_ww
#define hh mccinsert_hh
#define whalf mccinsert_whalf
#define hhalf mccinsert_hhalf
#define lwhalf mccinsert_lwhalf
#define lhhalf mccinsert_lhhalf
#define w1 mccinsert_w1
#define h1 mccinsert_h1
#define w2 mccinsert_w2
#define h2 mccinsert_h2
#define l mccinsert_l
#define R0 mccinsert_R0
#define Qc mccinsert_Qc
#define alpha mccinsert_alpha
#define m mccinsert_m
#define nslit mccinsert_nslit
#define d mccinsert_d
#define Qcx mccinsert_Qcx
#define Qcy mccinsert_Qcy
#define alphax mccinsert_alphax
#define alphay mccinsert_alphay
#define W mccinsert_W
#define mx mccinsert_mx
#define my mccinsert_my
#define nu mccinsert_nu
#define phase mccinsert_phase
#line 88 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
if (!w2) w2=w1;
  if (!h2) h2=h1;
  if (nslit <= 0 || W <=0)
  { fprintf(stderr,"Guide_channeled: %s: nslit and W must be positive\n", NAME_CURRENT_COMP);
    exit(-1); }
  w1c = (w1 + d)/(double)nslit;
  w2c = (w2 + d)/(double)nslit;
  ww = .5*(w2c - w1c);
  hh = .5*(h2 - h1);
  whalf = .5*(w1c - d);
  hhalf = .5*h1;
  lwhalf = l*whalf;
  lhhalf = l*hhalf;

  if (m)     { mx=my=m; }
  if (Qc)    { Qcx=Qcy=Qc; }
  if (alpha) { alphax=alphay=alpha; }

  if ((nslit > 1) && (w1 != w2))
  {
    fprintf(stderr,"WARNING: Guide_channeled: %s:"
    "This component does not work with multichannel focusing guide\n"
    "Use Guide_gravity for that.\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_channeled: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (mcgravitation) fprintf(stderr,"WARNING: Guide_channeled: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  if (nu != 0 || phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_channeled: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_channeled: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, nu, phase);
    }
}
#line 11024 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component afterinsert. */
  SIG_MESSAGE("afterinsert (Init)");
#define mccompcurname  afterinsert
#define mccompcurtype  Guide_channeled
#define mccompcurindex 5
#define w1c mccafterinsert_w1c
#define w2c mccafterinsert_w2c
#define ww mccafterinsert_ww
#define hh mccafterinsert_hh
#define whalf mccafterinsert_whalf
#define hhalf mccafterinsert_hhalf
#define lwhalf mccafterinsert_lwhalf
#define lhhalf mccafterinsert_lhhalf
#define w1 mccafterinsert_w1
#define h1 mccafterinsert_h1
#define w2 mccafterinsert_w2
#define h2 mccafterinsert_h2
#define l mccafterinsert_l
#define R0 mccafterinsert_R0
#define Qc mccafterinsert_Qc
#define alpha mccafterinsert_alpha
#define m mccafterinsert_m
#define nslit mccafterinsert_nslit
#define d mccafterinsert_d
#define Qcx mccafterinsert_Qcx
#define Qcy mccafterinsert_Qcy
#define alphax mccafterinsert_alphax
#define alphay mccafterinsert_alphay
#define W mccafterinsert_W
#define mx mccafterinsert_mx
#define my mccafterinsert_my
#define nu mccafterinsert_nu
#define phase mccafterinsert_phase
#line 88 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
if (!w2) w2=w1;
  if (!h2) h2=h1;
  if (nslit <= 0 || W <=0)
  { fprintf(stderr,"Guide_channeled: %s: nslit and W must be positive\n", NAME_CURRENT_COMP);
    exit(-1); }
  w1c = (w1 + d)/(double)nslit;
  w2c = (w2 + d)/(double)nslit;
  ww = .5*(w2c - w1c);
  hh = .5*(h2 - h1);
  whalf = .5*(w1c - d);
  hhalf = .5*h1;
  lwhalf = l*whalf;
  lhhalf = l*hhalf;

  if (m)     { mx=my=m; }
  if (Qc)    { Qcx=Qcy=Qc; }
  if (alpha) { alphax=alphay=alpha; }

  if ((nslit > 1) && (w1 != w2))
  {
    fprintf(stderr,"WARNING: Guide_channeled: %s:"
    "This component does not work with multichannel focusing guide\n"
    "Use Guide_gravity for that.\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_channeled: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (mcgravitation) fprintf(stderr,"WARNING: Guide_channeled: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  if (nu != 0 || phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_channeled: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_channeled: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, nu, phase);
    }
}
#line 11131 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Res1. */
  SIG_MESSAGE("Res1 (Init)");
#define mccompcurname  Res1
#define mccompcurtype  DiskChopper
#define mccompcurindex 6
#define Tg mccRes1_Tg
#define To mccRes1_To
#define delta_y mccRes1_delta_y
#define height mccRes1_height
#define omega mccRes1_omega
#define theta_0 mccRes1_theta_0
#define radius mccRes1_radius
#define yheight mccRes1_yheight
#define nu mccRes1_nu
#define nslit mccRes1_nslit
#define jitter mccRes1_jitter
#define delay mccRes1_delay
#define isfirst mccRes1_isfirst
#define n_pulse mccRes1_n_pulse
#define abs_out mccRes1_abs_out
#define phase mccRes1_phase
#define xwidth mccRes1_xwidth
#define verbose mccRes1_verbose
#line 67 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{
/* If slit height 'unset', assume full opening */
if (yheight == 0) {
        height=radius;
      } else {
        height=yheight;
      }
      delta_y = radius-height/2; /* radius at beam center */
      omega=2.0*PI*nu; /* rad/s */
      if (xwidth && !theta_0 && radius) theta_0 = 2*RAD2DEG*asin(xwidth/2/delta_y);

      if (nslit<=0 || theta_0 <= 0 || radius <=0)
      { fprintf(stderr,"DiskChopper: %s: nslit, theta_0 and radius must be > 0\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (nslit*theta_0 >= 360)
      { fprintf(stderr,"DiskChopper: %s: nslit * theta_0 exceeds 2PI\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (yheight && yheight>radius) {
        fprintf(stderr,"DiskChopper: %s: yheight must be < radius\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (isfirst && n_pulse <=0)
      { fprintf(stderr,"DiskChopper: %s: wrong First chopper pulse number (n_pulse=%g)\n", NAME_CURRENT_COMP, n_pulse);
        exit(-1); }
      if (!omega) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper frequency is 0!\n", NAME_CURRENT_COMP);
        omega = 1e-15; /* We should actually use machine epsilon here... */
      }
      if (!abs_out) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper will NOT absorb neutrons outside radius %g [m]\n", NAME_CURRENT_COMP, radius);
      }

      theta_0*=DEG2RAD;


      /* Calulate delay from phase and vice versa */
      if (phase) {
        if (delay) {
          fprintf(stderr,"DiskChopper: %s WARNING: delay AND phase specified. Using phase setting\n", NAME_CURRENT_COMP);
        }
        phase*=DEG2RAD;
        /* 'Delay' should always be a delay, taking rotation direction into account: */
        delay=phase/fabs(omega);
      } else {
        phase=delay*omega;  /* rad */
      }

      /* Time from opening of slit to next opening of slit */
      Tg=2.0*PI/fabs(omega)/nslit;

      /* How long can neutrons pass the Chopper at a single point */
      To=theta_0/fabs(omega);

      if (!xwidth) xwidth=2*delta_y*sin(theta_0/2);

      if (verbose && nu) {
        printf("DiskChopper: %s: frequency=%g [Hz] %g [rpm], time frame=%g [s] phase=%g [deg]\n",
          NAME_CURRENT_COMP, nu, nu*60, Tg, phase*RAD2DEG);
        printf("             %g slits, angle=%g [deg] height=%g [m], width=%g [m] at radius=%g [m]\n",
          nslit, theta_0*RAD2DEG, height, xwidth, delta_y);
      }
}
#line 11249 "let_new.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component Res1_counter. */
  SIG_MESSAGE("Res1_counter (Init)");
#define mccompcurname  Res1_counter
#define mccompcurtype  DiskChopper
#define mccompcurindex 7
#define Tg mccRes1_counter_Tg
#define To mccRes1_counter_To
#define delta_y mccRes1_counter_delta_y
#define height mccRes1_counter_height
#define omega mccRes1_counter_omega
#define theta_0 mccRes1_counter_theta_0
#define radius mccRes1_counter_radius
#define yheight mccRes1_counter_yheight
#define nu mccRes1_counter_nu
#define nslit mccRes1_counter_nslit
#define jitter mccRes1_counter_jitter
#define delay mccRes1_counter_delay
#define isfirst mccRes1_counter_isfirst
#define n_pulse mccRes1_counter_n_pulse
#define abs_out mccRes1_counter_abs_out
#define phase mccRes1_counter_phase
#define xwidth mccRes1_counter_xwidth
#define verbose mccRes1_counter_verbose
#line 67 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{
/* If slit height 'unset', assume full opening */
if (yheight == 0) {
        height=radius;
      } else {
        height=yheight;
      }
      delta_y = radius-height/2; /* radius at beam center */
      omega=2.0*PI*nu; /* rad/s */
      if (xwidth && !theta_0 && radius) theta_0 = 2*RAD2DEG*asin(xwidth/2/delta_y);

      if (nslit<=0 || theta_0 <= 0 || radius <=0)
      { fprintf(stderr,"DiskChopper: %s: nslit, theta_0 and radius must be > 0\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (nslit*theta_0 >= 360)
      { fprintf(stderr,"DiskChopper: %s: nslit * theta_0 exceeds 2PI\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (yheight && yheight>radius) {
        fprintf(stderr,"DiskChopper: %s: yheight must be < radius\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (isfirst && n_pulse <=0)
      { fprintf(stderr,"DiskChopper: %s: wrong First chopper pulse number (n_pulse=%g)\n", NAME_CURRENT_COMP, n_pulse);
        exit(-1); }
      if (!omega) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper frequency is 0!\n", NAME_CURRENT_COMP);
        omega = 1e-15; /* We should actually use machine epsilon here... */
      }
      if (!abs_out) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper will NOT absorb neutrons outside radius %g [m]\n", NAME_CURRENT_COMP, radius);
      }

      theta_0*=DEG2RAD;


      /* Calulate delay from phase and vice versa */
      if (phase) {
        if (delay) {
          fprintf(stderr,"DiskChopper: %s WARNING: delay AND phase specified. Using phase setting\n", NAME_CURRENT_COMP);
        }
        phase*=DEG2RAD;
        /* 'Delay' should always be a delay, taking rotation direction into account: */
        delay=phase/fabs(omega);
      } else {
        phase=delay*omega;  /* rad */
      }

      /* Time from opening of slit to next opening of slit */
      Tg=2.0*PI/fabs(omega)/nslit;

      /* How long can neutrons pass the Chopper at a single point */
      To=theta_0/fabs(omega);

      if (!xwidth) xwidth=2*delta_y*sin(theta_0/2);

      if (verbose && nu) {
        printf("DiskChopper: %s: frequency=%g [Hz] %g [rpm], time frame=%g [s] phase=%g [deg]\n",
          NAME_CURRENT_COMP, nu, nu*60, Tg, phase*RAD2DEG);
        printf("             %g slits, angle=%g [deg] height=%g [m], width=%g [m] at radius=%g [m]\n",
          nslit, theta_0*RAD2DEG, height, xwidth, delta_y);
      }
}
#line 11357 "let_new.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component betweenchop1and2. */
  SIG_MESSAGE("betweenchop1and2 (Init)");
#define mccompcurname  betweenchop1and2
#define mccompcurtype  Guide_channeled
#define mccompcurindex 8
#define w1c mccbetweenchop1and2_w1c
#define w2c mccbetweenchop1and2_w2c
#define ww mccbetweenchop1and2_ww
#define hh mccbetweenchop1and2_hh
#define whalf mccbetweenchop1and2_whalf
#define hhalf mccbetweenchop1and2_hhalf
#define lwhalf mccbetweenchop1and2_lwhalf
#define lhhalf mccbetweenchop1and2_lhhalf
#define w1 mccbetweenchop1and2_w1
#define h1 mccbetweenchop1and2_h1
#define w2 mccbetweenchop1and2_w2
#define h2 mccbetweenchop1and2_h2
#define l mccbetweenchop1and2_l
#define R0 mccbetweenchop1and2_R0
#define Qc mccbetweenchop1and2_Qc
#define alpha mccbetweenchop1and2_alpha
#define m mccbetweenchop1and2_m
#define nslit mccbetweenchop1and2_nslit
#define d mccbetweenchop1and2_d
#define Qcx mccbetweenchop1and2_Qcx
#define Qcy mccbetweenchop1and2_Qcy
#define alphax mccbetweenchop1and2_alphax
#define alphay mccbetweenchop1and2_alphay
#define W mccbetweenchop1and2_W
#define mx mccbetweenchop1and2_mx
#define my mccbetweenchop1and2_my
#define nu mccbetweenchop1and2_nu
#define phase mccbetweenchop1and2_phase
#line 88 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
if (!w2) w2=w1;
  if (!h2) h2=h1;
  if (nslit <= 0 || W <=0)
  { fprintf(stderr,"Guide_channeled: %s: nslit and W must be positive\n", NAME_CURRENT_COMP);
    exit(-1); }
  w1c = (w1 + d)/(double)nslit;
  w2c = (w2 + d)/(double)nslit;
  ww = .5*(w2c - w1c);
  hh = .5*(h2 - h1);
  whalf = .5*(w1c - d);
  hhalf = .5*h1;
  lwhalf = l*whalf;
  lhhalf = l*hhalf;

  if (m)     { mx=my=m; }
  if (Qc)    { Qcx=Qcy=Qc; }
  if (alpha) { alphax=alphay=alpha; }

  if ((nslit > 1) && (w1 != w2))
  {
    fprintf(stderr,"WARNING: Guide_channeled: %s:"
    "This component does not work with multichannel focusing guide\n"
    "Use Guide_gravity for that.\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_channeled: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (mcgravitation) fprintf(stderr,"WARNING: Guide_channeled: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  if (nu != 0 || phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_channeled: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_channeled: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, nu, phase);
    }
}
#line 11454 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component betweenchop2and3. */
  SIG_MESSAGE("betweenchop2and3 (Init)");
#define mccompcurname  betweenchop2and3
#define mccompcurtype  Guide_channeled
#define mccompcurindex 9
#define w1c mccbetweenchop2and3_w1c
#define w2c mccbetweenchop2and3_w2c
#define ww mccbetweenchop2and3_ww
#define hh mccbetweenchop2and3_hh
#define whalf mccbetweenchop2and3_whalf
#define hhalf mccbetweenchop2and3_hhalf
#define lwhalf mccbetweenchop2and3_lwhalf
#define lhhalf mccbetweenchop2and3_lhhalf
#define w1 mccbetweenchop2and3_w1
#define h1 mccbetweenchop2and3_h1
#define w2 mccbetweenchop2and3_w2
#define h2 mccbetweenchop2and3_h2
#define l mccbetweenchop2and3_l
#define R0 mccbetweenchop2and3_R0
#define Qc mccbetweenchop2and3_Qc
#define alpha mccbetweenchop2and3_alpha
#define m mccbetweenchop2and3_m
#define nslit mccbetweenchop2and3_nslit
#define d mccbetweenchop2and3_d
#define Qcx mccbetweenchop2and3_Qcx
#define Qcy mccbetweenchop2and3_Qcy
#define alphax mccbetweenchop2and3_alphax
#define alphay mccbetweenchop2and3_alphay
#define W mccbetweenchop2and3_W
#define mx mccbetweenchop2and3_mx
#define my mccbetweenchop2and3_my
#define nu mccbetweenchop2and3_nu
#define phase mccbetweenchop2and3_phase
#line 88 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
if (!w2) w2=w1;
  if (!h2) h2=h1;
  if (nslit <= 0 || W <=0)
  { fprintf(stderr,"Guide_channeled: %s: nslit and W must be positive\n", NAME_CURRENT_COMP);
    exit(-1); }
  w1c = (w1 + d)/(double)nslit;
  w2c = (w2 + d)/(double)nslit;
  ww = .5*(w2c - w1c);
  hh = .5*(h2 - h1);
  whalf = .5*(w1c - d);
  hhalf = .5*h1;
  lwhalf = l*whalf;
  lhhalf = l*hhalf;

  if (m)     { mx=my=m; }
  if (Qc)    { Qcx=Qcy=Qc; }
  if (alpha) { alphax=alphay=alpha; }

  if ((nslit > 1) && (w1 != w2))
  {
    fprintf(stderr,"WARNING: Guide_channeled: %s:"
    "This component does not work with multichannel focusing guide\n"
    "Use Guide_gravity for that.\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_channeled: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (mcgravitation) fprintf(stderr,"WARNING: Guide_channeled: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  if (nu != 0 || phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_channeled: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_channeled: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, nu, phase);
    }
}
#line 11561 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component PR. */
  SIG_MESSAGE("PR (Init)");
#define mccompcurname  PR
#define mccompcurtype  DiskChopper
#define mccompcurindex 10
#define Tg mccPR_Tg
#define To mccPR_To
#define delta_y mccPR_delta_y
#define height mccPR_height
#define omega mccPR_omega
#define theta_0 mccPR_theta_0
#define radius mccPR_radius
#define yheight mccPR_yheight
#define nu mccPR_nu
#define nslit mccPR_nslit
#define jitter mccPR_jitter
#define delay mccPR_delay
#define isfirst mccPR_isfirst
#define n_pulse mccPR_n_pulse
#define abs_out mccPR_abs_out
#define phase mccPR_phase
#define xwidth mccPR_xwidth
#define verbose mccPR_verbose
#line 67 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{
/* If slit height 'unset', assume full opening */
if (yheight == 0) {
        height=radius;
      } else {
        height=yheight;
      }
      delta_y = radius-height/2; /* radius at beam center */
      omega=2.0*PI*nu; /* rad/s */
      if (xwidth && !theta_0 && radius) theta_0 = 2*RAD2DEG*asin(xwidth/2/delta_y);

      if (nslit<=0 || theta_0 <= 0 || radius <=0)
      { fprintf(stderr,"DiskChopper: %s: nslit, theta_0 and radius must be > 0\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (nslit*theta_0 >= 360)
      { fprintf(stderr,"DiskChopper: %s: nslit * theta_0 exceeds 2PI\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (yheight && yheight>radius) {
        fprintf(stderr,"DiskChopper: %s: yheight must be < radius\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (isfirst && n_pulse <=0)
      { fprintf(stderr,"DiskChopper: %s: wrong First chopper pulse number (n_pulse=%g)\n", NAME_CURRENT_COMP, n_pulse);
        exit(-1); }
      if (!omega) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper frequency is 0!\n", NAME_CURRENT_COMP);
        omega = 1e-15; /* We should actually use machine epsilon here... */
      }
      if (!abs_out) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper will NOT absorb neutrons outside radius %g [m]\n", NAME_CURRENT_COMP, radius);
      }

      theta_0*=DEG2RAD;


      /* Calulate delay from phase and vice versa */
      if (phase) {
        if (delay) {
          fprintf(stderr,"DiskChopper: %s WARNING: delay AND phase specified. Using phase setting\n", NAME_CURRENT_COMP);
        }
        phase*=DEG2RAD;
        /* 'Delay' should always be a delay, taking rotation direction into account: */
        delay=phase/fabs(omega);
      } else {
        phase=delay*omega;  /* rad */
      }

      /* Time from opening of slit to next opening of slit */
      Tg=2.0*PI/fabs(omega)/nslit;

      /* How long can neutrons pass the Chopper at a single point */
      To=theta_0/fabs(omega);

      if (!xwidth) xwidth=2*delta_y*sin(theta_0/2);

      if (verbose && nu) {
        printf("DiskChopper: %s: frequency=%g [Hz] %g [rpm], time frame=%g [s] phase=%g [deg]\n",
          NAME_CURRENT_COMP, nu, nu*60, Tg, phase*RAD2DEG);
        printf("             %g slits, angle=%g [deg] height=%g [m], width=%g [m] at radius=%g [m]\n",
          nslit, theta_0*RAD2DEG, height, xwidth, delta_y);
      }
}
#line 11679 "let_new.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component betweenchop3and4. */
  SIG_MESSAGE("betweenchop3and4 (Init)");
#define mccompcurname  betweenchop3and4
#define mccompcurtype  Guide_channeled
#define mccompcurindex 11
#define w1c mccbetweenchop3and4_w1c
#define w2c mccbetweenchop3and4_w2c
#define ww mccbetweenchop3and4_ww
#define hh mccbetweenchop3and4_hh
#define whalf mccbetweenchop3and4_whalf
#define hhalf mccbetweenchop3and4_hhalf
#define lwhalf mccbetweenchop3and4_lwhalf
#define lhhalf mccbetweenchop3and4_lhhalf
#define w1 mccbetweenchop3and4_w1
#define h1 mccbetweenchop3and4_h1
#define w2 mccbetweenchop3and4_w2
#define h2 mccbetweenchop3and4_h2
#define l mccbetweenchop3and4_l
#define R0 mccbetweenchop3and4_R0
#define Qc mccbetweenchop3and4_Qc
#define alpha mccbetweenchop3and4_alpha
#define m mccbetweenchop3and4_m
#define nslit mccbetweenchop3and4_nslit
#define d mccbetweenchop3and4_d
#define Qcx mccbetweenchop3and4_Qcx
#define Qcy mccbetweenchop3and4_Qcy
#define alphax mccbetweenchop3and4_alphax
#define alphay mccbetweenchop3and4_alphay
#define W mccbetweenchop3and4_W
#define mx mccbetweenchop3and4_mx
#define my mccbetweenchop3and4_my
#define nu mccbetweenchop3and4_nu
#define phase mccbetweenchop3and4_phase
#line 88 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
if (!w2) w2=w1;
  if (!h2) h2=h1;
  if (nslit <= 0 || W <=0)
  { fprintf(stderr,"Guide_channeled: %s: nslit and W must be positive\n", NAME_CURRENT_COMP);
    exit(-1); }
  w1c = (w1 + d)/(double)nslit;
  w2c = (w2 + d)/(double)nslit;
  ww = .5*(w2c - w1c);
  hh = .5*(h2 - h1);
  whalf = .5*(w1c - d);
  hhalf = .5*h1;
  lwhalf = l*whalf;
  lhhalf = l*hhalf;

  if (m)     { mx=my=m; }
  if (Qc)    { Qcx=Qcy=Qc; }
  if (alpha) { alphax=alphay=alpha; }

  if ((nslit > 1) && (w1 != w2))
  {
    fprintf(stderr,"WARNING: Guide_channeled: %s:"
    "This component does not work with multichannel focusing guide\n"
    "Use Guide_gravity for that.\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_channeled: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (mcgravitation) fprintf(stderr,"WARNING: Guide_channeled: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  if (nu != 0 || phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_channeled: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_channeled: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, nu, phase);
    }
}
#line 11776 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component CR. */
  SIG_MESSAGE("CR (Init)");
#define mccompcurname  CR
#define mccompcurtype  DiskChopper
#define mccompcurindex 12
#define Tg mccCR_Tg
#define To mccCR_To
#define delta_y mccCR_delta_y
#define height mccCR_height
#define omega mccCR_omega
#define theta_0 mccCR_theta_0
#define radius mccCR_radius
#define yheight mccCR_yheight
#define nu mccCR_nu
#define nslit mccCR_nslit
#define jitter mccCR_jitter
#define delay mccCR_delay
#define isfirst mccCR_isfirst
#define n_pulse mccCR_n_pulse
#define abs_out mccCR_abs_out
#define phase mccCR_phase
#define xwidth mccCR_xwidth
#define verbose mccCR_verbose
#line 67 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{
/* If slit height 'unset', assume full opening */
if (yheight == 0) {
        height=radius;
      } else {
        height=yheight;
      }
      delta_y = radius-height/2; /* radius at beam center */
      omega=2.0*PI*nu; /* rad/s */
      if (xwidth && !theta_0 && radius) theta_0 = 2*RAD2DEG*asin(xwidth/2/delta_y);

      if (nslit<=0 || theta_0 <= 0 || radius <=0)
      { fprintf(stderr,"DiskChopper: %s: nslit, theta_0 and radius must be > 0\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (nslit*theta_0 >= 360)
      { fprintf(stderr,"DiskChopper: %s: nslit * theta_0 exceeds 2PI\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (yheight && yheight>radius) {
        fprintf(stderr,"DiskChopper: %s: yheight must be < radius\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (isfirst && n_pulse <=0)
      { fprintf(stderr,"DiskChopper: %s: wrong First chopper pulse number (n_pulse=%g)\n", NAME_CURRENT_COMP, n_pulse);
        exit(-1); }
      if (!omega) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper frequency is 0!\n", NAME_CURRENT_COMP);
        omega = 1e-15; /* We should actually use machine epsilon here... */
      }
      if (!abs_out) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper will NOT absorb neutrons outside radius %g [m]\n", NAME_CURRENT_COMP, radius);
      }

      theta_0*=DEG2RAD;


      /* Calulate delay from phase and vice versa */
      if (phase) {
        if (delay) {
          fprintf(stderr,"DiskChopper: %s WARNING: delay AND phase specified. Using phase setting\n", NAME_CURRENT_COMP);
        }
        phase*=DEG2RAD;
        /* 'Delay' should always be a delay, taking rotation direction into account: */
        delay=phase/fabs(omega);
      } else {
        phase=delay*omega;  /* rad */
      }

      /* Time from opening of slit to next opening of slit */
      Tg=2.0*PI/fabs(omega)/nslit;

      /* How long can neutrons pass the Chopper at a single point */
      To=theta_0/fabs(omega);

      if (!xwidth) xwidth=2*delta_y*sin(theta_0/2);

      if (verbose && nu) {
        printf("DiskChopper: %s: frequency=%g [Hz] %g [rpm], time frame=%g [s] phase=%g [deg]\n",
          NAME_CURRENT_COMP, nu, nu*60, Tg, phase*RAD2DEG);
        printf("             %g slits, angle=%g [deg] height=%g [m], width=%g [m] at radius=%g [m]\n",
          nslit, theta_0*RAD2DEG, height, xwidth, delta_y);
      }
}
#line 11894 "let_new.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component betweenchop4andatt. */
  SIG_MESSAGE("betweenchop4andatt (Init)");
#define mccompcurname  betweenchop4andatt
#define mccompcurtype  Guide_channeled
#define mccompcurindex 13
#define w1c mccbetweenchop4andatt_w1c
#define w2c mccbetweenchop4andatt_w2c
#define ww mccbetweenchop4andatt_ww
#define hh mccbetweenchop4andatt_hh
#define whalf mccbetweenchop4andatt_whalf
#define hhalf mccbetweenchop4andatt_hhalf
#define lwhalf mccbetweenchop4andatt_lwhalf
#define lhhalf mccbetweenchop4andatt_lhhalf
#define w1 mccbetweenchop4andatt_w1
#define h1 mccbetweenchop4andatt_h1
#define w2 mccbetweenchop4andatt_w2
#define h2 mccbetweenchop4andatt_h2
#define l mccbetweenchop4andatt_l
#define R0 mccbetweenchop4andatt_R0
#define Qc mccbetweenchop4andatt_Qc
#define alpha mccbetweenchop4andatt_alpha
#define m mccbetweenchop4andatt_m
#define nslit mccbetweenchop4andatt_nslit
#define d mccbetweenchop4andatt_d
#define Qcx mccbetweenchop4andatt_Qcx
#define Qcy mccbetweenchop4andatt_Qcy
#define alphax mccbetweenchop4andatt_alphax
#define alphay mccbetweenchop4andatt_alphay
#define W mccbetweenchop4andatt_W
#define mx mccbetweenchop4andatt_mx
#define my mccbetweenchop4andatt_my
#define nu mccbetweenchop4andatt_nu
#define phase mccbetweenchop4andatt_phase
#line 88 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
if (!w2) w2=w1;
  if (!h2) h2=h1;
  if (nslit <= 0 || W <=0)
  { fprintf(stderr,"Guide_channeled: %s: nslit and W must be positive\n", NAME_CURRENT_COMP);
    exit(-1); }
  w1c = (w1 + d)/(double)nslit;
  w2c = (w2 + d)/(double)nslit;
  ww = .5*(w2c - w1c);
  hh = .5*(h2 - h1);
  whalf = .5*(w1c - d);
  hhalf = .5*h1;
  lwhalf = l*whalf;
  lhhalf = l*hhalf;

  if (m)     { mx=my=m; }
  if (Qc)    { Qcx=Qcy=Qc; }
  if (alpha) { alphax=alphay=alpha; }

  if ((nslit > 1) && (w1 != w2))
  {
    fprintf(stderr,"WARNING: Guide_channeled: %s:"
    "This component does not work with multichannel focusing guide\n"
    "Use Guide_gravity for that.\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_channeled: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (mcgravitation) fprintf(stderr,"WARNING: Guide_channeled: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  if (nu != 0 || phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_channeled: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_channeled: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, nu, phase);
    }
}
#line 11991 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component guide_move. */
  SIG_MESSAGE("guide_move (Init)");
#define mccompcurname  guide_move
#define mccompcurtype  Guide_channeled
#define mccompcurindex 14
#define w1c mccguide_move_w1c
#define w2c mccguide_move_w2c
#define ww mccguide_move_ww
#define hh mccguide_move_hh
#define whalf mccguide_move_whalf
#define hhalf mccguide_move_hhalf
#define lwhalf mccguide_move_lwhalf
#define lhhalf mccguide_move_lhhalf
#define w1 mccguide_move_w1
#define h1 mccguide_move_h1
#define w2 mccguide_move_w2
#define h2 mccguide_move_h2
#define l mccguide_move_l
#define R0 mccguide_move_R0
#define Qc mccguide_move_Qc
#define alpha mccguide_move_alpha
#define m mccguide_move_m
#define nslit mccguide_move_nslit
#define d mccguide_move_d
#define Qcx mccguide_move_Qcx
#define Qcy mccguide_move_Qcy
#define alphax mccguide_move_alphax
#define alphay mccguide_move_alphay
#define W mccguide_move_W
#define mx mccguide_move_mx
#define my mccguide_move_my
#define nu mccguide_move_nu
#define phase mccguide_move_phase
#line 88 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
if (!w2) w2=w1;
  if (!h2) h2=h1;
  if (nslit <= 0 || W <=0)
  { fprintf(stderr,"Guide_channeled: %s: nslit and W must be positive\n", NAME_CURRENT_COMP);
    exit(-1); }
  w1c = (w1 + d)/(double)nslit;
  w2c = (w2 + d)/(double)nslit;
  ww = .5*(w2c - w1c);
  hh = .5*(h2 - h1);
  whalf = .5*(w1c - d);
  hhalf = .5*h1;
  lwhalf = l*whalf;
  lhhalf = l*hhalf;

  if (m)     { mx=my=m; }
  if (Qc)    { Qcx=Qcy=Qc; }
  if (alpha) { alphax=alphay=alpha; }

  if ((nslit > 1) && (w1 != w2))
  {
    fprintf(stderr,"WARNING: Guide_channeled: %s:"
    "This component does not work with multichannel focusing guide\n"
    "Use Guide_gravity for that.\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_channeled: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (mcgravitation) fprintf(stderr,"WARNING: Guide_channeled: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  if (nu != 0 || phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_channeled: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_channeled: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, nu, phase);
    }
}
#line 12098 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component funnel. */
  SIG_MESSAGE("funnel (Init)");
#define mccompcurname  funnel
#define mccompcurtype  Guide_channeled
#define mccompcurindex 15
#define w1c mccfunnel_w1c
#define w2c mccfunnel_w2c
#define ww mccfunnel_ww
#define hh mccfunnel_hh
#define whalf mccfunnel_whalf
#define hhalf mccfunnel_hhalf
#define lwhalf mccfunnel_lwhalf
#define lhhalf mccfunnel_lhhalf
#define w1 mccfunnel_w1
#define h1 mccfunnel_h1
#define w2 mccfunnel_w2
#define h2 mccfunnel_h2
#define l mccfunnel_l
#define R0 mccfunnel_R0
#define Qc mccfunnel_Qc
#define alpha mccfunnel_alpha
#define m mccfunnel_m
#define nslit mccfunnel_nslit
#define d mccfunnel_d
#define Qcx mccfunnel_Qcx
#define Qcy mccfunnel_Qcy
#define alphax mccfunnel_alphax
#define alphay mccfunnel_alphay
#define W mccfunnel_W
#define mx mccfunnel_mx
#define my mccfunnel_my
#define nu mccfunnel_nu
#define phase mccfunnel_phase
#line 88 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
if (!w2) w2=w1;
  if (!h2) h2=h1;
  if (nslit <= 0 || W <=0)
  { fprintf(stderr,"Guide_channeled: %s: nslit and W must be positive\n", NAME_CURRENT_COMP);
    exit(-1); }
  w1c = (w1 + d)/(double)nslit;
  w2c = (w2 + d)/(double)nslit;
  ww = .5*(w2c - w1c);
  hh = .5*(h2 - h1);
  whalf = .5*(w1c - d);
  hhalf = .5*h1;
  lwhalf = l*whalf;
  lhhalf = l*hhalf;

  if (m)     { mx=my=m; }
  if (Qc)    { Qcx=Qcy=Qc; }
  if (alpha) { alphax=alphay=alpha; }

  if ((nslit > 1) && (w1 != w2))
  {
    fprintf(stderr,"WARNING: Guide_channeled: %s:"
    "This component does not work with multichannel focusing guide\n"
    "Use Guide_gravity for that.\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_channeled: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (mcgravitation) fprintf(stderr,"WARNING: Guide_channeled: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  if (nu != 0 || phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_channeled: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_channeled: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, nu, phase);
    }
}
#line 12205 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component C5. */
  SIG_MESSAGE("C5 (Init)");
#define mccompcurname  C5
#define mccompcurtype  DiskChopper
#define mccompcurindex 16
#define Tg mccC5_Tg
#define To mccC5_To
#define delta_y mccC5_delta_y
#define height mccC5_height
#define omega mccC5_omega
#define theta_0 mccC5_theta_0
#define radius mccC5_radius
#define yheight mccC5_yheight
#define nu mccC5_nu
#define nslit mccC5_nslit
#define jitter mccC5_jitter
#define delay mccC5_delay
#define isfirst mccC5_isfirst
#define n_pulse mccC5_n_pulse
#define abs_out mccC5_abs_out
#define phase mccC5_phase
#define xwidth mccC5_xwidth
#define verbose mccC5_verbose
#line 67 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{
/* If slit height 'unset', assume full opening */
if (yheight == 0) {
        height=radius;
      } else {
        height=yheight;
      }
      delta_y = radius-height/2; /* radius at beam center */
      omega=2.0*PI*nu; /* rad/s */
      if (xwidth && !theta_0 && radius) theta_0 = 2*RAD2DEG*asin(xwidth/2/delta_y);

      if (nslit<=0 || theta_0 <= 0 || radius <=0)
      { fprintf(stderr,"DiskChopper: %s: nslit, theta_0 and radius must be > 0\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (nslit*theta_0 >= 360)
      { fprintf(stderr,"DiskChopper: %s: nslit * theta_0 exceeds 2PI\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (yheight && yheight>radius) {
        fprintf(stderr,"DiskChopper: %s: yheight must be < radius\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (isfirst && n_pulse <=0)
      { fprintf(stderr,"DiskChopper: %s: wrong First chopper pulse number (n_pulse=%g)\n", NAME_CURRENT_COMP, n_pulse);
        exit(-1); }
      if (!omega) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper frequency is 0!\n", NAME_CURRENT_COMP);
        omega = 1e-15; /* We should actually use machine epsilon here... */
      }
      if (!abs_out) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper will NOT absorb neutrons outside radius %g [m]\n", NAME_CURRENT_COMP, radius);
      }

      theta_0*=DEG2RAD;


      /* Calulate delay from phase and vice versa */
      if (phase) {
        if (delay) {
          fprintf(stderr,"DiskChopper: %s WARNING: delay AND phase specified. Using phase setting\n", NAME_CURRENT_COMP);
        }
        phase*=DEG2RAD;
        /* 'Delay' should always be a delay, taking rotation direction into account: */
        delay=phase/fabs(omega);
      } else {
        phase=delay*omega;  /* rad */
      }

      /* Time from opening of slit to next opening of slit */
      Tg=2.0*PI/fabs(omega)/nslit;

      /* How long can neutrons pass the Chopper at a single point */
      To=theta_0/fabs(omega);

      if (!xwidth) xwidth=2*delta_y*sin(theta_0/2);

      if (verbose && nu) {
        printf("DiskChopper: %s: frequency=%g [Hz] %g [rpm], time frame=%g [s] phase=%g [deg]\n",
          NAME_CURRENT_COMP, nu, nu*60, Tg, phase*RAD2DEG);
        printf("             %g slits, angle=%g [deg] height=%g [m], width=%g [m] at radius=%g [m]\n",
          nslit, theta_0*RAD2DEG, height, xwidth, delta_y);
      }
}
#line 12323 "let_new.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component C5_counter. */
  SIG_MESSAGE("C5_counter (Init)");
#define mccompcurname  C5_counter
#define mccompcurtype  DiskChopper
#define mccompcurindex 17
#define Tg mccC5_counter_Tg
#define To mccC5_counter_To
#define delta_y mccC5_counter_delta_y
#define height mccC5_counter_height
#define omega mccC5_counter_omega
#define theta_0 mccC5_counter_theta_0
#define radius mccC5_counter_radius
#define yheight mccC5_counter_yheight
#define nu mccC5_counter_nu
#define nslit mccC5_counter_nslit
#define jitter mccC5_counter_jitter
#define delay mccC5_counter_delay
#define isfirst mccC5_counter_isfirst
#define n_pulse mccC5_counter_n_pulse
#define abs_out mccC5_counter_abs_out
#define phase mccC5_counter_phase
#define xwidth mccC5_counter_xwidth
#define verbose mccC5_counter_verbose
#line 67 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{
/* If slit height 'unset', assume full opening */
if (yheight == 0) {
        height=radius;
      } else {
        height=yheight;
      }
      delta_y = radius-height/2; /* radius at beam center */
      omega=2.0*PI*nu; /* rad/s */
      if (xwidth && !theta_0 && radius) theta_0 = 2*RAD2DEG*asin(xwidth/2/delta_y);

      if (nslit<=0 || theta_0 <= 0 || radius <=0)
      { fprintf(stderr,"DiskChopper: %s: nslit, theta_0 and radius must be > 0\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (nslit*theta_0 >= 360)
      { fprintf(stderr,"DiskChopper: %s: nslit * theta_0 exceeds 2PI\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (yheight && yheight>radius) {
        fprintf(stderr,"DiskChopper: %s: yheight must be < radius\n", NAME_CURRENT_COMP);
        exit(-1); }
      if (isfirst && n_pulse <=0)
      { fprintf(stderr,"DiskChopper: %s: wrong First chopper pulse number (n_pulse=%g)\n", NAME_CURRENT_COMP, n_pulse);
        exit(-1); }
      if (!omega) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper frequency is 0!\n", NAME_CURRENT_COMP);
        omega = 1e-15; /* We should actually use machine epsilon here... */
      }
      if (!abs_out) {
        fprintf(stderr,"DiskChopper: %s WARNING: chopper will NOT absorb neutrons outside radius %g [m]\n", NAME_CURRENT_COMP, radius);
      }

      theta_0*=DEG2RAD;


      /* Calulate delay from phase and vice versa */
      if (phase) {
        if (delay) {
          fprintf(stderr,"DiskChopper: %s WARNING: delay AND phase specified. Using phase setting\n", NAME_CURRENT_COMP);
        }
        phase*=DEG2RAD;
        /* 'Delay' should always be a delay, taking rotation direction into account: */
        delay=phase/fabs(omega);
      } else {
        phase=delay*omega;  /* rad */
      }

      /* Time from opening of slit to next opening of slit */
      Tg=2.0*PI/fabs(omega)/nslit;

      /* How long can neutrons pass the Chopper at a single point */
      To=theta_0/fabs(omega);

      if (!xwidth) xwidth=2*delta_y*sin(theta_0/2);

      if (verbose && nu) {
        printf("DiskChopper: %s: frequency=%g [Hz] %g [rpm], time frame=%g [s] phase=%g [deg]\n",
          NAME_CURRENT_COMP, nu, nu*60, Tg, phase*RAD2DEG);
        printf("             %g slits, angle=%g [deg] height=%g [m], width=%g [m] at radius=%g [m]\n",
          nslit, theta_0*RAD2DEG, height, xwidth, delta_y);
      }
}
#line 12431 "let_new.c"
#undef verbose
#undef xwidth
#undef phase
#undef abs_out
#undef n_pulse
#undef isfirst
#undef delay
#undef jitter
#undef nslit
#undef nu
#undef yheight
#undef radius
#undef theta_0
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component endguide. */
  SIG_MESSAGE("endguide (Init)");
#define mccompcurname  endguide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 18
#define w1c mccendguide_w1c
#define w2c mccendguide_w2c
#define ww mccendguide_ww
#define hh mccendguide_hh
#define whalf mccendguide_whalf
#define hhalf mccendguide_hhalf
#define lwhalf mccendguide_lwhalf
#define lhhalf mccendguide_lhhalf
#define w1 mccendguide_w1
#define h1 mccendguide_h1
#define w2 mccendguide_w2
#define h2 mccendguide_h2
#define l mccendguide_l
#define R0 mccendguide_R0
#define Qc mccendguide_Qc
#define alpha mccendguide_alpha
#define m mccendguide_m
#define nslit mccendguide_nslit
#define d mccendguide_d
#define Qcx mccendguide_Qcx
#define Qcy mccendguide_Qcy
#define alphax mccendguide_alphax
#define alphay mccendguide_alphay
#define W mccendguide_W
#define mx mccendguide_mx
#define my mccendguide_my
#define nu mccendguide_nu
#define phase mccendguide_phase
#line 88 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
if (!w2) w2=w1;
  if (!h2) h2=h1;
  if (nslit <= 0 || W <=0)
  { fprintf(stderr,"Guide_channeled: %s: nslit and W must be positive\n", NAME_CURRENT_COMP);
    exit(-1); }
  w1c = (w1 + d)/(double)nslit;
  w2c = (w2 + d)/(double)nslit;
  ww = .5*(w2c - w1c);
  hh = .5*(h2 - h1);
  whalf = .5*(w1c - d);
  hhalf = .5*h1;
  lwhalf = l*whalf;
  lhhalf = l*hhalf;

  if (m)     { mx=my=m; }
  if (Qc)    { Qcx=Qcy=Qc; }
  if (alpha) { alphax=alphay=alpha; }

  if ((nslit > 1) && (w1 != w2))
  {
    fprintf(stderr,"WARNING: Guide_channeled: %s:"
    "This component does not work with multichannel focusing guide\n"
    "Use Guide_gravity for that.\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_channeled: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (mcgravitation) fprintf(stderr,"WARNING: Guide_channeled: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  if (nu != 0 || phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_channeled: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_channeled: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, nu, phase);
    }
}
#line 12528 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component snot. */
  SIG_MESSAGE("snot (Init)");
#define mccompcurname  snot
#define mccompcurtype  Guide_channeled
#define mccompcurindex 19
#define w1c mccsnot_w1c
#define w2c mccsnot_w2c
#define ww mccsnot_ww
#define hh mccsnot_hh
#define whalf mccsnot_whalf
#define hhalf mccsnot_hhalf
#define lwhalf mccsnot_lwhalf
#define lhhalf mccsnot_lhhalf
#define w1 mccsnot_w1
#define h1 mccsnot_h1
#define w2 mccsnot_w2
#define h2 mccsnot_h2
#define l mccsnot_l
#define R0 mccsnot_R0
#define Qc mccsnot_Qc
#define alpha mccsnot_alpha
#define m mccsnot_m
#define nslit mccsnot_nslit
#define d mccsnot_d
#define Qcx mccsnot_Qcx
#define Qcy mccsnot_Qcy
#define alphax mccsnot_alphax
#define alphay mccsnot_alphay
#define W mccsnot_W
#define mx mccsnot_mx
#define my mccsnot_my
#define nu mccsnot_nu
#define phase mccsnot_phase
#line 88 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
if (!w2) w2=w1;
  if (!h2) h2=h1;
  if (nslit <= 0 || W <=0)
  { fprintf(stderr,"Guide_channeled: %s: nslit and W must be positive\n", NAME_CURRENT_COMP);
    exit(-1); }
  w1c = (w1 + d)/(double)nslit;
  w2c = (w2 + d)/(double)nslit;
  ww = .5*(w2c - w1c);
  hh = .5*(h2 - h1);
  whalf = .5*(w1c - d);
  hhalf = .5*h1;
  lwhalf = l*whalf;
  lhhalf = l*hhalf;

  if (m)     { mx=my=m; }
  if (Qc)    { Qcx=Qcy=Qc; }
  if (alpha) { alphax=alphay=alpha; }

  if ((nslit > 1) && (w1 != w2))
  {
    fprintf(stderr,"WARNING: Guide_channeled: %s:"
    "This component does not work with multichannel focusing guide\n"
    "Use Guide_gravity for that.\n", NAME_CURRENT_COMP);
    exit(-1);
  }

  if (d*nslit > w1) exit(fprintf(stderr, "Guide_channeled: %s: absorbing walls fill input window. No space left for transmission (d*nslit > w1).\n", NAME_CURRENT_COMP));

  if (mcgravitation) fprintf(stderr,"WARNING: Guide_channeled: %s: "
    "This component produces wrong results with gravitation !\n"
    "Use Guide_gravity.\n",
    NAME_CURRENT_COMP);
  if (nu != 0 || phase != 0) {
      if (w1 != w2 || h1 != h2)
      exit(fprintf(stderr,"Guide_channeled: %s: rotating slit pack must be straight (w1=w2 and h1=h2).\n", NAME_CURRENT_COMP));
      printf("Guide_channeled: %s: Fermi Chopper mode: frequency=%g [Hz] phase=%g [deg]\n",
        NAME_CURRENT_COMP, nu, phase);
    }
}
#line 12635 "let_new.c"
#undef phase
#undef nu
#undef my
#undef mx
#undef W
#undef alphay
#undef alphax
#undef Qcy
#undef Qcx
#undef d
#undef nslit
#undef m
#undef alpha
#undef Qc
#undef R0
#undef l
#undef h2
#undef w2
#undef h1
#undef w1
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component detslit. */
  SIG_MESSAGE("detslit (Init)");
#define mccompcurname  detslit
#define mccompcurtype  Slit
#define mccompcurindex 20
#define xmin mccdetslit_xmin
#define xmax mccdetslit_xmax
#define ymin mccdetslit_ymin
#define ymax mccdetslit_ymax
#define radius mccdetslit_radius
#define xwidth mccdetslit_xwidth
#define yheight mccdetslit_yheight
#line 47 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Slit.comp"
{
if (xwidth > 0)  { xmax=xwidth/2;  xmin=-xmax; }
  if (yheight > 0) { ymax=yheight/2; ymin=-ymax; }
  if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Error: give geometry\n", NAME_CURRENT_COMP); exit(-1); }

}
#line 12688 "let_new.c"
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component mcplout. */
  SIG_MESSAGE("mcplout (Init)");
#define mccompcurname  mcplout
#define mccompcurtype  MCPL_output_horace
#define mccompcurindex 21
#define polarisationuse mccmcplout_polarisationuse
#define doubleprec mccmcplout_doubleprec
#define verbose mccmcplout_verbose
#define radius mccmcplout_radius
#define thickness mccmcplout_thickness
#define yheight mccmcplout_yheight
#define restore_neutron mccmcplout_restore_neutron
#define userflag mccmcplout_userflag
#define plate mccmcplout_plate
#define outputfile mccmcplout_outputfile
#define particle mccmcplout_particle
#define Particle mccmcplout_Particle
#define userflagenabled mccmcplout_userflagenabled
#define filename mccmcplout_filename
#define userflagcomment mccmcplout_userflagcomment
#define merge_mpi mccmcplout_merge_mpi
#define keep_mpi_unmerged mccmcplout_keep_mpi_unmerged
#line 82 "MCPL_output_horace.comp"
{
    char extension[128]="";
    char *myfilename;

#if defined (USE_MPI)
  /* In case of MPI, simply redefine the filename used by each node */
    MPI_MASTER(fprintf(stdout, "Message(%s): You are using MCPL_output with MPI, hence your will get %i filenames %s_node_#i as output.\n",NAME_CURRENT_COMP,mpi_node_count,filename); );
    sprintf(extension,"node_%i.mcpl",mpi_node_rank);
#else
    sprintf(extension,"mcpl");
#endif
    /*add output dir (if applicable) to the output filename and add extension if */
    myfilename=mcfull_file(filename,extension);

    char line[256];
    outputfile = mcpl_create_outfile(myfilename);
    /*reset filename to be whatever mcpl actually calls it. It may have added .mcpl*/
    snprintf(myfilename,strlen(myfilename)+5,"%s",mcpl_outfile_filename(outputfile));

    snprintf(line,255,"%s %s %s",MCCODE_NAME,MCCODE_VERSION,mcinstrument_name);
    mcpl_hdr_set_srcname(outputfile,line);
    mcpl_enable_universal_pdgcode(outputfile,2112);/*all particles are neutrons*/
    snprintf(line,255,"Output by COMPONENT: %s",NAME_CURRENT_COMP);
    mcpl_hdr_add_comment(outputfile,line);

    /*also add the instrument file and the command line as blobs*/
    FILE *fp;
    if( (fp=fopen(mcinstrument_source,"rb"))!=NULL){
        unsigned char *buffer;
        int size,status;
        /*find the file size by seeking to end, "tell" the position, and then go back again*/
        fseek(fp, 0L, SEEK_END);
        size = ftell(fp); // get current file pointer
        fseek(fp, 0L, SEEK_SET); // seek back to beginning of file
        if ( size && (buffer=malloc(size))!=NULL){
            if (size!=(fread(buffer,1,size,fp))){
	      fprintf(stderr,"Warning (%s): Source instrument file not read cleanly\n", NAME_CURRENT_COMP);
            }
            mcpl_hdr_add_data(outputfile, "mccode_instr_file", size, buffer);
            free(buffer);
        }
    }
    fclose(fp);

    int ii;
    char clr[2048],*clrp;
    char parval[CHAR_BUF_LENGTH];
    clrp=clr;
    clrp+=snprintf(clrp,2048,"%s",mcinstrument_exe);
    for (ii=0;ii<mcnumipar;ii++){
        if (mcinputtable[ii].par == NULL)
            strncpy(parval, (mcinputtable[ii].val ? mcinputtable[ii].val : ""), CHAR_BUF_LENGTH);
        else
            (*mcinputtypes[mcinputtable[ii].type].printer)(parval, mcinputtable[ii].par);
        clrp+=snprintf(clrp,2048-(clrp-clr)," %s=%s",mcinputtable[ii].name, parval);
    }
    *(clrp)='\0';
    mcpl_hdr_add_data(outputfile, "mccode_cmd_line" , strlen(clr), clr);

    if (polarisationuse) {
        mcpl_enable_polarisation(outputfile);
    }
    if (doubleprec){
        mcpl_enable_doubleprec(outputfile);
    }

#if defined (USE_MPI)
  MPI_MASTER(
#endif

    if (verbose==1) {
    printf("MCPL_output verbose mode: after generating the mcpl-file it will be reread and a summary printed.\n");
    }

#if defined (USE_MPI)
	    );
#endif

  /*pointer to the single particle storage area*/
  particle=&Particle;

  /*Add comments on what the orientation and position of this component is.*/
  /*Include the instrument file itself as a binary blob in the mcpl file*/

  userflagenabled=0;
  /*Have the option of including a user-flag like they do at Loki.*/
  if (strlen(userflagcomment)!=0){
      mcpl_enable_userflags(outputfile);
      userflagenabled=1;
      /*Don't add the comment if it's empty*/
      if(userflagcomment && strlen(userflagcomment)){
          snprintf(line,255,"userflags: %s",userflagcomment);
          mcpl_hdr_add_comment(outputfile,line);
      }
  }
   if (myfilename){
       free(myfilename);
   }

}
#line 12823 "let_new.c"
#undef keep_mpi_unmerged
#undef merge_mpi
#undef userflagcomment
#undef filename
#undef userflagenabled
#undef Particle
#undef particle
#undef outputfile
#undef plate
#undef userflag
#undef restore_neutron
#undef yheight
#undef thickness
#undef radius
#undef verbose
#undef doubleprec
#undef polarisationuse
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if(mcdotrace) mcdisplay();
    mcDEBUG_INSTR_END()
  }

} /* end init */

void mcraytrace(void) {
  /* Neutronics-specific defines */
#ifdef NEUTRONICS
extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;
#endif
  /* End of Neutronics-specific defines */
  /* Copy neutron state to local variables. */
  MCNUM mcnlx = mcnx;
  MCNUM mcnly = mcny;
  MCNUM mcnlz = mcnz;
  MCNUM mcnlvx = mcnvx;
  MCNUM mcnlvy = mcnvy;
  MCNUM mcnlvz = mcnvz;
  MCNUM mcnlt = mcnt;
  MCNUM mcnlsx = mcnsx;
  MCNUM mcnlsy = mcnsy;
  MCNUM mcnlsz = mcnsz;
  MCNUM mcnlp = mcnp;

  mcDEBUG_ENTER()
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define mcabsorb mcabsorbAll
  /* TRACE Component Origin [1] */
  mccoordschange(mcposrOrigin, mcrotrOrigin,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Origin (without coords transformations) */
  mcJumpTrace_Origin:
  SIG_MESSAGE("Origin (Trace)");
  mcDEBUG_COMP("Origin")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompOrigin
  STORE_NEUTRON(1,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[1]++;
  mcPCounter[1] += p;
  mcP2Counter[1] += p*p;
#define mccompcurname  Origin
#define mccompcurtype  Arm
#define mccompcurindex 1
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompOrigin:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(1,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Source [2] */
  mccoordschange(mcposrSource, mcrotrSource,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Source (without coords transformations) */
  mcJumpTrace_Source:
  SIG_MESSAGE("Source (Trace)");
  mcDEBUG_COMP("Source")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompSource
  STORE_NEUTRON(2,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[2]++;
  mcPCounter[2] += p;
  mcP2Counter[2] += p*p;
#define mccompcurname  Source
#define mccompcurtype  ViewModISISver1
#define mccompcurindex 2
{   /* Declarations of Source=ViewModISISver1() SETTING parameters. */
char* Face = mccSource_Face;
MCNUM E0 = mccSource_E0;
MCNUM E1 = mccSource_E1;
MCNUM tally = mccSource_tally;
MCNUM modPosition = mccSource_modPosition;
MCNUM modXsize = mccSource_modXsize;
MCNUM modZsize = mccSource_modZsize;
MCNUM xw = mccSource_xw;
MCNUM yh = mccSource_yh;
MCNUM dist = mccSource_dist;
#line 1067 "ViewModISISver1.comp"
{
  double v,r,E;
  double xf,yf,dx,dy;    /* mxp ->max var in param space */
  double Ival,Tval,Eval;

  Ncount++;

  x += TS.XAxis*(0.5-rand01());            /* Get point on shutter enterance */
  y += TS.ZAxis*(0.5-rand01());            /* Get point on shutter enterance */
    
  // if (modPosition) z=TS.rdumMid;    //   Goran          

  xf = xw*(0.5-rand01());          /* Choose focusing position uniformly */
  yf = yh*(0.5-rand01());          /* Choose focusing position uniformly */

  getPoint(&Tval,&Eval,&rtE0,&rtE1);
  
//  Ival=TS.Total*6.2415093e+12; // * 1.1879451;  /* ( of proton in 1uAmp) * (1-cos(30))*2*Pi  */
//  Ival=TS.Total*5.9294338e+12; // Number of proton in 1uAmp * 0.95 (0.95 because of muon target loss)
  Ival=TS.Total*6.2415093e+12; // Number of proton in 1uAmp

  dx = xf-x;
  dy = yf-y;
  r = sqrt(dx*dx+dy*dy+dist*dist);      // Actual distance to point 
  v = SE2V*sqrt(Eval);                  // Calculate the velocity 
  vz = v*fabs(dist)/r;
  vy = v*dy/r;
  vx = v*dx/r;


   t=Tval-(TS.rdumMid-TS.timeOffset)/vz; ///vz;
   
     if (modPosition)
     {
     t+=TS.rdumMid/vz;
     }

  p=angleArea*Ival/Nsim;

  //   testing only - Goran              
  //  if (!modPosition) p*=TS.timeOffset;


  if (!(Ncount % 100000))
    {
      fprintf(stderr,"FPos[%d]=> %g %g %g  \n",Ncount,x,y,z);
      fprintf(stderr,"FDir[%d]=> %g %g %g  \n",Ncount,vx,vy,vz);
      fprintf(stderr,"PlaneAxis %g %g \n",TS.XAxis,fullAngle);
      fprintf(stderr,"RMID %g \n",TS.rdumMid);
      fprintf(stderr,"TimeMid[%d]=> %g\n",Ncount,TS.rdumMid);
      fprintf(stderr,"TimeOffset[%d]=> %g\n",Ncount,TS.timeOffset);
      fprintf(stderr,"TimeTval[%d]=> %g\n",Ncount,Tval);
      fprintf(stderr,"TimeZero[%d]=> %g\n",Ncount,t);
    }  

}
#line 13116 "let_new.c"
}   /* End of Source=ViewModISISver1() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompSource:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(2,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component shutter [3] */
  mccoordschange(mcposrshutter, mcrotrshutter,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component shutter (without coords transformations) */
  mcJumpTrace_shutter:
  SIG_MESSAGE("shutter (Trace)");
  mcDEBUG_COMP("shutter")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompshutter
  STORE_NEUTRON(3,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[3]++;
  mcPCounter[3] += p;
  mcP2Counter[3] += p*p;
#define mccompcurname  shutter
#define mccompcurtype  Guide_channeled
#define mccompcurindex 3
#define w1c mccshutter_w1c
#define w2c mccshutter_w2c
#define ww mccshutter_ww
#define hh mccshutter_hh
#define whalf mccshutter_whalf
#define hhalf mccshutter_hhalf
#define lwhalf mccshutter_lwhalf
#define lhhalf mccshutter_lhhalf
{   /* Declarations of shutter=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccshutter_w1;
MCNUM h1 = mccshutter_h1;
MCNUM w2 = mccshutter_w2;
MCNUM h2 = mccshutter_h2;
MCNUM l = mccshutter_l;
MCNUM R0 = mccshutter_R0;
MCNUM Qc = mccshutter_Qc;
MCNUM alpha = mccshutter_alpha;
MCNUM m = mccshutter_m;
MCNUM nslit = mccshutter_nslit;
MCNUM d = mccshutter_d;
MCNUM Qcx = mccshutter_Qcx;
MCNUM Qcy = mccshutter_Qcy;
MCNUM alphax = mccshutter_alphax;
MCNUM alphay = mccshutter_alphay;
MCNUM W = mccshutter_W;
MCNUM mx = mccshutter_mx;
MCNUM my = mccshutter_my;
MCNUM nu = mccshutter_nu;
MCNUM phase = mccshutter_phase;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double angle=0;

  if (nu != 0 || phase != 0) { /* rotate neutron w/r to guide element */
    /* approximation of rotating straight Fermi Chopper */
    Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
    Rotation R;
    double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
    angle=fmod(360*nu*(t+dt)+phase, 360); /* in deg */
    /* modify angle so that Z0 guide side is always in front of incoming neutron */
    if (angle > 90 && angle < 270) { angle -= 180; }
    angle *= DEG2RAD;
    rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
    /* apply rotation to centered coordinates */
    Coords   RX = rot_apply(R, X);
    coords_get(RX, &x, &y, &z);
    z = z+l/2;
    /* rotate speed */
    X  = coords_set(vx,vy,vz);
    RX = rot_apply(R, X);
    coords_get(RX, &vx, &vy, &vz);
  }

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= w1/-2.0 || x >= w1/2.0 || y <= -hhalf || y >= hhalf)
    ABSORB;
  /* Shift origin to center of channel hit (absorb if hit dividing walls) */
  x += w1/2.0;
  edge = floor(x/w1c)*w1c;
  if(x - edge > w1c - d)
  {
    x -= w1/2.0; /* Re-adjust origin */
    ABSORB;
  }
  x -= (edge + (w1c - d)/2.0);
  hadj = edge + (w1c - d)/2.0 - w1/2.0;
  for(;;)
  {
    /* Compute the dot products of v and n for the four mirrors. */
    av = l*vx; bv = ww*vz;
    ah = l*vy; bh = hh*vz;
    vdotn_v1 = bv + av;         /* Left vertical */
    vdotn_v2 = bv - av;         /* Right vertical */
    vdotn_h1 = bh + ah;         /* Lower horizontal */
    vdotn_h2 = bh - ah;         /* Upper horizontal */
    /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
    cv1 = -whalf*l - z*ww; cv2 = x*l;
    ch1 = -hhalf*l - z*hh; ch2 = y*l;
    /* Compute intersection times. */
    t1 = (l - z)/vz;
    i = 0;
    if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
    {
      t1 = t2;
      i = 1;
    }
    if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
    {
      t1 = t2;
      i = 2;
    }
    if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
    {
      t1 = t2;
      i = 3;
    }
    if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
    {
      t1 = t2;
      i = 4;
    }
    if(i == 0)
      break;                    /* Neutron left guide. */
    PROP_DT(t1);
    switch(i)
    {
      case 1:                   /* Left vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
        dd = 2*vdotn_v1/nlen2;
        vx = vx - dd*l;
        vz = vz - dd*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        dd = 2*vdotn_v2/nlen2;
        vx = vx + dd*l;
        vz = vz - dd*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        dd = 2*vdotn_h1/nlen2;
        vy = vy - dd*l;
        vz = vz - dd*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        dd = 2*vdotn_h2/nlen2;
        vy = vy + dd*l;
        vz = vz - dd*hh;
        break;
    }
    /* Now compute reflectivity. */
    if((i <= 2 && mx == 0) || (i > 2 && my == 0))
    {
      x += hadj; /* Re-adjust origin */
      ABSORB;
    } else {
      double ref=1;
      if (i <= 2)
      {
        double par[] = {R0, Qcx, alphax, mx, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      } else {
        double par[] = {R0, Qcy, alphay, my, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      }
    }
    x += hadj; SCATTER; x -= hadj;
  } /* end for */
  x += hadj; /* Re-adjust origin */
  if (nu != 0 || phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }
}
#line 13415 "let_new.c"
}   /* End of shutter=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompshutter:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(3,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component insert [4] */
  mccoordschange(mcposrinsert, mcrotrinsert,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component insert (without coords transformations) */
  mcJumpTrace_insert:
  SIG_MESSAGE("insert (Trace)");
  mcDEBUG_COMP("insert")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompinsert
  STORE_NEUTRON(4,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[4]++;
  mcPCounter[4] += p;
  mcP2Counter[4] += p*p;
#define mccompcurname  insert
#define mccompcurtype  Guide_channeled
#define mccompcurindex 4
#define w1c mccinsert_w1c
#define w2c mccinsert_w2c
#define ww mccinsert_ww
#define hh mccinsert_hh
#define whalf mccinsert_whalf
#define hhalf mccinsert_hhalf
#define lwhalf mccinsert_lwhalf
#define lhhalf mccinsert_lhhalf
{   /* Declarations of insert=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccinsert_w1;
MCNUM h1 = mccinsert_h1;
MCNUM w2 = mccinsert_w2;
MCNUM h2 = mccinsert_h2;
MCNUM l = mccinsert_l;
MCNUM R0 = mccinsert_R0;
MCNUM Qc = mccinsert_Qc;
MCNUM alpha = mccinsert_alpha;
MCNUM m = mccinsert_m;
MCNUM nslit = mccinsert_nslit;
MCNUM d = mccinsert_d;
MCNUM Qcx = mccinsert_Qcx;
MCNUM Qcy = mccinsert_Qcy;
MCNUM alphax = mccinsert_alphax;
MCNUM alphay = mccinsert_alphay;
MCNUM W = mccinsert_W;
MCNUM mx = mccinsert_mx;
MCNUM my = mccinsert_my;
MCNUM nu = mccinsert_nu;
MCNUM phase = mccinsert_phase;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double angle=0;

  if (nu != 0 || phase != 0) { /* rotate neutron w/r to guide element */
    /* approximation of rotating straight Fermi Chopper */
    Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
    Rotation R;
    double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
    angle=fmod(360*nu*(t+dt)+phase, 360); /* in deg */
    /* modify angle so that Z0 guide side is always in front of incoming neutron */
    if (angle > 90 && angle < 270) { angle -= 180; }
    angle *= DEG2RAD;
    rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
    /* apply rotation to centered coordinates */
    Coords   RX = rot_apply(R, X);
    coords_get(RX, &x, &y, &z);
    z = z+l/2;
    /* rotate speed */
    X  = coords_set(vx,vy,vz);
    RX = rot_apply(R, X);
    coords_get(RX, &vx, &vy, &vz);
  }

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= w1/-2.0 || x >= w1/2.0 || y <= -hhalf || y >= hhalf)
    ABSORB;
  /* Shift origin to center of channel hit (absorb if hit dividing walls) */
  x += w1/2.0;
  edge = floor(x/w1c)*w1c;
  if(x - edge > w1c - d)
  {
    x -= w1/2.0; /* Re-adjust origin */
    ABSORB;
  }
  x -= (edge + (w1c - d)/2.0);
  hadj = edge + (w1c - d)/2.0 - w1/2.0;
  for(;;)
  {
    /* Compute the dot products of v and n for the four mirrors. */
    av = l*vx; bv = ww*vz;
    ah = l*vy; bh = hh*vz;
    vdotn_v1 = bv + av;         /* Left vertical */
    vdotn_v2 = bv - av;         /* Right vertical */
    vdotn_h1 = bh + ah;         /* Lower horizontal */
    vdotn_h2 = bh - ah;         /* Upper horizontal */
    /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
    cv1 = -whalf*l - z*ww; cv2 = x*l;
    ch1 = -hhalf*l - z*hh; ch2 = y*l;
    /* Compute intersection times. */
    t1 = (l - z)/vz;
    i = 0;
    if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
    {
      t1 = t2;
      i = 1;
    }
    if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
    {
      t1 = t2;
      i = 2;
    }
    if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
    {
      t1 = t2;
      i = 3;
    }
    if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
    {
      t1 = t2;
      i = 4;
    }
    if(i == 0)
      break;                    /* Neutron left guide. */
    PROP_DT(t1);
    switch(i)
    {
      case 1:                   /* Left vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
        dd = 2*vdotn_v1/nlen2;
        vx = vx - dd*l;
        vz = vz - dd*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        dd = 2*vdotn_v2/nlen2;
        vx = vx + dd*l;
        vz = vz - dd*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        dd = 2*vdotn_h1/nlen2;
        vy = vy - dd*l;
        vz = vz - dd*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        dd = 2*vdotn_h2/nlen2;
        vy = vy + dd*l;
        vz = vz - dd*hh;
        break;
    }
    /* Now compute reflectivity. */
    if((i <= 2 && mx == 0) || (i > 2 && my == 0))
    {
      x += hadj; /* Re-adjust origin */
      ABSORB;
    } else {
      double ref=1;
      if (i <= 2)
      {
        double par[] = {R0, Qcx, alphax, mx, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      } else {
        double par[] = {R0, Qcy, alphay, my, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      }
    }
    x += hadj; SCATTER; x -= hadj;
  } /* end for */
  x += hadj; /* Re-adjust origin */
  if (nu != 0 || phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }
}
#line 13722 "let_new.c"
}   /* End of insert=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompinsert:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(4,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component afterinsert [5] */
  mccoordschange(mcposrafterinsert, mcrotrafterinsert,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component afterinsert (without coords transformations) */
  mcJumpTrace_afterinsert:
  SIG_MESSAGE("afterinsert (Trace)");
  mcDEBUG_COMP("afterinsert")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompafterinsert
  STORE_NEUTRON(5,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[5]++;
  mcPCounter[5] += p;
  mcP2Counter[5] += p*p;
#define mccompcurname  afterinsert
#define mccompcurtype  Guide_channeled
#define mccompcurindex 5
#define w1c mccafterinsert_w1c
#define w2c mccafterinsert_w2c
#define ww mccafterinsert_ww
#define hh mccafterinsert_hh
#define whalf mccafterinsert_whalf
#define hhalf mccafterinsert_hhalf
#define lwhalf mccafterinsert_lwhalf
#define lhhalf mccafterinsert_lhhalf
{   /* Declarations of afterinsert=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccafterinsert_w1;
MCNUM h1 = mccafterinsert_h1;
MCNUM w2 = mccafterinsert_w2;
MCNUM h2 = mccafterinsert_h2;
MCNUM l = mccafterinsert_l;
MCNUM R0 = mccafterinsert_R0;
MCNUM Qc = mccafterinsert_Qc;
MCNUM alpha = mccafterinsert_alpha;
MCNUM m = mccafterinsert_m;
MCNUM nslit = mccafterinsert_nslit;
MCNUM d = mccafterinsert_d;
MCNUM Qcx = mccafterinsert_Qcx;
MCNUM Qcy = mccafterinsert_Qcy;
MCNUM alphax = mccafterinsert_alphax;
MCNUM alphay = mccafterinsert_alphay;
MCNUM W = mccafterinsert_W;
MCNUM mx = mccafterinsert_mx;
MCNUM my = mccafterinsert_my;
MCNUM nu = mccafterinsert_nu;
MCNUM phase = mccafterinsert_phase;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double angle=0;

  if (nu != 0 || phase != 0) { /* rotate neutron w/r to guide element */
    /* approximation of rotating straight Fermi Chopper */
    Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
    Rotation R;
    double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
    angle=fmod(360*nu*(t+dt)+phase, 360); /* in deg */
    /* modify angle so that Z0 guide side is always in front of incoming neutron */
    if (angle > 90 && angle < 270) { angle -= 180; }
    angle *= DEG2RAD;
    rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
    /* apply rotation to centered coordinates */
    Coords   RX = rot_apply(R, X);
    coords_get(RX, &x, &y, &z);
    z = z+l/2;
    /* rotate speed */
    X  = coords_set(vx,vy,vz);
    RX = rot_apply(R, X);
    coords_get(RX, &vx, &vy, &vz);
  }

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= w1/-2.0 || x >= w1/2.0 || y <= -hhalf || y >= hhalf)
    ABSORB;
  /* Shift origin to center of channel hit (absorb if hit dividing walls) */
  x += w1/2.0;
  edge = floor(x/w1c)*w1c;
  if(x - edge > w1c - d)
  {
    x -= w1/2.0; /* Re-adjust origin */
    ABSORB;
  }
  x -= (edge + (w1c - d)/2.0);
  hadj = edge + (w1c - d)/2.0 - w1/2.0;
  for(;;)
  {
    /* Compute the dot products of v and n for the four mirrors. */
    av = l*vx; bv = ww*vz;
    ah = l*vy; bh = hh*vz;
    vdotn_v1 = bv + av;         /* Left vertical */
    vdotn_v2 = bv - av;         /* Right vertical */
    vdotn_h1 = bh + ah;         /* Lower horizontal */
    vdotn_h2 = bh - ah;         /* Upper horizontal */
    /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
    cv1 = -whalf*l - z*ww; cv2 = x*l;
    ch1 = -hhalf*l - z*hh; ch2 = y*l;
    /* Compute intersection times. */
    t1 = (l - z)/vz;
    i = 0;
    if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
    {
      t1 = t2;
      i = 1;
    }
    if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
    {
      t1 = t2;
      i = 2;
    }
    if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
    {
      t1 = t2;
      i = 3;
    }
    if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
    {
      t1 = t2;
      i = 4;
    }
    if(i == 0)
      break;                    /* Neutron left guide. */
    PROP_DT(t1);
    switch(i)
    {
      case 1:                   /* Left vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
        dd = 2*vdotn_v1/nlen2;
        vx = vx - dd*l;
        vz = vz - dd*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        dd = 2*vdotn_v2/nlen2;
        vx = vx + dd*l;
        vz = vz - dd*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        dd = 2*vdotn_h1/nlen2;
        vy = vy - dd*l;
        vz = vz - dd*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        dd = 2*vdotn_h2/nlen2;
        vy = vy + dd*l;
        vz = vz - dd*hh;
        break;
    }
    /* Now compute reflectivity. */
    if((i <= 2 && mx == 0) || (i > 2 && my == 0))
    {
      x += hadj; /* Re-adjust origin */
      ABSORB;
    } else {
      double ref=1;
      if (i <= 2)
      {
        double par[] = {R0, Qcx, alphax, mx, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      } else {
        double par[] = {R0, Qcy, alphay, my, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      }
    }
    x += hadj; SCATTER; x -= hadj;
  } /* end for */
  x += hadj; /* Re-adjust origin */
  if (nu != 0 || phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }
}
#line 14029 "let_new.c"
}   /* End of afterinsert=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompafterinsert:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(5,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Res1 [6] */
  mccoordschange(mcposrRes1, mcrotrRes1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Res1 (without coords transformations) */
  mcJumpTrace_Res1:
  SIG_MESSAGE("Res1 (Trace)");
  mcDEBUG_COMP("Res1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompRes1
  STORE_NEUTRON(6,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[6]++;
  mcPCounter[6] += p;
  mcP2Counter[6] += p*p;
#define mccompcurname  Res1
#define mccompcurtype  DiskChopper
#define mccompcurindex 6
#define Tg mccRes1_Tg
#define To mccRes1_To
#define delta_y mccRes1_delta_y
#define height mccRes1_height
#define omega mccRes1_omega
{   /* Declarations of Res1=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccRes1_theta_0;
MCNUM radius = mccRes1_radius;
MCNUM yheight = mccRes1_yheight;
MCNUM nu = mccRes1_nu;
MCNUM nslit = mccRes1_nslit;
MCNUM jitter = mccRes1_jitter;
MCNUM delay = mccRes1_delay;
MCNUM isfirst = mccRes1_isfirst;
MCNUM n_pulse = mccRes1_n_pulse;
MCNUM abs_out = mccRes1_abs_out;
MCNUM phase = mccRes1_phase;
MCNUM xwidth = mccRes1_xwidth;
MCNUM verbose = mccRes1_verbose;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{
    double toff;
    double yprime;
    PROP_Z0;
    yprime = y+delta_y;

    /* Is neutron outside the vertical slit range and should we absorb? */
    if (abs_out && (x*x+yprime*yprime)>radius*radius) {
      ABSORB;
    }
    /* Does neutron hit inner solid part of chopper in case of yheight!=radius? */
    if ((x*x+yprime*yprime)<(radius-height)*(radius-height)) {
      ABSORB;
    }


    if (isfirst)
      {
        /* all events are put in the transmitted time frame */
        t=atan2(x,yprime)/omega + To*randpm1()/2.0 + delay + (jitter ? jitter*randnorm():0) + (n_pulse > 1 ? floor(n_pulse*rand01())*Tg : 0);
        /* correction: chopper slits transmission opening/full disk */
        p *= nslit*theta_0/2.0/PI;
      }
    else
      {
        toff=fabs(t-atan2(x,yprime)/omega - delay - (jitter ? jitter*randnorm():0));

        /* does neutron hit outside slit? */
        if (fmod(toff+To/2.0,Tg)>To) ABSORB;
      }
    SCATTER;

}
#line 14195 "let_new.c"
}   /* End of Res1=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompRes1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(6,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component Res1_counter [7] */
  mccoordschange(mcposrRes1_counter, mcrotrRes1_counter,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component Res1_counter (without coords transformations) */
  mcJumpTrace_Res1_counter:
  SIG_MESSAGE("Res1_counter (Trace)");
  mcDEBUG_COMP("Res1_counter")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompRes1_counter
  STORE_NEUTRON(7,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[7]++;
  mcPCounter[7] += p;
  mcP2Counter[7] += p*p;
#define mccompcurname  Res1_counter
#define mccompcurtype  DiskChopper
#define mccompcurindex 7
#define Tg mccRes1_counter_Tg
#define To mccRes1_counter_To
#define delta_y mccRes1_counter_delta_y
#define height mccRes1_counter_height
#define omega mccRes1_counter_omega
{   /* Declarations of Res1_counter=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccRes1_counter_theta_0;
MCNUM radius = mccRes1_counter_radius;
MCNUM yheight = mccRes1_counter_yheight;
MCNUM nu = mccRes1_counter_nu;
MCNUM nslit = mccRes1_counter_nslit;
MCNUM jitter = mccRes1_counter_jitter;
MCNUM delay = mccRes1_counter_delay;
MCNUM isfirst = mccRes1_counter_isfirst;
MCNUM n_pulse = mccRes1_counter_n_pulse;
MCNUM abs_out = mccRes1_counter_abs_out;
MCNUM phase = mccRes1_counter_phase;
MCNUM xwidth = mccRes1_counter_xwidth;
MCNUM verbose = mccRes1_counter_verbose;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{
    double toff;
    double yprime;
    PROP_Z0;
    yprime = y+delta_y;

    /* Is neutron outside the vertical slit range and should we absorb? */
    if (abs_out && (x*x+yprime*yprime)>radius*radius) {
      ABSORB;
    }
    /* Does neutron hit inner solid part of chopper in case of yheight!=radius? */
    if ((x*x+yprime*yprime)<(radius-height)*(radius-height)) {
      ABSORB;
    }


    if (isfirst)
      {
        /* all events are put in the transmitted time frame */
        t=atan2(x,yprime)/omega + To*randpm1()/2.0 + delay + (jitter ? jitter*randnorm():0) + (n_pulse > 1 ? floor(n_pulse*rand01())*Tg : 0);
        /* correction: chopper slits transmission opening/full disk */
        p *= nslit*theta_0/2.0/PI;
      }
    else
      {
        toff=fabs(t-atan2(x,yprime)/omega - delay - (jitter ? jitter*randnorm():0));

        /* does neutron hit outside slit? */
        if (fmod(toff+To/2.0,Tg)>To) ABSORB;
      }
    SCATTER;

}
#line 14358 "let_new.c"
}   /* End of Res1_counter=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompRes1_counter:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(7,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component betweenchop1and2 [8] */
  mccoordschange(mcposrbetweenchop1and2, mcrotrbetweenchop1and2,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component betweenchop1and2 (without coords transformations) */
  mcJumpTrace_betweenchop1and2:
  SIG_MESSAGE("betweenchop1and2 (Trace)");
  mcDEBUG_COMP("betweenchop1and2")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompbetweenchop1and2
  STORE_NEUTRON(8,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[8]++;
  mcPCounter[8] += p;
  mcP2Counter[8] += p*p;
#define mccompcurname  betweenchop1and2
#define mccompcurtype  Guide_channeled
#define mccompcurindex 8
#define w1c mccbetweenchop1and2_w1c
#define w2c mccbetweenchop1and2_w2c
#define ww mccbetweenchop1and2_ww
#define hh mccbetweenchop1and2_hh
#define whalf mccbetweenchop1and2_whalf
#define hhalf mccbetweenchop1and2_hhalf
#define lwhalf mccbetweenchop1and2_lwhalf
#define lhhalf mccbetweenchop1and2_lhhalf
{   /* Declarations of betweenchop1and2=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccbetweenchop1and2_w1;
MCNUM h1 = mccbetweenchop1and2_h1;
MCNUM w2 = mccbetweenchop1and2_w2;
MCNUM h2 = mccbetweenchop1and2_h2;
MCNUM l = mccbetweenchop1and2_l;
MCNUM R0 = mccbetweenchop1and2_R0;
MCNUM Qc = mccbetweenchop1and2_Qc;
MCNUM alpha = mccbetweenchop1and2_alpha;
MCNUM m = mccbetweenchop1and2_m;
MCNUM nslit = mccbetweenchop1and2_nslit;
MCNUM d = mccbetweenchop1and2_d;
MCNUM Qcx = mccbetweenchop1and2_Qcx;
MCNUM Qcy = mccbetweenchop1and2_Qcy;
MCNUM alphax = mccbetweenchop1and2_alphax;
MCNUM alphay = mccbetweenchop1and2_alphay;
MCNUM W = mccbetweenchop1and2_W;
MCNUM mx = mccbetweenchop1and2_mx;
MCNUM my = mccbetweenchop1and2_my;
MCNUM nu = mccbetweenchop1and2_nu;
MCNUM phase = mccbetweenchop1and2_phase;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double angle=0;

  if (nu != 0 || phase != 0) { /* rotate neutron w/r to guide element */
    /* approximation of rotating straight Fermi Chopper */
    Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
    Rotation R;
    double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
    angle=fmod(360*nu*(t+dt)+phase, 360); /* in deg */
    /* modify angle so that Z0 guide side is always in front of incoming neutron */
    if (angle > 90 && angle < 270) { angle -= 180; }
    angle *= DEG2RAD;
    rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
    /* apply rotation to centered coordinates */
    Coords   RX = rot_apply(R, X);
    coords_get(RX, &x, &y, &z);
    z = z+l/2;
    /* rotate speed */
    X  = coords_set(vx,vy,vz);
    RX = rot_apply(R, X);
    coords_get(RX, &vx, &vy, &vz);
  }

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= w1/-2.0 || x >= w1/2.0 || y <= -hhalf || y >= hhalf)
    ABSORB;
  /* Shift origin to center of channel hit (absorb if hit dividing walls) */
  x += w1/2.0;
  edge = floor(x/w1c)*w1c;
  if(x - edge > w1c - d)
  {
    x -= w1/2.0; /* Re-adjust origin */
    ABSORB;
  }
  x -= (edge + (w1c - d)/2.0);
  hadj = edge + (w1c - d)/2.0 - w1/2.0;
  for(;;)
  {
    /* Compute the dot products of v and n for the four mirrors. */
    av = l*vx; bv = ww*vz;
    ah = l*vy; bh = hh*vz;
    vdotn_v1 = bv + av;         /* Left vertical */
    vdotn_v2 = bv - av;         /* Right vertical */
    vdotn_h1 = bh + ah;         /* Lower horizontal */
    vdotn_h2 = bh - ah;         /* Upper horizontal */
    /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
    cv1 = -whalf*l - z*ww; cv2 = x*l;
    ch1 = -hhalf*l - z*hh; ch2 = y*l;
    /* Compute intersection times. */
    t1 = (l - z)/vz;
    i = 0;
    if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
    {
      t1 = t2;
      i = 1;
    }
    if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
    {
      t1 = t2;
      i = 2;
    }
    if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
    {
      t1 = t2;
      i = 3;
    }
    if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
    {
      t1 = t2;
      i = 4;
    }
    if(i == 0)
      break;                    /* Neutron left guide. */
    PROP_DT(t1);
    switch(i)
    {
      case 1:                   /* Left vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
        dd = 2*vdotn_v1/nlen2;
        vx = vx - dd*l;
        vz = vz - dd*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        dd = 2*vdotn_v2/nlen2;
        vx = vx + dd*l;
        vz = vz - dd*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        dd = 2*vdotn_h1/nlen2;
        vy = vy - dd*l;
        vz = vz - dd*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        dd = 2*vdotn_h2/nlen2;
        vy = vy + dd*l;
        vz = vz - dd*hh;
        break;
    }
    /* Now compute reflectivity. */
    if((i <= 2 && mx == 0) || (i > 2 && my == 0))
    {
      x += hadj; /* Re-adjust origin */
      ABSORB;
    } else {
      double ref=1;
      if (i <= 2)
      {
        double par[] = {R0, Qcx, alphax, mx, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      } else {
        double par[] = {R0, Qcy, alphay, my, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      }
    }
    x += hadj; SCATTER; x -= hadj;
  } /* end for */
  x += hadj; /* Re-adjust origin */
  if (nu != 0 || phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }
}
#line 14662 "let_new.c"
}   /* End of betweenchop1and2=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbetweenchop1and2:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(8,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component betweenchop2and3 [9] */
  mccoordschange(mcposrbetweenchop2and3, mcrotrbetweenchop2and3,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component betweenchop2and3 (without coords transformations) */
  mcJumpTrace_betweenchop2and3:
  SIG_MESSAGE("betweenchop2and3 (Trace)");
  mcDEBUG_COMP("betweenchop2and3")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompbetweenchop2and3
  STORE_NEUTRON(9,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[9]++;
  mcPCounter[9] += p;
  mcP2Counter[9] += p*p;
#define mccompcurname  betweenchop2and3
#define mccompcurtype  Guide_channeled
#define mccompcurindex 9
#define w1c mccbetweenchop2and3_w1c
#define w2c mccbetweenchop2and3_w2c
#define ww mccbetweenchop2and3_ww
#define hh mccbetweenchop2and3_hh
#define whalf mccbetweenchop2and3_whalf
#define hhalf mccbetweenchop2and3_hhalf
#define lwhalf mccbetweenchop2and3_lwhalf
#define lhhalf mccbetweenchop2and3_lhhalf
{   /* Declarations of betweenchop2and3=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccbetweenchop2and3_w1;
MCNUM h1 = mccbetweenchop2and3_h1;
MCNUM w2 = mccbetweenchop2and3_w2;
MCNUM h2 = mccbetweenchop2and3_h2;
MCNUM l = mccbetweenchop2and3_l;
MCNUM R0 = mccbetweenchop2and3_R0;
MCNUM Qc = mccbetweenchop2and3_Qc;
MCNUM alpha = mccbetweenchop2and3_alpha;
MCNUM m = mccbetweenchop2and3_m;
MCNUM nslit = mccbetweenchop2and3_nslit;
MCNUM d = mccbetweenchop2and3_d;
MCNUM Qcx = mccbetweenchop2and3_Qcx;
MCNUM Qcy = mccbetweenchop2and3_Qcy;
MCNUM alphax = mccbetweenchop2and3_alphax;
MCNUM alphay = mccbetweenchop2and3_alphay;
MCNUM W = mccbetweenchop2and3_W;
MCNUM mx = mccbetweenchop2and3_mx;
MCNUM my = mccbetweenchop2and3_my;
MCNUM nu = mccbetweenchop2and3_nu;
MCNUM phase = mccbetweenchop2and3_phase;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double angle=0;

  if (nu != 0 || phase != 0) { /* rotate neutron w/r to guide element */
    /* approximation of rotating straight Fermi Chopper */
    Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
    Rotation R;
    double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
    angle=fmod(360*nu*(t+dt)+phase, 360); /* in deg */
    /* modify angle so that Z0 guide side is always in front of incoming neutron */
    if (angle > 90 && angle < 270) { angle -= 180; }
    angle *= DEG2RAD;
    rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
    /* apply rotation to centered coordinates */
    Coords   RX = rot_apply(R, X);
    coords_get(RX, &x, &y, &z);
    z = z+l/2;
    /* rotate speed */
    X  = coords_set(vx,vy,vz);
    RX = rot_apply(R, X);
    coords_get(RX, &vx, &vy, &vz);
  }

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= w1/-2.0 || x >= w1/2.0 || y <= -hhalf || y >= hhalf)
    ABSORB;
  /* Shift origin to center of channel hit (absorb if hit dividing walls) */
  x += w1/2.0;
  edge = floor(x/w1c)*w1c;
  if(x - edge > w1c - d)
  {
    x -= w1/2.0; /* Re-adjust origin */
    ABSORB;
  }
  x -= (edge + (w1c - d)/2.0);
  hadj = edge + (w1c - d)/2.0 - w1/2.0;
  for(;;)
  {
    /* Compute the dot products of v and n for the four mirrors. */
    av = l*vx; bv = ww*vz;
    ah = l*vy; bh = hh*vz;
    vdotn_v1 = bv + av;         /* Left vertical */
    vdotn_v2 = bv - av;         /* Right vertical */
    vdotn_h1 = bh + ah;         /* Lower horizontal */
    vdotn_h2 = bh - ah;         /* Upper horizontal */
    /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
    cv1 = -whalf*l - z*ww; cv2 = x*l;
    ch1 = -hhalf*l - z*hh; ch2 = y*l;
    /* Compute intersection times. */
    t1 = (l - z)/vz;
    i = 0;
    if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
    {
      t1 = t2;
      i = 1;
    }
    if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
    {
      t1 = t2;
      i = 2;
    }
    if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
    {
      t1 = t2;
      i = 3;
    }
    if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
    {
      t1 = t2;
      i = 4;
    }
    if(i == 0)
      break;                    /* Neutron left guide. */
    PROP_DT(t1);
    switch(i)
    {
      case 1:                   /* Left vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
        dd = 2*vdotn_v1/nlen2;
        vx = vx - dd*l;
        vz = vz - dd*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        dd = 2*vdotn_v2/nlen2;
        vx = vx + dd*l;
        vz = vz - dd*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        dd = 2*vdotn_h1/nlen2;
        vy = vy - dd*l;
        vz = vz - dd*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        dd = 2*vdotn_h2/nlen2;
        vy = vy + dd*l;
        vz = vz - dd*hh;
        break;
    }
    /* Now compute reflectivity. */
    if((i <= 2 && mx == 0) || (i > 2 && my == 0))
    {
      x += hadj; /* Re-adjust origin */
      ABSORB;
    } else {
      double ref=1;
      if (i <= 2)
      {
        double par[] = {R0, Qcx, alphax, mx, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      } else {
        double par[] = {R0, Qcy, alphay, my, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      }
    }
    x += hadj; SCATTER; x -= hadj;
  } /* end for */
  x += hadj; /* Re-adjust origin */
  if (nu != 0 || phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }
}
#line 14969 "let_new.c"
}   /* End of betweenchop2and3=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbetweenchop2and3:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(9,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component PR [10] */
  mccoordschange(mcposrPR, mcrotrPR,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component PR (without coords transformations) */
  mcJumpTrace_PR:
  SIG_MESSAGE("PR (Trace)");
  mcDEBUG_COMP("PR")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompPR
  STORE_NEUTRON(10,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[10]++;
  mcPCounter[10] += p;
  mcP2Counter[10] += p*p;
#define mccompcurname  PR
#define mccompcurtype  DiskChopper
#define mccompcurindex 10
#define Tg mccPR_Tg
#define To mccPR_To
#define delta_y mccPR_delta_y
#define height mccPR_height
#define omega mccPR_omega
{   /* Declarations of PR=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccPR_theta_0;
MCNUM radius = mccPR_radius;
MCNUM yheight = mccPR_yheight;
MCNUM nu = mccPR_nu;
MCNUM nslit = mccPR_nslit;
MCNUM jitter = mccPR_jitter;
MCNUM delay = mccPR_delay;
MCNUM isfirst = mccPR_isfirst;
MCNUM n_pulse = mccPR_n_pulse;
MCNUM abs_out = mccPR_abs_out;
MCNUM phase = mccPR_phase;
MCNUM xwidth = mccPR_xwidth;
MCNUM verbose = mccPR_verbose;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{
    double toff;
    double yprime;
    PROP_Z0;
    yprime = y+delta_y;

    /* Is neutron outside the vertical slit range and should we absorb? */
    if (abs_out && (x*x+yprime*yprime)>radius*radius) {
      ABSORB;
    }
    /* Does neutron hit inner solid part of chopper in case of yheight!=radius? */
    if ((x*x+yprime*yprime)<(radius-height)*(radius-height)) {
      ABSORB;
    }


    if (isfirst)
      {
        /* all events are put in the transmitted time frame */
        t=atan2(x,yprime)/omega + To*randpm1()/2.0 + delay + (jitter ? jitter*randnorm():0) + (n_pulse > 1 ? floor(n_pulse*rand01())*Tg : 0);
        /* correction: chopper slits transmission opening/full disk */
        p *= nslit*theta_0/2.0/PI;
      }
    else
      {
        toff=fabs(t-atan2(x,yprime)/omega - delay - (jitter ? jitter*randnorm():0));

        /* does neutron hit outside slit? */
        if (fmod(toff+To/2.0,Tg)>To) ABSORB;
      }
    SCATTER;

}
#line 15135 "let_new.c"
}   /* End of PR=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompPR:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(10,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component betweenchop3and4 [11] */
  mccoordschange(mcposrbetweenchop3and4, mcrotrbetweenchop3and4,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component betweenchop3and4 (without coords transformations) */
  mcJumpTrace_betweenchop3and4:
  SIG_MESSAGE("betweenchop3and4 (Trace)");
  mcDEBUG_COMP("betweenchop3and4")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompbetweenchop3and4
  STORE_NEUTRON(11,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[11]++;
  mcPCounter[11] += p;
  mcP2Counter[11] += p*p;
#define mccompcurname  betweenchop3and4
#define mccompcurtype  Guide_channeled
#define mccompcurindex 11
#define w1c mccbetweenchop3and4_w1c
#define w2c mccbetweenchop3and4_w2c
#define ww mccbetweenchop3and4_ww
#define hh mccbetweenchop3and4_hh
#define whalf mccbetweenchop3and4_whalf
#define hhalf mccbetweenchop3and4_hhalf
#define lwhalf mccbetweenchop3and4_lwhalf
#define lhhalf mccbetweenchop3and4_lhhalf
{   /* Declarations of betweenchop3and4=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccbetweenchop3and4_w1;
MCNUM h1 = mccbetweenchop3and4_h1;
MCNUM w2 = mccbetweenchop3and4_w2;
MCNUM h2 = mccbetweenchop3and4_h2;
MCNUM l = mccbetweenchop3and4_l;
MCNUM R0 = mccbetweenchop3and4_R0;
MCNUM Qc = mccbetweenchop3and4_Qc;
MCNUM alpha = mccbetweenchop3and4_alpha;
MCNUM m = mccbetweenchop3and4_m;
MCNUM nslit = mccbetweenchop3and4_nslit;
MCNUM d = mccbetweenchop3and4_d;
MCNUM Qcx = mccbetweenchop3and4_Qcx;
MCNUM Qcy = mccbetweenchop3and4_Qcy;
MCNUM alphax = mccbetweenchop3and4_alphax;
MCNUM alphay = mccbetweenchop3and4_alphay;
MCNUM W = mccbetweenchop3and4_W;
MCNUM mx = mccbetweenchop3and4_mx;
MCNUM my = mccbetweenchop3and4_my;
MCNUM nu = mccbetweenchop3and4_nu;
MCNUM phase = mccbetweenchop3and4_phase;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double angle=0;

  if (nu != 0 || phase != 0) { /* rotate neutron w/r to guide element */
    /* approximation of rotating straight Fermi Chopper */
    Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
    Rotation R;
    double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
    angle=fmod(360*nu*(t+dt)+phase, 360); /* in deg */
    /* modify angle so that Z0 guide side is always in front of incoming neutron */
    if (angle > 90 && angle < 270) { angle -= 180; }
    angle *= DEG2RAD;
    rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
    /* apply rotation to centered coordinates */
    Coords   RX = rot_apply(R, X);
    coords_get(RX, &x, &y, &z);
    z = z+l/2;
    /* rotate speed */
    X  = coords_set(vx,vy,vz);
    RX = rot_apply(R, X);
    coords_get(RX, &vx, &vy, &vz);
  }

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= w1/-2.0 || x >= w1/2.0 || y <= -hhalf || y >= hhalf)
    ABSORB;
  /* Shift origin to center of channel hit (absorb if hit dividing walls) */
  x += w1/2.0;
  edge = floor(x/w1c)*w1c;
  if(x - edge > w1c - d)
  {
    x -= w1/2.0; /* Re-adjust origin */
    ABSORB;
  }
  x -= (edge + (w1c - d)/2.0);
  hadj = edge + (w1c - d)/2.0 - w1/2.0;
  for(;;)
  {
    /* Compute the dot products of v and n for the four mirrors. */
    av = l*vx; bv = ww*vz;
    ah = l*vy; bh = hh*vz;
    vdotn_v1 = bv + av;         /* Left vertical */
    vdotn_v2 = bv - av;         /* Right vertical */
    vdotn_h1 = bh + ah;         /* Lower horizontal */
    vdotn_h2 = bh - ah;         /* Upper horizontal */
    /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
    cv1 = -whalf*l - z*ww; cv2 = x*l;
    ch1 = -hhalf*l - z*hh; ch2 = y*l;
    /* Compute intersection times. */
    t1 = (l - z)/vz;
    i = 0;
    if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
    {
      t1 = t2;
      i = 1;
    }
    if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
    {
      t1 = t2;
      i = 2;
    }
    if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
    {
      t1 = t2;
      i = 3;
    }
    if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
    {
      t1 = t2;
      i = 4;
    }
    if(i == 0)
      break;                    /* Neutron left guide. */
    PROP_DT(t1);
    switch(i)
    {
      case 1:                   /* Left vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
        dd = 2*vdotn_v1/nlen2;
        vx = vx - dd*l;
        vz = vz - dd*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        dd = 2*vdotn_v2/nlen2;
        vx = vx + dd*l;
        vz = vz - dd*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        dd = 2*vdotn_h1/nlen2;
        vy = vy - dd*l;
        vz = vz - dd*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        dd = 2*vdotn_h2/nlen2;
        vy = vy + dd*l;
        vz = vz - dd*hh;
        break;
    }
    /* Now compute reflectivity. */
    if((i <= 2 && mx == 0) || (i > 2 && my == 0))
    {
      x += hadj; /* Re-adjust origin */
      ABSORB;
    } else {
      double ref=1;
      if (i <= 2)
      {
        double par[] = {R0, Qcx, alphax, mx, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      } else {
        double par[] = {R0, Qcy, alphay, my, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      }
    }
    x += hadj; SCATTER; x -= hadj;
  } /* end for */
  x += hadj; /* Re-adjust origin */
  if (nu != 0 || phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }
}
#line 15439 "let_new.c"
}   /* End of betweenchop3and4=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbetweenchop3and4:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(11,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component CR [12] */
  mccoordschange(mcposrCR, mcrotrCR,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component CR (without coords transformations) */
  mcJumpTrace_CR:
  SIG_MESSAGE("CR (Trace)");
  mcDEBUG_COMP("CR")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompCR
  STORE_NEUTRON(12,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[12]++;
  mcPCounter[12] += p;
  mcP2Counter[12] += p*p;
#define mccompcurname  CR
#define mccompcurtype  DiskChopper
#define mccompcurindex 12
#define Tg mccCR_Tg
#define To mccCR_To
#define delta_y mccCR_delta_y
#define height mccCR_height
#define omega mccCR_omega
{   /* Declarations of CR=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccCR_theta_0;
MCNUM radius = mccCR_radius;
MCNUM yheight = mccCR_yheight;
MCNUM nu = mccCR_nu;
MCNUM nslit = mccCR_nslit;
MCNUM jitter = mccCR_jitter;
MCNUM delay = mccCR_delay;
MCNUM isfirst = mccCR_isfirst;
MCNUM n_pulse = mccCR_n_pulse;
MCNUM abs_out = mccCR_abs_out;
MCNUM phase = mccCR_phase;
MCNUM xwidth = mccCR_xwidth;
MCNUM verbose = mccCR_verbose;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{
    double toff;
    double yprime;
    PROP_Z0;
    yprime = y+delta_y;

    /* Is neutron outside the vertical slit range and should we absorb? */
    if (abs_out && (x*x+yprime*yprime)>radius*radius) {
      ABSORB;
    }
    /* Does neutron hit inner solid part of chopper in case of yheight!=radius? */
    if ((x*x+yprime*yprime)<(radius-height)*(radius-height)) {
      ABSORB;
    }


    if (isfirst)
      {
        /* all events are put in the transmitted time frame */
        t=atan2(x,yprime)/omega + To*randpm1()/2.0 + delay + (jitter ? jitter*randnorm():0) + (n_pulse > 1 ? floor(n_pulse*rand01())*Tg : 0);
        /* correction: chopper slits transmission opening/full disk */
        p *= nslit*theta_0/2.0/PI;
      }
    else
      {
        toff=fabs(t-atan2(x,yprime)/omega - delay - (jitter ? jitter*randnorm():0));

        /* does neutron hit outside slit? */
        if (fmod(toff+To/2.0,Tg)>To) ABSORB;
      }
    SCATTER;

}
#line 15605 "let_new.c"
}   /* End of CR=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompCR:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(12,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component betweenchop4andatt [13] */
  mccoordschange(mcposrbetweenchop4andatt, mcrotrbetweenchop4andatt,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component betweenchop4andatt (without coords transformations) */
  mcJumpTrace_betweenchop4andatt:
  SIG_MESSAGE("betweenchop4andatt (Trace)");
  mcDEBUG_COMP("betweenchop4andatt")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompbetweenchop4andatt
  STORE_NEUTRON(13,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[13]++;
  mcPCounter[13] += p;
  mcP2Counter[13] += p*p;
#define mccompcurname  betweenchop4andatt
#define mccompcurtype  Guide_channeled
#define mccompcurindex 13
#define w1c mccbetweenchop4andatt_w1c
#define w2c mccbetweenchop4andatt_w2c
#define ww mccbetweenchop4andatt_ww
#define hh mccbetweenchop4andatt_hh
#define whalf mccbetweenchop4andatt_whalf
#define hhalf mccbetweenchop4andatt_hhalf
#define lwhalf mccbetweenchop4andatt_lwhalf
#define lhhalf mccbetweenchop4andatt_lhhalf
{   /* Declarations of betweenchop4andatt=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccbetweenchop4andatt_w1;
MCNUM h1 = mccbetweenchop4andatt_h1;
MCNUM w2 = mccbetweenchop4andatt_w2;
MCNUM h2 = mccbetweenchop4andatt_h2;
MCNUM l = mccbetweenchop4andatt_l;
MCNUM R0 = mccbetweenchop4andatt_R0;
MCNUM Qc = mccbetweenchop4andatt_Qc;
MCNUM alpha = mccbetweenchop4andatt_alpha;
MCNUM m = mccbetweenchop4andatt_m;
MCNUM nslit = mccbetweenchop4andatt_nslit;
MCNUM d = mccbetweenchop4andatt_d;
MCNUM Qcx = mccbetweenchop4andatt_Qcx;
MCNUM Qcy = mccbetweenchop4andatt_Qcy;
MCNUM alphax = mccbetweenchop4andatt_alphax;
MCNUM alphay = mccbetweenchop4andatt_alphay;
MCNUM W = mccbetweenchop4andatt_W;
MCNUM mx = mccbetweenchop4andatt_mx;
MCNUM my = mccbetweenchop4andatt_my;
MCNUM nu = mccbetweenchop4andatt_nu;
MCNUM phase = mccbetweenchop4andatt_phase;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double angle=0;

  if (nu != 0 || phase != 0) { /* rotate neutron w/r to guide element */
    /* approximation of rotating straight Fermi Chopper */
    Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
    Rotation R;
    double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
    angle=fmod(360*nu*(t+dt)+phase, 360); /* in deg */
    /* modify angle so that Z0 guide side is always in front of incoming neutron */
    if (angle > 90 && angle < 270) { angle -= 180; }
    angle *= DEG2RAD;
    rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
    /* apply rotation to centered coordinates */
    Coords   RX = rot_apply(R, X);
    coords_get(RX, &x, &y, &z);
    z = z+l/2;
    /* rotate speed */
    X  = coords_set(vx,vy,vz);
    RX = rot_apply(R, X);
    coords_get(RX, &vx, &vy, &vz);
  }

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= w1/-2.0 || x >= w1/2.0 || y <= -hhalf || y >= hhalf)
    ABSORB;
  /* Shift origin to center of channel hit (absorb if hit dividing walls) */
  x += w1/2.0;
  edge = floor(x/w1c)*w1c;
  if(x - edge > w1c - d)
  {
    x -= w1/2.0; /* Re-adjust origin */
    ABSORB;
  }
  x -= (edge + (w1c - d)/2.0);
  hadj = edge + (w1c - d)/2.0 - w1/2.0;
  for(;;)
  {
    /* Compute the dot products of v and n for the four mirrors. */
    av = l*vx; bv = ww*vz;
    ah = l*vy; bh = hh*vz;
    vdotn_v1 = bv + av;         /* Left vertical */
    vdotn_v2 = bv - av;         /* Right vertical */
    vdotn_h1 = bh + ah;         /* Lower horizontal */
    vdotn_h2 = bh - ah;         /* Upper horizontal */
    /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
    cv1 = -whalf*l - z*ww; cv2 = x*l;
    ch1 = -hhalf*l - z*hh; ch2 = y*l;
    /* Compute intersection times. */
    t1 = (l - z)/vz;
    i = 0;
    if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
    {
      t1 = t2;
      i = 1;
    }
    if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
    {
      t1 = t2;
      i = 2;
    }
    if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
    {
      t1 = t2;
      i = 3;
    }
    if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
    {
      t1 = t2;
      i = 4;
    }
    if(i == 0)
      break;                    /* Neutron left guide. */
    PROP_DT(t1);
    switch(i)
    {
      case 1:                   /* Left vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
        dd = 2*vdotn_v1/nlen2;
        vx = vx - dd*l;
        vz = vz - dd*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        dd = 2*vdotn_v2/nlen2;
        vx = vx + dd*l;
        vz = vz - dd*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        dd = 2*vdotn_h1/nlen2;
        vy = vy - dd*l;
        vz = vz - dd*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        dd = 2*vdotn_h2/nlen2;
        vy = vy + dd*l;
        vz = vz - dd*hh;
        break;
    }
    /* Now compute reflectivity. */
    if((i <= 2 && mx == 0) || (i > 2 && my == 0))
    {
      x += hadj; /* Re-adjust origin */
      ABSORB;
    } else {
      double ref=1;
      if (i <= 2)
      {
        double par[] = {R0, Qcx, alphax, mx, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      } else {
        double par[] = {R0, Qcy, alphay, my, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      }
    }
    x += hadj; SCATTER; x -= hadj;
  } /* end for */
  x += hadj; /* Re-adjust origin */
  if (nu != 0 || phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }
}
#line 15909 "let_new.c"
}   /* End of betweenchop4andatt=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompbetweenchop4andatt:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(13,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component guide_move [14] */
  mccoordschange(mcposrguide_move, mcrotrguide_move,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component guide_move (without coords transformations) */
  mcJumpTrace_guide_move:
  SIG_MESSAGE("guide_move (Trace)");
  mcDEBUG_COMP("guide_move")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompguide_move
  STORE_NEUTRON(14,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[14]++;
  mcPCounter[14] += p;
  mcP2Counter[14] += p*p;
#define mccompcurname  guide_move
#define mccompcurtype  Guide_channeled
#define mccompcurindex 14
#define w1c mccguide_move_w1c
#define w2c mccguide_move_w2c
#define ww mccguide_move_ww
#define hh mccguide_move_hh
#define whalf mccguide_move_whalf
#define hhalf mccguide_move_hhalf
#define lwhalf mccguide_move_lwhalf
#define lhhalf mccguide_move_lhhalf
{   /* Declarations of guide_move=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccguide_move_w1;
MCNUM h1 = mccguide_move_h1;
MCNUM w2 = mccguide_move_w2;
MCNUM h2 = mccguide_move_h2;
MCNUM l = mccguide_move_l;
MCNUM R0 = mccguide_move_R0;
MCNUM Qc = mccguide_move_Qc;
MCNUM alpha = mccguide_move_alpha;
MCNUM m = mccguide_move_m;
MCNUM nslit = mccguide_move_nslit;
MCNUM d = mccguide_move_d;
MCNUM Qcx = mccguide_move_Qcx;
MCNUM Qcy = mccguide_move_Qcy;
MCNUM alphax = mccguide_move_alphax;
MCNUM alphay = mccguide_move_alphay;
MCNUM W = mccguide_move_W;
MCNUM mx = mccguide_move_mx;
MCNUM my = mccguide_move_my;
MCNUM nu = mccguide_move_nu;
MCNUM phase = mccguide_move_phase;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double angle=0;

  if (nu != 0 || phase != 0) { /* rotate neutron w/r to guide element */
    /* approximation of rotating straight Fermi Chopper */
    Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
    Rotation R;
    double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
    angle=fmod(360*nu*(t+dt)+phase, 360); /* in deg */
    /* modify angle so that Z0 guide side is always in front of incoming neutron */
    if (angle > 90 && angle < 270) { angle -= 180; }
    angle *= DEG2RAD;
    rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
    /* apply rotation to centered coordinates */
    Coords   RX = rot_apply(R, X);
    coords_get(RX, &x, &y, &z);
    z = z+l/2;
    /* rotate speed */
    X  = coords_set(vx,vy,vz);
    RX = rot_apply(R, X);
    coords_get(RX, &vx, &vy, &vz);
  }

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= w1/-2.0 || x >= w1/2.0 || y <= -hhalf || y >= hhalf)
    ABSORB;
  /* Shift origin to center of channel hit (absorb if hit dividing walls) */
  x += w1/2.0;
  edge = floor(x/w1c)*w1c;
  if(x - edge > w1c - d)
  {
    x -= w1/2.0; /* Re-adjust origin */
    ABSORB;
  }
  x -= (edge + (w1c - d)/2.0);
  hadj = edge + (w1c - d)/2.0 - w1/2.0;
  for(;;)
  {
    /* Compute the dot products of v and n for the four mirrors. */
    av = l*vx; bv = ww*vz;
    ah = l*vy; bh = hh*vz;
    vdotn_v1 = bv + av;         /* Left vertical */
    vdotn_v2 = bv - av;         /* Right vertical */
    vdotn_h1 = bh + ah;         /* Lower horizontal */
    vdotn_h2 = bh - ah;         /* Upper horizontal */
    /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
    cv1 = -whalf*l - z*ww; cv2 = x*l;
    ch1 = -hhalf*l - z*hh; ch2 = y*l;
    /* Compute intersection times. */
    t1 = (l - z)/vz;
    i = 0;
    if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
    {
      t1 = t2;
      i = 1;
    }
    if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
    {
      t1 = t2;
      i = 2;
    }
    if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
    {
      t1 = t2;
      i = 3;
    }
    if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
    {
      t1 = t2;
      i = 4;
    }
    if(i == 0)
      break;                    /* Neutron left guide. */
    PROP_DT(t1);
    switch(i)
    {
      case 1:                   /* Left vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
        dd = 2*vdotn_v1/nlen2;
        vx = vx - dd*l;
        vz = vz - dd*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        dd = 2*vdotn_v2/nlen2;
        vx = vx + dd*l;
        vz = vz - dd*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        dd = 2*vdotn_h1/nlen2;
        vy = vy - dd*l;
        vz = vz - dd*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        dd = 2*vdotn_h2/nlen2;
        vy = vy + dd*l;
        vz = vz - dd*hh;
        break;
    }
    /* Now compute reflectivity. */
    if((i <= 2 && mx == 0) || (i > 2 && my == 0))
    {
      x += hadj; /* Re-adjust origin */
      ABSORB;
    } else {
      double ref=1;
      if (i <= 2)
      {
        double par[] = {R0, Qcx, alphax, mx, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      } else {
        double par[] = {R0, Qcy, alphay, my, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      }
    }
    x += hadj; SCATTER; x -= hadj;
  } /* end for */
  x += hadj; /* Re-adjust origin */
  if (nu != 0 || phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }
}
#line 16216 "let_new.c"
}   /* End of guide_move=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompguide_move:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(14,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component funnel [15] */
  mccoordschange(mcposrfunnel, mcrotrfunnel,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component funnel (without coords transformations) */
  mcJumpTrace_funnel:
  SIG_MESSAGE("funnel (Trace)");
  mcDEBUG_COMP("funnel")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompfunnel
  STORE_NEUTRON(15,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[15]++;
  mcPCounter[15] += p;
  mcP2Counter[15] += p*p;
#define mccompcurname  funnel
#define mccompcurtype  Guide_channeled
#define mccompcurindex 15
#define w1c mccfunnel_w1c
#define w2c mccfunnel_w2c
#define ww mccfunnel_ww
#define hh mccfunnel_hh
#define whalf mccfunnel_whalf
#define hhalf mccfunnel_hhalf
#define lwhalf mccfunnel_lwhalf
#define lhhalf mccfunnel_lhhalf
{   /* Declarations of funnel=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccfunnel_w1;
MCNUM h1 = mccfunnel_h1;
MCNUM w2 = mccfunnel_w2;
MCNUM h2 = mccfunnel_h2;
MCNUM l = mccfunnel_l;
MCNUM R0 = mccfunnel_R0;
MCNUM Qc = mccfunnel_Qc;
MCNUM alpha = mccfunnel_alpha;
MCNUM m = mccfunnel_m;
MCNUM nslit = mccfunnel_nslit;
MCNUM d = mccfunnel_d;
MCNUM Qcx = mccfunnel_Qcx;
MCNUM Qcy = mccfunnel_Qcy;
MCNUM alphax = mccfunnel_alphax;
MCNUM alphay = mccfunnel_alphay;
MCNUM W = mccfunnel_W;
MCNUM mx = mccfunnel_mx;
MCNUM my = mccfunnel_my;
MCNUM nu = mccfunnel_nu;
MCNUM phase = mccfunnel_phase;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double angle=0;

  if (nu != 0 || phase != 0) { /* rotate neutron w/r to guide element */
    /* approximation of rotating straight Fermi Chopper */
    Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
    Rotation R;
    double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
    angle=fmod(360*nu*(t+dt)+phase, 360); /* in deg */
    /* modify angle so that Z0 guide side is always in front of incoming neutron */
    if (angle > 90 && angle < 270) { angle -= 180; }
    angle *= DEG2RAD;
    rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
    /* apply rotation to centered coordinates */
    Coords   RX = rot_apply(R, X);
    coords_get(RX, &x, &y, &z);
    z = z+l/2;
    /* rotate speed */
    X  = coords_set(vx,vy,vz);
    RX = rot_apply(R, X);
    coords_get(RX, &vx, &vy, &vz);
  }

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= w1/-2.0 || x >= w1/2.0 || y <= -hhalf || y >= hhalf)
    ABSORB;
  /* Shift origin to center of channel hit (absorb if hit dividing walls) */
  x += w1/2.0;
  edge = floor(x/w1c)*w1c;
  if(x - edge > w1c - d)
  {
    x -= w1/2.0; /* Re-adjust origin */
    ABSORB;
  }
  x -= (edge + (w1c - d)/2.0);
  hadj = edge + (w1c - d)/2.0 - w1/2.0;
  for(;;)
  {
    /* Compute the dot products of v and n for the four mirrors. */
    av = l*vx; bv = ww*vz;
    ah = l*vy; bh = hh*vz;
    vdotn_v1 = bv + av;         /* Left vertical */
    vdotn_v2 = bv - av;         /* Right vertical */
    vdotn_h1 = bh + ah;         /* Lower horizontal */
    vdotn_h2 = bh - ah;         /* Upper horizontal */
    /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
    cv1 = -whalf*l - z*ww; cv2 = x*l;
    ch1 = -hhalf*l - z*hh; ch2 = y*l;
    /* Compute intersection times. */
    t1 = (l - z)/vz;
    i = 0;
    if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
    {
      t1 = t2;
      i = 1;
    }
    if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
    {
      t1 = t2;
      i = 2;
    }
    if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
    {
      t1 = t2;
      i = 3;
    }
    if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
    {
      t1 = t2;
      i = 4;
    }
    if(i == 0)
      break;                    /* Neutron left guide. */
    PROP_DT(t1);
    switch(i)
    {
      case 1:                   /* Left vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
        dd = 2*vdotn_v1/nlen2;
        vx = vx - dd*l;
        vz = vz - dd*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        dd = 2*vdotn_v2/nlen2;
        vx = vx + dd*l;
        vz = vz - dd*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        dd = 2*vdotn_h1/nlen2;
        vy = vy - dd*l;
        vz = vz - dd*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        dd = 2*vdotn_h2/nlen2;
        vy = vy + dd*l;
        vz = vz - dd*hh;
        break;
    }
    /* Now compute reflectivity. */
    if((i <= 2 && mx == 0) || (i > 2 && my == 0))
    {
      x += hadj; /* Re-adjust origin */
      ABSORB;
    } else {
      double ref=1;
      if (i <= 2)
      {
        double par[] = {R0, Qcx, alphax, mx, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      } else {
        double par[] = {R0, Qcy, alphay, my, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      }
    }
    x += hadj; SCATTER; x -= hadj;
  } /* end for */
  x += hadj; /* Re-adjust origin */
  if (nu != 0 || phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }
}
#line 16523 "let_new.c"
}   /* End of funnel=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompfunnel:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(15,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component C5 [16] */
  mccoordschange(mcposrC5, mcrotrC5,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component C5 (without coords transformations) */
  mcJumpTrace_C5:
  SIG_MESSAGE("C5 (Trace)");
  mcDEBUG_COMP("C5")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompC5
  STORE_NEUTRON(16,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[16]++;
  mcPCounter[16] += p;
  mcP2Counter[16] += p*p;
#define mccompcurname  C5
#define mccompcurtype  DiskChopper
#define mccompcurindex 16
#define Tg mccC5_Tg
#define To mccC5_To
#define delta_y mccC5_delta_y
#define height mccC5_height
#define omega mccC5_omega
{   /* Declarations of C5=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccC5_theta_0;
MCNUM radius = mccC5_radius;
MCNUM yheight = mccC5_yheight;
MCNUM nu = mccC5_nu;
MCNUM nslit = mccC5_nslit;
MCNUM jitter = mccC5_jitter;
MCNUM delay = mccC5_delay;
MCNUM isfirst = mccC5_isfirst;
MCNUM n_pulse = mccC5_n_pulse;
MCNUM abs_out = mccC5_abs_out;
MCNUM phase = mccC5_phase;
MCNUM xwidth = mccC5_xwidth;
MCNUM verbose = mccC5_verbose;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{
    double toff;
    double yprime;
    PROP_Z0;
    yprime = y+delta_y;

    /* Is neutron outside the vertical slit range and should we absorb? */
    if (abs_out && (x*x+yprime*yprime)>radius*radius) {
      ABSORB;
    }
    /* Does neutron hit inner solid part of chopper in case of yheight!=radius? */
    if ((x*x+yprime*yprime)<(radius-height)*(radius-height)) {
      ABSORB;
    }


    if (isfirst)
      {
        /* all events are put in the transmitted time frame */
        t=atan2(x,yprime)/omega + To*randpm1()/2.0 + delay + (jitter ? jitter*randnorm():0) + (n_pulse > 1 ? floor(n_pulse*rand01())*Tg : 0);
        /* correction: chopper slits transmission opening/full disk */
        p *= nslit*theta_0/2.0/PI;
      }
    else
      {
        toff=fabs(t-atan2(x,yprime)/omega - delay - (jitter ? jitter*randnorm():0));

        /* does neutron hit outside slit? */
        if (fmod(toff+To/2.0,Tg)>To) ABSORB;
      }
    SCATTER;

}
#line 16689 "let_new.c"
}   /* End of C5=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompC5:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(16,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component C5_counter [17] */
  mccoordschange(mcposrC5_counter, mcrotrC5_counter,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component C5_counter (without coords transformations) */
  mcJumpTrace_C5_counter:
  SIG_MESSAGE("C5_counter (Trace)");
  mcDEBUG_COMP("C5_counter")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompC5_counter
  STORE_NEUTRON(17,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[17]++;
  mcPCounter[17] += p;
  mcP2Counter[17] += p*p;
#define mccompcurname  C5_counter
#define mccompcurtype  DiskChopper
#define mccompcurindex 17
#define Tg mccC5_counter_Tg
#define To mccC5_counter_To
#define delta_y mccC5_counter_delta_y
#define height mccC5_counter_height
#define omega mccC5_counter_omega
{   /* Declarations of C5_counter=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccC5_counter_theta_0;
MCNUM radius = mccC5_counter_radius;
MCNUM yheight = mccC5_counter_yheight;
MCNUM nu = mccC5_counter_nu;
MCNUM nslit = mccC5_counter_nslit;
MCNUM jitter = mccC5_counter_jitter;
MCNUM delay = mccC5_counter_delay;
MCNUM isfirst = mccC5_counter_isfirst;
MCNUM n_pulse = mccC5_counter_n_pulse;
MCNUM abs_out = mccC5_counter_abs_out;
MCNUM phase = mccC5_counter_phase;
MCNUM xwidth = mccC5_counter_xwidth;
MCNUM verbose = mccC5_counter_verbose;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{
    double toff;
    double yprime;
    PROP_Z0;
    yprime = y+delta_y;

    /* Is neutron outside the vertical slit range and should we absorb? */
    if (abs_out && (x*x+yprime*yprime)>radius*radius) {
      ABSORB;
    }
    /* Does neutron hit inner solid part of chopper in case of yheight!=radius? */
    if ((x*x+yprime*yprime)<(radius-height)*(radius-height)) {
      ABSORB;
    }


    if (isfirst)
      {
        /* all events are put in the transmitted time frame */
        t=atan2(x,yprime)/omega + To*randpm1()/2.0 + delay + (jitter ? jitter*randnorm():0) + (n_pulse > 1 ? floor(n_pulse*rand01())*Tg : 0);
        /* correction: chopper slits transmission opening/full disk */
        p *= nslit*theta_0/2.0/PI;
      }
    else
      {
        toff=fabs(t-atan2(x,yprime)/omega - delay - (jitter ? jitter*randnorm():0));

        /* does neutron hit outside slit? */
        if (fmod(toff+To/2.0,Tg)>To) ABSORB;
      }
    SCATTER;

}
#line 16852 "let_new.c"
}   /* End of C5_counter=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompC5_counter:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(17,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component endguide [18] */
  mccoordschange(mcposrendguide, mcrotrendguide,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component endguide (without coords transformations) */
  mcJumpTrace_endguide:
  SIG_MESSAGE("endguide (Trace)");
  mcDEBUG_COMP("endguide")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompendguide
  STORE_NEUTRON(18,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[18]++;
  mcPCounter[18] += p;
  mcP2Counter[18] += p*p;
#define mccompcurname  endguide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 18
#define w1c mccendguide_w1c
#define w2c mccendguide_w2c
#define ww mccendguide_ww
#define hh mccendguide_hh
#define whalf mccendguide_whalf
#define hhalf mccendguide_hhalf
#define lwhalf mccendguide_lwhalf
#define lhhalf mccendguide_lhhalf
{   /* Declarations of endguide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccendguide_w1;
MCNUM h1 = mccendguide_h1;
MCNUM w2 = mccendguide_w2;
MCNUM h2 = mccendguide_h2;
MCNUM l = mccendguide_l;
MCNUM R0 = mccendguide_R0;
MCNUM Qc = mccendguide_Qc;
MCNUM alpha = mccendguide_alpha;
MCNUM m = mccendguide_m;
MCNUM nslit = mccendguide_nslit;
MCNUM d = mccendguide_d;
MCNUM Qcx = mccendguide_Qcx;
MCNUM Qcy = mccendguide_Qcy;
MCNUM alphax = mccendguide_alphax;
MCNUM alphay = mccendguide_alphay;
MCNUM W = mccendguide_W;
MCNUM mx = mccendguide_mx;
MCNUM my = mccendguide_my;
MCNUM nu = mccendguide_nu;
MCNUM phase = mccendguide_phase;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double angle=0;

  if (nu != 0 || phase != 0) { /* rotate neutron w/r to guide element */
    /* approximation of rotating straight Fermi Chopper */
    Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
    Rotation R;
    double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
    angle=fmod(360*nu*(t+dt)+phase, 360); /* in deg */
    /* modify angle so that Z0 guide side is always in front of incoming neutron */
    if (angle > 90 && angle < 270) { angle -= 180; }
    angle *= DEG2RAD;
    rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
    /* apply rotation to centered coordinates */
    Coords   RX = rot_apply(R, X);
    coords_get(RX, &x, &y, &z);
    z = z+l/2;
    /* rotate speed */
    X  = coords_set(vx,vy,vz);
    RX = rot_apply(R, X);
    coords_get(RX, &vx, &vy, &vz);
  }

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= w1/-2.0 || x >= w1/2.0 || y <= -hhalf || y >= hhalf)
    ABSORB;
  /* Shift origin to center of channel hit (absorb if hit dividing walls) */
  x += w1/2.0;
  edge = floor(x/w1c)*w1c;
  if(x - edge > w1c - d)
  {
    x -= w1/2.0; /* Re-adjust origin */
    ABSORB;
  }
  x -= (edge + (w1c - d)/2.0);
  hadj = edge + (w1c - d)/2.0 - w1/2.0;
  for(;;)
  {
    /* Compute the dot products of v and n for the four mirrors. */
    av = l*vx; bv = ww*vz;
    ah = l*vy; bh = hh*vz;
    vdotn_v1 = bv + av;         /* Left vertical */
    vdotn_v2 = bv - av;         /* Right vertical */
    vdotn_h1 = bh + ah;         /* Lower horizontal */
    vdotn_h2 = bh - ah;         /* Upper horizontal */
    /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
    cv1 = -whalf*l - z*ww; cv2 = x*l;
    ch1 = -hhalf*l - z*hh; ch2 = y*l;
    /* Compute intersection times. */
    t1 = (l - z)/vz;
    i = 0;
    if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
    {
      t1 = t2;
      i = 1;
    }
    if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
    {
      t1 = t2;
      i = 2;
    }
    if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
    {
      t1 = t2;
      i = 3;
    }
    if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
    {
      t1 = t2;
      i = 4;
    }
    if(i == 0)
      break;                    /* Neutron left guide. */
    PROP_DT(t1);
    switch(i)
    {
      case 1:                   /* Left vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
        dd = 2*vdotn_v1/nlen2;
        vx = vx - dd*l;
        vz = vz - dd*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        dd = 2*vdotn_v2/nlen2;
        vx = vx + dd*l;
        vz = vz - dd*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        dd = 2*vdotn_h1/nlen2;
        vy = vy - dd*l;
        vz = vz - dd*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        dd = 2*vdotn_h2/nlen2;
        vy = vy + dd*l;
        vz = vz - dd*hh;
        break;
    }
    /* Now compute reflectivity. */
    if((i <= 2 && mx == 0) || (i > 2 && my == 0))
    {
      x += hadj; /* Re-adjust origin */
      ABSORB;
    } else {
      double ref=1;
      if (i <= 2)
      {
        double par[] = {R0, Qcx, alphax, mx, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      } else {
        double par[] = {R0, Qcy, alphay, my, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      }
    }
    x += hadj; SCATTER; x -= hadj;
  } /* end for */
  x += hadj; /* Re-adjust origin */
  if (nu != 0 || phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }
}
#line 17156 "let_new.c"
}   /* End of endguide=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompendguide:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(18,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component snot [19] */
  mccoordschange(mcposrsnot, mcrotrsnot,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component snot (without coords transformations) */
  mcJumpTrace_snot:
  SIG_MESSAGE("snot (Trace)");
  mcDEBUG_COMP("snot")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompsnot
  STORE_NEUTRON(19,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[19]++;
  mcPCounter[19] += p;
  mcP2Counter[19] += p*p;
#define mccompcurname  snot
#define mccompcurtype  Guide_channeled
#define mccompcurindex 19
#define w1c mccsnot_w1c
#define w2c mccsnot_w2c
#define ww mccsnot_ww
#define hh mccsnot_hh
#define whalf mccsnot_whalf
#define hhalf mccsnot_hhalf
#define lwhalf mccsnot_lwhalf
#define lhhalf mccsnot_lhhalf
{   /* Declarations of snot=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccsnot_w1;
MCNUM h1 = mccsnot_h1;
MCNUM w2 = mccsnot_w2;
MCNUM h2 = mccsnot_h2;
MCNUM l = mccsnot_l;
MCNUM R0 = mccsnot_R0;
MCNUM Qc = mccsnot_Qc;
MCNUM alpha = mccsnot_alpha;
MCNUM m = mccsnot_m;
MCNUM nslit = mccsnot_nslit;
MCNUM d = mccsnot_d;
MCNUM Qcx = mccsnot_Qcx;
MCNUM Qcy = mccsnot_Qcy;
MCNUM alphax = mccsnot_alphax;
MCNUM alphay = mccsnot_alphay;
MCNUM W = mccsnot_W;
MCNUM mx = mccsnot_mx;
MCNUM my = mccsnot_my;
MCNUM nu = mccsnot_nu;
MCNUM phase = mccsnot_phase;
#line 130 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  double t1,t2;                                 /* Intersection times. */
  double av,ah,bv,bh,cv1,cv2,ch1,ch2,dd;        /* Intermediate values */
  double vdotn_v1,vdotn_v2,vdotn_h1,vdotn_h2;   /* Dot products. */
  int i;                                        /* Which mirror hit? */
  double q;                                     /* Q [1/AA] of reflection */
  double nlen2;                                 /* Vector lengths squared */
  double edge;
  double hadj;                                  /* Channel displacement */
  double angle=0;

  if (nu != 0 || phase != 0) { /* rotate neutron w/r to guide element */
    /* approximation of rotating straight Fermi Chopper */
    Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
    Rotation R;
    double dt=(-z+l/2)/vz; /* time shift to each center of slit package */
    angle=fmod(360*nu*(t+dt)+phase, 360); /* in deg */
    /* modify angle so that Z0 guide side is always in front of incoming neutron */
    if (angle > 90 && angle < 270) { angle -= 180; }
    angle *= DEG2RAD;
    rot_set_rotation(R, 0, -angle, 0); /* will rotate neutron instead of comp: negative side */
    /* apply rotation to centered coordinates */
    Coords   RX = rot_apply(R, X);
    coords_get(RX, &x, &y, &z);
    z = z+l/2;
    /* rotate speed */
    X  = coords_set(vx,vy,vz);
    RX = rot_apply(R, X);
    coords_get(RX, &vx, &vy, &vz);
  }

  /* Propagate neutron to guide entrance. */
  PROP_Z0;
  /* Scatter here to ensure that fully transmitted neutrons will not be
     absorbed in a GROUP construction, e.g. all neutrons - even the
     later absorbed ones are scattered at the guide entry. */
  SCATTER;
  if(x <= w1/-2.0 || x >= w1/2.0 || y <= -hhalf || y >= hhalf)
    ABSORB;
  /* Shift origin to center of channel hit (absorb if hit dividing walls) */
  x += w1/2.0;
  edge = floor(x/w1c)*w1c;
  if(x - edge > w1c - d)
  {
    x -= w1/2.0; /* Re-adjust origin */
    ABSORB;
  }
  x -= (edge + (w1c - d)/2.0);
  hadj = edge + (w1c - d)/2.0 - w1/2.0;
  for(;;)
  {
    /* Compute the dot products of v and n for the four mirrors. */
    av = l*vx; bv = ww*vz;
    ah = l*vy; bh = hh*vz;
    vdotn_v1 = bv + av;         /* Left vertical */
    vdotn_v2 = bv - av;         /* Right vertical */
    vdotn_h1 = bh + ah;         /* Lower horizontal */
    vdotn_h2 = bh - ah;         /* Upper horizontal */
    /* Compute the dot products of (O - r) and n as c1+c2 and c1-c2 */
    cv1 = -whalf*l - z*ww; cv2 = x*l;
    ch1 = -hhalf*l - z*hh; ch2 = y*l;
    /* Compute intersection times. */
    t1 = (l - z)/vz;
    i = 0;
    if(vdotn_v1 < 0 && (t2 = (cv1 - cv2)/vdotn_v1) < t1)
    {
      t1 = t2;
      i = 1;
    }
    if(vdotn_v2 < 0 && (t2 = (cv1 + cv2)/vdotn_v2) < t1)
    {
      t1 = t2;
      i = 2;
    }
    if(vdotn_h1 < 0 && (t2 = (ch1 - ch2)/vdotn_h1) < t1)
    {
      t1 = t2;
      i = 3;
    }
    if(vdotn_h2 < 0 && (t2 = (ch1 + ch2)/vdotn_h2) < t1)
    {
      t1 = t2;
      i = 4;
    }
    if(i == 0)
      break;                    /* Neutron left guide. */
    PROP_DT(t1);
    switch(i)
    {
      case 1:                   /* Left vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v1/sqrt(nlen2);
        dd = 2*vdotn_v1/nlen2;
        vx = vx - dd*l;
        vz = vz - dd*ww;
        break;
      case 2:                   /* Right vertical mirror */
        nlen2 = l*l + ww*ww;
        q = V2Q*(-2)*vdotn_v2/sqrt(nlen2);
        dd = 2*vdotn_v2/nlen2;
        vx = vx + dd*l;
        vz = vz - dd*ww;
        break;
      case 3:                   /* Lower horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h1/sqrt(nlen2);
        dd = 2*vdotn_h1/nlen2;
        vy = vy - dd*l;
        vz = vz - dd*hh;
        break;
      case 4:                   /* Upper horizontal mirror */
        nlen2 = l*l + hh*hh;
        q = V2Q*(-2)*vdotn_h2/sqrt(nlen2);
        dd = 2*vdotn_h2/nlen2;
        vy = vy + dd*l;
        vz = vz - dd*hh;
        break;
    }
    /* Now compute reflectivity. */
    if((i <= 2 && mx == 0) || (i > 2 && my == 0))
    {
      x += hadj; /* Re-adjust origin */
      ABSORB;
    } else {
      double ref=1;
      if (i <= 2)
      {
        double par[] = {R0, Qcx, alphax, mx, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      } else {
        double par[] = {R0, Qcy, alphay, my, W};
        StdReflecFunc(q, par, &ref);
        if (ref > 0)
          p *= ref;
        else {
          x += hadj; /* Re-adjust origin */
          ABSORB;                               /* Cutoff ~ 1E-10 */
        }
      }
    }
    x += hadj; SCATTER; x -= hadj;
  } /* end for */
  x += hadj; /* Re-adjust origin */
  if (nu != 0 || phase != 0) { /* rotate back neutron w/r to guide element */
      /* approximation of rotating straight Fermi Chopper */
      Coords   X = coords_set(x,y,z-l/2);  /* current coordinates of neutron in centered static frame */
      Rotation R;
      rot_set_rotation(R, 0, angle, 0); /* will rotate back neutron: positive side */
      /* apply rotation to centered coordinates */
      Coords   RX = rot_apply(R, X);
      coords_get(RX, &x, &y, &z);
      z = z+l/2;
      /* rotate speed */
      X  = coords_set(vx,vy,vz);
      RX = rot_apply(R, X);
      coords_get(RX, &vx, &vy, &vz);
    }
}
#line 17463 "let_new.c"
}   /* End of snot=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsnot:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(19,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component detslit [20] */
  mccoordschange(mcposrdetslit, mcrotrdetslit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component detslit (without coords transformations) */
  mcJumpTrace_detslit:
  SIG_MESSAGE("detslit (Trace)");
  mcDEBUG_COMP("detslit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompdetslit
  STORE_NEUTRON(20,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[20]++;
  mcPCounter[20] += p;
  mcP2Counter[20] += p*p;
#define mccompcurname  detslit
#define mccompcurtype  Slit
#define mccompcurindex 20
{   /* Declarations of detslit=Slit() SETTING parameters. */
MCNUM xmin = mccdetslit_xmin;
MCNUM xmax = mccdetslit_xmax;
MCNUM ymin = mccdetslit_ymin;
MCNUM ymax = mccdetslit_ymax;
MCNUM radius = mccdetslit_radius;
MCNUM xwidth = mccdetslit_xwidth;
MCNUM yheight = mccdetslit_yheight;
#line 56 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
    || ((radius != 0) && (x*x + y*y > radius*radius)))
      ABSORB;
    else
        SCATTER;
}
#line 17593 "let_new.c"
}   /* End of detslit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompdetslit:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(20,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component mcplout [21] */
  mccoordschange(mcposrmcplout, mcrotrmcplout,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component mcplout (without coords transformations) */
  mcJumpTrace_mcplout:
  SIG_MESSAGE("mcplout (Trace)");
  mcDEBUG_COMP("mcplout")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompmcplout
  STORE_NEUTRON(21,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[21]++;
  mcPCounter[21] += p;
  mcP2Counter[21] += p*p;
#define mccompcurname  mcplout
#define mccompcurtype  MCPL_output_horace
#define mccompcurindex 21
#define polarisationuse mccmcplout_polarisationuse
#define doubleprec mccmcplout_doubleprec
#define verbose mccmcplout_verbose
#define radius mccmcplout_radius
#define thickness mccmcplout_thickness
#define yheight mccmcplout_yheight
#define restore_neutron mccmcplout_restore_neutron
#define userflag mccmcplout_userflag
#define plate mccmcplout_plate
#define outputfile mccmcplout_outputfile
#define particle mccmcplout_particle
#define Particle mccmcplout_Particle
#define userflagenabled mccmcplout_userflagenabled
{   /* Declarations of mcplout=MCPL_output_horace() SETTING parameters. */
char* filename = mccmcplout_filename;
char* userflagcomment = mccmcplout_userflagcomment;
MCNUM merge_mpi = mccmcplout_merge_mpi;
MCNUM keep_mpi_unmerged = mccmcplout_keep_mpi_unmerged;
#line 184 "MCPL_output_horace.comp"
{

    double nrm;
    // t0: reaches z=-a, t1/t2: enters/leaves shell first time, 
    // t3/t4: enters/leaves shell second time, t5: reaches z=+a
    double t0, t1, t2, t3, t4, t5, dt;
    double diameter = 2 * radius;
    double frac, delta;
    Rotation Rt;
    Coords r, v, s, rg, vg, sg;
    double rx, ry, rz, rvx, rvy, rvz, rsx, rsy, rsz;
    int keep_neutron = 0;

    // Determines the fractional time this trajectory spends within the annulus
    // Which will be used by the rejection sampling algorithm
    if (!plate) {
        if (box_intersect(&t0, &t5, x, y, z, vx, vy, vz, diameter, yheight, diameter)) {
            if (cylinder_intersect(&t1, &t4, x, y, z, vx, vy, vz, radius, yheight)) {
                cylinder_intersect(&t2, &t3, x, y, z, vx, vy, vz, radius-thickness, yheight);
                frac = ((t2 - t1) + (t4 - t3)) / (t5 - t0);
                delta = rand01();
                if (delta <= frac) {
                    keep_neutron = 1;
                    dt = delta * (t5 - t0) + t1;
                    if (t2 > 0 && dt >= t2)
                        dt += t3 - t2;
                    PROP_DT(dt);
                    r = coords_set(x,y,z);
                    v = coords_set(vx,vy,vz);
                    s = coords_set(sx,sy,sz);
                    rot_transpose(ROT_R_CURRENT_COMP, Rt);
                    rg = rot_apply(Rt,r);
                    vg = rot_apply(Rt,v);
                    sg = rot_apply(Rt,s);
                    coords_get(rg, &rx, &ry, &rz);
                    coords_get(vg, &rvx, &rvy, &rvz);
                    coords_get(sg, &rsx, &rsy, &rsz);
                }
            }
        }
    }
    else {
        PROP_Z0;
        keep_neutron = 1;
        rx = x; ry = y; rz = z;
        rsx = sx; rsy = sy; rsz = sz;
        rvx = vx; rvy = vy; rvz = vz;
    }

    if (keep_neutron) {
        /*positions are in cm*/
        particle->position[0] = rx*100;
        particle->position[1] = ry*100;
        particle->position[2] = rz*100;

        if(polarisationuse){
            particle->polarisation[0] = rsx;
            particle->polarisation[1] = rsy;
            particle->polarisation[2] = rsz;
        }

        nrm = sqrt(rvx*rvx + rvy*rvy + rvz*rvz);
        /*ekin is in MeV*/
        particle->ekin = VS2E * nrm * nrm / 1e9;
        particle->direction[0] = rvx/nrm;
        particle->direction[1] = rvy/nrm;
        particle->direction[2] = rvz/nrm;
        /*time in ms:*/
        particle->time = t*1e3;
        /*weight in unspecified units:*/
        particle->weight = p;
        /*if specified also add the userflags*/
        if(userflagenabled){
            particle->userflags = (uint32_t) userflag;
        }

#if defined (USE_MPI)
        MPI_MASTER(
#endif
            if (verbose==3 && mcrun_num<10) {
                printf("id=%ld\tpdg=2112\tekin=%g MeV\tx=%g cm\ty=%g cm\tz=%g cm\tux=%g\tuy=%g\tuz=%g\tt=%g ms\tweight=%g\tpolx=%g\tpoly=%g\tpolz=%g\n",
                mcrun_num, particle->ekin, particle->position[0], particle->position[1], particle->position[2],
	            particle->direction[0], particle->direction[1], particle->direction[2], particle->time, particle->weight,
       	        particle->polarisation[0], particle->polarisation[1], particle->polarisation[2]);
            }
#if defined (USE_MPI)
        );
#endif

        mcpl_add_particle(outputfile,particle);

        SCATTER;
        if (restore_neutron) {
            RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
        }
    }
}
#line 17814 "let_new.c"
}   /* End of mcplout=MCPL_output_horace() SETTING parameter declarations. */
#undef userflagenabled
#undef Particle
#undef particle
#undef outputfile
#undef plate
#undef userflag
#undef restore_neutron
#undef yheight
#undef thickness
#undef radius
#undef verbose
#undef doubleprec
#undef polarisationuse
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompmcplout:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(21,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  mcabsorbAll:
  mcDEBUG_LEAVE()
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)
  /* Copy neutron state to global variables. */
  mcnx = mcnlx;
  mcny = mcnly;
  mcnz = mcnlz;
  mcnvx = mcnlvx;
  mcnvy = mcnlvy;
  mcnvz = mcnlvz;
  mcnt = mcnlt;
  mcnsx = mcnlsx;
  mcnsy = mcnlsy;
  mcnsz = mcnlsz;
  mcnp = mcnlp;

} /* end trace */

void mcsave(FILE *handle) {
  if (!handle) mcsiminfo_init(NULL);
  /* User component SAVE code. */

  /* User SAVE code for component 'mcplout'. */
  SIG_MESSAGE("mcplout (Save)");
#define mccompcurname  mcplout
#define mccompcurtype  MCPL_output_horace
#define mccompcurindex 21
#define polarisationuse mccmcplout_polarisationuse
#define doubleprec mccmcplout_doubleprec
#define verbose mccmcplout_verbose
#define radius mccmcplout_radius
#define thickness mccmcplout_thickness
#define yheight mccmcplout_yheight
#define restore_neutron mccmcplout_restore_neutron
#define userflag mccmcplout_userflag
#define plate mccmcplout_plate
#define outputfile mccmcplout_outputfile
#define particle mccmcplout_particle
#define Particle mccmcplout_Particle
#define userflagenabled mccmcplout_userflagenabled
{   /* Declarations of mcplout=MCPL_output_horace() SETTING parameters. */
char* filename = mccmcplout_filename;
char* userflagcomment = mccmcplout_userflagcomment;
MCNUM merge_mpi = mccmcplout_merge_mpi;
MCNUM keep_mpi_unmerged = mccmcplout_keep_mpi_unmerged;
#line 283 "MCPL_output_horace.comp"
{
#ifdef USE_MPI
  if (merge_mpi && mpi_node_count > 1) {
    mcpl_close_outfile(outputfile);
  } else {
    mcpl_closeandgzip_outfile(outputfile);
  }
#else
  mcpl_closeandgzip_outfile(outputfile);
#endif
}
#line 17940 "let_new.c"
}   /* End of mcplout=MCPL_output_horace() SETTING parameter declarations. */
#undef userflagenabled
#undef Particle
#undef particle
#undef outputfile
#undef plate
#undef userflag
#undef restore_neutron
#undef yheight
#undef thickness
#undef radius
#undef verbose
#undef doubleprec
#undef polarisationuse
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  if (!handle) mcsiminfo_close(); 
} /* end save */
void mcfinally(void) {
  /* User component FINALLY code. */
  mcsiminfo_init(NULL);
  mcsave(mcsiminfo_file); /* save data when simulation ends */

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] Origin\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] Origin=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] Source\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] Source=ViewModISISver1()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] shutter\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] shutter=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] insert\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] insert=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] afterinsert\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] afterinsert=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] Res1\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] Res1=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] Res1_counter\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] Res1_counter=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] betweenchop1and2\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] betweenchop1and2=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] betweenchop2and3\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] betweenchop2and3=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] PR\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] PR=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
    if (!mcNCounter[11]) fprintf(stderr, "Warning: No neutron could reach Component[11] betweenchop3and4\n");
    if (mcAbsorbProp[11]) fprintf(stderr, "Warning: %g events were removed in Component[11] betweenchop3and4=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[11]);
    if (!mcNCounter[12]) fprintf(stderr, "Warning: No neutron could reach Component[12] CR\n");
    if (mcAbsorbProp[12]) fprintf(stderr, "Warning: %g events were removed in Component[12] CR=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[12]);
    if (!mcNCounter[13]) fprintf(stderr, "Warning: No neutron could reach Component[13] betweenchop4andatt\n");
    if (mcAbsorbProp[13]) fprintf(stderr, "Warning: %g events were removed in Component[13] betweenchop4andatt=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[13]);
    if (!mcNCounter[14]) fprintf(stderr, "Warning: No neutron could reach Component[14] guide_move\n");
    if (mcAbsorbProp[14]) fprintf(stderr, "Warning: %g events were removed in Component[14] guide_move=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[14]);
    if (!mcNCounter[15]) fprintf(stderr, "Warning: No neutron could reach Component[15] funnel\n");
    if (mcAbsorbProp[15]) fprintf(stderr, "Warning: %g events were removed in Component[15] funnel=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[15]);
    if (!mcNCounter[16]) fprintf(stderr, "Warning: No neutron could reach Component[16] C5\n");
    if (mcAbsorbProp[16]) fprintf(stderr, "Warning: %g events were removed in Component[16] C5=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[16]);
    if (!mcNCounter[17]) fprintf(stderr, "Warning: No neutron could reach Component[17] C5_counter\n");
    if (mcAbsorbProp[17]) fprintf(stderr, "Warning: %g events were removed in Component[17] C5_counter=DiskChopper()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[17]);
    if (!mcNCounter[18]) fprintf(stderr, "Warning: No neutron could reach Component[18] endguide\n");
    if (mcAbsorbProp[18]) fprintf(stderr, "Warning: %g events were removed in Component[18] endguide=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[18]);
    if (!mcNCounter[19]) fprintf(stderr, "Warning: No neutron could reach Component[19] snot\n");
    if (mcAbsorbProp[19]) fprintf(stderr, "Warning: %g events were removed in Component[19] snot=Guide_channeled()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[19]);
    if (!mcNCounter[20]) fprintf(stderr, "Warning: No neutron could reach Component[20] detslit\n");
    if (mcAbsorbProp[20]) fprintf(stderr, "Warning: %g events were removed in Component[20] detslit=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[20]);
  /* User FINALLY code for component 'mcplout'. */
  SIG_MESSAGE("mcplout (Finally)");
#define mccompcurname  mcplout
#define mccompcurtype  MCPL_output_horace
#define mccompcurindex 21
#define polarisationuse mccmcplout_polarisationuse
#define doubleprec mccmcplout_doubleprec
#define verbose mccmcplout_verbose
#define radius mccmcplout_radius
#define thickness mccmcplout_thickness
#define yheight mccmcplout_yheight
#define restore_neutron mccmcplout_restore_neutron
#define userflag mccmcplout_userflag
#define plate mccmcplout_plate
#define outputfile mccmcplout_outputfile
#define particle mccmcplout_particle
#define Particle mccmcplout_Particle
#define userflagenabled mccmcplout_userflagenabled
{   /* Declarations of mcplout=MCPL_output_horace() SETTING parameters. */
char* filename = mccmcplout_filename;
char* userflagcomment = mccmcplout_userflagcomment;
MCNUM merge_mpi = mccmcplout_merge_mpi;
MCNUM keep_mpi_unmerged = mccmcplout_keep_mpi_unmerged;
#line 296 "MCPL_output_horace.comp"
{
#ifdef USE_MPI
 MPI_Barrier(MPI_COMM_WORLD);
 MPI_MASTER(
     /* Only attempt merge if requested and meaningful */
     if (merge_mpi && mpi_node_count > 1) {
        char **mpi_node_files;
        char *merge_outfilename;
        char extension[128]="mcpl";
        int j;
        mcpl_outfile_t merge_outfile;

        merge_outfilename=mcfull_file(filename,extension);

        mpi_node_files=(char **) calloc(mpi_node_count,sizeof(char *));
        for (j=0;j<mpi_node_count;j++){
            sprintf(extension,"node_%i.mcpl",j);
            mpi_node_files[j]=mcfull_file(filename,extension);
        }
        /*now do the merge through the call to mcpl_merge_files*/
        merge_outfile = mcpl_merge_files(merge_outfilename,mpi_node_count,(const char **) mpi_node_files);
        mcpl_closeandgzip_outfile(merge_outfile);

        /*remove the original unmerged files if wanted*/
        if(!keep_mpi_unmerged){
            int status=0;
            for (j=0;j<mpi_node_count;j++){
                status+=remove(mpi_node_files[j]);
            }
            if (status){
                fprintf(stderr,"Warning (%s): Could not remove one or more unmerged files.\n",NAME_CURRENT_COMP);
            }
        }

        /*free the string storage*/
        free(merge_outfilename);
        for (j=0;j<mpi_node_count;j++){
            free(mpi_node_files[j]);
        }
        free(mpi_node_files);
    }
  );
#endif
}
#line 18074 "let_new.c"
}   /* End of mcplout=MCPL_output_horace() SETTING parameter declarations. */
#undef userflagenabled
#undef Particle
#undef particle
#undef outputfile
#undef plate
#undef userflag
#undef restore_neutron
#undef yheight
#undef thickness
#undef radius
#undef verbose
#undef doubleprec
#undef polarisationuse
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[21]) fprintf(stderr, "Warning: No neutron could reach Component[21] mcplout\n");
    if (mcAbsorbProp[21]) fprintf(stderr, "Warning: %g events were removed in Component[21] mcplout=MCPL_output_horace()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[21]);
  mcsiminfo_close(); 
} /* end finally */
#define magnify mcdis_magnify
#define line mcdis_line
#define dashed_line mcdis_dashed_line
#define multiline mcdis_multiline
#define rectangle mcdis_rectangle
#define box mcdis_box
#define circle mcdis_circle
void mcdisplay(void) {
  printf("MCDISPLAY: start\n");
  /* Components MCDISPLAY code. */

  /* MCDISPLAY code for component 'Origin'. */
  SIG_MESSAGE("Origin (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Origin");
#define mccompcurname  Origin
#define mccompcurtype  Arm
#define mccompcurindex 1
#line 40 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  magnify("");
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 18122 "let_new.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Source'. */
  SIG_MESSAGE("Source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Source");
#define mccompcurname  Source
#define mccompcurtype  ViewModISISver1
#define mccompcurindex 2
{   /* Declarations of Source=ViewModISISver1() SETTING parameters. */
char* Face = mccSource_Face;
MCNUM E0 = mccSource_E0;
MCNUM E1 = mccSource_E1;
MCNUM tally = mccSource_tally;
MCNUM modPosition = mccSource_modPosition;
MCNUM modXsize = mccSource_modXsize;
MCNUM modZsize = mccSource_modZsize;
MCNUM xw = mccSource_xw;
MCNUM yh = mccSource_yh;
MCNUM dist = mccSource_dist;
#line 1125 "ViewModISISver1.comp"
{
  double cirp=0.0,cirq=0.3,pi=3.141592654;
  int pp=0; /* circle drawing parameter*/



  magnify("xy");
  multiline(5,-0.5*TS.XAxis,-0.5*TS.ZAxis,0.0,
	    0.5*TS.XAxis,-0.5*TS.ZAxis,0.0,
	    0.5*TS.XAxis,0.5*TS.ZAxis,0.0,
	    -0.5*TS.XAxis,0.5*TS.ZAxis,0.0,
	    -0.5*TS.XAxis,-0.5*TS.ZAxis,0.0);
  /* circle("xy",0.0,0.0,0.0,cos(cirp)); */

  /*line(0.5*sin(cirp),0.0,0.5*cos(cirp),0.5*sin(cirq),0.0,0.5*cos(cirq));*/

  /*line(-0.5,0.0,0.0,0.0,0.0,0.5);*/

  for (pp=0;pp<=20;pp=pp+2)
    {
      cirp= (pp*(pi/21.0))-(0.5*pi);
      cirq= ((pp+1)*(pi/21.0))-(0.5*pi);
      line(0.5*sin(cirp),0.0,0.5*cos(cirp),0.5*sin(cirq),0.0,0.5*cos(cirq));
    }

}
#line 18171 "let_new.c"
}   /* End of Source=ViewModISISver1() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'shutter'. */
  SIG_MESSAGE("shutter (McDisplay)");
  printf("MCDISPLAY: component %s\n", "shutter");
#define mccompcurname  shutter
#define mccompcurtype  Guide_channeled
#define mccompcurindex 3
#define w1c mccshutter_w1c
#define w2c mccshutter_w2c
#define ww mccshutter_ww
#define hh mccshutter_hh
#define whalf mccshutter_whalf
#define hhalf mccshutter_hhalf
#define lwhalf mccshutter_lwhalf
#define lhhalf mccshutter_lhhalf
{   /* Declarations of shutter=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccshutter_w1;
MCNUM h1 = mccshutter_h1;
MCNUM w2 = mccshutter_w2;
MCNUM h2 = mccshutter_h2;
MCNUM l = mccshutter_l;
MCNUM R0 = mccshutter_R0;
MCNUM Qc = mccshutter_Qc;
MCNUM alpha = mccshutter_alpha;
MCNUM m = mccshutter_m;
MCNUM nslit = mccshutter_nslit;
MCNUM d = mccshutter_d;
MCNUM Qcx = mccshutter_Qcx;
MCNUM Qcy = mccshutter_Qcy;
MCNUM alphax = mccshutter_alphax;
MCNUM alphay = mccshutter_alphay;
MCNUM W = mccshutter_W;
MCNUM mx = mccshutter_mx;
MCNUM my = mccshutter_my;
MCNUM nu = mccshutter_nu;
MCNUM phase = mccshutter_phase;
#line 296 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  int i;

  magnify("xy");
  for(i = 0; i < nslit; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)l,
              i*w2c - w2/2.0,  h2/2.0, (double)l,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)l,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)l,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)l, w2/2.0, -h2/2.0, (double)l);

  if (nu || phase) {
    double radius = sqrt(w1*w1+l*l);
    /* cylinder top/center/bottom  */
    circle("xz", 0,-h1/2,l/2,radius);
    circle("xz", 0,0    ,l/2,radius);
    circle("xz", 0, h1/2,l/2,radius);
  }
}
#line 18243 "let_new.c"
}   /* End of shutter=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'insert'. */
  SIG_MESSAGE("insert (McDisplay)");
  printf("MCDISPLAY: component %s\n", "insert");
#define mccompcurname  insert
#define mccompcurtype  Guide_channeled
#define mccompcurindex 4
#define w1c mccinsert_w1c
#define w2c mccinsert_w2c
#define ww mccinsert_ww
#define hh mccinsert_hh
#define whalf mccinsert_whalf
#define hhalf mccinsert_hhalf
#define lwhalf mccinsert_lwhalf
#define lhhalf mccinsert_lhhalf
{   /* Declarations of insert=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccinsert_w1;
MCNUM h1 = mccinsert_h1;
MCNUM w2 = mccinsert_w2;
MCNUM h2 = mccinsert_h2;
MCNUM l = mccinsert_l;
MCNUM R0 = mccinsert_R0;
MCNUM Qc = mccinsert_Qc;
MCNUM alpha = mccinsert_alpha;
MCNUM m = mccinsert_m;
MCNUM nslit = mccinsert_nslit;
MCNUM d = mccinsert_d;
MCNUM Qcx = mccinsert_Qcx;
MCNUM Qcy = mccinsert_Qcy;
MCNUM alphax = mccinsert_alphax;
MCNUM alphay = mccinsert_alphay;
MCNUM W = mccinsert_W;
MCNUM mx = mccinsert_mx;
MCNUM my = mccinsert_my;
MCNUM nu = mccinsert_nu;
MCNUM phase = mccinsert_phase;
#line 296 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  int i;

  magnify("xy");
  for(i = 0; i < nslit; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)l,
              i*w2c - w2/2.0,  h2/2.0, (double)l,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)l,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)l,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)l, w2/2.0, -h2/2.0, (double)l);

  if (nu || phase) {
    double radius = sqrt(w1*w1+l*l);
    /* cylinder top/center/bottom  */
    circle("xz", 0,-h1/2,l/2,radius);
    circle("xz", 0,0    ,l/2,radius);
    circle("xz", 0, h1/2,l/2,radius);
  }
}
#line 18323 "let_new.c"
}   /* End of insert=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'afterinsert'. */
  SIG_MESSAGE("afterinsert (McDisplay)");
  printf("MCDISPLAY: component %s\n", "afterinsert");
#define mccompcurname  afterinsert
#define mccompcurtype  Guide_channeled
#define mccompcurindex 5
#define w1c mccafterinsert_w1c
#define w2c mccafterinsert_w2c
#define ww mccafterinsert_ww
#define hh mccafterinsert_hh
#define whalf mccafterinsert_whalf
#define hhalf mccafterinsert_hhalf
#define lwhalf mccafterinsert_lwhalf
#define lhhalf mccafterinsert_lhhalf
{   /* Declarations of afterinsert=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccafterinsert_w1;
MCNUM h1 = mccafterinsert_h1;
MCNUM w2 = mccafterinsert_w2;
MCNUM h2 = mccafterinsert_h2;
MCNUM l = mccafterinsert_l;
MCNUM R0 = mccafterinsert_R0;
MCNUM Qc = mccafterinsert_Qc;
MCNUM alpha = mccafterinsert_alpha;
MCNUM m = mccafterinsert_m;
MCNUM nslit = mccafterinsert_nslit;
MCNUM d = mccafterinsert_d;
MCNUM Qcx = mccafterinsert_Qcx;
MCNUM Qcy = mccafterinsert_Qcy;
MCNUM alphax = mccafterinsert_alphax;
MCNUM alphay = mccafterinsert_alphay;
MCNUM W = mccafterinsert_W;
MCNUM mx = mccafterinsert_mx;
MCNUM my = mccafterinsert_my;
MCNUM nu = mccafterinsert_nu;
MCNUM phase = mccafterinsert_phase;
#line 296 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  int i;

  magnify("xy");
  for(i = 0; i < nslit; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)l,
              i*w2c - w2/2.0,  h2/2.0, (double)l,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)l,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)l,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)l, w2/2.0, -h2/2.0, (double)l);

  if (nu || phase) {
    double radius = sqrt(w1*w1+l*l);
    /* cylinder top/center/bottom  */
    circle("xz", 0,-h1/2,l/2,radius);
    circle("xz", 0,0    ,l/2,radius);
    circle("xz", 0, h1/2,l/2,radius);
  }
}
#line 18403 "let_new.c"
}   /* End of afterinsert=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Res1'. */
  SIG_MESSAGE("Res1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Res1");
#define mccompcurname  Res1
#define mccompcurtype  DiskChopper
#define mccompcurindex 6
#define Tg mccRes1_Tg
#define To mccRes1_To
#define delta_y mccRes1_delta_y
#define height mccRes1_height
#define omega mccRes1_omega
{   /* Declarations of Res1=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccRes1_theta_0;
MCNUM radius = mccRes1_radius;
MCNUM yheight = mccRes1_yheight;
MCNUM nu = mccRes1_nu;
MCNUM nslit = mccRes1_nslit;
MCNUM jitter = mccRes1_jitter;
MCNUM delay = mccRes1_delay;
MCNUM isfirst = mccRes1_isfirst;
MCNUM n_pulse = mccRes1_n_pulse;
MCNUM abs_out = mccRes1_abs_out;
MCNUM phase = mccRes1_phase;
MCNUM xwidth = mccRes1_xwidth;
MCNUM verbose = mccRes1_verbose;
#line 165 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{

  int j;
  /* Arrays for storing geometry of slit/beamstop */
  magnify("xy");
  circle("xy", 0, -delta_y, 0, radius);

  /* Drawing the slit(s) */
  for (j=0; j<nslit; j++) {
    /* Angular start/end of slit */
    double tmin = j*(2.0*PI/nslit) - theta_0/2.0 + phase;
    double tmax = tmin+theta_0;
    /* Draw lines for each slit. */

    line(
      radius*sin(tmin),          radius*cos(tmin)-delta_y,          0,
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0
      );
    line(
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0,
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0);
    line(
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0,
      radius*sin(tmax),          radius*cos(tmax)-delta_y,          0);
  }
}
#line 18469 "let_new.c"
}   /* End of Res1=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'Res1_counter'. */
  SIG_MESSAGE("Res1_counter (McDisplay)");
  printf("MCDISPLAY: component %s\n", "Res1_counter");
#define mccompcurname  Res1_counter
#define mccompcurtype  DiskChopper
#define mccompcurindex 7
#define Tg mccRes1_counter_Tg
#define To mccRes1_counter_To
#define delta_y mccRes1_counter_delta_y
#define height mccRes1_counter_height
#define omega mccRes1_counter_omega
{   /* Declarations of Res1_counter=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccRes1_counter_theta_0;
MCNUM radius = mccRes1_counter_radius;
MCNUM yheight = mccRes1_counter_yheight;
MCNUM nu = mccRes1_counter_nu;
MCNUM nslit = mccRes1_counter_nslit;
MCNUM jitter = mccRes1_counter_jitter;
MCNUM delay = mccRes1_counter_delay;
MCNUM isfirst = mccRes1_counter_isfirst;
MCNUM n_pulse = mccRes1_counter_n_pulse;
MCNUM abs_out = mccRes1_counter_abs_out;
MCNUM phase = mccRes1_counter_phase;
MCNUM xwidth = mccRes1_counter_xwidth;
MCNUM verbose = mccRes1_counter_verbose;
#line 165 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{

  int j;
  /* Arrays for storing geometry of slit/beamstop */
  magnify("xy");
  circle("xy", 0, -delta_y, 0, radius);

  /* Drawing the slit(s) */
  for (j=0; j<nslit; j++) {
    /* Angular start/end of slit */
    double tmin = j*(2.0*PI/nslit) - theta_0/2.0 + phase;
    double tmax = tmin+theta_0;
    /* Draw lines for each slit. */

    line(
      radius*sin(tmin),          radius*cos(tmin)-delta_y,          0,
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0
      );
    line(
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0,
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0);
    line(
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0,
      radius*sin(tmax),          radius*cos(tmax)-delta_y,          0);
  }
}
#line 18532 "let_new.c"
}   /* End of Res1_counter=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'betweenchop1and2'. */
  SIG_MESSAGE("betweenchop1and2 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "betweenchop1and2");
#define mccompcurname  betweenchop1and2
#define mccompcurtype  Guide_channeled
#define mccompcurindex 8
#define w1c mccbetweenchop1and2_w1c
#define w2c mccbetweenchop1and2_w2c
#define ww mccbetweenchop1and2_ww
#define hh mccbetweenchop1and2_hh
#define whalf mccbetweenchop1and2_whalf
#define hhalf mccbetweenchop1and2_hhalf
#define lwhalf mccbetweenchop1and2_lwhalf
#define lhhalf mccbetweenchop1and2_lhhalf
{   /* Declarations of betweenchop1and2=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccbetweenchop1and2_w1;
MCNUM h1 = mccbetweenchop1and2_h1;
MCNUM w2 = mccbetweenchop1and2_w2;
MCNUM h2 = mccbetweenchop1and2_h2;
MCNUM l = mccbetweenchop1and2_l;
MCNUM R0 = mccbetweenchop1and2_R0;
MCNUM Qc = mccbetweenchop1and2_Qc;
MCNUM alpha = mccbetweenchop1and2_alpha;
MCNUM m = mccbetweenchop1and2_m;
MCNUM nslit = mccbetweenchop1and2_nslit;
MCNUM d = mccbetweenchop1and2_d;
MCNUM Qcx = mccbetweenchop1and2_Qcx;
MCNUM Qcy = mccbetweenchop1and2_Qcy;
MCNUM alphax = mccbetweenchop1and2_alphax;
MCNUM alphay = mccbetweenchop1and2_alphay;
MCNUM W = mccbetweenchop1and2_W;
MCNUM mx = mccbetweenchop1and2_mx;
MCNUM my = mccbetweenchop1and2_my;
MCNUM nu = mccbetweenchop1and2_nu;
MCNUM phase = mccbetweenchop1and2_phase;
#line 296 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  int i;

  magnify("xy");
  for(i = 0; i < nslit; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)l,
              i*w2c - w2/2.0,  h2/2.0, (double)l,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)l,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)l,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)l, w2/2.0, -h2/2.0, (double)l);

  if (nu || phase) {
    double radius = sqrt(w1*w1+l*l);
    /* cylinder top/center/bottom  */
    circle("xz", 0,-h1/2,l/2,radius);
    circle("xz", 0,0    ,l/2,radius);
    circle("xz", 0, h1/2,l/2,radius);
  }
}
#line 18609 "let_new.c"
}   /* End of betweenchop1and2=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'betweenchop2and3'. */
  SIG_MESSAGE("betweenchop2and3 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "betweenchop2and3");
#define mccompcurname  betweenchop2and3
#define mccompcurtype  Guide_channeled
#define mccompcurindex 9
#define w1c mccbetweenchop2and3_w1c
#define w2c mccbetweenchop2and3_w2c
#define ww mccbetweenchop2and3_ww
#define hh mccbetweenchop2and3_hh
#define whalf mccbetweenchop2and3_whalf
#define hhalf mccbetweenchop2and3_hhalf
#define lwhalf mccbetweenchop2and3_lwhalf
#define lhhalf mccbetweenchop2and3_lhhalf
{   /* Declarations of betweenchop2and3=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccbetweenchop2and3_w1;
MCNUM h1 = mccbetweenchop2and3_h1;
MCNUM w2 = mccbetweenchop2and3_w2;
MCNUM h2 = mccbetweenchop2and3_h2;
MCNUM l = mccbetweenchop2and3_l;
MCNUM R0 = mccbetweenchop2and3_R0;
MCNUM Qc = mccbetweenchop2and3_Qc;
MCNUM alpha = mccbetweenchop2and3_alpha;
MCNUM m = mccbetweenchop2and3_m;
MCNUM nslit = mccbetweenchop2and3_nslit;
MCNUM d = mccbetweenchop2and3_d;
MCNUM Qcx = mccbetweenchop2and3_Qcx;
MCNUM Qcy = mccbetweenchop2and3_Qcy;
MCNUM alphax = mccbetweenchop2and3_alphax;
MCNUM alphay = mccbetweenchop2and3_alphay;
MCNUM W = mccbetweenchop2and3_W;
MCNUM mx = mccbetweenchop2and3_mx;
MCNUM my = mccbetweenchop2and3_my;
MCNUM nu = mccbetweenchop2and3_nu;
MCNUM phase = mccbetweenchop2and3_phase;
#line 296 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  int i;

  magnify("xy");
  for(i = 0; i < nslit; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)l,
              i*w2c - w2/2.0,  h2/2.0, (double)l,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)l,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)l,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)l, w2/2.0, -h2/2.0, (double)l);

  if (nu || phase) {
    double radius = sqrt(w1*w1+l*l);
    /* cylinder top/center/bottom  */
    circle("xz", 0,-h1/2,l/2,radius);
    circle("xz", 0,0    ,l/2,radius);
    circle("xz", 0, h1/2,l/2,radius);
  }
}
#line 18689 "let_new.c"
}   /* End of betweenchop2and3=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'PR'. */
  SIG_MESSAGE("PR (McDisplay)");
  printf("MCDISPLAY: component %s\n", "PR");
#define mccompcurname  PR
#define mccompcurtype  DiskChopper
#define mccompcurindex 10
#define Tg mccPR_Tg
#define To mccPR_To
#define delta_y mccPR_delta_y
#define height mccPR_height
#define omega mccPR_omega
{   /* Declarations of PR=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccPR_theta_0;
MCNUM radius = mccPR_radius;
MCNUM yheight = mccPR_yheight;
MCNUM nu = mccPR_nu;
MCNUM nslit = mccPR_nslit;
MCNUM jitter = mccPR_jitter;
MCNUM delay = mccPR_delay;
MCNUM isfirst = mccPR_isfirst;
MCNUM n_pulse = mccPR_n_pulse;
MCNUM abs_out = mccPR_abs_out;
MCNUM phase = mccPR_phase;
MCNUM xwidth = mccPR_xwidth;
MCNUM verbose = mccPR_verbose;
#line 165 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{

  int j;
  /* Arrays for storing geometry of slit/beamstop */
  magnify("xy");
  circle("xy", 0, -delta_y, 0, radius);

  /* Drawing the slit(s) */
  for (j=0; j<nslit; j++) {
    /* Angular start/end of slit */
    double tmin = j*(2.0*PI/nslit) - theta_0/2.0 + phase;
    double tmax = tmin+theta_0;
    /* Draw lines for each slit. */

    line(
      radius*sin(tmin),          radius*cos(tmin)-delta_y,          0,
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0
      );
    line(
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0,
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0);
    line(
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0,
      radius*sin(tmax),          radius*cos(tmax)-delta_y,          0);
  }
}
#line 18755 "let_new.c"
}   /* End of PR=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'betweenchop3and4'. */
  SIG_MESSAGE("betweenchop3and4 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "betweenchop3and4");
#define mccompcurname  betweenchop3and4
#define mccompcurtype  Guide_channeled
#define mccompcurindex 11
#define w1c mccbetweenchop3and4_w1c
#define w2c mccbetweenchop3and4_w2c
#define ww mccbetweenchop3and4_ww
#define hh mccbetweenchop3and4_hh
#define whalf mccbetweenchop3and4_whalf
#define hhalf mccbetweenchop3and4_hhalf
#define lwhalf mccbetweenchop3and4_lwhalf
#define lhhalf mccbetweenchop3and4_lhhalf
{   /* Declarations of betweenchop3and4=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccbetweenchop3and4_w1;
MCNUM h1 = mccbetweenchop3and4_h1;
MCNUM w2 = mccbetweenchop3and4_w2;
MCNUM h2 = mccbetweenchop3and4_h2;
MCNUM l = mccbetweenchop3and4_l;
MCNUM R0 = mccbetweenchop3and4_R0;
MCNUM Qc = mccbetweenchop3and4_Qc;
MCNUM alpha = mccbetweenchop3and4_alpha;
MCNUM m = mccbetweenchop3and4_m;
MCNUM nslit = mccbetweenchop3and4_nslit;
MCNUM d = mccbetweenchop3and4_d;
MCNUM Qcx = mccbetweenchop3and4_Qcx;
MCNUM Qcy = mccbetweenchop3and4_Qcy;
MCNUM alphax = mccbetweenchop3and4_alphax;
MCNUM alphay = mccbetweenchop3and4_alphay;
MCNUM W = mccbetweenchop3and4_W;
MCNUM mx = mccbetweenchop3and4_mx;
MCNUM my = mccbetweenchop3and4_my;
MCNUM nu = mccbetweenchop3and4_nu;
MCNUM phase = mccbetweenchop3and4_phase;
#line 296 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  int i;

  magnify("xy");
  for(i = 0; i < nslit; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)l,
              i*w2c - w2/2.0,  h2/2.0, (double)l,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)l,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)l,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)l, w2/2.0, -h2/2.0, (double)l);

  if (nu || phase) {
    double radius = sqrt(w1*w1+l*l);
    /* cylinder top/center/bottom  */
    circle("xz", 0,-h1/2,l/2,radius);
    circle("xz", 0,0    ,l/2,radius);
    circle("xz", 0, h1/2,l/2,radius);
  }
}
#line 18832 "let_new.c"
}   /* End of betweenchop3and4=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'CR'. */
  SIG_MESSAGE("CR (McDisplay)");
  printf("MCDISPLAY: component %s\n", "CR");
#define mccompcurname  CR
#define mccompcurtype  DiskChopper
#define mccompcurindex 12
#define Tg mccCR_Tg
#define To mccCR_To
#define delta_y mccCR_delta_y
#define height mccCR_height
#define omega mccCR_omega
{   /* Declarations of CR=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccCR_theta_0;
MCNUM radius = mccCR_radius;
MCNUM yheight = mccCR_yheight;
MCNUM nu = mccCR_nu;
MCNUM nslit = mccCR_nslit;
MCNUM jitter = mccCR_jitter;
MCNUM delay = mccCR_delay;
MCNUM isfirst = mccCR_isfirst;
MCNUM n_pulse = mccCR_n_pulse;
MCNUM abs_out = mccCR_abs_out;
MCNUM phase = mccCR_phase;
MCNUM xwidth = mccCR_xwidth;
MCNUM verbose = mccCR_verbose;
#line 165 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{

  int j;
  /* Arrays for storing geometry of slit/beamstop */
  magnify("xy");
  circle("xy", 0, -delta_y, 0, radius);

  /* Drawing the slit(s) */
  for (j=0; j<nslit; j++) {
    /* Angular start/end of slit */
    double tmin = j*(2.0*PI/nslit) - theta_0/2.0 + phase;
    double tmax = tmin+theta_0;
    /* Draw lines for each slit. */

    line(
      radius*sin(tmin),          radius*cos(tmin)-delta_y,          0,
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0
      );
    line(
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0,
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0);
    line(
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0,
      radius*sin(tmax),          radius*cos(tmax)-delta_y,          0);
  }
}
#line 18898 "let_new.c"
}   /* End of CR=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'betweenchop4andatt'. */
  SIG_MESSAGE("betweenchop4andatt (McDisplay)");
  printf("MCDISPLAY: component %s\n", "betweenchop4andatt");
#define mccompcurname  betweenchop4andatt
#define mccompcurtype  Guide_channeled
#define mccompcurindex 13
#define w1c mccbetweenchop4andatt_w1c
#define w2c mccbetweenchop4andatt_w2c
#define ww mccbetweenchop4andatt_ww
#define hh mccbetweenchop4andatt_hh
#define whalf mccbetweenchop4andatt_whalf
#define hhalf mccbetweenchop4andatt_hhalf
#define lwhalf mccbetweenchop4andatt_lwhalf
#define lhhalf mccbetweenchop4andatt_lhhalf
{   /* Declarations of betweenchop4andatt=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccbetweenchop4andatt_w1;
MCNUM h1 = mccbetweenchop4andatt_h1;
MCNUM w2 = mccbetweenchop4andatt_w2;
MCNUM h2 = mccbetweenchop4andatt_h2;
MCNUM l = mccbetweenchop4andatt_l;
MCNUM R0 = mccbetweenchop4andatt_R0;
MCNUM Qc = mccbetweenchop4andatt_Qc;
MCNUM alpha = mccbetweenchop4andatt_alpha;
MCNUM m = mccbetweenchop4andatt_m;
MCNUM nslit = mccbetweenchop4andatt_nslit;
MCNUM d = mccbetweenchop4andatt_d;
MCNUM Qcx = mccbetweenchop4andatt_Qcx;
MCNUM Qcy = mccbetweenchop4andatt_Qcy;
MCNUM alphax = mccbetweenchop4andatt_alphax;
MCNUM alphay = mccbetweenchop4andatt_alphay;
MCNUM W = mccbetweenchop4andatt_W;
MCNUM mx = mccbetweenchop4andatt_mx;
MCNUM my = mccbetweenchop4andatt_my;
MCNUM nu = mccbetweenchop4andatt_nu;
MCNUM phase = mccbetweenchop4andatt_phase;
#line 296 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  int i;

  magnify("xy");
  for(i = 0; i < nslit; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)l,
              i*w2c - w2/2.0,  h2/2.0, (double)l,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)l,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)l,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)l, w2/2.0, -h2/2.0, (double)l);

  if (nu || phase) {
    double radius = sqrt(w1*w1+l*l);
    /* cylinder top/center/bottom  */
    circle("xz", 0,-h1/2,l/2,radius);
    circle("xz", 0,0    ,l/2,radius);
    circle("xz", 0, h1/2,l/2,radius);
  }
}
#line 18975 "let_new.c"
}   /* End of betweenchop4andatt=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'guide_move'. */
  SIG_MESSAGE("guide_move (McDisplay)");
  printf("MCDISPLAY: component %s\n", "guide_move");
#define mccompcurname  guide_move
#define mccompcurtype  Guide_channeled
#define mccompcurindex 14
#define w1c mccguide_move_w1c
#define w2c mccguide_move_w2c
#define ww mccguide_move_ww
#define hh mccguide_move_hh
#define whalf mccguide_move_whalf
#define hhalf mccguide_move_hhalf
#define lwhalf mccguide_move_lwhalf
#define lhhalf mccguide_move_lhhalf
{   /* Declarations of guide_move=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccguide_move_w1;
MCNUM h1 = mccguide_move_h1;
MCNUM w2 = mccguide_move_w2;
MCNUM h2 = mccguide_move_h2;
MCNUM l = mccguide_move_l;
MCNUM R0 = mccguide_move_R0;
MCNUM Qc = mccguide_move_Qc;
MCNUM alpha = mccguide_move_alpha;
MCNUM m = mccguide_move_m;
MCNUM nslit = mccguide_move_nslit;
MCNUM d = mccguide_move_d;
MCNUM Qcx = mccguide_move_Qcx;
MCNUM Qcy = mccguide_move_Qcy;
MCNUM alphax = mccguide_move_alphax;
MCNUM alphay = mccguide_move_alphay;
MCNUM W = mccguide_move_W;
MCNUM mx = mccguide_move_mx;
MCNUM my = mccguide_move_my;
MCNUM nu = mccguide_move_nu;
MCNUM phase = mccguide_move_phase;
#line 296 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  int i;

  magnify("xy");
  for(i = 0; i < nslit; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)l,
              i*w2c - w2/2.0,  h2/2.0, (double)l,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)l,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)l,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)l, w2/2.0, -h2/2.0, (double)l);

  if (nu || phase) {
    double radius = sqrt(w1*w1+l*l);
    /* cylinder top/center/bottom  */
    circle("xz", 0,-h1/2,l/2,radius);
    circle("xz", 0,0    ,l/2,radius);
    circle("xz", 0, h1/2,l/2,radius);
  }
}
#line 19055 "let_new.c"
}   /* End of guide_move=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'funnel'. */
  SIG_MESSAGE("funnel (McDisplay)");
  printf("MCDISPLAY: component %s\n", "funnel");
#define mccompcurname  funnel
#define mccompcurtype  Guide_channeled
#define mccompcurindex 15
#define w1c mccfunnel_w1c
#define w2c mccfunnel_w2c
#define ww mccfunnel_ww
#define hh mccfunnel_hh
#define whalf mccfunnel_whalf
#define hhalf mccfunnel_hhalf
#define lwhalf mccfunnel_lwhalf
#define lhhalf mccfunnel_lhhalf
{   /* Declarations of funnel=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccfunnel_w1;
MCNUM h1 = mccfunnel_h1;
MCNUM w2 = mccfunnel_w2;
MCNUM h2 = mccfunnel_h2;
MCNUM l = mccfunnel_l;
MCNUM R0 = mccfunnel_R0;
MCNUM Qc = mccfunnel_Qc;
MCNUM alpha = mccfunnel_alpha;
MCNUM m = mccfunnel_m;
MCNUM nslit = mccfunnel_nslit;
MCNUM d = mccfunnel_d;
MCNUM Qcx = mccfunnel_Qcx;
MCNUM Qcy = mccfunnel_Qcy;
MCNUM alphax = mccfunnel_alphax;
MCNUM alphay = mccfunnel_alphay;
MCNUM W = mccfunnel_W;
MCNUM mx = mccfunnel_mx;
MCNUM my = mccfunnel_my;
MCNUM nu = mccfunnel_nu;
MCNUM phase = mccfunnel_phase;
#line 296 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  int i;

  magnify("xy");
  for(i = 0; i < nslit; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)l,
              i*w2c - w2/2.0,  h2/2.0, (double)l,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)l,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)l,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)l, w2/2.0, -h2/2.0, (double)l);

  if (nu || phase) {
    double radius = sqrt(w1*w1+l*l);
    /* cylinder top/center/bottom  */
    circle("xz", 0,-h1/2,l/2,radius);
    circle("xz", 0,0    ,l/2,radius);
    circle("xz", 0, h1/2,l/2,radius);
  }
}
#line 19135 "let_new.c"
}   /* End of funnel=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'C5'. */
  SIG_MESSAGE("C5 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "C5");
#define mccompcurname  C5
#define mccompcurtype  DiskChopper
#define mccompcurindex 16
#define Tg mccC5_Tg
#define To mccC5_To
#define delta_y mccC5_delta_y
#define height mccC5_height
#define omega mccC5_omega
{   /* Declarations of C5=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccC5_theta_0;
MCNUM radius = mccC5_radius;
MCNUM yheight = mccC5_yheight;
MCNUM nu = mccC5_nu;
MCNUM nslit = mccC5_nslit;
MCNUM jitter = mccC5_jitter;
MCNUM delay = mccC5_delay;
MCNUM isfirst = mccC5_isfirst;
MCNUM n_pulse = mccC5_n_pulse;
MCNUM abs_out = mccC5_abs_out;
MCNUM phase = mccC5_phase;
MCNUM xwidth = mccC5_xwidth;
MCNUM verbose = mccC5_verbose;
#line 165 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{

  int j;
  /* Arrays for storing geometry of slit/beamstop */
  magnify("xy");
  circle("xy", 0, -delta_y, 0, radius);

  /* Drawing the slit(s) */
  for (j=0; j<nslit; j++) {
    /* Angular start/end of slit */
    double tmin = j*(2.0*PI/nslit) - theta_0/2.0 + phase;
    double tmax = tmin+theta_0;
    /* Draw lines for each slit. */

    line(
      radius*sin(tmin),          radius*cos(tmin)-delta_y,          0,
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0
      );
    line(
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0,
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0);
    line(
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0,
      radius*sin(tmax),          radius*cos(tmax)-delta_y,          0);
  }
}
#line 19201 "let_new.c"
}   /* End of C5=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'C5_counter'. */
  SIG_MESSAGE("C5_counter (McDisplay)");
  printf("MCDISPLAY: component %s\n", "C5_counter");
#define mccompcurname  C5_counter
#define mccompcurtype  DiskChopper
#define mccompcurindex 17
#define Tg mccC5_counter_Tg
#define To mccC5_counter_To
#define delta_y mccC5_counter_delta_y
#define height mccC5_counter_height
#define omega mccC5_counter_omega
{   /* Declarations of C5_counter=DiskChopper() SETTING parameters. */
MCNUM theta_0 = mccC5_counter_theta_0;
MCNUM radius = mccC5_counter_radius;
MCNUM yheight = mccC5_counter_yheight;
MCNUM nu = mccC5_counter_nu;
MCNUM nslit = mccC5_counter_nslit;
MCNUM jitter = mccC5_counter_jitter;
MCNUM delay = mccC5_counter_delay;
MCNUM isfirst = mccC5_counter_isfirst;
MCNUM n_pulse = mccC5_counter_n_pulse;
MCNUM abs_out = mccC5_counter_abs_out;
MCNUM phase = mccC5_counter_phase;
MCNUM xwidth = mccC5_counter_xwidth;
MCNUM verbose = mccC5_counter_verbose;
#line 165 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\DiskChopper.comp"
{

  int j;
  /* Arrays for storing geometry of slit/beamstop */
  magnify("xy");
  circle("xy", 0, -delta_y, 0, radius);

  /* Drawing the slit(s) */
  for (j=0; j<nslit; j++) {
    /* Angular start/end of slit */
    double tmin = j*(2.0*PI/nslit) - theta_0/2.0 + phase;
    double tmax = tmin+theta_0;
    /* Draw lines for each slit. */

    line(
      radius*sin(tmin),          radius*cos(tmin)-delta_y,          0,
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0
      );
    line(
      (radius-height)*sin(tmin), (radius-height)*cos(tmin)-delta_y, 0,
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0);
    line(
      (radius-height)*sin(tmax), (radius-height)*cos(tmax)-delta_y, 0,
      radius*sin(tmax),          radius*cos(tmax)-delta_y,          0);
  }
}
#line 19264 "let_new.c"
}   /* End of C5_counter=DiskChopper() SETTING parameter declarations. */
#undef omega
#undef height
#undef delta_y
#undef To
#undef Tg
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'endguide'. */
  SIG_MESSAGE("endguide (McDisplay)");
  printf("MCDISPLAY: component %s\n", "endguide");
#define mccompcurname  endguide
#define mccompcurtype  Guide_channeled
#define mccompcurindex 18
#define w1c mccendguide_w1c
#define w2c mccendguide_w2c
#define ww mccendguide_ww
#define hh mccendguide_hh
#define whalf mccendguide_whalf
#define hhalf mccendguide_hhalf
#define lwhalf mccendguide_lwhalf
#define lhhalf mccendguide_lhhalf
{   /* Declarations of endguide=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccendguide_w1;
MCNUM h1 = mccendguide_h1;
MCNUM w2 = mccendguide_w2;
MCNUM h2 = mccendguide_h2;
MCNUM l = mccendguide_l;
MCNUM R0 = mccendguide_R0;
MCNUM Qc = mccendguide_Qc;
MCNUM alpha = mccendguide_alpha;
MCNUM m = mccendguide_m;
MCNUM nslit = mccendguide_nslit;
MCNUM d = mccendguide_d;
MCNUM Qcx = mccendguide_Qcx;
MCNUM Qcy = mccendguide_Qcy;
MCNUM alphax = mccendguide_alphax;
MCNUM alphay = mccendguide_alphay;
MCNUM W = mccendguide_W;
MCNUM mx = mccendguide_mx;
MCNUM my = mccendguide_my;
MCNUM nu = mccendguide_nu;
MCNUM phase = mccendguide_phase;
#line 296 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  int i;

  magnify("xy");
  for(i = 0; i < nslit; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)l,
              i*w2c - w2/2.0,  h2/2.0, (double)l,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)l,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)l,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)l, w2/2.0, -h2/2.0, (double)l);

  if (nu || phase) {
    double radius = sqrt(w1*w1+l*l);
    /* cylinder top/center/bottom  */
    circle("xz", 0,-h1/2,l/2,radius);
    circle("xz", 0,0    ,l/2,radius);
    circle("xz", 0, h1/2,l/2,radius);
  }
}
#line 19341 "let_new.c"
}   /* End of endguide=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'snot'. */
  SIG_MESSAGE("snot (McDisplay)");
  printf("MCDISPLAY: component %s\n", "snot");
#define mccompcurname  snot
#define mccompcurtype  Guide_channeled
#define mccompcurindex 19
#define w1c mccsnot_w1c
#define w2c mccsnot_w2c
#define ww mccsnot_ww
#define hh mccsnot_hh
#define whalf mccsnot_whalf
#define hhalf mccsnot_hhalf
#define lwhalf mccsnot_lwhalf
#define lhhalf mccsnot_lhhalf
{   /* Declarations of snot=Guide_channeled() SETTING parameters. */
MCNUM w1 = mccsnot_w1;
MCNUM h1 = mccsnot_h1;
MCNUM w2 = mccsnot_w2;
MCNUM h2 = mccsnot_h2;
MCNUM l = mccsnot_l;
MCNUM R0 = mccsnot_R0;
MCNUM Qc = mccsnot_Qc;
MCNUM alpha = mccsnot_alpha;
MCNUM m = mccsnot_m;
MCNUM nslit = mccsnot_nslit;
MCNUM d = mccsnot_d;
MCNUM Qcx = mccsnot_Qcx;
MCNUM Qcy = mccsnot_Qcy;
MCNUM alphax = mccsnot_alphax;
MCNUM alphay = mccsnot_alphay;
MCNUM W = mccsnot_W;
MCNUM mx = mccsnot_mx;
MCNUM my = mccsnot_my;
MCNUM nu = mccsnot_nu;
MCNUM phase = mccsnot_phase;
#line 296 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Guide_channeled.comp"
{
  int i;

  magnify("xy");
  for(i = 0; i < nslit; i++)
  {
    multiline(5,
              i*w1c - w1/2.0, -h1/2.0, 0.0,
              i*w2c - w2/2.0, -h2/2.0, (double)l,
              i*w2c - w2/2.0,  h2/2.0, (double)l,
              i*w1c - w1/2.0,  h1/2.0, 0.0,
              i*w1c - w1/2.0, -h1/2.0, 0.0);
    multiline(5,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0,
              (i+1)*w2c - d - w2/2.0, -h2/2.0, (double)l,
              (i+1)*w2c - d - w2/2.0,  h2/2.0, (double)l,
              (i+1)*w1c - d - w1/2.0,  h1/2.0, 0.0,
              (i+1)*w1c - d - w1/2.0, -h1/2.0, 0.0);
  }
  line(-w1/2.0, -h1/2.0, 0.0, w1/2.0, -h1/2.0, 0.0);
  line(-w2/2.0, -h2/2.0, (double)l, w2/2.0, -h2/2.0, (double)l);

  if (nu || phase) {
    double radius = sqrt(w1*w1+l*l);
    /* cylinder top/center/bottom  */
    circle("xz", 0,-h1/2,l/2,radius);
    circle("xz", 0,0    ,l/2,radius);
    circle("xz", 0, h1/2,l/2,radius);
  }
}
#line 19421 "let_new.c"
}   /* End of snot=Guide_channeled() SETTING parameter declarations. */
#undef lhhalf
#undef lwhalf
#undef hhalf
#undef whalf
#undef hh
#undef ww
#undef w2c
#undef w1c
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'detslit'. */
  SIG_MESSAGE("detslit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "detslit");
#define mccompcurname  detslit
#define mccompcurtype  Slit
#define mccompcurindex 20
{   /* Declarations of detslit=Slit() SETTING parameters. */
MCNUM xmin = mccdetslit_xmin;
MCNUM xmax = mccdetslit_xmax;
MCNUM ymin = mccdetslit_ymin;
MCNUM ymax = mccdetslit_ymax;
MCNUM radius = mccdetslit_radius;
MCNUM xwidth = mccdetslit_xwidth;
MCNUM yheight = mccdetslit_yheight;
#line 66 "C:\\mcstas-2.4.1\\lib\\tools\\Python\\mcrun\\..\\mccodelib\\..\\..\\..\\optics\\Slit.comp"
{
  magnify("xy");
  if (radius == 0) {
    double xw, yh;
    xw = (xmax - xmin)/2.0;
    yh = (ymax - ymin)/2.0;
    multiline(3, xmin-xw, (double)ymax, 0.0,
              (double)xmin, (double)ymax, 0.0,
              (double)xmin, ymax+yh, 0.0);
    multiline(3, xmax+xw, (double)ymax, 0.0,
              (double)xmax, (double)ymax, 0.0,
              (double)xmax, ymax+yh, 0.0);
    multiline(3, xmin-xw, (double)ymin, 0.0,
              (double)xmin, (double)ymin, 0.0,
              (double)xmin, ymin-yh, 0.0);
    multiline(3, xmax+xw, (double)ymin, 0.0,
              (double)xmax, (double)ymin, 0.0,
              (double)xmax, ymin-yh, 0.0);
  } else {
    circle("xy",0,0,0,radius);
  }
}
#line 19472 "let_new.c"
}   /* End of detslit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'mcplout'. */
  SIG_MESSAGE("mcplout (McDisplay)");
  printf("MCDISPLAY: component %s\n", "mcplout");
#define mccompcurname  mcplout
#define mccompcurtype  MCPL_output_horace
#define mccompcurindex 21
#define polarisationuse mccmcplout_polarisationuse
#define doubleprec mccmcplout_doubleprec
#define verbose mccmcplout_verbose
#define radius mccmcplout_radius
#define thickness mccmcplout_thickness
#define yheight mccmcplout_yheight
#define restore_neutron mccmcplout_restore_neutron
#define userflag mccmcplout_userflag
#define plate mccmcplout_plate
#define outputfile mccmcplout_outputfile
#define particle mccmcplout_particle
#define Particle mccmcplout_Particle
#define userflagenabled mccmcplout_userflagenabled
{   /* Declarations of mcplout=MCPL_output_horace() SETTING parameters. */
char* filename = mccmcplout_filename;
char* userflagcomment = mccmcplout_userflagcomment;
MCNUM merge_mpi = mccmcplout_merge_mpi;
MCNUM keep_mpi_unmerged = mccmcplout_keep_mpi_unmerged;
#line 342 "MCPL_output_horace.comp"
{
    double t,dt;
    int i;
    multiline(5, 0.2,0.2,0.0, -0.2,0.2,0.0, -0.2,-0.2,0.0, 0.2,-0.2,0.0, 0.2,0.2,0.0);
    /*M*/
    multiline(5,-0.085,-0.085,0.0, -0.085,0.085,0.0, -0.045,-0.085,0.0, -0.005,0.085,0.0, -0.005,-0.085,0.0);
    /*O*/
    dt=2*M_PI/32;
    t=0;
    for (i=0;i<32;i++){
        line(0.04*cos(t)+0.045,0.08*sin(t),0, 0.04*cos(t+dt)+0.045,0.08*sin(t+dt),0);
        t+=dt;
    }
}
#line 19517 "let_new.c"
}   /* End of mcplout=MCPL_output_horace() SETTING parameter declarations. */
#undef userflagenabled
#undef Particle
#undef particle
#undef outputfile
#undef plate
#undef userflag
#undef restore_neutron
#undef yheight
#undef thickness
#undef radius
#undef verbose
#undef doubleprec
#undef polarisationuse
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  printf("MCDISPLAY: end\n");
} /* end display */
#undef magnify
#undef line
#undef dashed_line
#undef multiline
#undef rectangle
#undef box
#undef circle
/* end of generated C code let_new.c */
