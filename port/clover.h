# ifndef QOP_CLOVER_DEFAULT_PRECISION
#  define QOP_CLOVER_DEFAULT_PRECISION 'D'
# endif

#ifndef MARK_21eedb1e5eff4d4ba0462475b43ae655
#define MARK_21eedb1e5eff4d4ba0462475b43ae655

# include <qop-clover.h>
# include <stdlib.h>
# include <string.h>
# include <qmp.h>
# include <sys/time.h>

# define q(x) qop_clover_##x
# define qf(x) qop_f3_clover_##x
# define qd(x) qop_d3_clover_##x
# define Q(x) QOP_CLOVER_##x
# define QF(x) QOP_F3_CLOVER_##x
# define QD(x) QOP_D3_CLOVER_##x

/* Cache size */
#define CACHE_LINE_SIZE 128
#define ALIGN(p) ((void *)((((ptrdiff_t)(p))+CACHE_LINE_SIZE-1) & \
                           ~(CACHE_LINE_SIZE-1)))


/* QCD types (qa0 controls these definitions) */
struct SUnF;
struct SUnD;
struct CloverF;
struct CloverD;
struct FermionF;
struct FermionD;
struct ProjectedFermionF;
struct ProjectedFermionD;
struct MxM_workspaceF;
struct MxM_workspaceD;

/* Internal types */
struct local {
  int lo[Q(DIM)];
  int hi[Q(DIM)];
  int dx[Q(DIM)];
};

/* structs neighbor and up_pack are defined by qa0 */
struct neighbor;
struct up_pack;
struct down_pack;

struct eo_lattice {
  struct Q(State) *state;                  /* back pointer to state */
  int              face_size;              /* 4-d size of the face */
  int              body_size;              /* 4-d size of the body */
  int              full_size;              /* face + body */
  int             *lx2v;                   /* 4-d layout 2 vector translation */
  int             *v2lx;                   /* only for init */
  int              Ls;                     /* Ls */
  struct local    *local;                  /* points to state.local */

  struct neighbor *neighbor;               /* neighbor data (body,face) */
  int              send_up_size[Q(DIM)];   /* 4-d send size in each up-dir */
  struct up_pack  *up_pack[Q(DIM)];        /* 4-d (U,f) for up-face packing */
  int              send_down_size[Q(DIM)]; /* 4-d send size in each down-dir */
  struct down_pack *down_pack[Q(DIM)];      /* 4-d (f) for down-face packing */
  int              receive_up_size[Q(DIM)]; /* 4-d (U,f) up receive size */
  int              receive_down_size[Q(DIM)]; /* 4-d (f) down receive size */

  int              real_size;              /* 0, 4 or 8 */ 
  int              h_valid;                /* is .handle valid? */
  QMP_msghandle_t  handle;                 /* global send&receive handle */
  QMP_msghandle_t  th[4*Q(DIM)];           /* transitody handles */
  int              th_count;               /* number of valid th[] */
  QMP_msgmem_t     mh[4*Q(DIM)];           /* memory handles for th[] */
  int              mh_count;               /* number of valid mh[] */
  QMP_mem_t       *mem[4*Q(DIM)];          /* memory for mh[] */
  int              mem_count;              /* number of valid mem[] */
  void            *send_up_buf[Q(DIM)];    /* pf up-bufs */
  void            *send_down_buf[Q(DIM)];  /* pf down-bufs */
  void            *receive_buf[2*Q(DIM)];  /* pf receive bufs (up[], down[]) */
  int              total_send;             /* bytes to send */
  int              total_receive;          /* bytes to receive */
};

struct Q(State) {
  const char        *version;         /* to get version string into app */
  int                used;            /* gc ref counter */
  int                saved;           /* gc internal counter */
  size_t             allocated;       /* currently allocated bytes */
  size_t             max_allocated;   /* maximum allocation */

  int                error_latched;   /* if 0, allow error recording */
  int                fatal_error;     /* if 0, allow reseting latch */
  const char        *error;           /* error string */

  int                real_size;       /* 0, 4 or 8 */ 
  struct eo_lattice  even;            /* even sublattice */
  struct eo_lattice  odd;             /* odd sublattice */

  struct timeval     t0, t1;          /* for timing */
  double             time_sec;        /* seconds in the last routine */
  long long          flops;           /* FLOP in the last routine */
  long long          sent;            /* bytes sent in the last routine */
  long long          received;        /* bytes received in the last routine */

  int                volume;          /* 4-d volume */
  int                lattice[Q(DIM)]; /* 4-d lattice size */
  struct local       local;           /* 4-d local sublattice */
  int                node[Q(DIM)];    /* local node address */
  int                network[Q(DIM)]; /* the network geometry */
  int                master_p;        /* are we the master? */
  int               *lx2v;            /* Sublattice 1-d -> 4-d translation */
  int               *v2lx;            /* Only for init */
};

typedef enum {
    CG_SUCCESS,
    CG_MAXITER,
    CG_EIGCONV,
    CG_ZEROMODE,
    CG_NOEMEM
} CG_STATUS;



/* debug printing */
extern int QDP_this_node;
extern int QDP_is_initialized(void);

#define XXX_DEBUG

#if defined(XXX_DEBUG)
#  define LOG_PRINT(...)      do { \
    if (!QDP_is_initialized() || 0 == QDP_this_node) \
        printf(__VA_ARGS__); \
} while (0)
#  define LOG_FPRINT(fp, ...) do { \
    if (!QDP_is_initialized() || 0 == QDP_this_node) \
        fprintf(fp, __VA_ARGS__);\
} while (0)
#  define LOG_ECHO(fp,...)    do { \
    printf("%s[%d]: ", __func__, QDP_this_node); \
    printf(__VA_ARGS__); \
} while (0)
#  define LOG_FECHO(fp,...)   do { \
    fprintf(fp, "%s[%d]: ", __func__, QDP_this_node); \
    fprintf(fp, __VA_ARGS__); \
} while (0)
#else
#  define LOG_PRINT(...)        do {} while (0)
#  define LOG_FPRINT(fp, ...)   do {} while (0)
#  define LOG_ECHO(fp,...)      do {} while (0)
#  define LOG_FECHO(fp,...)     do {} while (0)
#endif

/* Deflator state */
#include <deflator-la.h>

#if defined(HAVE_LAPACK)
#elif defined(HAVE_GSL)
#  include <gsl/gsl_vector.h>
#  include <gsl/gsl_matrix.h>
#  include <gsl/gsl_eigen.h>
#else
#  error "no linear algebra library"
#endif

#define DEFLATOR_VEC_SIZE(pstate) ((pstate)->even.full_size)
struct Q(Deflator) {
    struct Q(State) *state;

    /* XXX other pieces of the deflation state */
    int                 dim;     /* size of problem vectors */

    int                 vmax;
    int                 vsize;
    int                 nev;
    int                 umax;
    int                 usize;
    int                 frozen;

    /* eig current state */
    double              eps;
    double              resid_norm_sq_min;
    
    latmat_c            V;
    doublecomplex       *T;

    /* incr_eig current state */
    latmat_c            U;
    doublecomplex       *H;
    doublecomplex       *C;


    long int            lwork;
    doublecomplex       *zwork;
    doublecomplex       *hevecs2;
    double              *hevals;
#if defined(HAVE_LAPACK)
    doublecomplex       *hevecs1;
    doublecomplex       *tau;
    double              *rwork;
    
    double              *debug_hevals;
    long int            debug_lwork;
    doublecomplex       *debug_zwork;
    double              *debug_rwork;

#elif defined(HAVE_GSL)
    doublecomplex       *zwork2;
    gsl_matrix_complex  *gsl_T_full;
    gsl_matrix_complex  *gsl_hevecs1;
    gsl_vector          *gsl_hevals1;
    gsl_eigen_hermv_workspace *gsl_wkspace1;
    gsl_matrix_complex  *gsl_T_m1;
    gsl_matrix_complex  *gsl_hevecs2;
    gsl_vector          *gsl_hevals2;
    gsl_eigen_hermv_workspace *gsl_wkspace2;
    gsl_matrix_complex  *gsl_T_proj;
    gsl_matrix_complex  *gsl_hevecs3;
    gsl_eigen_hermv_workspace *gsl_wkspace3;
    gsl_matrix_complex  *gsl_QR;
    gsl_matrix_complex  *gsl_Q_unpack;
    gsl_matrix_complex  *gsl_tmp_MxS;
    gsl_vector_complex  *gsl_tau;
    size_t              *hevals_select1;
    size_t              *hevals_select2;

    gsl_vector          *debug_gsl_hevals;
    gsl_eigen_herm_workspace *debug_gsl_wkspace;

#else
#  error "no linear algebra library"
#endif

    latmat_c            tmp_V;
    latvec_c            work_c_1;
    latvec_c            work_c_2;
    latvec_c            work_c_3;

};

/* layout translation */
void q(l2v)(int x[Q(DIM)], const struct local *local, int p);
int q(v2l)(const int x[Q(DIM)], const struct local *local);

/* Implementation functions */
int q(set_error)(struct Q(State) *state, int fatal, const char *error);

int q(setup_comm)(struct Q(State) *state, int real_size);
int q(free_comm)(struct Q(State) *state);

void *q(malloc)(struct Q(State) *state, size_t bytes);
void *q(allocate_aligned)(struct Q(State) *state,
                          size_t *size, void **aligned_ptr,
                          size_t hdr_size, size_t bulk_size);
void q(free)(struct Q(State) *state, void *ptr, size_t bytes);
void q(cleanup_state)(struct Q(State) *state);

/* Backend controled structure sizes */
int q(sizeof_neighbor)(int volume);
int q(sizeof_up_pack)(int volume);
int q(sizeof_down_pack)(int volume);

/* qa0 level data access routines */
int q(get_down_pack_f)(const struct down_pack *up, int p);
int q(get_up_pack_f)(const struct up_pack *up, int p);
void q(put_down_pack)(struct down_pack *down, int p, int f);
void q(get_down_pack)(int *f, const struct down_pack *up, int p);
void q(put_up_pack)(struct up_pack *up, int p, int f, int u);
void q(get_up_pack)(int *f, int *u, const struct up_pack *up, int p);
void q(put_neighbor)(struct neighbor *n, int p,
                     int m,
                     const int f_up[Q(DIM)], int u_up,
                     const int f_down[Q(DIM)], const int u_down[Q(DIM)]);
void q(get_neighbor)(int *m, int *f_up, int *u_up,
                     int *f_down, int *u_down,
                     const struct neighbor *n, int p);
void q(fix_neighbor_f_up)(struct neighbor *n, int p, int f_up, int d);
void q(fix_neighbor_f_down)(struct neighbor *n, int p, int f_down, int d);

/* mixed precision operations */
/* Fd = Fd + Ff */
unsigned int q(f_d_eq_dpf)(struct FermionD *dst,
                           int size,
                           const struct FermionD *src_d,
                           const struct FermionF *src_f);
/* Ff = Fd - Fd */
unsigned int q(f_f_eq_dmd_norm2)(struct FermionF *dst,
                                 double *local_norm,
                                 int size,
                                 const struct FermionD *src_a,
                                 const struct FermionD *src_b);

/* converting gauge and clover from double down to float */
void q(g_f_eq_d)(struct SUnF *dst,
                 int size,
                 const struct SUnD *src);
void q(c_f_eq_d)(struct CloverF *dst,
                 int size,
                 const struct CloverD *src);

/* the mixed solver */
int q(mixed_cg)(struct Q(State)             *state,
                const char                  *name,
                struct QD(Fermion)          *psi,
                int                         *out_iterations,
                double                      *out_epsilon,
                const struct QD(Fermion)    *psi_0,
                const struct QD(Gauge)      *gauge,
                const struct QD(Fermion)    *eta,
                struct Q(Deflator)          *deflator,
                int                          f_iter,
                double                       f_epsilon,
                int                          max_iterations,
                double                       min_epsilon,
                unsigned int                 options);

/* handling eig deflator */
int Q(create_deflator)(
        struct Q(Deflator) **deflator_ptr,
        struct Q(State) *s,
        int vmax, int nev,
        double eps, int umax);
void Q(free_deflator)(struct Q(Deflator) **deflator_ptr);
void Q(deflator_reset)(struct Q(Deflator) *deflator);
void Q(deflator_stop)(struct Q(Deflator) *deflator);
void Q(deflator_resume)(struct Q(Deflator) *deflator);
int Q(deflator_is_stopped)(struct Q(Deflator) *deflator);
int q(df_preamble)(struct Q(State)           *state,
                   struct Q(Deflator)        *deflator,
                   struct FermionF           *psi_e,
                   struct FermionF           *rho_e,
                   double                    *rho_norm2,
                   struct FermionF           *chi_e, /* const ! */
                   struct MxM_workspaceF     *ws);
int q(df_update0)(struct Q(State)          *state,
                  struct Q(Deflator)       *deflator,
                  double                    a1,
                  double                    b1,
                  double                    a0,
                  double                    b0,
                  double                    r,
                  struct FermionF          *rho);
int q(df_update1)(struct Q(State)          *state,
                   struct Q(Deflator)       *deflator,
                   double                    a1,
                   double                    b1,
                   double                    a0,
                   double                    b0,
                   double                    r,
                   struct FermionF          *rho,
                   struct FermionF          *A_rho);
int q(df_postamble)(struct Q(State)           *state,
                    struct Q(Deflator)        *deflator,
                    struct MxM_workspaceF     *ws);

/* Timing */
#define BEGIN_TIMING(s) do { gettimeofday(&((s)->t0), NULL); } while (0)
#define END_TIMING(s, f, snd, rcv) do { \
    gettimeofday(&((s)->t1), NULL); \
    (s)->time_sec = ((s)->t1.tv_sec - (s)->t0.tv_sec) \
      + 1e-6 * ((s)->t1.tv_usec - (s)->t0.tv_usec); \
    (s)->flops = (f); (s)->sent = (snd); (s)->received = (rcv); } while (0)

/* Argument checking */
#define DECLARE_STATE struct Q(State) *state = NULL
#define CHECK_ARG0(n) do { if ((n) == 0) return 1;      \
    state = (n)->state; } while (0)
#define CHECK_ARGn(n,f) do { if ((n) == 0)                              \
      return q(set_error)(state, 0, f "(): NULL argument");             \
    if ((n)->state != state)                                            \
      return q(set_error)(state, 0, f "(): geometry mismatch"); } while (0)
#define CHECK_POINTER(n,f) do { if ((n) == 0)                           \
      return q(set_error)(state, 0, f "(): NULL argument"); } while (0)

#endif /* !defined(MARK_21eedb1e5eff4d4ba0462475b43ae655) */

# undef qx
# undef QX
# undef REAL
# undef SUn
# undef Clover
# undef Fermion
# undef ProjectedFermion
# undef MxM_workspace
# if QOP_CLOVER_DEFAULT_PRECISION=='D'
#  define qx(x) qop_d3_clover_##x
#  define QX(x) QOP_D3_CLOVER_##x
#  define REAL double
#  define SUn SUnD
#  define Clover CloverD
#  define Fermion FermionD
#  define ProjectedFermion ProjectedFermionD
#  define MxM_workspace MxM_workspaceD
# endif
# if QOP_CLOVER_DEFAULT_PRECISION=='F'
#  define qx(x) qop_f3_clover_##x
#  define QX(x) QOP_F3_CLOVER_##x
#  define REAL float
#  define SUn SUnF
#  define Clover CloverF
#  define Fermion FermionF
#  define ProjectedFermion ProjectedFermionF
#  define MxM_workspace MxM_workspaceF
# endif

/* CLOVER types */
struct QX(Fermion) {
  struct Q(State) *state;
  size_t size;
  struct Fermion *even;  
  struct Fermion *odd;
};

struct QX(HalfFermion) {
  struct Q(State) *state;
  size_t size;
  struct Fermion *even;  
};

struct QX(Gauge) {
  struct Q(State) *state;
  size_t size;
  struct SUn *g_data;
  struct Clover *ce_data;
  struct Clover *co_data;
  struct Clover *cox_data;
};

/* allocation routines */
void *qx(allocate_eo)(struct Q(State) *state,
                      size_t *size, void **aligned_ptr,
                      size_t hdr_size, int even_count, int odd_count);
void *qx(step_even)(struct Q(State) *state, void *aligned_ptr);
void *qx(step_odd)(struct Q(State) *state, void *aligned_ptr);

/* data interface routines types */
void qx(x_import)(struct eo_lattice *eo,
                  struct Fermion *data, 
                  double (*reader)(const int pos[Q(DIM)],
                                   int color,
                                   int dirac, 
                                   int re_im,
                                   void *env),
                  void *env);
void qx(x_export)(struct eo_lattice *eo,
                  const struct Fermion *data, 
                  void (*writer)(const int pos[Q(DIM)],
                                 int color,
                                 int dirac, 
                                 int re_im,
                                 double value,
                                 void *env),
                  void *env);

/* Projections */
typedef unsigned int (*qx(Up_project))(struct ProjectedFermion *r,
                                       int size,
                                       const struct up_pack *link,
                                       const struct SUn *U,
                                       const struct Fermion *f);
typedef unsigned int (*qx(Down_project))(struct ProjectedFermion *r,
                                         int size,
                                         const struct down_pack *link,
                                         const struct Fermion *f);
unsigned int qx(proj_g0plus)(struct ProjectedFermion *r,
                             int size,
                             const struct down_pack *link,
                             const struct Fermion *f);
unsigned int qx(proj_g1plus)(struct ProjectedFermion *r,
                             int size,
                             const struct down_pack *link,
                             const struct Fermion *f);
unsigned int qx(proj_g2plus)(struct ProjectedFermion *r,
                             int size,
                             const struct down_pack *link,
                             const struct Fermion *f);
unsigned int qx(proj_g3plus)(struct ProjectedFermion *r,
                             int size,
                             const struct down_pack *link,
                             const struct Fermion *f);
unsigned int qx(proj_g0minus)(struct ProjectedFermion *r,
                              int size,
                              const struct down_pack *link,
                              const struct Fermion *f);
unsigned int qx(proj_g1minus)(struct ProjectedFermion *r,
                              int size,
                              const struct down_pack *link,
                              const struct Fermion *f);
unsigned int qx(proj_g2minus)(struct ProjectedFermion *r,
                              int size,
                              const struct down_pack *link,
                              const struct Fermion *f);
unsigned int qx(proj_g3minus)(struct ProjectedFermion *r,
                              int size,
                              const struct down_pack *link,
                              const struct Fermion *f);
unsigned int qx(proj_Ucg0plus)(struct ProjectedFermion *r,
                               int size,
                               const struct up_pack *link,
                               const struct SUn *U,
                               const struct Fermion *f);
unsigned int qx(proj_Ucg1plus)(struct ProjectedFermion *r,
                               int size,
                               const struct up_pack *link,
                               const struct SUn *U,
                               const struct Fermion *f);
unsigned int qx(proj_Ucg2plus)(struct ProjectedFermion *r,
                               int size,
                               const struct up_pack *link,
                               const struct SUn *U,
                               const struct Fermion *f);
unsigned int qx(proj_Ucg3plus)(struct ProjectedFermion *r,
                               int size,
                               const struct up_pack *link,
                               const struct SUn *U,
                               const struct Fermion *f);
unsigned int qx(proj_Ucg0minus)(struct ProjectedFermion *r,
                                int size,
                                const struct up_pack *link,
                                const struct SUn *U,
                                const struct Fermion *f);
unsigned int qx(proj_Ucg1minus)(struct ProjectedFermion *r,
                                int size,
                                const struct up_pack *link,
                                const struct SUn *U,
                                const struct Fermion *f);
unsigned int qx(proj_Ucg2minus)(struct ProjectedFermion *r,
                                int size,
                                const struct up_pack *link,
                                const struct SUn *U,
                                const struct Fermion *f);
unsigned int qx(proj_Ucg3minus)(struct ProjectedFermion *r,
                                int size,
                                const struct up_pack *link,
                                const struct SUn *U,
                                const struct Fermion *f);

/* projection tables */
/*  normal projection */
extern qx(Up_project) qx(up_project_n)[Q(DIM)];
extern qx(Down_project) qx(down_project_n)[Q(DIM)];
/*  conjugated projection */
extern qx(Up_project) qx(up_project_x)[Q(DIM)];
extern qx(Down_project) qx(down_project_x)[Q(DIM)];
/* compute projections on the boundary and fill the send buffers */
void qx(boundary)(struct eo_lattice *xy,
                  const qx(Up_project) up_proj[],
                  const qx(Down_project) down_proj[],
                  const struct SUn *U,
                  const struct Fermion *src_y,
                  long long *flops);

/* Backend controled structure sizes */
int qx(sizeof_fermion)(int volume);
int qx(sizeof_gauge)(int volume);
int qx(sizeof_clover)(int volume);
int qx(sizeof_vfermion)(int volume, int count);

/* qa0 level data access routines */
void qx(put_fermion)(struct Fermion *data, int pos, const double r[]);
void qx(get_fermion)(double r[], const struct Fermion *data, int pos);
void qx(put_gauge)(struct SUn *ptr, int pos, const double r[]);
void qx(put_clover_lo)(struct Clover *ptr, int pos, const double r[]);
void qx(put_clover_hi)(struct Clover *ptr, int pos, const double r[]);

/* Linear algebra on fermions */
void qx(f_zero)(struct Fermion *dst, 
                int size);
void qx(f_copy)(struct Fermion *dst, 
                int size,
                const struct Fermion *src);
unsigned int qx(f_dot)(double *v_r, double *v_i,
                       int size,
                       const struct Fermion *a,
                       const struct Fermion *b);
unsigned int qx(f_add3)(struct Fermion *r,
                        int size,
                        const struct Fermion *a,
                        double s,
                        const struct Fermion *b);
unsigned int qx(f_add2)(struct Fermion *r,
                        int size,
                        double s,
                        const struct Fermion *b);
unsigned int qx(f_cadd2)(struct Fermion *r,
                         int size,
                         double zr, double zi,
                         const struct Fermion *b);
unsigned int qx(f_add2_norm)(struct Fermion *r,
                             double *local_norm,
                             int size,
                             double s,
                             const struct Fermion *b);
unsigned int qx(f_rmul1)(struct Fermion *r,
                         int size,
                         double s);
unsigned int qx(f_add2x)(struct Fermion *r,
                         int size,
                         double s,
                         const struct Fermion *b);
unsigned int qx(f_norm)(double *s,
                        int size,
                        const struct Fermion *a);
unsigned int qx(f_diff_norm)(double *s,
                             int size,
                             const struct Fermion *a,
                             const struct Fermion *b);

/* algebra for arrays of fermions */

unsigned int qx(fv_zero)(struct vFermion *dst, 
                         int size, int len);

/* fv[fv_begin + (0 .. len-1)] = gv[gv_begin + (0 .. len-1)]
*/
unsigned int qx(fv_copy)(
        int size, int len,
        struct vFermion *fv, int fv_size, int fv_begin,
        const struct vFermion *gv, int gv_size, int gv_begin
        );
/*
 * set fv[idx] = x
*/
unsigned int qx(fv_put)(
        int size,
        struct vFermion *fv, int fv_size, int fv_idx,
        const struct Fermion *x
        );

/*
 * read x = fv[idx]
*/
unsigned int qx(fv_get)(
        int size,
        struct Fermion *x,
        const struct vFermion *fv, int fv_size, int fv_idx
        );

/*
*   g = fv[fv_begin + (0 .. f_vlen-1)] . v
*   v is a complex vector [fv_len] indexed as [re:0/im:1 + 2 * i]
*/
unsigned int qx(fv_dot_zv)(
        int size,
        struct Fermion *g,
        const struct vFermion *fv, int fv_size, int fv_begin, int fv_len,
        const double *v
        );

/*
*   gv[gv_begin + (0 .. gv_len-1)] = fv[fv_begin + (0 .. f_len - 1)] . m
*   m is a complex matrix [fv_len*gv_len] indexed as [re:0/im:1 + 2 * (row + ldm * col) ]
*/
unsigned int qx(fv_dot_zm)(
        int size,
        struct vFermion *gv, int gv_row_size, int gv_begin, int gv_len,
        const struct vFermion *fv, int fv_row_size, int fv_begin, int fv_len,
        const double *m, int ldm
        );

/*  XXX this includes global reduction
 *  c[i] = herm(fv[fv_begin+i]) * g 
 *      for all i = (0 .. fv_len-1)
 *  c is complex vector as [re:0/im:1 + 2 * i]
 */
unsigned int qx(fvH_dot_f)(
        int size,
        double *c,
        const struct vFermion *fv, int fv_size, int fv_begin, int fv_len,
        const struct Fermion *g
        );

/* Local part of the above */
unsigned int qx(do_fvH_dot_f)(
        int size,
        double *c,
        const struct vFermion *fv, int fv_size, int fv_begin, int fv_len,
        const struct Fermion *g);

/* XXX this includes global reduction
 * c[i,j] = herm(fv[fv_begin + i]) . g[gv_begin+j] 
 *      for all i = (0 .. fv_len-1), 
 *              j = (0 .. gv_len-1),
 * c is a complex matrix as [re:0/im:1 + 2 * (i + ldc * j)]
 */
unsigned int qx(fvH_dot_fv)(
        int size,
        double *c, int ldc,
        const struct vFermion *fv, int fv_size, int fv_begin, int fv_len,
        const struct vFermion *gv, int gv_size, int gv_begin, int gv_len
        );

/* local part of the above */
unsigned int qx(do_fvH_dot_fv)(
        int size,
        double *c, int ldc,
        const struct vFermion *fv, int fv_size, int fv_begin, int fv_len,
        const struct vFermion *gv, int gv_size, int gv_begin, int gv_len);

/* basic matrices */
unsigned int qx(op_norm2)(double *global_norm,
                          const struct QX(Fermion) *psi,
                          struct Q(State) *state);
unsigned int qx(do_A)(struct Fermion *r_x,
                      int size,
                      const struct Clover *C,
                      const struct Fermion *s_x);

/* basic A+B, A, B, and their combinations  */
unsigned int qx(do_ApB)(struct Fermion *r_x,
                        int start, int size,
                        const struct neighbor *neighbor,
                        const struct SUn *U,
                        const struct Clover *C,
                        const struct Fermion *s_x,
                        const struct Fermion *s_y,
                        void *rb[]);
unsigned int qx(do_AxpBx)(struct Fermion *r_x,
                          int start, int size,
                          const struct neighbor *neighbor,
                          const struct SUn *U,
                          const struct Clover *C,
                          const struct Fermion *s_x,
                          const struct Fermion *s_y,
                          void *rb[]);
unsigned int qx(do_CmB)(struct Fermion *r_x,
                        int start, int size,
                        const struct neighbor *neighbor,
                        const struct SUn *U,
                        const struct Fermion *s_x,
                        const struct Fermion *s_y,
                        void *rb[]);
unsigned int qx(do_AmB)(struct Fermion *r_x,
                        int start, int size,
                        const struct neighbor *neighbor,
                        const struct SUn *U,
                        const struct Clover *C,
                        const struct Fermion *s_x,
                        const struct Fermion *s_y,
                        void *rb[]);
unsigned int qx(do_AmB_norm)(struct Fermion *r_x,
                             double *local_norm,
                             int start, int size,
                             const struct neighbor *neighbor,
                             const struct SUn *U,
                             const struct Clover *C,
                             const struct Fermion *s_x,
                             const struct Fermion *s_y,
                             void *rb[]);
unsigned int qx(do_AxmBx)(struct Fermion *r_x,
                          int start, int size,
                          const struct neighbor *neighbor,
                          const struct SUn *U,
                          const struct Clover *C,
                          const struct Fermion *s_x,
                          const struct Fermion *s_y,
                          void *rb[]);
unsigned int qx(do_AB)(struct Fermion *r_x,
                       int start, int size,
                       const struct neighbor *neighbor,
                       const struct SUn *U,
                       const struct Clover *C,
                       const struct Fermion *s_y,
                       void *rb[]);
unsigned int qx(do_AxBx)(struct Fermion *r_x,
                         int start, int size,
                         const struct neighbor *neighbor,
                         const struct SUn *U,
                         const struct Clover *C,
                         const struct Fermion *s_y,
                         void *rb[]);

/* even/odd level routines */
void qx(op_CmB)(struct Fermion *res_x,
                struct eo_lattice *r_x,
                const struct SUn *g_data,
                const struct Fermion *a_x,
                const struct Fermion *b_y,
                long long *flops,
                long long *sent,
                long long *received);
void qx(op_A)(struct Fermion *r_x,
              struct eo_lattice *xy,
              const struct Clover *cx_data,
              const struct Fermion *s_x,
              long long *flops);
void qx(op_AB)(struct Fermion *r_x,
               struct eo_lattice *xy,
               const struct SUn *g_data,
               const struct Clover *cx_data,
               const struct Fermion *s_y,
               long long *flops,
               long long *sent,
               long long *received);
void qx(op_AmB)(struct Fermion *r_x,
                struct eo_lattice *xy,
                const struct SUn *g_data,
                const struct Clover *cx_data,
                const struct Fermion *s_x,
                const struct Fermion *s_y,
                long long *flops,
                long long *sent,
                long long *received);
void qx(op_AmB_norm)(struct Fermion *r_x,
                     double *global_norm,
                     struct eo_lattice *xy,
                     const struct SUn *g_data,
                     const struct Clover *cx_data,
                     const struct Fermion *s_x,
                     const struct Fermion *s_y,
                     long long *flops,
                     long long *sent,
                     long long *received);
void qx(op_AxBx)(struct Fermion *r_x,
                 struct eo_lattice *xy,
                 const struct SUn *g_data,
                 const struct Clover *cx_data,
                 const struct Fermion *s_y,
                 long long *flops,
                 long long *sent,
                 long long *received);
void qx(op_AxmBx)(struct Fermion *r_x,
                  struct eo_lattice *xy,
                  const struct SUn *g_data,
                  const struct Clover *cx_data,
                  const struct Fermion *s_x,
                  const struct Fermion *s_y,
                  long long *flops,
                  long long *sent,
                  long long *received);
void qx(op_ApB)(struct Fermion *r_x,
                struct eo_lattice *xy,
                const struct SUn *U,
                const struct Clover *C,
                const struct Fermion *a_x,
                const struct Fermion *a_y,
                long long *flops,
                long long *sent,
                long long *received);
void qx(op_AxpBx)(struct Fermion *r_x,
                  struct eo_lattice *xy,
                  const struct SUn *U,
                  const struct Clover *C,
                  const struct Fermion *a_x,
                  const struct Fermion *a_y,
                  long long *flops,
                  long long *sent,
                  long long *received);
void qx(op_even_M)(struct Fermion *r_x,
                   struct Q(State) *state,
                   const struct QX(Gauge) *gauge,
                   const struct Fermion *a_x,
                   long long *flops,
                   long long *sent,
                   long long *received,
                   struct Fermion *tmp_y);
void qx(op_even_Mn)(struct Fermion *r_x,
                    double *global_norm,
                    struct Q(State) *state,
                    const struct QX(Gauge) *gauge,
                    const struct Fermion *a_x,
                    long long *flops,
                    long long *sent,
                    long long *received,
                    struct Fermion *tmp_y);
void qx(op_even_Mx)(struct Fermion *r_x,
                    struct Q(State) *state,
                    const struct QX(Gauge) *gauge,
                    const struct Fermion *a_x,
                    long long *flops,
                    long long *sent,
                    long long *received,
                    struct Fermion *tmp_y);

/* logging */
void qx(zprint)(struct Q(State) *state,
                const char *source,
                const char *fmt,
                ...);
/* parts of the CG solver */
void qx(cg_precondition)(struct Fermion *chi_e,
                         struct Q(State) *state,
                         const struct QX(Gauge) *gauge,
                         const struct Fermion *eta_e,
                         const struct Fermion *eta_o,
                         long long *flops,
                         long long *sent,
                         long long *received,
                         struct Fermion *t0_e,
                         struct Fermion *t0_o);
void qx(cg_inflate)(struct Fermion *psi_o,
                    struct Q(State) *state,
                    const struct QX(Gauge) *gauge,
                    const struct Fermion *eta_o,
                    const struct Fermion *psi_e,
                    long long *flops,
                    long long *sent,
                    long long *received,
                    struct Fermion *t_o);
double qx(cg_dirac_error)(const struct Fermion *psi_e,
                          const struct Fermion *psi_o,
                          struct Q(State) *state,
                          const struct QX(Gauge) *gauge,
                          const struct Fermion *eta_e,
                          const struct Fermion *eta_o,
                          long long *flops,
                          long long *sent,
                          long long *received,
                          struct Fermion *t0_e,
                          struct Fermion *t0_o);
void qx(cg_log)(double cg_res, const char *source, int iter,
                const struct Fermion *xi_e,
                struct Q(State) *state,
                const struct QX(Gauge) *gauge,
                const struct Fermion *chi_e,
                long long *flops,
                long long *sent,
                long long *received,
                unsigned int options,
                struct Fermion *t0_e,
                struct Fermion *t1_e,
                struct Fermion *t0_o,
                struct Fermion *t1_o);

struct MxM_workspace {
    struct Q(State)        *state;
    const struct QX(Gauge) *gauge;
    struct Fermion         *tmp_e;
    struct Fermion         *tmp_o;
    long long              *flops;
    long long              *sent;
    long long              *received;
};

void qx(cg_operator)(struct Fermion            *res_e,
                     const struct Fermion      *psi_e,
                     struct MxM_workspace      *ws);

CG_STATUS qx(cg_solver)(struct Fermion *psi_e,
                        const char *source,
                        int *out_iter,
                        double *out_epsilon,
                        struct Q(State) *state,
                        const struct QX(Gauge) *gauge,
                        const struct Fermion *chi_e,
                        struct Q(Deflator) *deflator,
                        int max_iter,
                        double epsilon,
                        unsigned options,
                        long long *flops,
                        long long *sent,
                        long long *received,
                        struct Fermion *rho_e,
                        struct Fermion *pi_e,
                        struct Fermion *zeta_e,
                        struct Fermion *t0_e,
                        struct Fermion *t1_e,
                        struct Fermion *t0_o,
                        struct Fermion *t1_o);

/*
 *  compute x <- x + alpha p
 *          p <- r + beta p
 */
unsigned int qx(cg_xp)(struct Fermion *x,
                       struct Fermion *p,
                       int size,
                       double alpha,
                       double beta,
                       const struct Fermion *r);

