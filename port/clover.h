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

/* Deflator state */
struct Q(Deflator) {
    struct Q(State) *state;
    int nev;
    int vsize;
    double eps;
    /* XXX other pieces of the deflation state */
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
# if QOP_CLOVER_DEFAULT_PRECISION=='D'
#  define qx(x) qop_d3_clover_##x
#  define QX(x) QOP_D3_CLOVER_##x
#  define REAL double
#  define SUn SUnD
#  define Clover CloverD
#  define Fermion FermionD
#  define ProjectedFermion ProjectedFermionD
# endif
# if QOP_CLOVER_DEFAULT_PRECISION=='F'
#  define qx(x) qop_f3_clover_##x
#  define QX(x) QOP_F3_CLOVER_##x
#  define REAL float
#  define SUn SUnF
#  define Clover CloverF
#  define Fermion FermionF
#  define ProjectedFermion ProjectedFermionF
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
unsigned int qx(f_add2_norm)(struct Fermion *r,
                             double *local_norm,
                             int size,
                             double s,
                             const struct Fermion *b);
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
int qx(cg_solver)(struct Fermion *psi_e,
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

