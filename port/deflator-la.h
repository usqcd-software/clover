#ifndef DEFLATOR_LA_H_hbq6DwyKpy5jzcxzbemt
#define DEFLATOR_LA_H_hbq6DwyKpy5jzcxzbemt

#include <assert.h>
#define NOT_IMPLEMENTED     assert(NULL == "implement me first!")


#if defined(HAVE_LAPACK)
#  include <f2c_types.h>
#  include <lapack.h>
#elif defined(HAVE_GSL)
#  define CHECK_GSL_STATUS(func) do { int status__ = func; assert(0 == status__); } while(0)
#  include <gsl/gsl_vector.h>
#  include <gsl/gsl_matrix.h>
    typedef struct {
        double r, i;
    } doublecomplex;
#else
#  error "no linear algebra library"
#endif 


/***************************/
/*  latvector definitions  */
/***************************/

struct FermionF;
typedef struct {
    int dim;
    struct FermionF *f;
} latvec_c;

#define latvec_c_null(p) { (p)->f = NULL; }
#define latvec_c_is_null(p) (NULL == (p)->f)

latvec_c latvec_c_alloc(struct Q(State) *state, int dim);
latvec_c latvec_c_view(int dim, struct FermionF *f);
void latvec_c_copy(latvec_c x, latvec_c y); /* y <- x */
void latvec_c_zero(latvec_c x); /* x <- 0 */
void latvec_c_free(struct Q(State) *state, latvec_c *v);

doublecomplex lat_c_dotu(latvec_c x, latvec_c y);
void lat_c_scal_d(double alpha, latvec_c x);
void lat_c_axpy_d(double alpha, latvec_c x, latvec_c y);
double lat_c_nrm2(latvec_c x);

struct FermionD;
typedef struct {
    int dim;
    struct FermionD *f;
} latvec_z;

#define latvec_z_null(p) { (p)->f = NULL; }
#define latvec_z_is_null(p) (NULL == (p)->f)

latvec_z latvec_z_alloc(struct Q(State) *state, int dim);
void latvec_z_free(struct Q(State) *state, latvec_z *v);
latvec_z latvec_z_view(int dim, struct FermionD *f);
#if 0
void latvec_z_copy(latvec_z x, latvec_z y); /* y <- x */

/* copy + conversion: y <- x */
void latvec_cz_copy(latvec_c x, latvec_z y);
void latvec_zc_copy(latvec_z x, latvec_c y);


doublecomplex lat_z_dotu(latvec_z x, latvec_z y);
doublecomplex lat_cz_dotu(latvec_c x, latvec_z y);
void lat_c_scal(doublecomplex alpha, latvec_c x);
void lat_c_axpy(doublecomplex alpha, latvec_c x, latvec_c y);
void lat_cz_axpy(doublecomplex alpha, latvec_c x, latvec_z y);
void lat_z_axpy(doublecomplex alpha, latvec_z x, latvec_z y);
void lat_cz_axpy_d(double alpha, latvec_c x, latvec_z y);
void lat_z_axpy_d(double alpha, latvec_z x, latvec_z y);
double lat_z_nrm2(latvec_z x);
#endif


/*******************************/
/*  latmatrix definitions      */
/*******************************/
struct vFermion;
typedef struct {
    int dim;
    int size;
    int begin;
    int len;
    struct vFermion *fv;
} latmat_c;

#define latmat_c_null(p) { (p)->fv = NULL; }
#define latmat_c_is_null(p) (NULL == (p)->fv)

latmat_c latmat_c_alloc(struct Q(State) *state, int dim, int ncol);
void latmat_c_free(struct Q(State) *state, latmat_c *m);
latmat_c latmat_c_view(int dim, int size, struct vFermion *fv);
void latmat_c_copy(latmat_c m1, latmat_c m2); /* m2 <- m1 */
/* create a submatrix of subset of columns
   only a 'view' is created, no allocation is performed
   do not try to 'free' submatrix: it may result in memory error */
latmat_c latmat_c_submat_col(latmat_c m, int col, int ncol);
void latmat_c_insert_col(latmat_c m, int col, latvec_c v);
void latmat_c_get_col(latmat_c m, int col, latvec_c v);


/* use fortran/BLAS conventions for matrices as arrays of column vectors:
   row index runs fastest; function semantics is chosen close to BLAS, with
   known dimension (lattice vector size) omitted
 */
/* C <- A^\dag * B, A:lat*m, B:lat*n, C:m*n */
void lat_lmH_dot_lm(int m, int n,
                    latmat_c a, 
                    latmat_c b, 
                    doublecomplex *c, int ldc);
/* y <- A^\dag x, A:lat*m, x:lat, y:m */
void lat_lmH_dot_lv(int m,  
                    latmat_c a, 
                    latvec_c x, 
                    doublecomplex *y);
/* C <- A * B, A:lat*k, B:k*n, C:lat*n */
void lat_lm_dot_zm(int n, int k, 
                   latmat_c a,
                   doublecomplex *b, int ldb, 
                   latmat_c c);
/* y <- A * x, A:lat*n, x:n, y:lat */
void lat_lm_dot_zv(int n, 
                   latmat_c a, 
                   doublecomplex *x,
                   latvec_c y);


/* lin. operator tie-back */
void latvec_c_linop(latvec_c y, latvec_c x, latvec_c aux);
//void latvec_z_linop();

#endif/*DEFLATOR_LA_H_hbq6DwyKpy5jzcxzbemt*/
