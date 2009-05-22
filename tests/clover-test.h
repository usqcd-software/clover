#define QDP_Precision 'D'
#define QDP_Nc 3
#include <qdp.h>
#include <qmp.h>

#define NDIM  4
#define NELEMS(x) (sizeof (x) / sizeof ((x)[0]))

extern int lattice[NDIM];

#define printf0 if (QDP_this_node == 0) printf

extern void create_Mvector(QDP_ColorMatrix *U[], int size);
extern void destroy_Mvector(QDP_ColorMatrix *U[], int size);

extern QLA_Real plaquette(QDP_ColorMatrix *U[]);
extern int read_gauge(QDP_ColorMatrix *U[], const char *name);
