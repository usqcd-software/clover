#include "clover-test.h"

int
read_gauge(QDP_ColorMatrix *link[], const char *name)
{
  QDP_String *id = QDP_string_create();
  QDP_Reader *reader;

  reader = QDP_open_read(id, (char *)name);
  if (reader == 0) {
    printf0("open_read(%s) failed\n", name);
    return 1;
  }
  printf0("---------------------------------------------------------------\n");
  printf0("file metadata(%s):\n", name);
  printf0("%s\n", QDP_string_ptr(id));
  printf0("---------------------------------------------------------------\n");
  if( QDP_vread_M(reader, id, link, NDIM) ) {
    printf0("vread_M(%s) failed\n", name);
    return 1;
  }
  printf0("---------------------------------------------------------------\n");
  printf0("record metadata(%s):\n", name);
  printf0("%s\n", QDP_string_ptr(id));
  printf0("---------------------------------------------------------------\n");
  QDP_close_read(reader);
  QDP_string_destroy(id);
  return 0;  
}

void
coord_gauge(QDP_ColorMatrix *U[])
{
    int d, sites, i, p[NDIM], a, b, j;
    QLA_Real v;
    QLA_Complex w;
    QLA_ColorMatrix *qu;

    sites = QDP_numsites(QDP_this_node);
    for (d = 0; d < NDIM; d++) {
        qu = QDP_expose_M(U[d]);
        for (i = 0; i < sites; i++) {
            QDP_get_coords(p, QDP_this_node, i);
            v = 0.0;
            for (j = 0; j < NDIM; j++)
                v = v * 10 + p[j];
            v = v * 10 + d;
            for (a = 0; a < QDP_Nc; a++) {
                for (b = 0; b < QDP_Nc; b++) {
                    QLA_Real z = v * 1000 + a * 100 + b * 10;
                    QLA_c_eq_r_plus_ir(w, z, z + 1);
                }
            }
        }
        QDP_reset_M(U[d]);
    }
}

void
unit_gauge(QDP_ColorMatrix *U[])
{
    int d;
    QLA_Complex w;

    QLA_c_eq_r(w, 1.0);
    for (d = 0; d < NDIM; d++) {
        QDP_M_eq_c(U[d], &w, QDP_all);
    }
}
