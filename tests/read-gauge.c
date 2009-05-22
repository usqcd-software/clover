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
