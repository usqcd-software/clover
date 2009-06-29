TOP = ..
V = 0
CONFIG = ../target

QMP_CFLAGS  = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)
QMP_LDFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --ldflags)
QMP_LIBS    = $(shell $(QMP_TOP:%=%/bin/qmp-config) --libs)

QDP_CFLAGS  = $(shell $(QDP_TOP:%=%/bin/qdp-config) --cflags)
QDP_LDFLAGS = $(shell $(QDP_TOP:%=%/bin/qdp-config) --ldflags)
QDP_LIBS    = $(shell $(QDP_TOP:%=%/bin/qdp-config) --libs)

TCFLAGS  = $(CFLAGS) $(QMP_CFLAGS) $(QDP_CFLAGS)
TLDFLAGS = $(LDFLAGS) $(QMP_LDFLAGS) $(QDP_LDFLAGS)
TCLIBS   = $(QDP_LIBS) $(QMP_LIBS) $(CLIBS)

include $(TOP)/config/$(CONFIG)

.PHONY: all clean realclean dist

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null " "
   C=
endif

programs = get-gauge \
           check-gamma

utils = gauge.c \
        fermion.c \
        clover.c

headers = clover-test.h

sources = $(programs:%=%.c) $(utils)

objects = $(sources:%.c=%.o)

###############################################################################
all: $(programs)

clean:
	$E RM objects
	$C$(RM) $(objects)

realclean: clean
	$E RM programs
	$C$(RM) $(programs)

dist: realclean

###############################################################################

$(programs):%: %.o
	$E LD $@
	$C$(CC) $(TLDFLAGS) -o $@ $^ $(TCLIBS)

$(programs): $(utils:%.c=%.o)

$(objects): %.o: %.c
	$E CC $@
	$C$(CC) $(TCFLAGS) -c -o $@ $<

$(objects): $(headers)
