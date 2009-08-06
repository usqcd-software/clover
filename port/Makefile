#
CONFIG     = macosx
LIBRARY    = libqop-clover.a
LIMPDIR    = .

QMP_CFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)
lib        = $(LIMPDIR)/$(LIBRARY)
HDR        = qop-clover.h

include ../config/$(CONFIG)

V          = 0

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

.PHONY: all clean realclean library dist

x.sources = version.c \
            error.c \
            perf.c \
            step-even.c \
            step-odd.c \
            alloc-aligned.c \
            alloc.c \
            dealloc.c \
            init.c \
            fini.c \
            comm-setup.c \
            comm-free.c

p.sources = fermion-dot.c \
            fermion-madd.c \
            fermion-norm2.c \
            fermion-alloc.c \
            fermion-free.c \
            fermion-import.c \
            fermion-export.c \
            hfermion-dot.c \
            hfermion-madd.c \
            hfermion-norm2.c  \
            hfermion-alloc.c \
            hfermion-free.c \
            hfermion-import.c \
            hfermion-export.c \
            alloc-eo.c \
            x-import.c \
            x-export.c \
            gauge-import.c \
            gauge-free.c \
            op-dirac.c \
            op-dirac-conj.c \
            op-apb.c \
            op-axpbx.c \
            op-boundary.c \
            boundary-n.c \
            boundary-x.c \
            solver-dirac.c \
            op-print.c \
            cg-dirac-error.c \
            cg-inflate.c \
	    cg-precondition.c \
            cg-solver.c \
            cg-log.c \
            op-m.c \
            op-mn.c \
            op-mx.c \
            op-a1b.c \
            op-a1xbx.c \
            op-amb.c \
            op-axmbx.c \
            op-ambn.c \
            op-cmb.c \
            op-a.c \

xxxx.p.sources =  \
            op-ax.c \
	    op-a1.c \
            op-a1x.c \
            op-b.c \
            op-bx.c \
	    op-b1.c \
            op-b1x.c \
            op-f.c \
            op-fx.c \
            op-a1xbx.c \
            op-1mf.c \
            op-1mfx.c \
            op-1mfxn.c \
            op-1mba1f.c \
            op-1mba1fn.c \
            op-a1xbxfx.c \
            op-ba1f.c \
            op-ba1.c \
            op-mxn.c \
            op-precond.c \
            op-precond-conj.c \
            scg-solver.c \
            solver-mxm.c \
            solver-shifted.c


x.objects = $(x.sources:%.c=$(LIMPDIR)/%.o)

p.objects = $(p.sources:%.c=$(LIMPDIR)/%f.o) \
            $(p.sources:%.c=$(LIMPDIR)/%d.o)

objects = $(x.objects) $(p.objects)

headers = clover.h \
          $(HDR)


all: $(lib)
ifeq ("$(LIMPDIR)", ".")
	$E $(HDR) in place
else
	$E CP $(HDR)
	$C cp $(HDR) $(LIMPDIR)/
endif

clean:
	$E RM objects
	$C$(RM) $(objects)

realclean: clean
ifneq ("$(LIMPDIR)", ".")
	$C$(RM) $(LIMPDIR)/$(HDR)
endif


dist: clean

$(lib): library

library: $(objects)
	$E AR $(lib)
	$C$(AR) cr $(lib) $^
	$C$(RANLIB) $(lib)

$(x.objects): $(LIMPDIR)/%.o: %.c
	$E CC $<
	$C$(CC) $(CFLAGS) -I. -c -o $@ $<

$(p.sources:%.c=$(LIMPDIR)/%f.o): $(LIMPDIR)/%f.o: %.c
	$E CC-f $<
	$C$(CC) $(CFLAGS) -I. -DQOP_CLOVER_DEFAULT_PRECISION=\'F\' -c -o $@ $<

$(p.sources:%.c=$(LIMPDIR)/%d.o): $(LIMPDIR)/%d.o: %.c
	$E CC-d $<
	$C$(CC) $(CFLAGS) -I. -DQOP_CLOVER_DEFAULT_PRECISION=\'D\' -c -o $@ $<

$(objects): $(headers)
