# XXX This file still needs editing.
# common makefile for standard targets
TOP        = ..
V          = 0
LIBRARY    = libqop-clover.a
QMP_CFLAGS = $(shell $(QMP_TOP:%=%/bin/qmp-config) --cflags)

include $(TOP)/config/$(CONFIG)

.PHONY: all clean realclean dist

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

headers = ../../port/mdwf.h

i.sources = sizeof-down-pack.c \
            sizeof-neighbor.c \
            sizeof-up-pack.c \
            put-neighbor.c \
            get-neighbor.c \
            put-up.c \
            get-up.c \
            get-up-f.c \
            fix-up.c \
            put-down.c \
            get-down.c \
            fix-down.c \
            get-down-f.c \

xxxx.i.sources = put-ab.c \
            put-abi.c \
            put-abi-z.c \
            sizeof-ab-table.c \
            sizeof-abi-table.c \

x.sources = f-norm \
            f-add3 \
            f-dot \
            get-fermion \
            put-fermion \
            sizeof-fermion \
            put-gauge \
            sizeof-gauge \
            put-clover-lo \
            put-clover-hi \
            sizeof-clover \
            proj-minus \
            proj-plus \
            proj-u-minus \
            proj-u-plus \
            do-A \
            do-A-times-B \
            do-A-conj-times-B-conj \
            do-A-plus-B \
            do-A-minus-B \
            do-A-minus-B-norm \
            do-C-minus-B \
            do-A-conj-plus-B-conj \
            do-A-conj-minus-B-conj \
            cg-xp \
            f-add2 \
            f-add2-norm \
            f-copy \
            f-diff-norm \

xxxx.x.sources = \
            sizeof-pfermion \
            fv-zero \
            fv-copy \
            fv-get \
            fv-put \
            f-zero \
            f-add2x \
            scg-madd \
            scg-xp \
            do-A \
            do-A-conj \
            do-A-inv \
            do-A-inv-conj \
            do-F \
            do-F-conj \
            do-B-A-inv \
            do-B-A-inv-F \
            do-1-sub-B-A-inv-F \
            do-1-sub-F \
            do-1-sub-F-conj \
            do-1-sub-F-conj-norm \
            do-1-sub-B-A-inv-F-norm \
            do-A-inv-conj-B-conj \
            do-A-inv-conj-B-conj-F-conj \

sources = $(i.sources) \
          $(x.sources:%=%f.c) \
          $(x.sources:%=%d.c)

objects = $(sources:%.c=%.o)

all: $(objects)
	$E AR $@/$(LIMP)
	$C$(AR) cr $(LIBRARY) $^
	$C$(RANLIB) $(LIBRARY)
	$C$(MAKE) CONFIG='$(CONFIG)' V='$V' LIBRARY='$(LIBRARY)' \
		LIMPDIR=../$(LIMP) -C ../port $@

dist clean:
	$E RM $(LIMP)/objects
	$C$(RM) $(objects)
	$C$(MAKE) CONFIG='$(CONFIG)' V='$V' LIBRARY='$(LIBRARY)' \
		LIMPDIR=../$(LIMP) -C ../port $@

realclean: clean
	$E RM $(LIMP)/sources
	$C$(RM) $(sources)
	$C$(MAKE) CONFIG='$(CONFIG)' V='$V' LIBRARY='$(LIBRARY)' \
		LIMPDIR=../$(LIMP) -C ../port $@
	$C$(RM) $(LIBRARY)

$(sources:%.c=%.o): %.o: %.c
	$E CC $<
	$C$(CC) $(CFLAGS) -I../port -c -o $@ $<

