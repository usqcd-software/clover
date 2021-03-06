TOP     = ..
V       = 0
CONFIG  = macosx
TARGETS = c99-32 c99-64 cee-32 cee-64 bgl/xlc bgq/xlc dry
TARGET  =
DEST    =
QA0     = qa0

include $(TOP)/config/$(CONFIG)

.PHONY: all clean realclean dist

ifeq ("$V", "0")
   E=@echo "  "
   C=@
else
   E=@echo > /dev/null
   C=
endif

DEST=$(TOP)/$(shell echo $(TARGET) | sed -e sx/x-xg)

headers = types.qa0 \
          defs-A.qa0 \
          defs-B.qa0 \
          defs-spin.qa0 \

sources = sizeof-down-pack.qa0 \
          sizeof-up-pack.qa0 \
          sizeof-neighbor.qa0 \
          put-neighbor.qa0 \
          get-neighbor.qa0 \
          put-up.qa0 \
          get-up.qa0 \
          get-up-f.qa0 \
          fix-up.qa0 \
          put-down.qa0 \
          get-down.qa0 \
          get-down-f.qa0 \
          fix-down.qa0 \
          f-d-eq-d-plus-f.qa0 \
          f-f-eq-d-minus-d.qa0 \
          g-f-eq-d.qa0 \
          c-f-eq-d.qa0 \

x.sources = f-dot.qa0 \
            f-add3.qa0 \
            f-norm.qa0 \
            f-zero.qa0 \
            f-cadd2.qa0 \
            f-rmul.qa0 \
            get-fermion.qa0 \
            put-fermion.qa0 \
            sizeof-fermion.qa0 \
            put-gauge.qa0 \
            sizeof-gauge.qa0 \
            put-clover-lo.qa0 \
            put-clover-hi.qa0 \
            sizeof-clover.qa0 \
            proj-minus.qa0 \
            proj-plus.qa0 \
            proj-u-minus.qa0 \
            proj-u-plus.qa0 \
            do-A.qa0 \
            do-A-times-B.qa0 \
            do-A-conj-times-B-conj.qa0 \
            do-A-plus-B.qa0 \
            do-A-minus-B.qa0 \
            do-A-minus-B-norm.qa0 \
            do-C-minus-B.qa0 \
            do-A-conj-plus-B-conj.qa0 \
            do-A-conj-minus-B-conj.qa0 \
            cg-xp.qa0 \
            f-add2.qa0 \
            f-add2-norm.qa0 \
            f-copy.qa0 \
            f-diff-norm.qa0 \
            sizeof-vfermion.qa0 \
            fv-zero.qa0 \
            fv-copy.qa0 \
            put-vfermion.qa0 \
            get-vfermion.qa0 \
            fv-mul-zv.qa0 \
            fv-mul-zm.qa0 \
            do-fvH-dot-f.qa0 \
            do-fvH-dot-fv.qa0

# Who include defs-A.qa0
A.sources = do-A.qa0 \
            do-A-times-B.qa0 \
            do-A-conj-times-B-conj.qa0 \
            do-A-plus-B.qa0 \
            do-A-minus-B.qa0 \
            do-A-minus-B-norm.qa0 \
            do-A-conj-plus-B-conj.qa0 \
            do-A-conj-minus-B-conj.qa0 \

# Who include defs-F.qa0 (they also need defs-spin.qa0)
B.sources = do-A-times-B.qa0 \
            do-A-conj-times-B-conj.qa0 \
            do-A-plus-B.qa0 \
            do-A-minus-B.qa0 \
            do-A-minus-B-norm.qa0 \
            do-C-minus-B.qa0 \
            do-A-conj-plus-B-conj.qa0 \
            do-A-conj-minus-B-conj.qa0 \

# Who include defs-spin.qa0
spin.sources = $(B.sources) \
               proj-plus.qa0 \
               proj-u-plus.qa0 \
               proj-minus.qa0 \
               proj-u-minus.qa0

# Convert qa0 to C names
dst.sources = $(sources:%.qa0=$(DEST)/%.c) \
              $(x.sources:%.qa0=$(DEST)/%d.c) \
              $(x.sources:%.qa0=$(DEST)/%f.c)

c.A.sources = $(A.sources:%.qa0=$(DEST)/%d.c) \
              $(A.sources:%.qa0=$(DEST)/%f.c)

c.B.sources = $(B.sources:%.qa0=$(DEST)/%d.c) \
              $(B.sources:%.qa0=$(DEST)/%f.c)

c.spin.sources = $(spin.sources:%.qa0=$(DEST)/%d.c) \
                 $(spin.sources:%.qa0=$(DEST)/%f.c)

.PHONY: do.all

dist all:
	@for t in $(TARGETS) ; do \
	  echo making target $$t; \
	  $(MAKE) "V=$(V)" CONFIG="$(CONFIG)" TARGET="$$t" do.all ; \
	done

do.all: $(dst.sources)

clean:

.PHONY: do.realclean

realclean: clean
	@for t in $(TARGETS); do \
	  echo making realclean $$t; \
	  $(MAKE) "V=$(V)" CONFIG="$(CONFIG)" TARGET="$$t" do.realclean ; \
	done

do.realclean:
	$E RM $(DEST)/sources
	$C $(RM) $(dst.sources)

$(sources:%.qa0=$(DEST)/%.c): $(DEST)/%.c: %.qa0
	$E QA0 $@
	$C Nc=$(NC) $(QA0) $(TARGET) double $@ $<

$(x.sources:%.qa0=$(DEST)/%d.c): $(DEST)/%d.c: %.qa0
	$E QA0-d $@
	$C Nc=$(NC) $(QA0) $(TARGET) double $@ $<

$(x.sources:%.qa0=$(DEST)/%f.c): $(DEST)/%f.c: %.qa0
	$E QA0-f $@
	$C Nc=$(NC) $(QA0) $(TARGET) float $@ $<

$(dst.sources): types.qa0
$(c.A.sources): defs-A.qa0
$(c.B.sources): defs-B.qa0
$(c.spin.sources): defs-spin.qa0

