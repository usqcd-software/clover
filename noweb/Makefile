NOWEAVE  = noweave -filter btdefn -index -latex -delay
NOTANGLE = notangle -filter btdefn
NOINDEX  = noindex
CPIF     = cpif
MZSCHEME = mzscheme
MZC      = mzc
INDENT   = indent
LATEX    = latex
DVIPDFM  = dvipdfm
#TFLAGS   = -L'\#line %L "%F"%N'

noweb = clover.nw \
        qcd-defs.nw

sources = ../port/qop-clover.h

f.tex = $(noweb:%.nw=%.tex)
f.dvi = $(noweb:%.nw=%.dvi)
f.pdf = $(noweb:%.nw=../doc/%$(NC).pdf)

.PHONY: all clean realclean 

all: pdf sources

clean:
	$(RM) $(noweb:%.nw=%.aux) $(noweb:%.nw=%.dvi) $(noweb:%.nw=%.log) \
              $(noweb:%.nw=%.nwi) $(noweb:%.nw=%.out) $(noweb:%.nw=%.tex) \
              $(noweb:%.nw=%.toc)

realclean: clean
	$(RM) $(f.pdf) $(sources)

sources: $(sources)

pdf: $(f.pdf)

###
$(f.pdf): ../doc/%$(NC).pdf: %.dvi
	$(DVIPDFM) -o $@ $<

$(f.dvi): %.dvi: %.tex
	$(LATEX) $<
	$(NOINDEX) $(<:%.tex=%)
	$(LATEX) $<
	$(NOINDEX) $(<:%.tex=%)
	$(LATEX) $<

$(f.tex): %.tex: %.nw
	$(NOWEAVE) $< > $@

$(sources): clover.nw
	$(NOTANGLE) $(TFLAGS) -R'File [[$@]]' -c $< | $(CPIF) $@
