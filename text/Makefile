FIGSGNU=$(shell find figures/*.gnu -type f)
FIGSEPS=$(subst gnu,eps,$(FIGSGNU))

all: manus.pdf clean

manus.pdf: fig manus.tex manus.bib
	latexmk -ps -pdfps -use-make manus.tex

diff: diff_manus.tex
	latex  diff_manus.tex &&\
	latex  diff_manus.tex &&\
	bibtex diff_manus     &&\
	latex  diff_manus.tex &&\
	latex  diff_manus.tex &&\
	dvips  diff_manus.dvi &&\
	ps2pdf diff_manus.ps

diff_manus.tex: manus.tex manus.bib $(FIGSEPS)
	latexdiff --config="PICTUREENV=(?:picture|DIFnomarkup|comment|longtable|table)[\w\d*@]*" foo.tex manus.tex > diff_manus.tex

manus.bib: manus.tex $(shell find bib/*.bib -type f)
	cat bib/*.bib >$@

bib:
	rm -f manus.bib manus.bbl manus.blg manus.aux
	make manus.bib

clean:
	rm -f manus.{ps,dvi,aux,log,out,blg,bbl,toc,tdo,tex.bak,bib,synctex.gz}
	rm -f diff_manus.{ps,dvi,aux,log,out,blg,bbl,toc,tdo}
	latexmk -c

fig:
	cd figures &&\
	gnuplot fig_alignscore.gnu && \
	gnuplot fig_cost.gnu && \
	gnuplot fig_time.gnu

figures/%.eps: figures/%.gnu
	cd figures && gnuplot ../$<

cover: cover_letter.pdf

cover_letter.pdf: cover_letter.tex
	latex cover_letter
	latex cover_letter
	dvips cover_letter
	ps2pdf cover_letter.ps

.PHONY: all bib clean cover submit
