FILE=Manual
Manual.pdf: ${FILE}.tex ${FILE}.aux
	pdflatex ${FILE}.tex
	bibtex ${FILE}.aux
	pdflatex ${FILE}.tex
	pdflatex ${FILE}.tex
open: 
	open ${FILE}.pdf
clean:
	rm -f ${FILE}.{ps,log,aux,out,dvi,bbl,blg}
#epj_coll.bbl: epj_coll.bib
#	latex ${@:.bbl=}
#	bibtex ${@:.bbl=}
