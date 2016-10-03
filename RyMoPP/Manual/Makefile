FILE=Manual
Manual.pdf: ${FILE}.tex ${FILE}.aux
	pdflatex ${FILE}
	bibtex ${FILE}
	pdflatex ${FILE}
	pdflatex ${FILE}
open: 
	open ${FILE}.pdf
clean:
	rm -f ${FILE}.{ps,log,aux,out,dvi,bbl,blg}
#epj_coll.bbl: epj_coll.bib
#	latex ${@:.bbl=}
#	bibtex ${@:.bbl=}
