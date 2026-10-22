all: reading_history.pdf flatvertex.pdf

reading_history.pdf: reading_history.tex dragon_folds_c.png dragon_hasse.png
	pdflatex -interaction=nonstopmode reading_history.tex >/dev/null
	pdflatex -interaction=nonstopmode reading_history.tex >/dev/null

flatvertex.pdf: flatvertex.tex
	pdflatex -interaction=nonstopmode flatvertex.tex >/dev/null
	pdflatex -interaction=nonstopmode flatvertex.tex >/dev/null

clean:
	rm -f *.aux *.log *.out *.toc
