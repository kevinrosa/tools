################### xelatex -ini -jobname="preamble" "&xelatex" mylatexformat.ltx """preamble.tex"""
xelatex -ini -jobname="preamble-memoir" "&xelatex preamble-memoir.tex\dump"
xelatex -ini -jobname="preamble-beamer" "&xelatex preamble-beamer.tex\dump"
pdflatex -ini -jobname="preamble-ametsoc" "&pdflatex preamble-ametsoc.tex\dump"
