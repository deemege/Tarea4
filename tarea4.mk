
Resultados_hw4.pdf : Resultados_hw4.tex cuerda.png cuerdaPerturbada.png
    pdflatex Resultados_hw4.tex


cuerda.png : cuerda.txt
    Plots.py

cuerdaPerturbada.png : cuerdaPerturbada.txt
    Plots.py
