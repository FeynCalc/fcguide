# FeynCalc guide

[FeynCalc](http://www.feyncalc.org/) is a Mathematica package for symbolic evaluation of Feynman diagrams and algebraic calculations in quantum field theory and elementary particle physics.

This guide is a hands-on, practical introduction to automatized 
quantum field theory calculations with FeynCalc and Mathematica.

# License

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">
<img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" />
</a>

This work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">
Creative Commons Attribution-ShareAlike 4.0 International License
</a>.

Copyright (C) 1990-2015 Rolf Mertig,

Copyright (C) 1997-2015 Frederik Orellana,

Copyright (C) 2014-2015 Vladyslav Shtabovenko.

# Compiling the LaTeX source

## Prerequisites

* Mathematica TeX style files and fonts - get them here: http://feyncalc.org/download/texmf.tar.gz
* A recent LaTeX installation (see http://latex-project.org/ftp.html) with 

## Mathematica styles and fonts

### MacOS X:

Run the following as root - modifying the path to match your Mathematica and TeX installation.

```
	curl -o texmf.tar.gz http://feyncalc.org/download/texmf.tar.gz
	tar -xvzf texmf.tar.gz
	sudo tar -cvzf /usr/local/texlive.orig.tar.gz /usr/local/texlive
	sudo rsync --recursive texmf/ /usr/local/texlive/texmf-local/
	sudo cp -i texmf/dvips/config/wolfram.map /usr/local/texlive/*/texmf-var/fonts/map/dvips/updmap/
	sudo mktexlsr
	updmap-sys --enable Map=wolfram.map
```

### Linux:

```
	curl -o texmf.tar.gz http://feyncalc.org/download/texmf.tar.gz
	tar -xvzf texmf.tar.gz
	cp -a texmf ~/.texmf-var
	mkdir ~/.texmf-var/fonts/map
	cp -a ~/.texmf-var/pdftex ~/.texmf-var/fonts/map/
	cp -a ~/.texmf-var/dvips ~/.texmf-var/fonts/map/
	echo "Map wolfram.map" > ~/.updmap.cfg
	echo "Map pdfwolfram.map" >> ~/.updmap.cfg
	updmap --cnffile ~/.updmap.cfg
	rm ~/.updmap.cfg
```

* On Fedora 20 please use `~/.texlive2013/texmf-var` instead of `~/.texmf-var`

### Windows:

--- not tested ---
The texmf directory with Mathematica TeX style files should be put in the search path of the TeX installation.

## LaTeX invocation

### Using PDFLaTex

	pdflatex FCGuide.tex

### Using Latexmk

	latexmk -pdfdvi FCGuide.tex 

### Using Latex
	latex FCGuide.tex
	dvips FCGuide.dvi
	ps2pdf FCGuide.ps

