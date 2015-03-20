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

* A Mathematica installation with TeX style files and fonts (<= 9.x)
* A recent LaTeX installation (see http://latex-project.org/ftp.html) with 

## Mathematica styles and fonts

### MacOS X:

Run the following as root - modifying the path to match your Mathematica and TeX installation.

	tar -cvzf /usr/local/texlive.orig.tar.gz /usr/local/texlive
	rsync --recursive /Applications/Mathematica.app/SystemFiles/IncludeFiles/TeX/texmf/ /usr/local/texlive/texmf-local/
	cp -i /Applications/Mathematica.app/SystemFiles/IncludeFiles/TeX/texmf/dvips/config/wolfram.map 2014/texmf-var/fonts/map/dvips/updmap/
	mktexlsr
	updmap-sys --enable Map=wolfram.map

### Linux:

In the following we assume that your Mathematica installation is located in 
`/opt/Mathematica/9.0/`. If not, please modify the path accordingly.

	cp -a /opt/Mathematica/9.0/SystemFiles/IncludeFiles/TeX/texmf ~/.texmf-var
	mkdir ~/.texmf-var/fonts/map
	cp -a ~/.texmf-var/pdftex ~/.texmf-var/fonts/map/
	cp -a ~/.texmf-var/dvips ~/.texmf-var/fonts/map/
	echo "Map wolfram.map" > ~/.updmap.cfg
	echo "Map pdfwolfram.map" >> ~/.updmap.cfg
	updmap --cnffile ~/.updmap.cfg
	rm ~/.updmap.cfg

* On Fedora 20 please use `~/.texlive2013/texmf-var` instead of `~/.texmf-var`

### Windows:

--- not tested ---
The texmf directory of Mathematica should be put in the search path of the TeX installation.

## LaTeX invocation

### Using PDFLaTex

	pdflatex FCGuide.tex

### Using Latexmk

	latexmk -pdfdvi FCGuide.tex 

### Using Latex
	latex FCGuide.tex
	dvips FCGuide.dvi
	ps2pdf FCGuide.ps

