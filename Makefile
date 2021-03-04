## This is Ali's SIRmodel

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

Sources += codes/ali_test_sir.Rmd

Ignore += codes/ali_test_sir.html
## codes/ali_test_sir.html: codes/ali_test_sir.Rmd

subdirs += note codes
note/SIR_notes.pdf:

Sources += smoothing.tex
Ignore += smoothing.pdf
smoothing.pdf: smoothing.tex
	$(pandocs)

## Ali et all manuscript
Sources += SIR_manuscript.tex
Ignore += SIR_manuscript.pdf

## tex2dvi -p automatically handles all of the repeated runs needed for updating references, citations, etc.
## sudo apt-get install texinfo ...
SIR_manuscript.pdf: SIR_manuscript.tex SIRlibrary.bib pix/R0contour_random.pdf pix/R0contour_TTI.pdf
	texi2dvi -p SIR_manuscript.tex

pix/R0contour_random.pdf: codes/sir_plot.R codes/params.R codes/SIRfunctions.R
	cd codes; R CMD BATCH --vanilla sir_plot.R

pix/R0contour_TTI.pdf: codes/sir_plot.R codes/params.R codes/SIRfunctions.R
	cd codes; R CMD BATCH --vanilla sir_plot.R

alldirs += $(subdirs)

######################################################################

## Cribbing
%.tex: ali/%.tex
	/bin/cp $< $@

codes/%.Rmd: ali/codes/%.Rmd
	/bin/cp $< $@

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

## Want to chain and make makestuff if it doesn't exist
## Compress this ¶ to choose default makestuff route
Makefile: makestuff/Makefile
makestuff/Makefile:
	git clone $(msrepo)/makestuff
	ls makestuff/Makefile

-include makestuff/os.mk

-include makestuff/texi.mk
-include makestuff/makeR.mk
-include makestuff/pandoc.mk

-include makestuff/git.mk
-include makestuff/visual.mk

