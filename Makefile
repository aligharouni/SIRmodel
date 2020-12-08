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

note/SIR_notes.pdf:
note/SIR_notes.pdf:

Sources += smoothing.tex
Ignore += smoothing.pdf
smoothing.pdf: smoothing.tex
	$(pandocs)

## Ali et all manuscript
Sources += SIR_manuscript.tex
Ignore += SIR_manuscript.pdf
SIR_manuscript.pdf: SIR_manuscript.tex SIRlibrary.bib pix/R0contour_random.pdf pix/R0contour_TTI.pdf

pix/R0contour_random.pdf: codes/ali_test_sir.Rmd
	cd codes; echo "rmarkdown::render('ali_test_sir.Rmd')" | R --slave

######################################################################

### Makestuff

Sources += Makefile

## Sources += content.mk
## include content.mk

Ignore += makestuff
msrepo = https://github.com/dushoff

## Want to chain and make makestuff if it doesn't exist
## Compress this ¶ to choose default makestuff route
Makefile: makestuff/Makefile
makestuff/Makefile:
	git clone $(msrepo)/makestuff
	ls makestuff/Makefile

-include makestuff/os.mk

-include makestuff/texdeps.mk
-include makestuff/makeR.mk
-include makestuff/pandoc.mk

-include makestuff/git.mk
-include makestuff/visual.mk
