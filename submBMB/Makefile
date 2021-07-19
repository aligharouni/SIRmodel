## This is Ali's SIRmodel submission to BMB

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

subdirs += note codes

hotdirs += codes

alldirs += $(subdirs)

Ignore += $(alldirs)

######################################################################

## Manuscript

Sources += $(wildcard *.tex *.bib *.R)

##SIR_manuscript.pdf: ../SIR_manuscript.tex ../codes/modeldefs.tex

cover.pdf: cover.tex


## Figures
## There should be a way to create all eps at once (I feel!)?
Ignore += sir_comp.eps
sir_comp.eps: ../pix/sir_comp.pdf 
	pdftops -eps $< $@
	
Ignore += R0contour_random.eps
R0contour_random.eps: ../codes/R0contour_random.pdf 
	pdftops -eps $< $@

Ignore += R0contour_TTI.eps
R0contour_TTI.eps: ../codes/R0contour_TTI.pdf 
	pdftops -eps $< $@
	
Ignore += R0contour_random2.eps
R0contour_random2.eps: ../codes/R0contour_random2.pdf 
	pdftops -eps $< $@
	
Ignore += R0contour_TTI2.eps
R0contour_TTI2.eps: ../codes/R0contour_TTI2.pdf 
	pdftops -eps $< $@

	
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
-include makestuff/pipeR.mk
-include makestuff/hotcold.mk

-include makestuff/git.mk
-include makestuff/visual.mk

