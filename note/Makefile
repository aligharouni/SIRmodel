## This is the note subdirectory of Ali's SIRmodel
## I made a Makefile here because the tex files are designed to be run from here.

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

SIR_notes.pdf: SIR_notes.tex

minutes.html: minutes.Rmd
######################################################################

### Makestuff

Sources += Makefile

Makefile: makestuff/Makefile
makestuff/Makefile:
	cd .. && $(MAKE) $@
	ln -s ../makestuff .
	ls $@

-include makestuff/os.mk

-include makestuff/texi.mk
-include makestuff/makeR.mk
-include makestuff/pandoc.mk

-include makestuff/git.mk
-include makestuff/visual.mk
