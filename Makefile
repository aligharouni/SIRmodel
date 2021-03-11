## This is Ali's SIRmodel (makestuff branch)

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

## Ali et all manuscript

Sources += $(wildcard *.tex *.bib)

SIR_manuscript.pdf: SIR_manuscript.tex

## compartmental flowchart in ipe
## sudo apt-get install ipe
Ignore += pix/sir_comp.pdf
pix/sir_comp.pdf: pix/sir_comp.ipe
	ipetoipe -pdf $< $@

######################################################################

## Old and busted notes
## smoothing.pdf: smoothing.tex

######################################################################

## Copying and cribbing

Ignore += ali
ali:
	git clone https://github.com/aligharouni/SIRmodel.git $@

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

