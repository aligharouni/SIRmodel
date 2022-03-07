## This is Ali's SIRmodel

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt manuscript_BMB.tex response1_BMB.tex"

######################################################################

subdirs += note codes

hotdirs += codes

alldirs += $(subdirs)

Ignore += $(alldirs)

######################################################################

## Manuscript

Sources += $(wildcard *.bib *.R) README.md
Sources += abstract.tex cover.tex manuscript_BMB.tex response1_BMB.tex SIR_manuscript.tex

authors.Rout: authors.R
	$(pipeR)

SIR_manuscript.pdf: SIR_manuscript.tex codes/modeldefs.tex

cover.pdf: cover.tex

## compartmental flowchart in ipe
## sudo apt-get install ipe
Ignore += pix/sir_comp.pdf
pix/sir_comp.pdf: pix/sir_comp.ipe
	ipetoipe -pdf $< $@

## Missing files .tex and/or .txt, probably doesn't matter 2021 Sep 13 (Mon)
cohort.pdf: cohort.tex 

######################################################################
## submission to BMB

allrev: submitted manuscript_BMB.pdf manuscript_BMB_olddiff.pdf manuscript_BMB_olddiff.pdf response1_BMB.pdf

manuscript_BMB.pdf: manuscript_BMB.tex codes/modeldefs.tex

submitted: manuscript_BMB.tex.b16a28171.oldfile ;

## manuscript_BMB_olddiff.pdf: manuscript_BMB.tex
manuscript_BMB_olddiff.pdf: manuscript_BMB.tex.b16a28171.oldfile

##First round of responses to reviewers' comments
response1_BMB.pdf: response1_BMB.tex

######################################################################

## Old and busted notes
## smoothing.pdf: smoothing.tex

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

