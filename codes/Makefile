## This is the codes subdirectory of Ali's SIRmodel

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

autopipeR = defined
Sources += $(wildcard *.R)

params.Rout: params.R 

SIRfunctions.Rout: SIRfunctions.R

sir_plot.Rout: sir_plot.R SIRfunctions.rda params.rda

Ignore += R0contour_TTI.pdf R0contour_random.pdf
R0contour_TTI.pdf R0contour_random.pdf: sir_plot.Rout ;

Ignore += targeted_test_df.csv random_test_df.csv
targeted_test_df.csv random_test_df.csv: sir_plot.Rout ;

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
	cd .. && $(MAKE) makestuff
	ln -s ../makestuff .

-include makestuff/os.mk

-include makestuff/pipeR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
