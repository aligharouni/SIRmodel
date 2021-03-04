## This is the notes subdirectory of Ali's SIRmodel

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

Sources += $(wildcard *.R)

## Cribbing rule

%.R: 
	$(CP) ../ali/note/$*.R .

######################################################################

## Try to rebuild scripts

autopipeR = defined

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
