
line 23.  We use "reproduction number" rather than "reproductive number" in most of the paper.  It is slightly odd to use one term in the intro and the other in the rest of the paper.

# Methods

lines 83-84.  It is better to use Roman font (not math italic) for subscripts that are not variables themselves.  I would say it is late to make this change, but I wonder if it is the copy editor's fault, since we do use Roman font subscripts in the flow chart in Fig. 1.  It would be good to be consistent.

- AG: yes, it is the copy editor's fault, we got it right in the submitted version. Somehow they didn't recognize \_ and they changed it to _. 

caption to Fig. 1: "sensitivity (probability of positive tests)".  This should be "sensitivity (probability that an infected individual tests positive)"

Table 1 header:  "Value" would be better as "Default Value" or "Assumed value".  I assume the meaning is the value used in simulations.

line 121.  The word waiting should probably be in quotes ("waiting") to be consistent with elsewhere and with "confirmed" in the same sentence.

lines 144-145.  Something is missing in this phrase and the result doesn't make sense:

   "how much these processes in reducing spread from low levels, 
   and is in turn given by"

I suggest:

   "how much these processes in reduce spread, 
   and is in turn given by"

line 179.  "maximum capacity of ρ"
would be better as "maximum testing intensity ρ"

Fig 3 caption: when referring to panels, we should have (a) and (b) not a and b.

# Results

# Discussion

line 278.  "other" --> "others"

line 287.  S_n is waiting and S_u is untested.  The subscripts need to be exchanged.

line 306.  "than typical" --> "than is typical"

line 319.  "Future steps" --> "Future work" or "Future research"

line 364.   Do we want everything public that is on this repo?  Is it worth careting an SIR_testing_model repo and putting only the source code there and not our revisions, responses to referees etc.  I would think that a simpler repo with only what is promised in the paper and no unnecessary clutter would be better.

# Appendix

line 406.  "summ" --> "sum"

line 408.  end line with : rather than .

line 412.  "only determined" --> "determined only"

line 413-414.  The typesetting of this equation is very unfortunate. Perhaps they can use a smaller font for this equation so it fits on one line.  I note that they have done that in eq A20 further down.

AG: I put a comment on this line.

line 415.  "G11 does so" --> "so does G11"

line 416.  There should be a comma before "which", i.e., "eigenvalue, which is"

line 435-436.  If Delta change either increase or decrease with rho, how does it follow that R0 necessarily decreases with rho??

AG: added "R0<0 or R0>0".

line 439.  delete "the one"

line 487 and the rest of this long equation: w_I is typeset with underscores rather than subscripts.  This needs to be fixed everywhere.  Actually, this problem persists for the rest of the appendix.  You need to go through carefullly and fix it many many times.  I'm puzzled by how this could happen.  It is as if they changed $w_I$ to $w\_I$ etc.

AG: It is indeed w\_{IS}, and since they don't recognize \_, it reads as w_S*I which is incorrect. I have to change it to w_{IS} thoroughly.  

line 570.  w_IS has underscore rather than subscript.  And what is w_{IS}??  shouldn't this say  w_I or w_S ?   or w_X  for X\in\{S,I}

AG: w_{IS}=w_I/w_S. It is defined in line 483. Somehow, they couldn't recognize w\_{IS}, so I have to change it to w_{IS}.

line 590.  It is too late now, but it would have been better to define sigma in eq (1) in this way and have a model that is biologically well-posed.  If it really doesn't affect anything, then why didn't we do it?  Maybe it isn't too late and we can simply replace eq (1) appropriately.  I had forgotten this problem.  Was there discussion about this point that led to having this "oops" appendix rather than fixing the problem in the main body of the ms?

