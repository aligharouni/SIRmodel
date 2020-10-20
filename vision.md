# Vision statement

An exploration of simple mathematical/computational models of COVID-19 to gain insight into the effects of different levels of testing and isolation on $R_t$/growth rates. The goal is insight.

What's out there already, that we want to differentiate from ("literature review of models of test/trace/isolate" should be on the to-do list)

* explicit models of TTI (trace/test/isolate) based on network or agent-based models
* models of repeated random testing of isolated populations (e.g. LTC/university/sports bubble) (e.g. Bergstrom, Bergstrom, and Hao(?))

What should or could be in this paper:

## Part 1

* analysis of the basic SIR model with testing: analytical results as far as possible, supplemented by numerical results as necessary/for illustration. Heuristic explanations/insights of what we learn from this.
* Maybe?? order-of-magnitude/back-of-the-envelope calculations of the magnitudes of testing & isolation necessary. (However, the model may be too simplified for even back-of-the-envelope calculations to be meaningful.)
* Maybe?? comparison with strength-and-speed framework (since T/T/I is essentially a speed-based intervention)

## Part 2

* extension to a SEPAIR model (i.e. including latent/exposed, presymptomatic, asymptomatic)
    * analytical results will probably be out of reach
	* but it _may_ be possible to extend some of the insights gained in part 1
	* in any case, some numerical results (the space of relevant/interesting parameters will be larger: figure out a sensible way to show results from this parameter space. Are there useful dimensionless parameters that could help frame this?)
* more detailed quantitative exploration: in particular, *either* Latin hypercube *or* (nearly the same!) a sample over estimated ('prior') distributions of the parameters. Goal: answer the general question "how likely is TTI to be able to control COVID?"

## Unresolved/open questions

* We will have to think a lot more about the relationship/mapping between testing policy (including contact tracing intensity) and the weighting (W_j) parameters.
* At present the model is ill-posed (i.e. it crashes due to the test-rate-scaling component). Is there a simple way to change this so that it's biologically well-posed (i.e. all trajectories that start non-negative remain non-negative), and has the properties we want?
