The file "opinionsInNetworksMisinformationGit.r" contains the code necessary to recreate the figures in the paper:

Goldstone, R. L., Dubova, M., Aiyappa, R., & Edinger, A. (under review). The Spread of Beliefs in Partially Modularized Communities.

The software has been tested on versions 4.3.0 of R, 1.4.3 of Igraph, 1.2.0.1 of Infotheory, 1.1.2 of Dplyr, and 3.4.2 of ggplot2.

Lines 46-58 set the default values of the model parameters.  These will be used unless they are overruled by other explicitly passed values, as happens with the parameter sweeps near the end of the code.

It is easy to create your own two-parameter sweeps to observe how the behavior of the simulated population changes depending on important factors in the model.  An example of a two-parameter sweep is:

ExpTest<-parameterSweep("numR", c(5,10,20,40,80),"evidenceI",c(0,0.25,0.75,1),reps=50)
plotResults(ExpTest)

This creates and then plots a 5 X 4 set of simulations, each of which is the result of 50 replications.  Five levels of number of rounds of opinion exchange (5, 10, 20, 40, or 80 rounds) is factorially combined with five levels of evidence integration (0 = purely summing agents, 0.25, 0.75, and 1 = purely averaging agents).  The parameters that can be simulated include:
"probO"=Probability of Connecting to Outgroup Member
"weightO"=Weight to Outgroup Members
"numO"=Number of Opinions
"numC"=Number of Communities
"evidenceI"=Evidence Integration
"utilityD"=Utility Dropoff Rate
"weightI"=Weight to Ingroup Members
"weightM"=Weight of Agent to Itself
"gam"=gamma
"probI"=Probability of Connecting to Ingroup Member
"numR"=Number of Rounds

See the above paper for more information on these parameters.
