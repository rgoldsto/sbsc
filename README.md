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

The "Simulations" folder contains plots of 9 classes of simulations:

  evidenceI X utilityD: a 4 (evidenceI = 0, 0.25, 0.75, 1) X 5 (utilityD = 0.02, 0.04, 0.1, 0.2, 0.4) parameter sweep

  gamma X evidenceIntegration

  numCommunities X numOpinions

  numCommunities X evidenceIntegration

  numOpinions X gamma

  numRounds X evidenceI

  probO X gamma

  probOut X evidenceIntegration

  weightO X evidenceI

For each of these simulation classes, we calculate 5 measures for the networks:
  Average Performance: The average quality of opinions of the agents in a population

  Percentage Best: The percentage of agents who have the best (highest ranked) opinion

  Entropy: The entropy of the opinions in the population.  As agents adopt a more diverse set of opinions (e.g. it becomes harder to predict what opinion any particular agent has based on the overall distribution of opinions in the population), entropy increases.

  Mutual information: Measures the dependency between two categorical variables via I(X;Y) = H(Y)-H(Y|X) = H(X)-H(X|Y) = H(X,Y)-H(X|Y)-H(Y|X).  In the present case, the variables are the communities to which agents belong and the opinions possessed by agents.  Note that if all agents in a population have the same opinion, the entropy will low and the mutual information will also be low because knowing an agent's community doesn't provide any (additional) information in determining their opinion and vice versa.

  Assortativity: the average probability of two network neighbors having the same opinion. Assortativity is not always the same as mutual information because if every agent in the population has the same opinion, assortativity is high, but mutual information is low.

