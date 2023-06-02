#R code for the Spread of Beliefs in Segregated Communities Model
# Programmed by Robert Goldstone (https://pc.cogs.indiana.edu), 2022


library(igraph)
library(infotheo)
library(dplyr)
library(ggplot2)

setwd( "/Users/rgoldsto/word_proc_docs/socialCAS/opinionsInNetworks") # This line should be set to your own preferred working directory

#Next steps:
# Could make people more stubborn to change opinion if they have held an opinion for a long time - previous steps opinions influence "myself"
# Could allow for edges to change based on homophily

numAgents<-200
numOpinions<-10 
numCommunities<-4
evidenceIntegration<-0 # 0=summed evidence (popularity biased), 1=averaged
utilityDropoff<-0.1 # how quickly the utility for different opinions drops off as they get worse. If <<1 then all opinions are equally good
weightIngroup<-1 # weight given to evidence from in-group members
weightOutgroup <-0.1 #weight given to out-group members
weightMyself<-1 # weight given to my own evidence
gamma<-1 #determinism in softmax decision rule. Lower number means more stochasticity
probIn<-0.12 #probability of creating edge within group
probOut<-0.004 #probability of creating edge between groups
numRounds<-20 # rounds of opinion exchange 

edgeProbabilities<-matrix(probOut, nrow=numCommunities, ncol=numCommunities)
diag(edgeProbabilities)<-probIn #set diagonal, in group probabilities

# generate a network according to a stochastic block model
g <- sample_sbm(numAgents, pref.matrix=edgeProbabilities, block.sizes=rep(numAgents/numCommunities,numCommunities))
com <- spinglass.community(g, spins=numCommunities) #one way of determining communities
minC <- rep(-Inf, vcount(g))
maxC <- rep(Inf, vcount(g))
minC[1] <- maxC[1] <- 0
co <- layout_with_fr(g, minx=minC, maxx=maxC,
                     miny=minC, maxy=maxC)
V(g)$community<-com$membership
V(g)$opinion<-sample(1:numOpinions, numAgents, replace=T) #give each agent an opinion
V(g)$color <- V(g)$opinion
plot(g,layout=co,vertex.size=4)



exchangeOpinions<-function (g,numA=numAgents,numO=numOpinions,numC=numCommunities,evidenceI=evidenceIntegration,
   utilityD=utilityDropoff,weightI=weightIngroup,weightO=weightOutgroup,weightM=weightMyself,
   gam=gamma)
{
  newOpinions<-c()
  for (i in V(g)) #for each agent, change opinion
  {
    evidence<-rep(0,numO)
    numPicking<-rep(0,numO) #number of neighbors picking each option
    neighbors<-neighbors(g,i)
    for (j in neighbors)
    {
      neighborOpinion<-vertex_attr(g, "opinion",j)
      utility<-exp(-neighborOpinion*utilityD)
      #numPicking[neighborOpinion]<-numPicking[neighborOpinion]+1 #if really want to count all opinions equally even if weighted less for outgroup
      if (vertex_attr(g, "community",i) == vertex_attr(g, "community",j)) #same group?
      {
        evidence[neighborOpinion]<-evidence[neighborOpinion]+weightI*utility
        numPicking[neighborOpinion]<- numPicking[neighborOpinion]+weightIngroup #Note, adding 
      }
      else
      {
        evidence[neighborOpinion]<-evidence[neighborOpinion]+weightO*utility
        numPicking[neighborOpinion]<- numPicking[neighborOpinion]+weightOutgroup        
        #If weightOut is low  then somebody in another group having opinion X could lead one NOT
        #to adopt X to the extent that evidence is divided by the numPicking - the denominator is increased a bit
        #by the additional opinion X, but the numerator isn't increased much because it is an outgroup opinion
      }
    }
    
    myOpinion<-vertex_attr(g, "opinion",i) #add in my own opinion
    evidence[myOpinion]<-evidence[myOpinion]+weightM*exp(-myOpinion*utilityD)
    numPicking[myOpinion]<-numPicking[myOpinion]+1
    
    for (j in 1:numO)
    {
      normalizer<-(evidenceI*(numPicking[j]-1)+1)
      if (normalizer==0)
      {
        evidence[j]=0 # could be >0 if it is possible to create a new opinion different from any neighbor
      }
      else
      {
        evidence[j]<-evidence[j]/normalizer
        if (evidence[j]!=0){
          evidence[j]<-exp(evidence[j]*gam) }
      }
    }
    newOpinions[i]<-sample(1:numO,1,prob=evidence,replace=TRUE) #Luce Choice Rule
    #cat("vertex ",i,": ",evidence,", ",vertex_attr(g, "opinion",i),"->",newOpinions[i],"\n")
  }
  V(g)$opinion<-newOpinions
  return(g)
}

oneRun<-function(numA=numAgents,numO=numOpinions,numC=numCommunities,evidenceI=evidenceIntegration,
                 utilityD=utilityDropoff,weightI=weightIngroup,weightO=weightOutgroup,weightM=weightMyself,
                 gam=gamma,probI=probIn,probO=probOut,numR=numRounds,var1="",var2="",val1=0,val2=0) #assign default values to whatever isn't explicitly passed as swept parameter
  {
    edgeProbabilities<-matrix(probO, nrow=numC, ncol=numC) #create entirely new graph each run
    diag(edgeProbabilities)<-probI #set diagonal, in group probabilities
    g <- sample_sbm(numA, pref.matrix=edgeProbabilities, block.sizes=rep(numA/numC,numC))
    while (!(is.connected(g)))
    {
      g <- sample_sbm(numA, pref.matrix=edgeProbabilities, block.sizes=rep(numA/numC,numC)) #keep on making graphs until connected
    }
    com <- spinglass.community(g, spins=numC) #one way of determining communities
    V(g)$community<-com$membership
    V(g)$opinion<-sample(1:numO, numA, replace=T) #give each agent a random opinion
    minC <- rep(-Inf, vcount(g))
    maxC <- rep(Inf, vcount(g))
    minC[1] <- maxC[1] <- 0
    co <- layout_with_fr(g, minx=minC, maxx=maxC,
                         miny=minC, maxy=maxC)
    for (k in 1:numR)
    {
      g<-exchangeOpinions(g,numA,numO,numC,evidenceI,utilityD,weightI,weightO,weightM,gam)
    }
    #V(g)$color <- scales::dscale(as.factor(11-V(g)$opinion),diverging_pal) #Good to see what opinions look like on continuum
    V(g)$color<-V(g)$opinion
    plot(g,layout=co,vertex.size=4,vertex.label=NA,main=paste(var1,"=",val1,"  ",var2,"=",val2)) # plot after opinion exchanges
    return(g) #for now, just returning opinions, not graph
  }
  
parameterSweep <- function (p1,p1vals,p2,p2vals,reps)
  {#need to change probOut to generalize in two places: 150 and 170
    results<-data.frame()
    for (par1 in p1vals) # ,0.006,.018,0.054
    {
      for (par2 in p2vals) # ,0.3,0.5,0.7,0.9
      {
        cat(p1,"=",par1,"  ",p2,"=",par2,"\n")
        for (r in 1:reps)
        { 
          #assortativity = average probability of a network neighbor having the same opinion
          #Assortativity is not always the same as mutual information because if everybody has same opinion, assortativity is high, but mutual information is low
          args<-list()
          args[[p1]]<-par1
          args[[p2]]<-par2
          args[["var1"]]<-p1
          args[["var2"]]<-p2
          args[["val1"]]<-par1
          args[["val2"]]<-par2
          network<-do.call(oneRun,args) #must convert string to actual variable name for passing to function
          A <- as.matrix(get.adjacency(network))
          assort <- mean(rowSums(outer(V(network)$opinion, V(network)$opinion, `==`) * A) / rowSums(A))
          newRow=data.frame(var1=par1,var2=par2,
                            pctBest=(sum(V(network)$opinion==1)/length(V(network))),aveValue=6-mean(V(network)$opinion), #note: reverse Average value so high value is good
                            entropy=entropy(V(network)$opinion),mutualInformation=mutinformation(V(network)$community,V(network)$opinion),
                            assortativity=assort)
          names(newRow)[names(newRow) == "var1"] <- p1
          names(newRow)[names(newRow) == "var2"] <- p2
          results<-rbind(results,newRow)
        }
      }
    }
    return(results)
  }

plotResults<-function(exp)
  { #assumes first variable is X axix, second is for different lines
    conversions=list("probO"="Probability of Connecting to Outgroup Member",
                     "weightO"="Weight to Outgroup Members","numO"="Number of Opinions","numC"="Number of\n Communities","evidenceI"="Evidence Integration",
                     "utilityD"="Utility Dropoff","weightI"="Weight to Ingroup Members","weightM"="Weight of Agent to Itself","gam"="gamma","probI"="Probability of Connecting to Ingroup Member",
                     "numR"="Number of Rounds")
    myBreaks<-unlist(unique(exp[2]))
    xBreaks<-unlist(unique(exp[1]))
    xRange<-max(xBreaks)-min(xBreaks)
    xAxisName=conversions[[names(exp)[1]]]
    legendTitle=conversions[[names(exp)[2]]]
    exp$rank<-match(exp[,2],myBreaks)
    p<-ggplot(exp,aes_string(x=names(exp)[1],y="aveValue",group="rank",color="rank")) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=xRange/16)+ 
      stat_summary(fun=mean,geom="point")+
      scale_x_continuous(labels=xBreaks,breaks=xBreaks)+
      scale_color_gradient(name=legendTitle,labels=myBreaks,guide = guide_legend()) +xlab(xAxisName) + ylab("Average Performance") +
      theme_bw() + theme(legend.title.align=0.5,panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    q<-ggplot(exp,aes_string(x=names(exp)[1],y="pctBest",group="rank",color="rank")) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=xRange/16)+ 
      stat_summary(fun=mean,geom="point")+
      scale_color_gradient(name=legendTitle,labels=myBreaks,guide = guide_legend()) +xlab(xAxisName) + ylab("Percentage of Agents Who Have the Best Opinion")+
      theme_bw() + theme(legend.title.align=0.5,panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    r<-ggplot(exp,aes_string(x=names(exp)[1],y="entropy",group="rank",color="rank")) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=xRange/16)+
      stat_summary(fun=mean,geom="point")+
      scale_color_gradient(name=legendTitle,labels=myBreaks,guide = guide_legend())+xlab(xAxisName) + ylab("Entropy")+
      theme_bw() + theme(legend.title.align=0.5,panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    s<-ggplot(exp,aes_string(x=names(exp)[1],y="mutualInformation",group="rank",color="rank")) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=xRange/16)+
      stat_summary(fun=mean,geom="point")+
      scale_color_gradient(name=legendTitle,labels=myBreaks,guide = guide_legend())+xlab(xAxisName) + ylab("Mutual Information")+
      theme_bw() + theme(legend.title.align=0.5,panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    t<-ggplot(exp,aes_string(x=names(exp)[1],y="assortativity",group="rank",color="rank")) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar",width=xRange/16)+
      stat_summary(fun=mean,geom="point")+
      scale_color_gradient(name=legendTitle,labels=myBreaks,guide = guide_legend())+xlab(xAxisName) + ylab("Assortativity")+
      theme_bw() + theme(legend.title.align=0.5,panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
     #old style for mapping continuous gradient to color:
    #p<-ggplot(exp,aes_string(x=names(exp)[1],y="aveValue",group=names(exp)[2],color=names(exp)[2])) + stat_summary(fun=mean,geom="line")+
     #    stat_summary(fun.data=mean_cl_normal,geom="errorbar")+ 
    #     xlab(xAxisName) + ylab("Average Performance")+ labs(color=legendTitle)
    #q<-ggplot(exp,aes_string(x=names(exp)[1],y="pctBest",group=names(exp)[2],color=names(exp)[2])) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar")+ xlab(xAxisName) + ylab("Percentage of Agents Who Have the Best Opinion")+ labs(color=legendTitle) 
    #r<-ggplot(exp,aes_string(x=names(exp)[1],y="entropy",group=names(exp)[2],color=names(exp)[2])) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar")+ xlab(xAxisName) + ylab("Entropy")+ labs(color=legendTitle) 
    #s<-ggplot(exp,aes_string(x=names(exp)[1],y="mutualInformation",group=names(exp)[2],color=names(exp)[2])) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar")+ xlab(xAxisName) + ylab("Mutual Information")+ labs(color=legendTitle) 
    #t<-ggplot(exp,aes_string(x=names(exp)[1],y="assortativity",group=names(exp)[2],color=names(exp)[2])) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar")+ xlab(xAxisName) + ylab("Assortativity")+ labs(color=legendTitle) 
    ggsave(paste(names(exp[1]),"X",names(exp[2]),"aveValue.jpg"),p,device="jpeg",width=1000,scale=3,height=700,units = "px")
    ggsave(paste(names(exp[1]),"X",names(exp[2]),"pctBest.jpg"),q,device="jpeg",width=1000,scale=3,height=700,units = "px")
    ggsave(paste(names(exp[1]),"X",names(exp[2]),"entropy.jpg"),r,device="jpeg",width=1000,scale=3,height=700,units = "px")
    ggsave(paste(names(exp[1]),"X",names(exp[2]),"mutualInformation.jpg"),s,device="jpeg",width=1000,scale=3,height=700,units = "px")
    ggsave(paste(names(exp[1]),"X",names(exp[2]),"assortativity.jpg"),t,device="jpeg",width=1000,scale=3,height=700,units = "px")
}


Exp1<-parameterSweep("probO", c(0.002,0.004,.008,0.016),"weightO",c(0.1,0.3,0.5,0.7,0.9),reps=50)
#Note: 3- average value so up is good
ggplot(Exp1,aes(x=probOut,y=3-aveValue,group=weightOut,color=weightOut)) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar")+ xlab("Probability of Connecting to Outgroup Member") + ylab("Average Performance")
#Both connecting to outgroup members is good and weighting their opinion, and combination is particularly good.  If don't weight others, doesn't matter if you connect to them
ggplot(Exp1,aes(x=probOut,y=pctBest,group=weightOut,color=weightOut)) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar")+ xlab("Probability of Connecting to Outgroup Member") + ylab("Percentage of Agents Who Have the Best Opinion")
ggplot(Exp1,aes(x=probOut,y=entropy,group=weightOut,color=weightOut)) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar")+ xlab("Probability of Connecting to Outgroup Member") + ylab("Entropy")
#High Entropy is mostly a bad thing because when people are in agreement, they are usually in agreement on the right opinion
ggplot(Exp1,aes(x=probOut,y=mutualInformation,group=weightOut,color=weightOut)) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar")+ xlab("Probability of Connecting to Outgroup Member") + ylab("Mutual Information")
#High mutual information is bad because it means that there are pockets of people who have fring opinions which are likely to be wrong
ggplot(Exp1,aes(x=probOut,y=assortativity,group=weightOut,color=weightOut)) + stat_summary(fun=mean,geom="line")+stat_summary(fun.data=mean_cl_normal,geom="errorbar")+ xlab("Probability of Connecting to Outgroup Member") + ylab("Assortativity")
#If don't connect to outgroup members, then weight you give to them doesn't matter.  If you are connected to many outgroup members but you don't weight them, then assortativity will be low- there'll be lots of times when you have opposite opinions of others
#But if you're connected to a lot of outgroup members and you care what they think, then assortativity will be high - all of your neighbors will have your opinion

Exp0<-parameterSweep("numR", c(5,10,20,40,80),"evidenceI",c(0,0.25,0.75,1),reps=50)
plotResults(ExpTest)

exp1<-parameterSweep("probO", c(0.002,0.004,.008,0.016),"weightO",c(0.1,0.3,0.5,0.7,0.9),reps=100)
plotResults(exp1)

numAgents=240 # and probOut=0.004
#exp2a<-parameterSweep("numO", c(12),"numC",c(12),reps=3) # numAgents is divisible evenly by all values
exp2<-parameterSweep("numO", c(2,5,8,12),"numC",c(2,5,8,12),reps=100) # numAgents is divisible evenly by all values
#Notice below - aes_string allows strings to be put in as variable names so that "variable" variables can be used
plotResults(exp2)
# With more opinions, hard to find better ones. More communities also uniformly hurts because found good opinions stay locked in their community. Superadditive: Lots of opinions with lots of groups particularly bad - no rallying around one good opinion
#mutualInformation is nonmonotonic with num communities. It is low with many communities (more than one community with same opinion) or few (one opinion takes over everywhere)
# assortativity is worse with more opinions and communities because people will disagree
#Take home: particularly hard to find good/best opinions as there are more messages and more communities. REALLY need to reduce number of communities when there are a lot of opinions

numAgents=240
exp9<-parameterSweep("numC",c(2,5,8,12),"evidenceI",c(0,0.25,0.75,1),reps=100)
#exp9a<-parameterSweep("numC",c(12),"evidenceI",c(1),reps=4)
plotResults(exp9)
#strong interaction. For small number of communities, better to sum opinions.  For larger number of communities, better to average opinions
# The interaction is cross over in both directions, so having greater number of communities HELPS performance of groups
# if the individuals will be averaging, and fewer communities helps for summation.
#Why? small number of communities means each community is bigger, so summation can have amplifying influence
#If large number of groups, then averaging is better because sums won't sample outside 

#probability out group connection X evidence Integration
numAgents=200
exp8<-parameterSweep("probO", c(0.002,0.004,.008,0.016),"evidenceI",c(0,0.25,0.75,1),reps=100)
#exp8a<-parameterSweep("probO", c(0.002),"evidenceI",c(1),reps=4)
plotResults(exp8)
#Averaging is bad as outgroup probability increases.  Too much heterogeneity in opinions when averaging - no single opinion catches on very widely
#when summing, one or two opinions catch on (and usually best ones)
#The problem with combining high outgroup connection with averaging is that BOTH promote diversity -- too much diversity of opinion
#See also Campbell, Izquierdo, and Goldstone on maintaining an optimal level of agent diversity - problems with too much or too little diversity

# numOpinions X gamma
exp7<-parameterSweep("numO", c(2,5,8,12),"gam",c(0.4,0.8,1,10),reps=100)
#exp7a<-parameterSweep("numO", c(12),"gam",c(0.4),reps=4)
plotResults(exp7)
#default evidence integration=0=summed
#Do worse with more options, obviously.  Gamma=1 is best, even better than more deterministic gamma=10.  That's because it explores the opinion space more
#If few opinions, gamma=0.8>gamma=10, but if many opinions gamma=10>gamma=0.8.  I might have thought you need more random exploration
#with more opinions, but the problem is apparently the network doesn't settle down.  If lots of possibilities, need MORE deterministic search

exp5<-parameterSweep("evidenceI",c(0,0.25,0.75,1),"gam",c(0.4,0.8,1,10),reps=100)
#exp5a<-parameterSweep("evidenceI",c(0),"gam",c(10),reps=4)
plotResults(exp5)
#EvidenceIntegration (four values) X gamma - if harder to choose better values, is it better if popularity based?
#0=summed evidence (popularity); 1=averaged. generally, summing evidence is good because each person has some sensitivity for telling truth
#If people have very little sensitivity (gamma is low) then averaging is better because summing can inflate noise/wrong answers
#implication - going with group is dangerous when people can't reliably choose better answer
#when gamma is high, then using average opinion makes people not choose best answer often - they'll choose something in local network which may not be "tried and true"
#explore-exploit: gamma can be too high (if summing evidence) - opinion catches on too quickly and can't be shaken lose by better opinion
#Intermediate values of integration between summing and averaging perform very differently depending on gamma
#If deterministic, then intermediate is best strategy.  if stochastic, then worse!
#For averaging strategy, generally each sub-group has some orange, but it never catches on
#For summing, some sub-groups are all orange and others are none. Hence high assortativity
#If gamma is high, then halfway between summing and averaging allows best answer to often catch on everywhere
#If gamma is low, then halfway is a hodgepodge of opinions. Noise is too high to let best catch on
#If gamma is low, then summing typically allows best to catch on in a couple of groups
#So, overall, there is a tradeoff between information catching across a network and premature convergence on not-great options
#most impressive case of this is that for summing, it is better to have gamma=1 not gamma=10.  Less sensitivity avoids things catching on too fast, particularly if based on sum of neighbors' iffy opinions


#evidenceIntegration<-0 # 0=summed evidence (popularity biased), 1=averaged
#utilityDropoff<-0.1 # how quickly the utility for different opinions drops off as they get worse. If <<1 then all opinions are equally good
numAgents=200 # and probOut=0.004
exp3<-parameterSweep("evidenceI",c(0,0.25,0.75,1),"utilityD",c(0.02,0.04,0.1,0.2,0.4),reps=100)
#exp3a<-parameterSweep("evidenceI",c(0,0.25,0.75,1),"utilityD",c(0.02,0.04,0.1,0.2,0.4),reps=2)
plotResults(exp3)
# a high utility dropoff makes a BIGGER difference between first and second best opinion, but decreases difference between 8th and 9th
#Summing generally better than averaging - summing amplifies some real signal
#

exp4<-parameterSweep("probO", c(0.002,0.004,.008,0.016),"gam",c(0.4,0.8,1,10),reps=100)
plotResults(exp4)
#Best gamma is 1, not 5.  Too little exploration if gamma is 5.  Too much randomness is very bad though, so lower gammas are really bad
#If more deterministic, then increasing Outgroup is relatively GOOD because it adds diversity (and their opinions can be trusted?)
#If more stochastic, then increasing Outgroup is BAD because you may be influenced by a wrong opinion



exp6<-parameterSweep("weightO",c(0.1,0.3,0.5,0.7,0.9),"evidenceI",c(0,0.25,0.75,1),reps=100)
plotResults(exp6)
# weightOutgroup X evidenceIntegration (two values) - if biased against others' opinions is it bad to be popularity based?
# people do modestly better when they give greater weight to outgroup members - use more information
# Intermediate levels of evidence Integration are worse than extremes!  Not clear why.
# Best if summing evidence, not averaging.

numAgents=200

#to see ordering of the colors:
x = 1:10
y = rep(1,10)
plot(x,y, pch=20, cex=10, col=rev(diverging_pal(10)), xlim=c(0.5,10.5)) #another good palette is 
plot(x,y, pch=20, cex=10, col=categorical_pal(8), xlim=c(0.5,10.5)) #this is the default palette
text(x,y-0.04, labels=c(1:10))
#ordering: orange, light blue, green yellow...

