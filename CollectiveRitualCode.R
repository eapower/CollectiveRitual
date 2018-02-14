## Code to accompany ``Social Support Networks and Religiosity in Rural South India'' 
## Includes code for: 
## (1) generating the social support networks
## (2) running the exponential graph models predicting supportive relationships
## (3) calculating cohesiveness of the networks of co-participants & evaluating how probable the observed measures are
## (4) predicting the likelihood that individuals will have supportive relationships with people of other religions
## (5) visualizing the networks


#################################
#################################
## 1. GENERATING THE NETWORKS  ##
#################################
#################################

require(igraph)

## Read in the files that include details of each individual, including age, gender, caste & religion, years of education, household wealth, religious participation, etc.
Ten <- read.csv("TenMetadata.csv",header=TRUE,as.is=TRUE)

## Read in the edge list, which includes columns for Ego, Alter, and then each of the 12 support types
elTen <- read.csv("TenSocialEdgeList.csv")

## generate a network from the edgelist, but this lacks the associated metadata
snaTen <- graph.data.frame(elTen)
IndivIDTen <- V(snaTen)$name
IndivIDTen <- data.frame(IndivIDTen)

## append the individual metadata to each individual included in the network
attTen <- merge(IndivIDTen,Ten,by.x="IndivIDTen",by.y="IndivID",sort=FALSE,all.x=TRUE)
colnames(attTen)[1] <- "IndivID"

## generate the full network
snaTenSupFull <- graph.data.frame(d = elTen, vertices = attTen, directed=TRUE)

## reduce the network down to include only those who completed the survey in the village (many people outside the village are named, but the analyses here require that each individual can have both incoming and outgoing ties)
snaTenSup <- delete.vertices(snaTenSupFull,V(snaTenSupFull)[degree(snaTenSupFull,mode="out")==0])

## limit the edge list to only Hindu respondents
elTenH<-merge(elTen,Ten[,c(1,8)],by.x="Ego",by.y="IndivID",all.x=TRUE)
elTenH<-subset(elTenH,Religion==0)
elTenH<-elTenH[,-16]

## generate networks only for Hindu respondents
snaTenH <- graph.data.frame(elTenH)
IndivIDTenH <- V(snaTenH)$name
IndivIDTenH <- data.frame(IndivIDTenH)
attTenH <- merge(IndivIDTenH,Ten,by.x="IndivIDTenH",by.y="IndivID",sort=FALSE,all.x=TRUE)
colnames(attTenH)[1] <- "IndivID"

snaTenHSupFull <- graph.data.frame(d = elTenH, vertices = attTenH, directed=TRUE)
snaTenHSup <- delete.vertices(snaTenHSupFull,V(snaTenHSupFull)[degree(snaTenHSupFull,mode="out")==0])

######################
## Kinship Networks ##
######################

## This file is an edgelist, including a line for each pair of individuals who are related to one another through a close kinship tie, here including: spousal, parent/child, sibling
closekinEL <- read.csv("CloseKinEdgeList.csv",header=TRUE)

require(car)

## This code generates a full list of all of the dyads in the support network, and then populates it with the data from the kinship edgelist, recording whether each dyad is related or not.
kinship<-function(g){
  kin <- V(g)$name
  kinel <- expand.grid(kin,kin)
  colnames(kinel) <- c("Ego","Alter")
  kinel$Edge <- paste(kinel[,1],kinel[,2],sep="_")
  kin1 <- merge(kinel,closekinEL[,c(3,4)],by="Edge",all.x=TRUE)
  kin1$Kin <- recode(kin1$Kin,"NA=0")
  kin1$Edge<- NULL
  knet <- graph.data.frame(kin1)
  return(knet)
}

knetTenH<-kinship(snaTenHSup)


########################
## Household Distance ##
########################

## These two files record the distance (as the bird flies) between every pair of households in the two villages
TenHHouseholdDist<-read.csv("TenHouseholdDist.csv",as.is=TRUE)

## Extract the household ID for each individual
TenHi<-V(snaTenHSup)$name
TenHh<-V(snaTenHSup)$HouseID
TenHv<-cbind(TenHi,TenHh)
TenHv<-as.data.frame(TenHv)
colnames(TenHv)<-c("IndivID","HouseID")
TenHv<-TenHv[with(TenHv, order(IndivID)), ]

## Make an empty matrix with names as the House ID
TenHdistancemat <-matrix(nrow=length(TenHv$IndivID),ncol=length(TenHv$IndivID),dimnames=list(as.vector(TenHv$HouseID),as.vector(TenHv$HouseID)))

## Populate the matrix with the values from the Distance file. This will take a moment.
for(i in 1:nrow(TenHHouseholdDist)){TenHdistancemat[rownames(TenHdistancemat)==TenHHouseholdDist[i,1],colnames(TenHdistancemat)==TenHHouseholdDist[i,2]] <- as.numeric(TenHHouseholdDist[i,3])}

## Rename points so that they now correspond properly to IndivID not House ID
dimnames(TenHdistancemat)<-list(as.vector(TenHv$IndivID),as.vector(TenHv$IndivID))

## Replace NAs, meaning individuals who are in the same house, with 0s
TenHdistancemat[is.na(TenHdistancemat)]<-0


########################
## PORTING TO STATNET ##
########################

detach(package:igraph)
require(intergraph)
require(network)
require(ergm)

Net_snaTenHSup <- asNetwork(snaTenHSup)
Net_kinTenH <- asNetwork(knetTenH)


##########################################
##########################################
## 2. EXPONTENTIAL RANDOM GRAPH MODELS  ##
##########################################
##########################################

## Note that the GWESP alpha values used here are the ones that result in the lowest AIC & BIC values. To establish these values, I ran models starting with GWESP alpha = 0.1 ((gwesp(alpha=0.1,T)) and increased by steps of 0.1 until AIC and BIC values reached their lowest point.
## Models were also evaluated by looking at the MCMC output (mcmc.diagnostics(modelname)) and goodness of fit (gof(modelname ~idegree + odegree + esp + distance)) to ensure that the generated networks are reasonable approximations of the observed network

## vows$Procession[vows$X2012.Agnichati == 1 | vows$X2012.Pookkuzhi == 1 | vows$X2012.Paalkudam == 1 | vows$X2012.Vel == 1 | vows$X2012.Suriya.Kavadi == 1 | vows$X2012.Paravai.Kavadi == 1 | vows$X2012.Shakthi.Karagam == 1 | vows$X2012.Small.vel == 1 | vows$X2012.Made.Agnichati == 1 | vows$X2012.Malai == 1]<- 1
## vows$Mulaipari[vows$X2012.Mulaipari == 1 | vows$X2012.Tavazumpillai == 1 | vows$X2012.Karumpuli.Sempuli == 1 | vows$X2012.Shakthi.Karagam == 1 | vows$X2012.Grew.Mulaipari == 1]<-1


## For MariAtt: nodematch("MariAtt",diff=TRUE,keep=2) ## WeeklyPlusWorship vs MariAtt (the former includes 14 people who visit the temple on their own each week, but don't regularly attend the monthly puja (mostly SC), and one person who doesn't go each week, but attends the puja)

## For MariFest: nodemix("mulvowNone2012", base = c(2, 4, 5, 6, 8))

## Also tried: edgecov(Net_CoFasting, attrname = "CoFasting") ;  edgecov(Net_CoNormal, attrname = "CoNormal") + edgecov(Net_CoBig, attrname = "CoBig") ;  edgecov(Net_CoFest, attrname = "CoFest") ;  edgecov(Net_CoWeight, attrname = "CoWeight")

## modelTenSupHindu <- readRDS("modelTenSupHindu.rds")
## modelTenSupHindu_coworship<-readRDS("modelTenSupHindu_coworship.rds")
## modelTenSupHindu_cofest<-readRDS("modelTenSupHindu_cofest.rds")
## modelTenSupHindu_all<-readRDS("modelTenSupHindu_all.rds")


modelTenSupHindu <- ergm(Net_snaTenHSup ~ edges + nodecov("Age") + nodematch("Gender") + nodefactor("Gender") + edgecov(Net_kinTenH,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + absdiff("EduYears") + nodeifactor("EverCommMember") + edgecov(TenHdistancemat/10) + mutual + gwesp(0.4,T) + gwdsp(0.4,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
summary(modelTenSupHindu)
saveRDS(modelTenSupHindu,"modelTenSupHindu.rds")

modelTenSupHindu_coworship <- ergm(Net_snaTenHSup ~ edges + nodecov("Age") + nodematch("Gender") + nodefactor("Gender") + edgecov(Net_kinTenH,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + absdiff("EduYears") + nodeifactor("EverCommMember") + edgecov(TenHdistancemat/10) + nodematch("MariAtt",diff=TRUE,keep=2) + mutual + gwesp(0.4,T) + gwdsp(0.4,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
summary(modelTenSupHindu_coworship)
saveRDS(modelTenSupHindu_coworship,"modelTenSupHindu_coworship.rds")

modelTenSupHindu_cofest <- ergm(Net_snaTenHSup ~ edges + nodecov("Age") + nodematch("Gender") + nodefactor("Gender") + edgecov(Net_kinTenH,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + absdiff("EduYears") + nodeifactor("EverCommMember") + edgecov(TenHdistancemat/10) + nodematch("Procession",diff=TRUE,keep=2) + nodematch("Mulaipari",diff=TRUE,keep=2) + mutual + gwesp(0.4,T) + gwdsp(0.4,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
summary(modelTenSupHindu_cofest)
saveRDS(modelTenSupHindu_cofest,"modelTenSupHindu_cofest.rds")

modelTenSupHindu_all <- ergm(Net_snaTenHSup ~ edges + nodecov("Age") + nodematch("Gender") + nodefactor("Gender") + edgecov(Net_kinTenH,attrname="Kin") + nodematch("Caste") + nodefactor("Caste") + absdiff("EduYears") + nodeifactor("EverCommMember") + edgecov(TenHdistancemat/10) + nodematch("MariAtt", diff = TRUE, keep = 2) + nodematch("Procession",diff=TRUE,keep=2) + nodematch("Mulaipari",diff=TRUE,keep=2) + nodeifactor("WeeklyPlusWorship") + nodeifactor("Possession") + nodeicov("WeightedReligPart") + mutual + gwesp(0.4,T) + gwdsp(0.4,T), control = control.ergm(MCMC.burnin=15000,MCMC.samplesize=50000,MCMC.interval=1000),verbose=FALSE)
summary(modelTenSupHindu_all)
saveRDS(modelTenSupHindu_all,"modelTenSupHindu_all.rds")

require(texreg)

texreg(list(modelTenSupHindu,modelTenSupHindu_coworship,modelTenSupHindu_cofest,modelTenSupHindu_all),digits=3)

##################################
## Reporting & Plotting Results ##
##################################

or <- exp( modelTenSupHindu_coworship$coef ) 
ste <- sqrt( diag( modelTenSupHindu_coworship$covar ) ) 
lci <- exp( modelTenSupHindu_coworship$coef-1.96*ste ) 
uci <- exp( modelTenSupHindu_coworship$coef+1.96*ste ) 
oddsratios <- rbind( round( lci,digits = 4 ),round( or,digits = 4 ),round( uci,digits = 4 ) ) 
oddsratios <- t( oddsratios ) 
colnames( oddsratios ) <- c( "Lower","OR","Upper" )
oddsratios
teststat <- modelTenSupHindu_coworship$coef/ste 
teststats <- rbind( round( teststat,digits = 4 )) 
teststats <- t( teststats ) 
colnames( teststats ) <- c("Wald")
teststats

TenResults<-cbind(summary(modelTenSupHindu_coworship)$coefs[1],summary(modelTenSupHindu_coworship)$coefs[2],or,summary(modelTenSupHindu_coworship)$coefs[4])

require(xtable)

xtable(TenResults,digits=c(3,3,3,3,4))


## Plotting ERGM with Co-Wor
TenHCoWorcoefs <- modelTenSupHindu_coworship$coef

estoprob <- function(b) {
  exp(b)/(1+exp(b))
}

TenHCoWorpredergm<-read.csv("TenHPredictERGMCollRit.csv",header=TRUE,as.is=TRUE)
TenHCoWorpredergm<-TenHCoWorpredergm[-c(19:23),]

pred.vect = rep(0,88)
for (i in 2:89) {
  pred.vect[i-1]=estoprob(sum(as.numeric(TenHCoWorpredergm[,i])*TenHCoWorcoefs))
}

require(RColorBrewer)
colors <- brewer.pal(5,"Set1")
plot(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)],ylim=c(0,1),col=colors[1],lwd=2,pch=1,ylab="Probability of Tie",xlab="Model Specifications",main="Probability of a Support Tie Between Two Hindu Women in Tenpatti",xaxt="n")
points(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)+1],ylim=c(0,1),col=colors[2],lwd=2,pch=6)

abline(0, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.1, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.2, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.3, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.4, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.5, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.6, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.7, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.8, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.9, 0, lwd=1,col="lightgrey",lty="dotted")
abline(1, 0, lwd=1,col="lightgrey",lty="dotted")

labloc=c(1,2,3,4,5,6,7,8,9,10,11)
lab=c("Diff caste","Same caste","+Mutual","+GWESP=1","GWESP=2","GWESP=3","Close kin","+Mutual","+GWESP=1","GWESP=2","GWESP=3")
axis(1,at=labloc,labels=FALSE)
text(x=labloc,y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[4]/3),labels=lab,cex=0.6,pos=3,adj=2,xpd=TRUE)
legend("topleft",c("No Co-Attendance","Co-Attendance"),cex=0.6,col=c(colors[c(1,2)]),pch=c(1,6))


## Plotting ERGM with Co-Fest
TenHCoFestpredergm<-read.csv("TenHPredictERGMCollRit.csv",header=TRUE,as.is=TRUE)
TenHCoFestpredergm<-TenHCoFestpredergm[-c(18,21:23),]

TenHCoFestcoefs <- modelTenSupHindu_cofest$coef

pred.vect = rep(0,88)
for (i in 2:89) {
  pred.vect[i-1]=estoprob(sum(as.numeric(TenHCoFestpredergm[,i])*TenHCoFestcoefs))
}


require(RColorBrewer)
colors <- brewer.pal(5,"Set1")
plot(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)],ylim=c(0,1),col=colors[1],lwd=2,pch=1,ylab="Probability of Tie",xlab="Model Specifications",main="Probability of a Support Tie Between Two Hindu Women in Tenpatti",xaxt="n")
points(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)+2],ylim=c(0,1),col=colors[2],lwd=2,pch=2)
points(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)+3],ylim=c(0,1),col=colors[3],lwd=2,pch=3)
points(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)+4],ylim=c(0,1),col=colors[4],lwd=2,pch=4)

abline(0, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.1, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.2, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.3, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.4, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.5, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.6, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.7, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.8, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.9, 0, lwd=1,col="lightgrey",lty="dotted")
abline(1, 0, lwd=1,col="lightgrey",lty="dotted")

labloc=c(1,2,3,4,5,6,7,8,9,10,11)
lab=c("Diff caste","Same caste","+Mutual","+GWESP=1","GWESP=2","GWESP=3","Close kin","+Mutual","+GWESP=1","GWESP=2","GWESP=3")
axis(1,at=labloc,labels=FALSE)
text(x=labloc,y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[4]/3),labels=lab,cex=0.6,pos=3,adj=2,xpd=TRUE)
legend("topleft",c("No Co-Participation","Vow Procession","Mulaipari","Both"),cex=0.6,col=c(colors[c(1:4)]),pch=c(1:4))


or <- exp( modelTenSupHindu_cofest$coef ) 
ste <- sqrt( diag( modelTenSupHindu_cofest$covar ) ) 
lci <- exp( modelTenSupHindu_cofest$coef-1.96*ste ) 
uci <- exp( modelTenSupHindu_cofest$coef+1.96*ste ) 
oddsratios <- rbind( round( lci,digits = 4 ),round( or,digits = 4 ),round( uci,digits = 4 ) ) 
oddsratios <- t( oddsratios ) 
colnames( oddsratios ) <- c( "Lower","OR","Upper" )
oddsratios
teststat <- modelTenSupHindu_cofest$coef/ste 
teststats <- rbind( round( teststat,digits = 4 )) 
teststats <- t( teststats ) 
colnames( teststats ) <- c("Wald")
teststats

TenResults<-cbind(summary(modelTenSupHindu_cofest)$coefs[1],summary(modelTenSupHindu_cofest)$coefs[2],or,summary(modelTenSupHindu_cofest)$coefs[4])
xtable(TenResults,digits=c(3,3,3,3,4))


## Plotting ERGM with All
TenHCoAllpredergm<-read.csv("TenHPredictERGMCollRit.csv",header=TRUE,as.is=TRUE)

TenHCoAllcoefs <- modelTenSupHindu_all$coef

pred.vect = rep(0,88)
for (i in 2:89) {
  pred.vect[i-1]=estoprob(sum(as.numeric(TenHCoAllpredergm[,i])*TenHCoAllcoefs))
}


require(RColorBrewer)
colors <- brewer.pal(8,"Set1")
plot(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)],ylim=c(0,1),col=colors[1],lwd=2,pch=1,ylab="Probability of Tie",xlab="Model Specifications",main="Probability of a Support Tie Between Two Hindu Women in Tenpatti",xaxt="n")
points(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)+1],ylim=c(0,1),col=colors[2],lwd=2,pch=2)
points(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)+2],ylim=c(0,1),col=colors[3],lwd=2,pch=3)
points(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)+3],ylim=c(0,1),col=colors[4],lwd=2,pch=4)
points(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)+4],ylim=c(0,1),col=colors[5],lwd=2,pch=5)
points(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)+5],ylim=c(0,1),col=colors[6],lwd=2,pch=6)
points(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)+6],ylim=c(0,1),col=colors[7],lwd=2,pch=7)
points(pred.vect[c(1,9,17,25,33,41,49,57,65,73,81)+7],ylim=c(0,1),col=colors[8],lwd=2,pch=8)


abline(0, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.1, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.2, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.3, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.4, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.5, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.6, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.7, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.8, 0, lwd=1,col="lightgrey",lty="dotted")
abline(0.9, 0, lwd=1,col="lightgrey",lty="dotted")
abline(1, 0, lwd=1,col="lightgrey",lty="dotted")

labloc=c(1,2,3,4,5,6,7,8,9,10,11)
lab=c("Diff caste","Same caste","+Mutual","+GWESP=1","GWESP=2","GWESP=3","Close kin","+Mutual","+GWESP=1","GWESP=2","GWESP=3")
axis(1,at=labloc,labels=FALSE)
text(x=labloc,y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[4]/3),labels=lab,cex=0.6,pos=3,adj=2,xpd=TRUE)
legend("topleft",c("No Co-Participation","Pujai","Vow Procession","Mulaipari","Vow&Mulaipari","Pujai&Vow","Pujai&Mulaipari","All"),cex=0.6,col=c(colors[c(1:8)]),pch=c(1:8))

or <- exp( modelTenSupHindu_all$coef ) 
ste <- sqrt( diag( modelTenSupHindu_all$covar ) ) 
lci <- exp( modelTenSupHindu_all$coef-1.96*ste ) 
uci <- exp( modelTenSupHindu_all$coef+1.96*ste ) 
oddsratios <- rbind( round( lci,digits = 4 ),round( or,digits = 4 ),round( uci,digits = 4 ) ) 
oddsratios <- t( oddsratios ) 
colnames( oddsratios ) <- c( "Lower","OR","Upper" )
oddsratios
teststat <- modelTenSupHindu_all$coef/ste 
teststats <- rbind( round( teststat,digits = 4 )) 
teststats <- t( teststats ) 
colnames( teststats ) <- c("Wald")
teststats

TenResults<-cbind(summary(modelTenSupHindu_all)$coefs[1],summary(modelTenSupHindu_all)$coefs[2],or,summary(modelTenSupHindu_all)$coefs[4])
xtable(TenResults,digits=c(3,3,3,3,4))



#############################################################################
#############################################################################
## 3. COLLECTIVE WORSHIP EXCESS EDGES, DENSITY, TRANSITIVITY, RECIPROCITY  ##
#############################################################################
#############################################################################

## Calculating the 'excess edges,' density, transitivity, and reciprocity for each of the collective worship subgraphs
## To see if increases in any of these measures of cohesion are greater than we would expect, I additionally calculate these same measures on 10,000 hypothetical subgraphs with the same number of nodes/individuals, pulled from the larger set of potential participants
## The comparisons are as follows:
## (1) The subgraph of Hindu residents to the whole village
## (2) The subgraph of monthly puja participants to all Hindu residents
## (3) The subgraph of annual festival participants to all Hindu residents
## (4) The subgraph of vow procession participants to all Hindu residents (not reported in the paper)
## (5) The subgraph of vow procession participants to all annual festival participants (reported in the paper)

## creating additional networks that include only those individuals who participated in each collective ritual
snaTenHSupMariFest2012<-delete.vertices(snaTenHSup,V(snaTenHSup)[get.vertex.attribute(snaTenHSup,name="Fest2012")==0])
snaTenHSupMariProcession<-delete.vertices(snaTenHSup,V(snaTenHSup)[get.vertex.attribute(snaTenHSup,name="Procession")==0])
snaTenHSupMariMulaipari<-delete.vertices(snaTenHSup,V(snaTenHSup)[get.vertex.attribute(snaTenHSup,name="Mulaipari")==0])

snaTenHSupMariAtt<-delete.vertices(snaTenHSup,V(snaTenHSup)[get.vertex.attribute(snaTenHSup,name="MariAtt")==0])


#########################
## All Hindu Residents ##
#########################

## First, seeing if the subgraph of Hindu residents is more cohesive than the whole village

## The relevant comparison here is to the whole village

## Creating a matrix of all dyads where I'm in a sense weighting the existence of each edge by the overall likelihood of that edge, given the outdegree of i and the indegree of j (over the total number of edges)
g3<-get.adjacency(snaTenSup)
for (i in 1:length(V(snaTenSup)$name)){
  for (j in 1:length(V(snaTenSup)$name)){
    g3[i,j]=are.connected(snaTenSup,i,j)*1 - ((degree(snaTenSup,i,mode="out")*degree(snaTenSup,j,mode="in"))/ecount(snaTenSup))
  }
}

## Excess edges measure
sum(g3[V(snaTenSup)$Caste!="RCYaathavar",V(snaTenSup)$Caste!="RCYaathavar"])

## Now creating the comparisons
excessedgesHindu = rep(0,10000)
for (i in 1:10000){
  V(snaTenSup)$HinduRandom = sample(V(snaTenSup)$Religion,362)
  excessedgesHindu[i] = sum(g3[V(snaTenSup)$HinduRandom==0,V(snaTenSup)$HinduRandom==0])
}

TenIndivIDs<-V(snaTenSup)$name

denHindu<-rep(0,10000)
transHindu<-rep(0,10000)
recipHindu<-rep(0,10000)

for (i in 1:10000){
  denHindu[i] = graph.density(induced.subgraph(snaTenSup,vids=sample(TenIndivIDs,248)))
  transHindu[i] = transitivity(induced.subgraph(snaTenSup,vids=sample(TenIndivIDs,248)))
  recipHindu[i] = reciprocity(induced.subgraph(snaTenSup,vids=sample(TenIndivIDs,248)))
}

HinduSim<-cbind(excessedgesHindu,denHindu,transHindu,recipHindu)
HinduSim<-as.data.frame(HinduSim)

1-ecdf(HinduSim$excessedgesHindu)(sum(g3[V(snaTenSup)$Religion==0,V(snaTenSup)$Religion==0]))
1-ecdf(HinduSim$denHindu)(graph.density(snaTenHSup))
1-ecdf(HinduSim$transHindu)(transitivity(snaTenHSup))
1-ecdf(HinduSim$recipHindu)(reciprocity(snaTenHSup))

AllHindu<-cbind(sum(g3[V(snaTenSup)$Religion==0,V(snaTenSup)$Religion==0]),1-ecdf(HinduSim$excessedgesHindu)(sum(g3[V(snaTenSup)$Religion==0,V(snaTenSup)$Religion==0])),graph.density(snaTenHSup),1-ecdf(HinduSim$denHindu)(graph.density(snaTenHSup)),transitivity(snaTenHSup),1-ecdf(HinduSim$transHindu)(transitivity(snaTenHSup)),reciprocity(snaTenHSup),1-ecdf(HinduSim$recipHindu)(reciprocity(snaTenHSup)))


##################
## Monthly Puja ##
##################

## The relevant comparison here is to the Hindu residents

g<-get.adjacency(snaTenHSup)

for (i in 1:length(V(snaTenHSup)$name)){
  for (j in 1:length(V(snaTenHSup)$name)){
    g[i,j]=are.connected(snaTenHSup,i,j)*1 - ((degree(snaTenHSup,i,mode="out")*degree(snaTenHSup,j,mode="in"))/ecount(snaTenHSup))
  }
}


sum(g[V(snaTenHSup)$MariAtt==1,V(snaTenHSup)$MariAtt==1])

excessedgesMariAtt = rep(0,10000)
for (i in 1:10000){
  V(snaTenHSup)$MariAttRandom = sample(V(snaTenHSup)$MariAtt,248)
  excessedgesMariAtt[i] = sum(g[V(snaTenHSup)$MariAttRandom==1,V(snaTenHSup)$MariAttRandom==1])
}

TenHIndivIDs<-V(snaTenHSup)$name

denMariAtt<-rep(0,10000)
transMariAtt<-rep(0,10000)
recipMariAtt<-rep(0,10000)

for (i in 1:10000){
  denMariAtt[i] = graph.density(induced.subgraph(snaTenHSup,vids=sample(TenHIndivIDs,123)))
  transMariAtt[i] = transitivity(induced.subgraph(snaTenHSup,vids=sample(TenHIndivIDs,123)))
  recipMariAtt[i] = reciprocity(induced.subgraph(snaTenHSup,vids=sample(TenHIndivIDs,123)))
}

MariAttSim<-cbind(excessedgesMariAtt,denMariAtt,transMariAtt,recipMariAtt)
MariAttSim<-as.data.frame(MariAttSim)

1-ecdf(MariAttSim$excessedgesMariAtt)(sum(g[V(snaTenHSup)$MariAtt==1,V(snaTenHSup)$MariAtt==1]))
1-ecdf(MariAttSim$denMariAtt)(graph.density(snaTenHSupMariAtt))
1-ecdf(MariAttSim$transMariAtt)(transitivity(snaTenHSupMariAtt))
1-ecdf(MariAttSim$recipMariAtt)(reciprocity(snaTenHSupMariAtt))

MariAtt<-cbind(sum(g[V(snaTenHSup)$MariAtt==1,V(snaTenHSup)$MariAtt==1]),1-ecdf(MariAttSim$excessedgesMariAtt)(sum(g[V(snaTenHSup)$MariAtt==1,V(snaTenHSup)$MariAtt==1])),graph.density(snaTenHSupMariAtt),1-ecdf(MariAttSim$denMariAtt)(graph.density(snaTenHSupMariAtt)),transitivity(snaTenHSupMariAtt),1-ecdf(MariAttSim$transMariAtt)(transitivity(snaTenHSupMariAtt)),reciprocity(snaTenHSupMariAtt),1-ecdf(MariAttSim$recipMariAtt)(reciprocity(snaTenHSupMariAtt)))


#####################
## Annual Festival ##
#####################

## As above, the relevant comparison here is to the Hindu residents

sum(g[V(snaTenHSup)$Fest2012==1,V(snaTenHSup)$Fest2012==1])

excessedgesMariFest = rep(0,10000)
for (i in 1:10000){
  V(snaTenHSup)$MariFestRandom = sample(V(snaTenHSup)$Fest2012,248)
  excessedgesMariFest[i] = sum(g[V(snaTenHSup)$MariFestRandom==1,V(snaTenHSup)$MariFestRandom==1])
}

denMariFest<-rep(0,10000)
transMariFest<-rep(0,10000)
recipMariFest<-rep(0,10000)

for (i in 1:10000){
  denMariFest[i] = graph.density(induced.subgraph(snaTenHSup,vids=sample(TenHIndivIDs,64)))
  transMariFest[i] = transitivity(induced.subgraph(snaTenHSup,vids=sample(TenHIndivIDs,64)))
  recipMariFest[i] = reciprocity(induced.subgraph(snaTenHSup,vids=sample(TenHIndivIDs,64)))
}

MariFestSim<-cbind(excessedgesMariFest,denMariFest,transMariFest,recipMariFest)
MariFestSim<-as.data.frame(MariFestSim)

1-ecdf(MariFestSim$excessedgesMariFest)(sum(g[V(snaTenHSup)$Fest2012==1,V(snaTenHSup)$Fest2012==1]))
1-ecdf(MariFestSim$denMariFest)(graph.density(snaTenHSupMariFest2012))
1-ecdf(MariFestSim$transMariFest)(transitivity(snaTenHSupMariFest2012))
1-ecdf(MariFestSim$recipMariFest)(reciprocity(snaTenHSupMariFest2012))

MariFest<-cbind(sum(g[V(snaTenHSup)$Fest2012==1,V(snaTenHSup)$Fest2012==1]),1-ecdf(MariFestSim$excessedgesMariFest)(sum(g[V(snaTenHSup)$Fest2012==1,V(snaTenHSup)$Fest2012==1])),graph.density(snaTenHSupMariFest2012),1-ecdf(MariFestSim$denMariFest)(graph.density(snaTenHSupMariFest2012)),transitivity(snaTenHSupMariFest2012),1-ecdf(MariFestSim$transMariFest)(transitivity(snaTenHSupMariFest2012)),reciprocity(snaTenHSupMariFest2012),1-ecdf(MariFestSim$recipMariFest)(reciprocity(snaTenHSupMariFest2012)))


####################
## Vow Procession ##
####################

## First, comparing to the Hindu residents (not reported in the paper)

sum(g[V(snaTenHSup)$Procession==1,V(snaTenHSup)$Procession==1])


excessedgesMariProcession = rep(0,10000)
for (i in 1:10000){
  V(snaTenHSup)$MariProcessionRandom = sample(V(snaTenHSup)$Procession,248)
  excessedgesMariProcession[i] = sum(g[V(snaTenHSup)$MariProcessionRandom==1,V(snaTenHSup)$MariProcessionRandom==1])
}


denMariProcession<-rep(0,10000)
transMariProcession<-rep(0,10000)
recipMariProcession<-rep(0,10000)


for (i in 1:10000){
  denMariProcession[i] = graph.density(induced.subgraph(snaTenHSup,vids=sample(TenHIndivIDs,28)))
  transMariProcession[i] = transitivity(induced.subgraph(snaTenHSup,vids=sample(TenHIndivIDs,28)))
  recipMariProcession[i] = reciprocity(induced.subgraph(snaTenHSup,vids=sample(TenHIndivIDs,28)))
}

MariProcessionSim<-cbind(excessedgesMariProcession,denMariProcession,transMariProcession,recipMariProcession)
MariProcessionSim<-as.data.frame(MariProcessionSim)

1-ecdf(MariProcessionSim$excessedgesMariProcession)(sum(g[V(snaTenHSup)$Procession==1,V(snaTenHSup)$Procession==1]))
1-ecdf(MariProcessionSim$denMariProcession)(graph.density(snaTenHSupMariProcession))
1-ecdf(MariProcessionSim$transMariProcession)(transitivity(snaTenHSupMariProcession))
1-ecdf(MariProcessionSim$recipMariProcession)(reciprocity(snaTenHSupMariProcession))

MariProcession<-cbind(sum(g[V(snaTenHSup)$Procession==1,V(snaTenHSup)$Procession==1]),1-ecdf(MariProcessionSim$excessedgesMariProcession)(sum(g[V(snaTenHSup)$Procession==1,V(snaTenHSup)$Procession==1])),graph.density(snaTenHSupMariProcession),1-ecdf(MariProcessionSim$denMariProcession)(graph.density(snaTenHSupMariProcession)),transitivity(snaTenHSupMariProcession),1-ecdf(MariProcessionSim$transMariProcession)(transitivity(snaTenHSupMariProcession)),reciprocity(snaTenHSupMariProcession),1-ecdf(MariProcessionSim$recipMariProcession)(reciprocity(snaTenHSupMariProcession)))


### Now, the comparison is to all festival participants (reported in the paper)

g2<-get.adjacency(snaTenHSupMariFest2012)

for (i in 1:length(V(snaTenHSupMariFest2012)$name)){
  for (j in 1:length(V(snaTenHSupMariFest2012)$name)){
    g2[i,j]=are.connected(snaTenHSupMariFest2012,i,j)*1 - ((degree(snaTenHSupMariFest2012,i,mode="out")*degree(snaTenHSupMariFest2012,j,mode="in"))/ecount(snaTenHSupMariFest2012))
  }
}

sum(g2[V(snaTenHSupMariFest2012)$Procession==1,V(snaTenHSupMariFest2012)$Procession==1])


excessedgesMariProcession1 = rep(0,10000)
for (i in 1:10000){
  V(snaTenHSupMariFest2012)$MariProcessionRandom1 = sample(V(snaTenHSupMariFest2012)$Procession,64)
  excessedgesMariProcession1[i] = sum(g2[V(snaTenHSupMariFest2012)$MariProcessionRandom1==1,V(snaTenHSupMariFest2012)$MariProcessionRandom1==1])
}

denMariProcession1<-rep(0,10000)
transMariProcession1<-rep(0,10000)
recipMariProcession1<-rep(0,10000)

TenFestIndivIDs<-V(snaTenHSupMariFest2012)$name

for (i in 1:10000){
  denMariProcession1[i] = graph.density(induced.subgraph(snaTenHSupMariFest2012,vids=sample(TenFestIndivIDs,28)))
  transMariProcession1[i] = transitivity(induced.subgraph(snaTenHSupMariFest2012,vids=sample(TenFestIndivIDs,28)))
  recipMariProcession1[i] = reciprocity(induced.subgraph(snaTenHSupMariFest2012,vids=sample(TenFestIndivIDs,28)))
}

MariProcessionSim1<-cbind(excessedgesMariProcession1,denMariProcession1,transMariProcession1,recipMariProcession1)
MariProcessionSim1<-as.data.frame(MariProcessionSim1)

1-ecdf(MariProcessionSim1$excessedgesMariProcession1)(sum(g2[V(snaTenHSupMariFest2012)$Procession==1,V(snaTenHSupMariFest2012)$Procession==1]))
1-ecdf(MariProcessionSim1$denMariProcession1)(graph.density(snaTenHSupMariProcession))
1-ecdf(MariProcessionSim1$transMariProcession1)(transitivity(snaTenHSupMariProcession))
1-ecdf(MariProcessionSim1$recipMariProcession1)(reciprocity(snaTenHSupMariProcession))

MariProcession1<-cbind(sum(g2[V(snaTenHSupMariFest2012)$Procession==1,V(snaTenHSupMariFest2012)$Procession==1]),1-ecdf(MariProcessionSim1$excessedgesMariProcession1)(sum(g2[V(snaTenHSupMariFest2012)$Procession==1,V(snaTenHSupMariFest2012)$Procession==1])),graph.density(snaTenHSupMariProcession),1-ecdf(MariProcessionSim1$denMariProcession1)(graph.density(snaTenHSupMariProcession)),transitivity(snaTenHSupMariProcession),1-ecdf(MariProcessionSim1$transMariProcession1)(transitivity(snaTenHSupMariProcession)),reciprocity(snaTenHSupMariProcession),1-ecdf(MariProcessionSim1$recipMariProcession1)(reciprocity(snaTenHSupMariProcession)))

## Combining results (these are what is reported in table 2 in the paper)
rbind(AllHindu,MariAtt,MariFest,MariProcession1)

## To account for multiple comparisons, we can also adjust the p-values to reflect the multiple comparisons. Here, I use the Bonferroni method
p.adjust(AllHindu[,c(2,4,6,8)],method="bonferroni")
p.adjust(MariAtt[,c(2,4,6,8)],method="bonferroni")
p.adjust(MariFest[,c(2,4,6,8)],method="bonferroni")
p.adjust(MariProcession[,c(2,4,6,8)],method="bonferroni")
p.adjust(MariProcession1[,c(2,4,6,8)],method="bonferroni")

##########################
##########################
## 4. RELIGIOUS ALTERS  ##
##########################
##########################

## Code used to count the number of a person's ties that are to people of other religions
## Religions are coded as 0 = Hindu, 1 = Catholic, 2 = Muslim, 3 = Protestant (but, no Tenpatti Hindu residents named Protestants as alters, so that level is not included)

get_other <- function(graph, attribute) {

  mat <- get.adjacency(graph)

  attr_levels = get.vertex.attribute(graph, attribute)
  attr_counts<-data.frame(0,0,0)

  num_levels = length(unique(attr_levels))

  for (ego in 1:nrow(mat)) {

    # initialize actor-specific variables
    alter_attr_counts = rep(0, num_levels)
    num_alters_this_ego = 0

    for (alter in 1:ncol(mat)) {

      # only examine alters that are actually tied to ego
      if (mat[ego, alter] == 1) {

        num_alters_this_ego = num_alters_this_ego + 1

        # get the alter's level on the attribute
        alter_attr = get.vertex.attribute(graph, attribute, alter)

        # increment the count of alters with this level
        # of the attribute by 1
        alter_attr_counts[alter_attr + 1] =
          alter_attr_counts[alter_attr + 1] + 1
      }
    }

    # append new ego to list
    attr_counts<-rbind(attr_counts,alter_attr_counts)
  }
  return(attr_counts)
}

## Note: this will take some time
otherRelig<-get_other(snaTenHSupFull,"Religion")

otherRelig<-cbind(V(snaTenHSupFull)$name, V(snaTenHSupFull)$Religion,otherRelig[-1,],degree(snaTenHSupFull,mode="out"))
names(otherRelig)<-c("IndivID","OwnRelig","AlterHindu","AlterCatholic","AlterMuslim","OutDegree")
otherRelig<-subset(otherRelig,OutDegree>0)

## Adding in other covariates that might influence the number of religious others
## Where possible, centering and scaling the variables
otherRelig<-merge(otherRelig,Ten[,c(1,3:7,14:17)],by="IndivID")
otherRelig$n_other<-otherRelig$AlterCatholic+otherRelig$AlterMuslim
otherRelig$is_SC<-otherRelig$Caste=="Pallar" | otherRelig$Caste == "Arundhathiyar"
otherRelig$is_male<-otherRelig$Gender=="Male"
otherRelig$eduyears <- (otherRelig$EduYears-mean(otherRelig$EduYears))/sd(otherRelig$EduYears)
otherRelig$age <- (otherRelig$Age-mean(otherRelig$Age))/sd(otherRelig$Age)
otherRelig$Age2 <- otherRelig$Age^2
otherRelig$age2 <- (otherRelig$Age2-mean(otherRelig$Age2))/sd(otherRelig$Age2)
otherRelig$wealth <- (otherRelig$HouseholdWealthCalc-mean(otherRelig$HouseholdWealthCalc))/sd(otherRelig$HouseholdWealthCalc)

otherRelig <- otherRelig[,c(1,6,12:20,22,23)]
names(otherRelig)<-c("indiv","n_alters","is_participant_fest","is_participant_puja","is_vow","is_mulai","n_other","is_SC","is_male","eduyears","age","age2","wealth")
otherRelig$indiv <- factor(otherRelig$indiv)
otherRelig$is_SC<-as.numeric(otherRelig$is_SC)
otherRelig$is_male<-as.numeric(otherRelig$is_male)

## some strange aspect of R data formatting throws an error with this current formatting; as a work around...
write.csv("otherRelig.csv")
otherRelig<-read.csv("otherRelig.csv",as.is=TRUE)
otherRelig$indiv<-factor(otherRelig$indiv)

require(rethinking)

## Binomial predicting religious others, and including participation in the monthly puja
m_wor<-map2stan(
  alist(
    n_other ~ dbinom(n_alters, p),
    logit(p) <- alpha + tau*z[indiv] + beta_wor*is_participant_puja + beta_edu*eduyears + beta_age*age + beta_age2*age2 + beta_male*is_male + beta_SC*is_SC+ beta_wealth*wealth,
    z[indiv] ~ dnorm(0, 1),
    alpha ~ dnorm(0, 1),
    tau ~ dcauchy(0, 1),
    beta_wor ~ dnorm(0, 1),
    beta_edu ~ dnorm(0, 1),
    beta_age ~ dnorm(0, 1),
    beta_age2 ~ dnorm(0, 1),
    beta_male ~ dnorm(0, 1),
    beta_SC ~ dnorm(0, 1),
    beta_wealth ~ dnorm(0, 1)
  ),
  data=otherRelig, chains=4, iter=2000, warmup=1000, control=list(max_treedepth=20),constraints=list(tau="lower=0")
)

## here including participation in the annual festival
m_fest<-map2stan(
  alist(
    n_other ~ dbinom(n_alters, p),
    logit(p) <- alpha + tau*z[indiv] + beta_fest*is_participant_fest + beta_edu*eduyears + beta_age*age + beta_age2*age2 + beta_male*is_male + beta_SC*is_SC+ beta_wealth*wealth,
    z[indiv] ~ dnorm(0, 1),
    alpha ~ dnorm(0, 1),
    tau ~ dcauchy(0, 1),
    beta_fest ~ dnorm(0, 1),
    beta_edu ~ dnorm(0, 1),
    beta_age ~ dnorm(0, 1),
    beta_age2 ~ dnorm(0, 1),
    beta_male ~ dnorm(0, 1),
    beta_SC ~ dnorm(0, 1),
    beta_wealth ~ dnorm(0, 1)
  ),
  data=otherRelig, chains=4, iter=2000, warmup=1000, control=list(max_treedepth=20),constraints=list(tau="lower=0")
)

## breaking out festival participation into the vow procession and the mulaipari procession
m_fest_levels<-map2stan(
  alist(
    n_other ~ dbinom(n_alters, p),
    logit(p) <- alpha + tau*z[indiv] + beta_mulai*is_mulai+beta_vow*is_vow + beta_edu*eduyears + beta_age*age + beta_age2*age2 + beta_male*is_male + beta_SC*is_SC+ beta_wealth*wealth,
    z[indiv] ~ dnorm(0, 1),
    alpha ~ dnorm(0, 1),
    tau ~ dcauchy(0, 1),
    beta_mulai ~ dnorm(0, 1),
    beta_vow ~ dnorm(0, 1),
    beta_edu ~ dnorm(0, 1),
    beta_age ~ dnorm(0, 1),
    beta_age2 ~ dnorm(0, 1),
    beta_male ~ dnorm(0, 1),
    beta_SC ~ dnorm(0, 1),
    beta_wealth ~ dnorm(0, 1)
  ),
  data=otherRelig, chains=4, iter=2000, warmup=1000, control=list(max_treedepth=20),constraints=list(tau="lower=0")
)

precis(m_wor,prob=0.95,digits=3)
precis(m_fest,prob=0.95,digits=3)
precis(m_fest_levels,prob=0.95,digits=3)



## Creating summary plots of the models
b <- precis(m_wor, prob=0.65)
b2 <- precis(m_wor, prob=0.95)
points <- cbind(b[,1:4], b2[,3:4])
par(mar=c(5,7,4,1), mfrow=c(1,2), lend=3)
plot(points[,1], c(9:1), xlim=c(-2.75, 1.25), ylim=c(0.5, 9.5), xlab=("Estimate"),ylab="", yaxt="n", main="Monthly Worship", cex.axis=0.8, type="n")
abline(v=0, lty=1, lwd=0.75)
for (i in 9:1) {
  abline(h=10-i, lty=3, lwd=0.5)
  lines(c(points[10-i,3], points[10-i,4]), c(i, i), lwd=5)
  lines(c(points[10-i,5], points[10-i,6]), c(i, i), lwd=2)
}
points(points[,1], c(9:1), pch=3, cex=1.5)
axis(2, at = c(9:1), labels = c("Alpha", "Tau", "Collective Ritual", "Education", "Age", "Age squared", "Gender (F=0)", "Caste (SC=1)", "Wealth"), cex.axis=0.9, tick = FALSE, line=FALSE, las=TRUE)
b <- precis(m_fest, prob=0.65)
b2 <- precis(m_fest, prob=0.95)
points <- cbind(b[,1:4], b2[,3:4])
par(mar=c(5,1,4,7), lend=3)
plot(points[,1], c(9:1), xlim=c(-2.75, 1.25), ylim=c(0.5, 9.5), xlab=("Estimate"),ylab="", yaxt="n", main="Annual Festival", cex.axis=0.8, type="n")
abline(v=0, lty=1, lwd=0.75)
for (i in 9:1) {
  abline(h=10-i, lty=3, lwd=0.5)
  lines(c(points[10-i,3], points[10-i,4]), c(i, i), lwd=5)
  lines(c(points[10-i,5], points[10-i,6]), c(i, i), lwd=2)
}
points(points[,1], c(9:1), pch=3, cex=1.5)




##############################
##############################
## 5. PLOTTING THE NETWORK  ##
##############################
##############################

castecol = get.vertex.attribute(snaTenHSup,"Caste")
unique(castecol)
require(RColorBrewer)
colors <- brewer.pal(12,"Paired")
castecol[castecol == "Pallar"] = colors[1]
castecol[castecol == "HinduYaathavar"] = colors[8]
castecol[castecol == "Agamudaiyaan"] = colors[5]
castecol[castecol == "Arundhathiyar"] = colors[2]
castecol[castecol == "Aasaari"] = colors[11]
castecol[castecol == "Naayakkar"] = "#666666"
castecol[castecol == "Kallar"] = colors[6]
castecol[castecol == "Kulaalar"] = colors[4]
castecol[castecol == "HinduVellalar"] = "#666666"

require(scales)
evcent <- rescale(evcent(snaTenHSup,directed=TRUE)$vector,to=c(2,8))
snlayout <- layout.fruchterman.reingold(snaTenHSup)

maricol = V(snaTenHSup)$MariAtt
maricol[maricol == 0] = "white"
maricol[maricol == 1] = "dodgerblue"

festcol = V(snaTenHSup)$PMB
festcol[festcol == "None"] = "white"
festcol[festcol == "Mulaipari"] = "dodgerblue"
festcol[festcol == "Procession"] = "red"
festcol[festcol == "Both"] = "mediumorchid"



plot.igraph(snaTenHSup,edge.width=E(snaTenHSup)$SupSum/3,edge.arrow.size=.1,vertex.size=evcent,vertex.color=castecol,vertex.label=NA,layout=snlayout,edge.curved=TRUE)
legend("topright",c("Acari","Akamutaiyar","Aruntatiyar","Yatavar","Kallar","Kulalar","Pallar","Rare"),fill=c(colors[c(11,5,2,8,6,4,1)],"#666666"))

plot.igraph(snaTenHSup,edge.width=E(snaTenHSup)$SupSum/3,edge.arrow.size=.1,vertex.size=evcent,vertex.color=maricol,vertex.label=NA,layout=snlayout,edge.curved=TRUE)
legend("topright",c("No","Yes"),fill=c("white","dodgerblue"))

plot.igraph(snaTenHSup,edge.width=E(snaTenHSup)$SupSum/3,edge.arrow.size=.1,vertex.size=evcent,vertex.color=festcol,vertex.label=NA,layout=snlayout,edge.curved=TRUE)
legend("topright",c("No","Mulaipari","Vow Procession","Both"),fill=c("white","dodgerblue","red","mediumorchid"))



