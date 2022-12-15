
# Code used for final assignment of From Problem to Analysis
# By Robin Vloeberghs


# Raw data comes from: Rouault, Seow, Gillan and Fleming. (2018) Biological Psychiatry
# Psychiatric symptom dimensions are associated with dissociable shifts in metacognition but not task performance.
# https://github.com/marionrouault/webstudy


# 1. ("ME_phase2_excludqnadata_all.mat") questionnaire data
# 2. ("ME_phase2_excludanalyseddat_all.mat") task performance data




##########
## CLEAR ALL
rm(list=ls())

##########
## LOADING LIBRARIES/TOOL
# loading tools
library(ggplot2) # for plotting graphs
library(gridExtra) # for ploting graphs
library(lme4) # for linear regression functions
library(plyr) # for collapse-and-mean functions like ddply
library(psych)
library(GPArotation)
library(paran)
library(reshape)
library(polycor)
library(nFactors)
library(R.matlab)
library(reshape)
library(doBy)
library(rstudioapi)
library(effects)
library(tidyr)
library(car)
library(jtools)
library(ppcor)
library(ggResidpanel)

options(scipen = 999) # to display variable quantities in decimals (not in scientific notation format)

##########
## SET DIRECTORY
curdir <- dirname(getSourceEditorContext()$path)
setwd(curdir)

##########
## LOADING DATA
qnData = readMat("ME_phase2_excludqnadata_all.mat") # load questionnaire data
taskData = readMat("ME_phase2_excludanalyseddat_all.mat") # load task performance data
confModel = read.csv('fitted_parameters_rouault.csv') # load ConfModel params

##########
## CREATE EMPTY OBJECTS
# create objects for variables from task performance data
id<-matrix(0,length(taskData$analyseddata),1) # subject id
age<-matrix(0,length(taskData$analyseddata),1)
gender<-matrix(0,length(taskData$analyseddata),1)
accuracy<-matrix(0,length(taskData$analyseddata),1) # accuracy
confMean<-matrix(0,length(taskData$analyseddata),1) # mean confidence

# create objects for variables from task questionnaire data
qnid<-matrix(0,length(qnData$allqna),1) # subject id
zung<-matrix(0,length(qnData$allqna),1)
anxiety<-matrix(0,length(qnData$allqna),1)
ocir<-matrix(0,length(qnData$allqna),1)
leb<-matrix(0,length(qnData$allqna),1)
iq<-matrix(0,length(qnData$allqna),1)
bis<-matrix(0,length(qnData$allqna),1)
schizo<-matrix(0,length(qnData$allqna),1)
eat<-matrix(0,length(qnData$allqna),1)
apathy<-matrix(0,length(qnData$allqna),1)
alcohol<-matrix(0,length(qnData$allqna),1)

##########
## EXTRACTING DATA
# extracting data from allqna data file
# loop over for all subjects
for (i in 1:length(qnData$allqna)) 
{
  qnid[i] = qnData$allqna[[i]][[1]][,,1]$id
  zung[i] = qnData$allqna[[i]][[1]][,,1]$zung[,,1]$score #first brackets is subject number
  anxiety[i] = qnData$allqna[[i]][[1]][,,1]$anxiety[,,1]$score
  ocir[i] = qnData$allqna[[i]][[1]][,,1]$ocir[,,1]$score
  leb[i] = qnData$allqna[[i]][[1]][,,1]$leb[,,1]$score
  iq[i] = qnData$allqna[[i]][[1]][,,1]$iq[,,1]$score
  bis[i] = qnData$allqna[[i]][[1]][,,1]$bis[,,1]$score[,,1]$total
  schizo[i] = qnData$allqna[[i]][[1]][,,1]$schizo[,,1]$score[,,1]$total
  eat[i] = qnData$allqna[[i]][[1]][,,1]$eat[,,1]$score[,,1]$total
  apathy[i] = qnData$allqna[[i]][[1]][,,1]$apathy[,,1]$score
  alcohol[i] = qnData$allqna[[i]][[1]][,,1]$alcohol[,,1]$score
}

# extracting data from analysed data
# loop over for all subjects
# taskData$analyseddata[[i]][[1]][,,1]$data[,5] is RT

for (i in 1:length(taskData$analyseddata))
{
  id[i]=taskData$analyseddata[[i]][[1]][,,1]$data[1,4]
  age[i] =taskData$analyseddata[[i]][[1]][,,1]$data[1,2]
  gender[i]=taskData$analyseddata[[i]][[1]][,,1]$data[1,3]
  confMean[i] = mean(taskData$analyseddata[[i]][[1]][,,1]$data[,9])
  accuracy[i] = mean(taskData$analyseddata[[i]][[1]][,,1]$data[,6])
}

alpha = confModel$alpha
beta = confModel$beta

##########
## MERGING DATA

# create dataframe to store questionnaire data
qnFrame = data.frame(qnid, anxiety, eat, apathy, alcohol, zung, ocir, leb, iq, bis, schizo)
# create dataframe to store task performance data
taskFrame = data.frame(id, age, gender, confMean, accuracy, alpha, beta)
# merge all data together into one data frame
allData = merge(taskFrame, qnFrame,by.x=c("id"), by.y=c("qnid"))


##########
## SCALING DATA
allData$age.sc = scale(allData$age)
allData$confMean.sc = scale(allData$confMean)
allData$accuracy.sc = scale(allData$accuracy)

# scaling the questionnaire scores (total score per subject)
allData$zung.sc = scale(log(allData$zung))
allData$anxiety.sc = scale(log(allData$anxiety))
allData$ocir.sc = scale(log(allData$ocir+1))
allData$leb.sc = scale(log(allData$leb+1))
allData$iq.sc = scale(allData$iq)
allData$schizo.sc = scale(log(allData$schizo+1))
allData$bis.sc = scale(log(allData$bis))
allData$eat.sc = scale(log(allData$eat+1))
allData$apathy.sc = scale(log(allData$apathy))
allData$alcohol.sc = scale(log(allData$alcohol+1))


####################
##  FACTOR ANALYSIS 
# LOAD ALL QUESTIONNAIRE (individual questions) DATA
# create objects
qnIndivid<-matrix(0,length(qnData$allqna),1)
zungAll<-matrix(0,length(qnData$allqna),length(qnData$allqna[[1]][[1]][,,1]$zung[,,1]$raw))
anxietyAll<-matrix(0,length(qnData$allqna),length(qnData$allqna[[1]][[1]][,,1]$anxiety[,,1]$raw))
ocirAll<-matrix(0,length(qnData$allqna),length(qnData$allqna[[1]][[1]][,,1]$ocir[,,1]$raw))
lebAll<-matrix(0,length(qnData$allqna),length(qnData$allqna[[1]][[1]][,,1]$leb[,,1]$raw[,,1]$avg))
bisAll<-matrix(0,length(qnData$allqna),length(qnData$allqna[[1]][[1]][,,1]$bis[,,1]$raw))
schizoAll<-matrix(0,length(qnData$allqna),length(qnData$allqna[[1]][[1]][,,1]$schizo[,,1]$raw))
eatAll<-matrix(0,length(qnData$allqna),length(qnData$allqna[[1]][[1]][,,1]$eat[,,1]$raw))
apathyAll<-matrix(0,length(qnData$allqna),length(qnData$allqna[[1]][[1]][,,1]$apathy[,,1]$raw))
alcoholAll<-matrix(0,length(qnData$allqna),length(qnData$allqna[[1]][[1]][,,1]$alcohol[,,1]$raw))

# extracting data from allqna
for (i in 1:length(qnData$allqna))
{
  qnIndivid[i,]=qnData$allqna[[i]][[1]][,,1]$id
  zungAll[i,] = qnData$allqna[[i]][[1]][,,1]$zung[,,1]$raw #first brackets is subject number
  anxietyAll[i,] = t(qnData$allqna[[i]][[1]][,,1]$anxiety[,,1]$raw)
  ocirAll[i,] = qnData$allqna[[i]][[1]][,,1]$ocir[,,1]$raw
  lebAll[i,] = (qnData$allqna[[i]][[1]][,,1]$leb[,,1]$raw[,,1]$avg)
  bisAll[i,] = qnData$allqna[[i]][[1]][,,1]$bis[,,1]$raw
  schizoAll[i,] = qnData$allqna[[i]][[1]][,,1]$schizo[,,1]$raw
  eatAll[i,]=qnData$allqna[[i]][[1]][,,1]$eat[,,1]$raw
  apathyAll[i,]=qnData$allqna[[i]][[1]][,,1]$apathy[,,1]$raw
  alcoholAll[i,]=qnData$allqna[[i]][[1]][,,1]$alcohol[,,1]$raw
}

qns = data.frame("qnid"=qnIndivid,"zung"=zungAll, "anxiety"=anxietyAll,"ocir"= ocirAll, "leb" =lebAll,"bis"= bisAll,"schizo"= schizoAll, 'alcohol'=alcoholAll,'eat'=eatAll,'apathy'=apathyAll)



# DO FACTOR ANALYSIS ON RAW QUESTIONAIRRE SCORES

# Produce correlation matrix using hetcor to account for both continuous and binary correlations
# first column is id so ignore
het.mat <- hetcor(qns[,2:length(qns)])$cor

# fa uses correlation matrix as input
# 3 factors because of semantic interpretation: see https://elifesciences.org/articles/11305
# and because of Cattell-Nelson-Gorsuch test (see below)
fa <- fa(r = het.mat, nfactors = 3, n.obs = nrow(qns), rotate = "oblimin", fm="ml", scores="regression")


# Factor scores represent estimates of common part of the variables
fa.scores <- factor.scores(x=qns[,2:length(qns)], f=fa)
scores = data.frame("id"=qns$qnid, fa.scores$scores)


# perform Cattell-Nelson-Gorsuch test for number of factors to retain based on eigenvalues of all factors
holdCNG <- nCng(fa$values, cor = TRUE, model = "factors", details = TRUE)

# join factor scores with main data matrix
colnames(scores) <- c("id", "AD", "Compul", "SW")
factorData <- merge(allData, scores,by.x=c("id"), by.y=c("id")) 


### What to report: https://www.spss-tutorials.com/apa-reporting-factor-analysis/
# Factor Loadings & Communalities: fa -> first table

n_items <- c(20,20,18,24,30,43,10,26,18)
label <- c("Depression","Generalized Anxiety","OCD","Social Anxiety","Impulsivity","Schizotypy","Alcoholism","Eating Disorders","Apathy")
labels_all <- rep(label,n_items)

communality <- fa$communality # is h2 is first table of fa,  the sum of squared factor loadings for that item
uniquenesses <- unname(fa$uniquenesses)

# Standardized loadings (pattern matrix) based upon correlation matrix
loadings <- as.data.frame(unname(unclass(print(fa$loadings,cutoff=0))))
colnames(loadings) <- c("ML1","ML2","ML3")
#colnames(loadings) <- c("Anxious-Depression","Compulsive Behavior and Intrusive Thought","Social Withdrawal")


df_plot <- cbind.data.frame(communality,uniquenesses,labels_all,loadings)
df_plot_reordered <- df_plot[order(df_plot[,'labels_all']), ] #data follows alphabetical order, consistent with legend
position <- c(1:nrow(loadings))
df_plot_reordered <- cbind(position,df_plot_reordered)


df_plot_long <- gather(df_plot_reordered, "factor", "loadings", 5:7)
df_plot_long$factor <- factor(df_plot_long$factor, levels = c("ML1", "ML2","ML3"),
                              labels = c("Factor 'Anxious-Depression'","Factor 'Compulsive Behavior and Intrusive Thought'","Factor 'Social Withdrawal'"))

colorBlindGrey <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#88CCEE")
ggplot(data=df_plot_long, aes(x=position, y=loadings,fill=labels_all)) +
  geom_bar(stat="identity", color="black" )+
  theme_minimal() +
  labs(x="Individual questionnaire items (209 items)", y="Factor loadings",fill='Questionnaires') +
  facet_wrap(~factor,ncol=1) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 25)) +
  scale_fill_manual(values=c(colorBlindGrey))

  
ggplot(data=df_plot_reordered, aes(x=position, y=communality,fill=labels_all)) +
  geom_bar(stat="identity", color="black" )+
  theme_minimal() +
  labs(x="Individual questionnaire items (209 items)", y="Communalities",fill='Questionnaires') +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_fill_manual(values=c(colorBlindGrey))



# Eigen values + Total Variance Explained
# Negative values not that uncommon (see https://stats.oarc.ucla.edu/sas/output/factor-analysis/)
eigenvalues <- fa$values
plot(eigenvalues)
component <- 1:length(fa$values)

df_eigenvalue <- cbind.data.frame(component,eigenvalues)

ggplot(data=df_eigenvalue[1:10,], aes(x=component, y=eigenvalues)) +
  geom_point(color="blue")+
  geom_line(group = 1) +
  theme_classic() +
  labs(x="Factor number", y="Eigenvalue") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

fa$Vaccounted

# Factor Correlations
fa$Phi


model_conf <- lm(confMean.sc ~ AD + Compul + SW, factorData) 
summary(model_conf)

model_alpha <- lm(alpha ~ AD + Compul + SW + accuracy.sc, factorData) 
summary(model_alpha)
vif(model_alpha)
effect_plot(model_alpha, pred = AD, interval = TRUE, plot.points = TRUE, x.label=c("Anxiety-Depression"), y.label=c("Alpha"))


model_beta <- lm(beta ~ AD + Compul + SW + accuracy.sc, factorData) 
summary(model_beta)
vif(model_beta)
resid_panel(model_beta)

plot(effect(c('AD'), model_beta))
plot(effect(c('Compul'), model_beta))

effect_plot(model_beta, pred = AD, interval = TRUE, plot.points = TRUE, x.label=c("Anxiety-Depression"), y.label=c("Beta"))
effect_plot(model_beta, pred = Compul, interval = TRUE, plot.points = TRUE, x.label=c("Compulsive Behavior and Intrusive Thoughts"), y.label=c("Beta"))


# partial correlation: beta and AD given Compul and SW
# df = n - 2 - k with k being the number of variables we conditioned upon
pcor.test(factorData$beta,factorData$AD,factorData[,c("Compul","SW","accuracy.sc")],method = c("spearman"))

# partial correlation: beta and Compul given AD and SW
pcor.test(factorData$beta,factorData$Compul,factorData[,c("AD","SW","accuracy.sc")],method = c("spearman"))


pcor.test(factorData$alpha,factorData$AD,factorData[,c("Compul","SW","accuracy.sc")],method = c("spearman"))
pcor.test(factorData$alpha,factorData$Compul,factorData[,c("AD","SW","accuracy.sc")],method = c("spearman"))




plot(factorData$alpha)
plot(factorData$beta)
cor.test(factorData$alpha,factorData$beta)



model_AD <- lm(AD ~ beta * alpha, factorData) 
summary(model_AD)
plot(effect(c('beta:alpha'), model_AD))
plot(effect(c('alpha'), model_AD))
vif(model_AD)


model_Compul <- lm(Compul ~ beta * alpha, factorData) 
summary(model_Compul)

plot(effect(c('beta'), model_Compul)) 
plot(effect(c('beta:alpha'), model_Compul))
vif(model_Compul)




