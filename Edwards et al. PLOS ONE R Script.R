#remove any lists you saved earlier
rm(list=ls())

#read in your data
SERT<-read.csv("Edwards et al. PLOS ONE GLMM Data.csv",header=T)

#find out how much data R read, to check for errors
dim(SERT)

#look at what the data looks like, to check for errors
head(SERT)
tail(SERT)
summary(SERT)
str(SERT)

#Change some variables to factor or numeric
SERT$BirdID<-as.factor(SERT$Field.Ring)
SERT$exp<-as.integer(SERT$Exploratory.score)
SERT$terr<-as.factor(SERT$Territory)
SERT$observ<-as.factor(SERT$TesterID)
SERT$Status<-as.factor(SERT$Sub.Dom)
SERT$Age<-as.numeric(SERT$Cent.Age)#mean centered and divided by 2*StDev
SERT$Age2<-as.numeric(SERT$Age2)#squared for quadratic effect
SERT$bold<-as.integer(SERT$Bold.score)
SERT$Exp.Test.Number<-as.integer(SERT$Exploration.Test.Number)
SERT$Bold.Test.Number<-as.integer(SERT$Boldness.Test.Number)
SERT$Hap1f<-as.factor(SERT$Hap1)#For overdominant model
SERT$Hap2f<-as.factor(SERT$Hap2)
SERT$Hap3f<-as.factor(SERT$Hap3)
SERT$Hap4f<-as.factor(SERT$Hap4)
SERT$Hap5f<-as.factor(SERT$Hap5)
SERT$Hap1I<-as.integer(SERT$Hap1)#For additive model
SERT$Hap2I<-as.integer(SERT$Hap2)
SERT$Hap3I<-as.integer(SERT$Hap3)
SERT$Hap4I<-as.integer(SERT$Hap4)
SERT$Hap5I<-as.integer(SERT$Hap5)

#normality
hist(SERT$exp)
qqnorm(SERT$exp)
shapiro.test(SERT$exp)#Not normally distributed- use poisson dist

table(SERT$BTO) # count variable

hist(SERT$bold)
qqnorm(SERT$bold)
shapiro.test(SERT$bold)#Not normally distributed- use poisson dist


#For variance
library(lattice)

#SNPs
bwplot(SERT$exp~SERT$SNP.147, xlab="SNP 147", ylab="Exploration score",
       ylim=c(-5, 155))#Not Homogenous
bwplot(SERT$exp~SERT$SNP.209, xlab="SNP 209", ylab="Exploration score",
       ylim=c(-5, 155))#Not Homogenous
bwplot(SERT$exp~SERT$SNP.446, xlab="SNP 446", ylab="Exploration score",
       ylim=c(-5, 205))#Homogenous
bwplot(SERT$exp~SERT$SNP.467, xlab="SNP 467", ylab="Exploration score",
       ylim=c(-5, 205))#Homogenous

bwplot(SERT$bold~SERT$SNP.147, xlab="SNP 147", ylab="Boldness score",
       ylim=c(-5, 200))#Homogenous-no G
bwplot(SERT$bold~SERT$SNP.209, xlab="SNP 209", ylab="Boldness score",
       ylim=c(-5, 200))#Homogenous-No G
bwplot(SERT$bold~SERT$SNP.446, xlab="SNP 446", ylab="Boldness score",
       ylim=c(-5, 205))#Homogenous
bwplot(SERT$bold~SERT$SNP.467, xlab="SNP 467", ylab="Boldness score",
       ylim=c(-5, 205))#Not Homogenous

#Haplotyopes
bwplot(SERT$exp~SERT$Hap1f, xlab="Hap1", ylab="Exploration score")#Homogenous
bwplot(SERT$exp~SERT$Hap2f, xlab="Hap2", ylab="Exploration score")#Homogenous
bwplot(SERT$exp~SERT$Hap3f, xlab="Haplotype 3", ylab="Exploration score",
       scales = list(y=list(limits=c(-5,155)), tck=c(1,0)))#
bwplot(SERT$exp~SERT$Hap4f, xlab="Hap4", ylab="Exploration score")#Homogenous
bwplot(SERT$exp~SERT$Hap5f, xlab="Hap5", ylab="Exploration score")#Homogenous

#Exploration with age

xyplot(exp~Age, data=SERT)

#Test for homogeneity in the varaince
fligner.test(SERT$exp, SERT$SNP.147)#Not Homogenous
fligner.test(SERT$exp, SERT$Exp.Test.Number)# Not Homogenous
fligner.test(SERT$exp, SERT$Bold.Test.Number)#Homogenous

#Check for Collinearity amongst the predictors (test no and age)
source("http://www.highstat.com/Book2/HighstatLibV6.R")

corvif(SERT[,c(25, 29)])#<4 so keep them in the model


##############Power test##############################
#u=degrees of freedom of the numerator (df of the model p-1)
#v= degrees of freedom in the demoninator (n-p), n being the sample size and p the predictors
#in the model. Continuous counts as one, categorical predictors=(no.levels-1).
#(In our case sex(2-1), status (2-1), age (1),age2(1), test no (1), bird ID+ObserverID=2)-1=6
#Korsten: (SNP(3-1), sex(2-1), brood ID (1))-1=3

library(pwr)

pwr.f2.test(u =6, v =NULL, f2 =0.24, sig.level =0.05, power =0.8)

###################GLMM###############################
library(lme4)

##random intercept model, we treat the SNP as a factor to account for dominance
#or overdominance of alleles. 

#Exploration and SNP
Modintercept<-glmer(exp~SNP.446+Exp.Test.Number+Age+Age2+Status+Sex+(1|BTO)+(1|observ), 
                    data=SERT, family=poisson(link=log), 
                    control=glmerControl(optimizer="bobyqa"))
    
summary(Modintercept)

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  (rdf <- nrow(model@frame)-model.df)
  rp <- residuals(model)
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE,log.p=TRUE)
  c(chisq=Pearson.chisq,ratio=prat,p=exp(pval),logp=pval)
}

overdisp_fun(Modintercept)# Poisson model badly overdipsersed

#To correct for overdispersion, fit individual level random effects
X<-factor(1:length(SERT$exp))
Modisp<-update(Modintercept, .~.+(1|X))
overdisp_fun(Modisp)
summary(Modisp)
confint(Modisp, method="Wald")


#Boldness and SNP
Modintercept2<-glmer(bold~SNP.467+Bold.Test.Number+Status+Sex+Age+Age2+
                       (1|BTO)+(1|observ), 
                     data=SERT, family=poisson(link=log),
                     control=glmerControl(optimizer="bobyqa"))

summary(Modintercept2)

#Boldness# 
overdisp_fun(Modintercept2)# Poisson model badly overdipsersed
Y<-factor(1:length(SERT$bold))
ModispBold<-update(Modintercept2, .~.+(1|Y))
overdisp_fun(ModispBold)
summary(ModispBold)
confint(ModispBold, method="Wald")

#Null Model
#Exploration
ModelNull<-glmer(exp~Exp.Test.Number+Age+Age2+Status+Sex+(1|BTO)+(1|observ), data=SERT, 
                 family=poisson(link=log),
                 control=glmerControl(optimizer="bobyqa"))

ModNullDisp<-update(ModelNull,.~.+(1|X))#Account for overdispersion in the null model

#Boldness
ModelNull2<-glmer(bold~Bold.Test.Number+Age+Age2+Status+Sex+(1|BTO)+(1|observ), 
                  family=poisson(link=log), data=SERT,
                  control=glmerControl(optimizer="bobyqa"))

ModNullDisp2<-update(ModelNull2,.~.+(1|Y))#Account for overdispersion in the null model

#Compare the two models using LRT
anova(Modisp,ModNullDisp)

#Bold
anova(ModispBold,ModNullDisp2)


########random intercept model, we treat the SNP as continuous in an additive model #################

library(lme4)
#Exploration and Allele

ModAllele<-glmer(exp~Allele467+Exp.Test.Number+Age+Age2+Status+Sex+
                   (1|BirdID)+(1|observ), 
                 data=SERT, family=poisson(link=log),
                 control=glmerControl(optimizer="bobyqa"))
summary(ModAllele)


overdisp_fun(ModAllele)# Poisson model badly overdipsersed
X<-factor(1:length(SERT$exp))
ModAlldisp<-update(ModAllele, .~.+(1|X))
overdisp_fun(ModAlldisp)
summary(ModAlldisp)
confint(ModAlldisp, method="Wald")


#Boldness and Allele
ModAllele2<-glmer(bold~Allele467+Bold.Test.Number+Age+Age2+Status+Sex+
                    (1|BirdID)+(1|observ), 
                  data=SERT, family=poisson(link=log),
                  control=glmerControl(optimizer="bobyqa"))
summary(ModAllele2)

overdisp_fun(ModAllele2)# Poisson model badly overdipsersed
Y<-factor(1:length(SERT$bold))
ModAlldisp2<-update(ModAllele2, .~.+(1|Y))
overdisp_fun(ModAlldisp2)
summary(ModAlldisp2)
confint(ModAlldisp2, method="Wald")

#Null Model
#Exploration
ModAlleleNull<-glmer(exp~Exp.Test.Number+Age+Age2+Status+Sex+
                       (1|BirdID)+(1|observ), data=SERT, 
                     family=poisson(link=log),
                     control=glmerControl(optimizer="bobyqa"))

ModNullAllDisp<-update(ModAlleleNull,.~.+(1|X))#Account for overdispersion 

#Boldness
ModAlleleNull2<-glmer(bold~Bold.Test.Number+Age+Age2+Status+Sex+(1|BirdID)+(1|observ), 
                      family=poisson(link=log), data=SERT,
                      control=glmerControl(optimizer="bobyqa"))

ModNullAllDisp2<-update(ModAlleleNull2,.~.+(1|Y))#Account for overdispersion 

#Compare the two models using LRT
anova(ModAlldisp,ModNullAllDisp)

anova(ModAlldisp2,ModNullAllDisp2)


########################Check for correlation with haplotypes##############################################
###################code as integer for addtive and factor for overdominant###################

#Exploration and Haplotype
library(lme4)

ModHap<-glmer(exp~Hap5f+Exp.Test.Number+Age+Age2+Status+Sex+
                (1|BirdID)+(1|observ), 
              data=SERT, family=poisson(link=log),
              control=glmerControl(optimizer="bobyqa"))
summary(ModHap)

overdisp_fun(ModHap)# Poisson model badly overdipsersed
X<-factor(1:length(SERT$exp))
ModAlldisp<-update(ModHap, .~.+(1|X))
overdisp_fun(ModAlldisp)
summary(ModAlldisp)
confint(ModAlldisp, method="Wald")

#Boldness and Hap
ModHap2<-glmer(bold~Hap4f+Bold.Test.Number+Age+Age2+Status+Sex+
                 (1|BirdID)+(1|observ), 
               data=SERT, family=poisson(link=log),
               control=glmerControl(optimizer="bobyqa"))
summary(ModHap2)

overdisp_fun(ModHap2)# Poisson model badly overdipsersed
Y<-factor(1:length(SERT$bold))
ModAlldisp2<-update(ModHap2, .~.+(1|Y))
overdisp_fun(ModAlldisp2)
summary(ModAlldisp2)
confint(ModAlldisp2, method="Wald")

#Null Model
#Exploration
ModHapNull<-glmer(exp~Exp.Test.Number+Age+Age2+Status+Sex
                  +(1|BirdID)+(1|observ), data=SERT, 
                  family=poisson(link=log),
                  control=glmerControl(optimizer="bobyqa"))

ModHapNullDisp<-update(ModHapNull,.~.+(1|X))#Account for overdispersion 

#Boldness
ModHapNull2<-glmer(bold~Bold.Test.Number+Age+Age2+Status+Sex+(1|BirdID)+(1|observ), 
                   family=poisson(link=log), data=SERT,
                   control=glmerControl(optimizer="bobyqa"))

ModHapNullDisp2<-update(ModHapNull2,.~.+(1|Y))#Account for overdispersion 


#Compare the two models using LRT
anova(ModHapNullDisp, ModAlldisp)

anova(ModHapNullDisp2, ModAlldisp2)

##########Adjust for multiple testing (false discovery control)##########################
#############Create csv with pvalues for each model type and adjust with p.adjust########
###e.g. 

p.adjust(PExpSNP$P.Value, method="BH")

###################################HardyWeinberg################

library(HardyWeinberg)

#147
x <- c(60, 2, 23)

names(x) <- c("TT","GG","GT")

HWExact(x, verbose=TRUE)

#209
head(HW2)
y <- c(60, 2, 23)

names(y) <- c("TT","CC","CT")

HWExact(y, verbose=TRUE)


#446
head(HW3)
z <- c(15, 36, 34)

names(z) <- c("TT","AA","AT")

HWExact(z, verbose=TRUE)

#467
head(HW4)
u <- c(35, 7, 43)

names(u) <- c("TT","CC","CT")

HWExact(u, verbose=TRUE)


