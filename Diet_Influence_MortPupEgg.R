#This is the script that we used to examine autogenous egg production and larval development
#under different nutritional conditions for Cx. pipiens.
#K. Bell and M. Fritz
#03-14-21

##### Loading Libraries for analysis #########

x <- c("gridExtra", "brms", "performance", "MASS", "ggpubr", "rstan", "rstanarm", "reshape") 
lapply(x, FUN = function(X) {do.call("library", list(X))})


#### Set Working Directory ####
#setwd("~/Google Drive/Kate_Bell_Research/Nutrition_Exp./CSVs/")#Kate uses this
#setwd("C:/Users/megan/Desktop/Culex_Nutrition_Exp")#Megan uses this on Windows 10, R v.4.0.2
setwd("~/Desktop/Culex_Nutrition_Exp")#Megan uses this on Ubuntu 18.04, R v.3.6.3

#### Custom Functions ####
#this function builds a table of summary statistics
SumStats <- function(data, col1, col2){
  SumMeans <- tapply(data[,col1], data[,col2], mean)
  SumSD <- tapply(data[,col1], data[,col2], sd)
  SumLen <- tapply(data[,col1], data[,col2], length)
  SumTab <- rbind(SumMeans, SumSD, SumLen)
  print(SumTab)
}

###### READ IN FOLLICLE DEVELOPMENT DATA ######
dataR1<-read.csv("Egg_Development_Rep1.csv",header=T)#Did rep2, but not using because growth chamber malfunctioned
dataR3<-read.csv("Egg_Development_Rep3.csv",header=T)
dataR4<-read.csv("Egg_Development_Rep4.csv",header=T)
dataR5<-read.csv("Egg_Development_Rep5.csv",header=T)

#fix column differences & problematic header
dataR5<-dataR5[,-c(5)]
names(dataR5)[1] <- "Strain"

#### Subset to just get molestus ####
data_molR1<-subset(dataR1,dataR1$Strain=="mol")
data_molR3<-subset(dataR3,dataR3$Strain=="mol")
data_molR4<-subset(dataR4,dataR4$Strain=="mol")
data_molR5<-subset(dataR5,dataR5$Strain=="mol")

#### Combine R1,3,4,5 ####
data_comb1 <- rbind(data_molR1,data_molR3,data_molR4,data_molR5)
data_comb1$Comb_Eggs_Per_Fem[110:245] <- data_comb1$Num_Dev_Foll[110:245] # Fix differences in columns
rep <- c(rep(1,51),rep(3,58),rep(4,76),rep(5,60)) # Add column with rep information
data_comb1$rep <- as.factor(rep) # Make sure rep is a factor

### Fix Names ###
data_comb1$Treat <- gsub("Hi$", "High", data_comb1$Treat)#$regex for "ends with"
data_comb1$Treat <- gsub("Med$", "Medium", data_comb1$Treat)
data_comb1$Treat <- gsub("X.low", "Xlow", data_comb1$Treat)
data_comb1$Treat <- factor(data_comb1$Treat)
str(data_comb1)#sanity check

keep <- data.frame(data_comb1$Comb_Eggs_Per_Fem,data_comb1$Treat,data_comb1$rep)
colnames(keep) <- c("Eggs","Treat","Rep")
keep$Treat <- factor(keep$Treat, ordered = F, levels = c("High", "Medium", "Low", "Xlow"))
levels(keep$Treat)#sanity check

##### GET SUMMARY DATA FOR FOLLICLE NUMBER #####
SumStats(subset(keep, Rep == 1),1,2)#call to SumStats function above
SumStats(subset(keep, Rep == 3),1,2)
SumStats(subset(keep, Rep == 4),1,2)
SumStats(subset(keep, Rep == 5),1,2)
SumStats(keep,1,2)#This gives overall summary


##### Model Construction Follicle Number #####
kruskal.test(keep$Eggs ~ keep$Treat) #quick look using non-parametric test.

# Bayesian model in BRMS
# Windows version of R/rstan throws the following warning after running these models, but this does not impact model output according to rstan forums
# Warning message: In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :  '-E' not found

fitEggs3 <- brm(Eggs ~ Treat + Rep, data = keep, family ='poisson', chains = 4, iter = 6000, warmup = 2000, thin = 2, seed = 12345)
summary(fitEggs3)
parnames(fitEggs3)
pp_check(fitEggs3, type = "bars", nsamples = 10) ## Make plot of simulated data vs. real data

fitEggs4 <- brm(Eggs ~ Treat + Rep, data = keep, family ='zero_inflated_poisson', chains = 4, iter = 6000, warmup = 2000, thin = 2, seed = 12345)
summary(fitEggs4)
parnames(fitEggs4)
pp_check(fitEggs4, type = "bars", nsamples = 10) ## Make plot of simulated data vs. real data

#model comparison - fitEggs4 better fit with lower looic and higher elpd
looic(fitEggs3)
looic(fitEggs4)

#visualizing posterior analysis for fitEggs4
pdf("Follicles_ppd_1.pdf")
grid.arrange(stan_trace(fitEggs4$fit,pars=c("b_Intercept","b_TreatMedium","b_TreatLow","b_TreatXlow"), ncol=1),stan_dens(fitEggs4$fit,pars=c("b_Intercept","b_TreatMedium","b_TreatLow","b_TreatXlow"), separate_chains=TRUE,ncol=1),ncol=2)
dev.off()

pdf("Follicles_ppd_2.pdf")
grid.arrange(stan_trace(fitEggs4$fit,pars=c("b_Rep3","b_Rep4","b_Rep5"), ncol=1),stan_dens(fitEggs4$fit,pars=c("b_Rep3","b_Rep4","b_Rep5"), separate_chains=TRUE,ncol=1),ncol=2)
dev.off()

#### Eggs as binary trait ####

#Quick check of zeros
zero_keep <- subset(keep, Eggs == 0)
summary(zero_keep)

#recoding to make egg development binary
keep$binary <- ifelse(keep$Eggs == 0, 1, 0)

Rep1 <- with(subset(keep, Rep==1), data.frame(t(table(binary, Treat))))
Rep3 <- with(subset(keep, Rep==3), data.frame(t(table(binary, Treat))))
Rep4 <- with(subset(keep, Rep==4), data.frame(t(table(binary, Treat))))
Rep5 <- with(subset(keep, Rep==5), data.frame(t(table(binary, Treat))))

#since no zeros for rep 1, adding extra dataframe with zeros
Treat <- as.character(c("High", "Medium","Low", "Xlow"))
binary <- c(1,1,1,1)
Freq <- c(0,0,0,0)
zeros_Rep1 <- data.frame(cbind(Treat, binary, Freq))

Tot_binary <- rbind(zeros_Rep1, Rep1, Rep3, Rep4, Rep5)
Tot_binary$Rep <- rep(c(1,3,4,5), each = 8, times = 1)
Tot_binary2 <- reshape(Tot_binary, idvar = c("Treat", "Rep"), timevar = "binary", direction = "wide")
str(Tot_binary2)# sanity check

#making frequencies numeric
Tot_binary2$Freq.1 <- as.numeric(Tot_binary2$Freq.1)#these are actually the zeros
Tot_binary2$Freq.0 <- as.numeric(Tot_binary2$Freq.0)
Tot_binary2$Num_trials <- Tot_binary2$Freq.1 + Tot_binary2$Freq.0 #this gives the number of trials

#here are the models
fitEggs5 <- brm(Freq.1 | trials(Num_trials) ~ Treat, data = Tot_binary2, family ='binomial', chains = 4, iter = 6000, warmup = 2000, thin = 2, seed = 12345)
summary(fitEggs5)
parnames(fitEggs5)
pp_check(fitEggs5, type = "bars", nsamples = 10) 

fitEggs6 <- brm(binary ~ Treat, data = keep, family ='bernoulli', chains = 4, iter = 6000, warmup = 2000, thin = 2, seed = 12345)
summary(fitEggs6)
parnames(fitEggs6)
pp_check(fitEggs6, type = "bars", nsamples = 10) ## Make plot of simulated data vs. real data

#model comparison - models give comparable estimates, but fitEggs5 has better fit with lower looic and higher elpd
looic(fitEggs5)
looic(fitEggs6)


########## READ IN PUPATION DATA ##########
dataR1P<-read.csv("Pupation_Rep1.csv",header=T)#Did rep2, but not using because growth chamber malfunctioned
dataR3P<-read.csv("Pupation_Rep3.csv",header=T)
dataR4P<-read.csv("Pupation_Rep4.csv",header=T)
dataR5P<-read.csv("Pupation_Rep5.csv",header=T)

#### Correct the disparity in days due to different data recorders ####
add<-function(x){x+2}
dataR3P$Days_To_Pupation<-lapply(dataR3P$Days_To_Pupation,add)
dataR4P$Days_To_Pupation<-lapply(dataR4P$Days_To_Pupation,add)
dataR3P$Days_To_Pupation<-as.numeric(dataR3P$Days_To_Pupation)
dataR4P$Days_To_Pupation<-as.numeric(dataR4P$Days_To_Pupation)

#### Set up data frame to run model ####
names(dataR5P)[1] <- "Strain"#fixing colname issue
data_comb<-rbind(dataR1P,dataR3P[,-c(6)],dataR4P,dataR5P)


data_comb$Food_Treatment <- gsub("Hi$", "High", data_comb$Food_Treatment)
data_comb$Food_Treatment <- gsub("Med$", "Medium", data_comb$Food_Treatment)
data_comb$Food_Treatment <- gsub("X.low", "Xlow", data_comb$Food_Treatment)
data_comb$Food_Treatment <- factor(data_comb$Food_Treatment)

data_comb$Food_Treatment_ord <- factor(data_comb$Food_Treatment, ordered = F, 
                                       levels = c("High", "Medium", "Low", "Xlow"))

dat<-data_comb[rep(1:nrow(data_comb), data_comb[["Num_Pupae"]]), ]

dat$Rep<-as.factor(dat$Rep)
str(dat)#sanity check

##### GET SUMMARY DATA FOR TIME TO PUPATION #####
#molestus first
SumStats(subset(dat, Rep == 1 & Strain =="mol"),4,7)#call to SumStats function above
SumStats(subset(dat, Rep == 3 & Strain == "mol"),4,7)
SumStats(subset(dat, Rep == 4 & Strain == "mol"),4,7)
SumStats(subset(dat, Rep == 5 & Strain == "mol"),4,7)
SumStats(subset(dat, Strain == "mol"),4,7)#This gives overall molestus summary

#now pipiens
SumStats(subset(dat, Rep == 1 & Strain =="pip"),4,7)
SumStats(subset(dat, Rep == 3 & Strain == "pip"),4,7)
SumStats(subset(dat, Rep == 4 & Strain == "pip"),4,7)
SumStats(subset(dat, Rep == 5 & Strain == "pip"),4,7)
SumStats(subset(dat, Strain == "pip"),4,7)#This gives overall pipiens summary



##### Model Construction Time To Pupation #####

kruskal.test(dat$Days_To_Pupation ~ dat$Food_Treatment_ord)#quick check with non-parametric test.


#Bayesian Model Comparison

dat.brm1 <- brm(Days_To_Pupation | trunc(lb = 6) ~ Food_Treatment_ord * Strain + Rep, data=dat, family = poisson(link = "log"), 
                   chains=4, iter=6000, warmup=2000, thin=2, sample_prior = TRUE, seed = 12345)

summary(dat.brm1)
pp_check(dat.brm1, type = "bars", nsamples = 10) ## Make plot of simulated data vs. real data

#updating dataset to try and get a better model fit.
subtract<-function(x){x-6} #updating dataset so we don't need the truncated poisson
dat_update<-dat ## Move data to new data frame
dat_update$Days_To_Pupation<-lapply(dat_update$Days_To_Pupation,subtract)
dat_update$Days_To_Pupation<-as.numeric(dat_update$Days_To_Pupation)
summary(dat_update$Days_To_Pupation) ## Make sure days changed as expected

dat.brm2 <-brm(Days_To_Pupation ~ Food_Treatment_ord * Strain + Rep, data=dat_update, family = poisson(link = "log"),
               chains=4, iter=6000, warmup=2000, thin=2, sample_prior = TRUE, seed = 12345) ## Run model with poisson

parnames(dat.brm2)
summary(dat.brm2 )
pp_check(dat.brm2, type = "bars", nsamples = 10) ## Make plot of simulated data vs. real data

#model comparison
looic(dat.brm1)
looic(dat.brm2)

pdf("Pupation_ppd1.pdf")
bayesplot::color_scheme_set("viridis")
plot(dat.brm2)
plot(dat.brm2, regex_pars = c("Days_To_Pupation", "Food_Treatment_ord"))
plot(dat.brm2, plotfun = "combo", regex_pars = "Food_Treatment_ord") 
dev.off()



####### COMBINED BOXPLOT #######
pdf("Combined_BoxPlot.pdf",width=8,height=8)
EN<-ggboxplot(keep, x = "Treat", y = "Eggs",
              fill = "Treat",palette = c("#4d4d4d","#878787","#bababa","#e0e0e0"),
              xlab="Treatment", ylab="Number of Elongated Follicles") + rremove("legend") +theme(plot.title = element_text(hjust = 0.5,face="bold"))

Pup<-ggboxplot(dat, x = "Food_Treatment_ord", y = "Days_To_Pupation",ylim=c(0,25),
               fill = "Strain",palette = c("grey8","white"),
               xlab="Treatment", ylab="Days Until Pupation") +theme(plot.title = element_text(hjust = 0.5,face="bold"),legend.position="right")

ggarrange(Pup, EN,
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()



############## BAYESIAN GLM FOR MORTALITY ##############

#### Read in Data ####
dataR1M<-read.csv("Mortality_Rep1.csv",header=T)#Did rep2, but not using because growth chamber malfunctioned
dataR3M<-read.csv("Mortality_Rep3.csv",header=T)
dataR4M<-read.csv("Mortality_Rep4.csv",header=T)
dataR5M<-read.csv("Mortality_Rep5.csv",header=T)
names(dataR5M)[1] <- "Strain"

#### Correct the disparity in days due to different data recorders####
add<-function(x){x+2}

dataR3M$Days_to_Mort<-lapply(dataR3M$Days_to_Mort,add)
dataR4M$Days_to_Mort<-lapply(dataR4M$Days_to_Mort,add)
dataR3M$Days_to_Mort<-as.numeric(dataR3M$Days_to_Mort)
dataR4M$Days_to_Mort<-as.numeric(dataR4M$Days_to_Mort)

str(dataR3M)#sanity check
str(dataR4M)


################ Prep Data for Analysis ##################
rep1<-dataR1M[dataR1M$Days_to_Mort==(max(dataR1M$Days_to_Mort)),c(1:4,6)]
rep3<-dataR3M[dataR3M$Days_to_Mort==(max(dataR3M$Days_to_Mort)),c(1:4,6)]
rep4<-dataR4M[dataR4M$Days_to_Mort==(max(dataR4M$Days_to_Mort)),c(1:4,6)]
rep5<-dataR5M[dataR5M$Days_to_Mort==(max(dataR5M$Days_to_Mort)),c(1:4,6)]

rep1p<-dataR1P[dataR1P$Days_To_Pupation==(max(dataR1P$Days_To_Pupation)),c(1:4,6)]#these from pupation datasets loaded above
rep3p<-dataR3P[dataR3P$Days_To_Pupation==(max(dataR3P$Days_To_Pupation)),c(1:4,7)]
rep4p<-dataR4P[dataR4P$Days_To_Pupation==(max(dataR4P$Days_To_Pupation)),c(1:4,6)]
rep5p<-dataR5P[dataR5P$Days_To_Pupation==(max(dataR5P$Days_To_Pupation)),c(1:4,6)]

datR1<-data.frame(rep1$Strain,rep1$Treatment,rep1$Rep,stringsAsFactors = FALSE)
datR1$Alive<-rep1p$Cum_Pupae
datR1$Dead<-rep1$Cumulative_dead

datR3<-data.frame(rep3$Strain,rep3$Treatment,rep3$Rep,stringsAsFactors = FALSE)
datR3$Alive<-rep3p$Cum_Pup  
datR3$Dead<-rep3$Cumulative_dead

datR4<-data.frame(rep4$Strain,rep4$Treatment,rep4$Rep,stringsAsFactors = FALSE)
datR4$Alive<-rep4p$Cum_Pupae
datR4$Dead<-rep4$Cumulative_dead

datR5<-data.frame(rep5$Strain,rep5$Treatment,rep5$Rep,stringsAsFactors = FALSE)
datR5$Alive<-rep5p$Cum_Pupae
datR5$Dead<-rep5$Cumulative_dead

colnames(datR1)<-c("Strain","Treatment","Rep","Alive","Dead")
colnames(datR3)<-c("Strain","Treatment","Rep", "Alive","Dead")
colnames(datR4)<-c("Strain","Treatment","Rep","Alive","Dead")
colnames(datR5)<-c("Strain","Treatment","Rep","Alive","Dead")
data_comb2<-rbind(datR1,datR3,datR4,datR5)

### Fix Names ###
data_comb2$Treatment <- gsub("Hi$", "High", data_comb2$Treatment)
data_comb2$Treatment <- gsub("Med$", "Medium", data_comb2$Treatment)
data_comb2$Treatment <- gsub("X.low", "Xlow", data_comb2$Treatment)
data_comb2$Treatment <- factor(data_comb2$Treatment)

data_comb2$Treatment <- factor(data_comb2$Treatment, ordered = F, 
                        levels = c("High", "Medium", "Low", "Xlow"))

data_comb2$Rep = factor(data_comb2$Rep)
str(data_comb2)#sanity check


##### Binomial Model Construction #####
Mort_out <- brm(Alive | trials (Alive + Dead) ~ Treatment + Rep + Strain + Treatment*Strain, data = data_comb2,
         family ="binomial", chains=4, iter=6000, warmup=2000, thin=2, seed = 12345)

summary(Mort_out)
parnames(Mort_out)
pp_check(Mort_out, type = "bars", nsamples = 10)


pdf("Mortality1_TP.pdf")
grid.arrange(stan_trace(Mort_out$fit,pars=c("b_Intercept","b_TreatmentMedium","b_TreatmentLow","b_TreatmentXlow"), ncol=1),
             stan_dens(Mort_out$fit,pars=c("b_Intercept","b_TreatmentMedium","b_TreatmentLow","b_TreatmentXlow"), separate_chains=TRUE, ncol=1), ncol=2)
dev.off()

pdf("Mortality2_TP.pdf")
grid.arrange(stan_trace(Mort_out$fit,pars=c("b_Rep3","b_Rep4","b_Rep5","b_Strainpip"), ncol=1),
             stan_dens(Mort_out$fit,pars=c("b_Rep3","b_Rep4","b_Rep5","b_Strainpip"), separate_chains=TRUE, ncol=1), ncol=2)
dev.off()

pdf("Mortality3_TP.pdf")
grid.arrange(stan_trace(Mort_out$fit,pars=c("b_TreatmentMedium:Strainpip","b_TreatmentLow:Strainpip","b_TreatmentXlow:Strainpip"), ncol=1),
             stan_dens(Mort_out$fit,pars=c("b_TreatmentMedium:Strainpip","b_TreatmentLow:Strainpip","b_TreatmentXlow:Strainpip"), separate_chains=TRUE, ncol=1),ncol=2)
dev.off()


###### Sex-specific influence Diet ######

#reading in data
dataSS <-read.csv("Pip_Mol_pupTime_bySex.csv", header=T)
str(dataSS)

#expanding dataset to get one row per observation
expanded <- untable(dataSS[,c(1:5)], num=dataSS[,6])
identical(nrow(expanded),sum(dataSS$num_adults))#sanity check
tail(dataSS)
tail(expanded, n = 10)

str(expanded)
expanded$start_date <- as.factor(expanded$start_date)

#making sex-by-form for plotting
expanded$sexByform <- paste0(expanded$sex,"_",expanded$form)
levels(expanded$diet)
levels(expanded$diet) <- c("High", "Extra-Low")


###### PLOTS - Sex*Diet*Form ######

#Exploratory
ggplot(expanded) +
  aes(x = sexByform, y = pup_time_days, color = factor(sexByform))+
  scale_color_manual(values=c("red", "blue", "darkred","darkblue")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "grey30", fill=NA, size=1.5)) +
  geom_violin(lwd = 1.2, trim=FALSE) + #geom_jitter(shape=21, fill = NA, colour = c("grey60"), position=position_jitter(0.2)) +
  facet_wrap(~diet)

#For Pub
png("Pupation_Time_SexByDietByForm.png", unit = "px", height = 500, width = 1000)
ggplot(expanded) +
  aes(x = sexByform, y = pup_time_days) + 
  xlab("Sex and Form") + ylab("Pupation Time in Days") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "grey30", fill=NA, size=1.5), strip.text = element_text(size=15),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x=element_text(angle=45, hjust=1, size = 11, face = "bold"),
        axis.text.y=element_text(size = 11, face = "bold")) +
  scale_x_discrete("",labels=c("female_mol"="Female Mol", "female_pip"="Female Pip", "male_mol"="Male Mol", "male_pip"="Male Pip")) +
  geom_boxplot(lwd=1.1, fatten = T, fill = "grey80") + geom_jitter(shape=21, fill = NA, colour = c("grey30"), position=position_jitter(0.2)) +
  facet_wrap(~diet)
dev.off()

#getting pip sumstats by sex & diet
SumStats(subset(dataSS, form == "pip" & sex == "male" & diet == "xlo"),6,4)
SumStats(subset(dataSS, form == "pip" & sex == "female" & diet == "xlo"),6,4)
SumStats(subset(dataSS, form == "pip" & sex == "male" & diet == "hi"),6,4)
SumStats(subset(dataSS, form == "pip" & sex == "female" & diet == "hi"),6,4)

#getting mol sumstats by sex & diet
SumStats(subset(dataSS, form == "mol" & sex == "male" & diet == "xlo"),6,4)
SumStats(subset(dataSS, form == "mol" & sex == "female" & diet == "xlo"),6,4)
SumStats(subset(dataSS, form == "mol" & sex == "male" & diet == "hi"),6,4)
SumStats(subset(dataSS, form == "mol" & sex == "female" & diet == "hi"),6,4)



#Bayesian Model Comparison

exp.brm1 <- brm(pup_time_days | trunc(lb = 6) ~ diet * sex * form + start_date, data=expanded, family = poisson(link = "log"), 
                chains=4, iter=6000, warmup=2000, thin=2, sample_prior = TRUE, seed = 12345)

summary(exp.brm1)
pp_check(exp.brm1, type = "bars", nsamples = 10) ## Make plot of simulated data vs. real data

#updating dataset to try and get a better model fit.
subtract<-function(x){x-6} #updating dataset so we don't need the truncated poisson
exp_update<-expanded ## Move data to new data frame
exp_update$pup_time_days<-lapply(exp_update$pup_time_days,subtract)
exp_update$pup_time_days<-as.numeric(exp_update$pup_time_days)
summary(exp_update$pup_time_days) ## Make sure days changed as expected

exp.brm2 <-brm(pup_time_days ~ diet * sex * form + start_date, data=exp_update, family = poisson(link = "log"),
               chains=4, iter=6000, warmup=2000, thin=2, sample_prior = TRUE, seed = 12345) ## Run model with poisson

parnames(exp.brm2)
summary(exp.brm2 )
pp_check(exp.brm2, type = "bars", nsamples = 10) ## Make plot of simulated data vs. real data

#model comparison
looic(exp.brm1)
looic(exp.brm2)
