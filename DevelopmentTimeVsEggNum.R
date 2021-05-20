#This is the script that we used to examine tradeoffs between larval development time and egg production for Cx. pipiens.
#M. Fritz
#05-14-21

x <- c("gridExtra", "brms", "performance", "MASS", "ggpubr", "rstan", "rstanarm", "reshape") 
lapply(x, FUN = function(X) {do.call("library", list(X))})


#### Set Working Directory ####
#setwd("~/Google Drive/Kate_Bell_Research/Nutrition_Exp./CSVs/")#Kate uses this
setwd("~/Desktop/Culex_Nutrition_Exp")#Megan uses this on Ubuntu 18.04, R v.3.6.3

#### Read in Data ####
data<-read.csv("EarlyvsLateEggNumber.csv",header=T) 


#### Quick Check with cor.test ####

#low
cor.test(data$Day[data$Diet=="low"], data$Number.Eggs[data$Diet=="low"], method=c("pearson")) #no corr

#xtra low
cor.test(data$Day[data$Diet=="x.low"], data$Number.Eggs[data$Diet=="x.low"], method=c("pearson")) #no corr


#prepping data for analysis
data_low <- subset(data, Diet == "low")
data_xlow <- subset(data, Diet == "x.low")


#### Bayesian Models ####

#low first - using both poisson and gaussian
dat.brm.low1 <- brm(Number.Eggs ~ Day, data=data_low, family = poisson(link = "log"), 
                chains=4, iter=8000, warmup=2000, thin=2, sample_prior = TRUE, seed = 12345)

summary(dat.brm.low1)


dat.brm.low2 <- brm(Number.Eggs ~ Day, data=data_low, family = gaussian, 
                    chains=4, iter=8000, warmup=2000, thin=2, sample_prior = TRUE, seed = 12345)

summary(dat.brm.low2)

looic(dat.brm.low1)
looic(dat.brm.low2) #indicates gaussian is better fit

#xlow - just gaussian

dat.brm.xlow1 <- brm(Number.Eggs ~ Day, data=data_xlow, family = gaussian, 
                    chains=4, iter=8000, warmup=2000, thin=2, sample_prior = TRUE, seed = 12345)

summary(dat.brm.xlow1)


#### Plot for Pub ####

pdf("DevelopmentTime_EggNumber_abline.pdf")
par(mfrow=c(2,1))
plot(data$Day[data$Diet=="low"],data$Number.Eggs[data$Diet=="low"],col="Black",pch=16,main="Low Diet",xlab="Day Pupae Picked",ylab="Number of Developed Follicles")
abline(lm(data$Number.Eggs[data$Diet=="low"]~data$Day[data$Diet=="low"]), col="red") # regression line (y~x) 
plot(data$Day[data$Diet=="x.low"],data$Number.Eggs[data$Diet=="x.low"],col="Black",pch=16,main="Extra Low Diet",xlab="Day Pupae Picked",ylab="Number of Developed Follicles")
abline(lm(data$Number.Eggs[data$Diet=="x.low"]~data$Day[data$Diet=="x.low"]), col="red") # regression line (y~x) 
dev.off()

