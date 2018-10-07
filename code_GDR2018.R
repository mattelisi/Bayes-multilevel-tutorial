# This script contains the code for the examples presented in the method
# course on multilevel modelling at the GDR Vision 2018, held in Paris
# 
# Matteo Lisi, 2018

# clear workspace
rm(list=ls())
setwd("~/sync/multilevel_tutorial_GDR2018/")

# load useful libraries
library(lme4)
library(ggplot2)
library(mlisi)
nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))

###################################################
# sleep study example (frequentist approach)
###################################################

str(sleepstudy)
ggplot(sleepstudy, aes(x=Days,y=Reaction))+geom_point()+geom_smooth(method="lm",color="black",lwd=1)+facet_wrap(~Subject,ncol=9)+nice_theme
m.1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
summary(m.1)

# boostrap confidence intervals
#CI_fixef <- confint(m.1, method="boot", nsim=500, oldNames=F)
#save(CI_fixef, file="sleep1.RData") # save result to external file for loading into slides
load("sleepCI.RData")
print(CI_fixef, digits=2)

# analyze shrinking
# use lmList function to compute individual fits
m.single <- coef(lmList(Reaction ~ Days | Subject, sleepstudy))
par.mixed <- as.matrix(ranef(m.1)$Subject) + repmat(t(as.matrix(fixef(m.1))),18,1)

# plot parameters from mixed model and individual fit
plot(m.single[,1], m.single[,2], xlab="Intercept",ylab="Slope",pch=19,cex.lab=1.2,col="dark grey",xlim=c(190,310),ylim=c(-5,25))

# draw ellipse illustrating covariance of random effects
vcov_m.1 <- matrix(as.vector(VarCorr(m.1)$Subject),ncol=2)
mixtools::ellipse(c(mean(par.mixed[,1]), mean(par.mixed[,2])), sigma=vcov_m.1,alpha=0.05, col="grey", lty=2)

points(mean(m.single[,1]), mean(m.single[,2]),pch=19,col="dark grey",cex=2)
points(mean(par.mixed[,1]), mean(par.mixed[,2]),pch=21,col="black",cex=2,lwd=2)
text(m.single[,1], m.single[,2],labels=rownames(m.single),pos=1,cex=0.8)
points(par.mixed[,1], par.mixed[,2])
arrows(m.single[,1], m.single[,2],par.mixed[,1], par.mixed[,2], length=0.1)
legend("bottomright",c("single-subject fits","multilevel model"),pch=c(19,21),col=c("dark grey", "black"),bty="n")


###################################################
# Metropolis algorithm example
###################################################

str(sleepstudy)
with(sleepstudy[sleepstudy$Subject=="308",],plot( Days,Reaction,pch=19,main="308"))

x <- sleepstudy$Days[sleepstudy$Subject=="308"]
y <- sleepstudy$Reaction[sleepstudy$Subject=="308"]

loglik <- function(par){
  pred <- par[1] + par[2]*x
  return(sum(dnorm(y, mean = pred, sd = exp(par[3]), log = T)))
}

# priors
par(mfrow=c(1,3))
curve(dnorm(x, mean=250, sd=180),from=0, to=1000, xlab="Intercept",ylab="prior density", col="blue")
curve(dnorm(x, mean=20, sd=20),from=-50, to=50, xlab="Slope",ylab="prior density", col="blue")
x_pSD <- seq(-1,6, length.out = 500)
y_pSD <- dnorm(x_pSD , mean=4,sd=1)
plot(exp(x_pSD),y_pSD, type="l", xlab=expression(sigma[epsilon]),ylab="prior density", col="blue")

logprior <- function(par){
  intercept_prior <- dnorm(par[1], mean=250, sd=180, log=T)
  slope_prior <- dnorm(par[2], mean=20, sd=20, log=T)
  sd_prior <- dnorm(par[3],mean=4, sd=1, log=T)
  return(intercept_prior+slope_prior+sd_prior)
}

logposterior <- function(par){
  return (loglik(par) + logprior(par))
}

# So, the aim of this algorithm is to jump around in parameter space, but in a way that the probability to be at a point is proportional to the function we sample from (this is usually called the target function). In our case this is the posterior defined above.

# metropolis
run_metropolis_MCMC <- function(startvalue, iterations){
  chain <- array(dim = c(iterations+1,3))
  chain[1,] <- startvalue
  
  for (i in 1:iterations){
    
    # draw a random proposal
    proposal <- proposalfunction(chain[i,])
    
    # ratio of posterior density between new and old values
    a <- exp(logposterior(proposal) - logposterior(chain[i,]))
    
    # accept/reject
    if (runif(1) < a){
      chain[i+1,] <- proposal
    }else{
      chain[i+1,] <- chain[i,]
    }
  }
  return(chain)
}

proposalfunction <- function(par){
  return(rnorm(3, mean = par, sd= c(15,5,0.2)))
}
startvalue <- c(250,20,5)

set.seed(1)
chain <- run_metropolis_MCMC(startvalue, 20000)
# #save(chain, file="chain_ex.RData") # save for slides
burnIn <- 5000
acceptance <- 1-mean(duplicated(chain[-(1:burnIn),]))

LSfit <- lm(y~x)
interceptLS <- coef(LSfit)[1]
slopeLS <- coef(LSfit)[2]
sigmaLS <- summary(LSfit)$sigma

# plot distributions and chain values
par(mfrow = c(2,3))
hist(chain[-(1:burnIn),1],main="Intercept",border="white",col="dark grey", breaks=20)
abline(v = interceptLS, col="red",lwd=2)
hist(chain[-(1:burnIn),2],main="Slope",border="white",col="dark grey", breaks=20)
abline(v = slopeLS , col="red",lwd=2)
hist(exp(chain[-(1:burnIn),3]),main=expression(sigma[epsilon]),border="white",col="dark grey", breaks=20)
abline(v = sigmaLS, col="red" ,lwd=2)

plot(chain[-(1:burnIn),1], type = "l", main = "Chain values")
abline(h = interceptLS, col="red",lwd=2)
plot(chain[-(1:burnIn),2], type = "l" , main = "Chain values")
abline(h = slopeLS, col="red",lwd=2)
plot(exp(chain[-(1:burnIn),3]), type = "l" , main = "Chain values")
abline(h = sigmaLS, col="red",lwd=2)

# remove initial 'burn in' samples
burnIn <- 5000
slope_samples <- chain[-(1:burnIn),2]
alpha <- 0.05

# mean of posterior distribution
mean(slope_samples)

# 95% Bayesian credible interval
round(c(quantile(slope_samples, probs = alpha/2),quantile(slope_samples, probs = 1-alpha/2)),digits=2)


###################################################
# Sleepstudy example using stan (Bayesian approach)
###################################################

# format data as list
d_stan <- list(Subject=as.numeric(factor(sleepstudy$Subject, labels=1:length(unique(sleepstudy$Subject)))), Days=sleepstudy$Days, RT = sleepstudy$Reaction/1000, N=nrow(sleepstudy), J=length(unique(sleepstudy$Subject)) )

library(rstan)
options(mc.cores = parallel::detectCores()) # indicate stan to use multiple cores if available
sleep_model <- stan(file = "sleep_model_v2.stan", data = d_stan, iter = 500, chains = 4, control = list(max_treedepth = 18))

# saveRDS(sleep_model,file="sleep_fit_ppc.rds")
# sleep_model <- readRDS("sleep_fit.rds")

# examine sampling results
print(sleep_model, pars = c("beta","sigma_e","sigma_u"), probs = c(0.025, 0.975),digits=3)

# examine posterior distribution
# chain convergence
sleep_model@model_pars
traceplot(sleep_model, pars = c("beta","sigma_e","sigma_u"), inc_warmup = F)
pairs(sleep_model, pars = c("beta","sigma_e","sigma_u"))

## Use the L matrices (Choleski factors) to compute the correlation matrices
L_u <- extract(sleep_model, pars = "L_u")$L_u
cor_u <- apply(L_u, 1, function(x) tcrossprod(x)[1, 2])
print(signif(quantile(cor_u, probs = c(0.025, 0.5, 0.975))))

# compute posterior predictive distribution
# this is the implied distribution of observations according to the model
# weighted by the posterior probability of parameters
# I use the simulated data to compute 95% confidence intervals on the 
# predicted values of observations for new subjects (random samples from 
# the population distribution)
nd <- {}
y_rep <- extract(sleep_model, pars = "y_rep")$y_rep
for(i in 1:nrow(y_rep)){
  nd <- rbind(nd, data.frame(Days=sleepstudy$Days, Subject=sleepstudy$Subject,  Reaction=y_rep[i,]))
}
ndag <- aggregate(Reaction ~ Days, nd, mean)
alpha <- 0.05
ndag$ci_lb <- aggregate(Reaction ~ Days, nd, function(x) quantile(x, probs = alpha/2))$Reaction
ndag$ci_ub <- aggregate(Reaction ~ Days, nd, function(x) quantile(x, probs = 1-alpha/2))$Reaction
# saveRDS(ndag,file="sleep_ppc_CI.rds") # again, that's for slides

plot(sleepstudy$Days, sleepstudy$Reaction/1000, xlab="Days",ylab="Reaction",ylim=c(0.17,0.5),type="n")
polygon(c(ndag$Days,rev(ndag$Days)), c(ndag$ci_ub,rev(ndag$ci_lb)), col="skyblue",border=NA)
points(jitter(sleepstudy$Days,factor=0.4), sleepstudy$Reaction/1000, pch=19,col=rgb(0,0,0,0.5),cex=1.2)

#########
# you can fit the same model using RStanArm
library(rstanarm)

# define some priors
prInter = normal(location=0.3, scale=1, autoscale=F)
prCoeff <- normal(location=c(0.01), scale =  c(0.1), autoscale=F)
prSigma <- cauchy(location=c(0), scale =  c(0.5), autoscale=F)
sleepstudy$Reaction_s <- sleepstudy$Reaction/1000

# model is defined with lmer-style syntax but using MCMC sampling
sleep_model_2 <- stan_lmer(Reaction_s ~ Days + (Days|Subject), data = sleepstudy, cores = 2, iter = 2000, chain=2, prior = prCoeff, prior_intercept = prInter, prior_covariance=prVarCov)

fixef(sleep_model_2)
ranef(sleep_model_2)

#################################################################################
### mcmc diagnostics: examples on how sampling can go wrong
#################################################################################

dat <- list(y=c(-1,1))
stan_code <- '
data{
real y[2];
}
parameters{
real mu;
real<lower=0> sigma;
}
model{
y ~ normal(mu, sigma);
}
'
m_ <- stan(model_code=stan_code, data = dat, iter = 4000, chains = 2)
traceplot(m_, pars = c("mu","sigma"), inc_warmup = F)
saveRDS(m_,file="wildchain1.rds")

# with weakly informative priors
dat <- list(y=c(-1,1))
stan_code <- "
data{
real y[2];
}
parameters{
real mu;
real<lower=0> sigma;
}
model{
mu ~ normal(1, 10);
sigma ~ cauchy(0, 1);
y ~ normal(mu, sigma);
}
"
m_ <- stan(model_code=stan_code, data = dat, iter = 4000, chains = 2)
traceplot(m_, pars = c("mu","sigma"), inc_warmup = F)
saveRDS(m_,file="wildchain2.rds")

# unidentifiable parameters
dat <- list(y=rnorm(100, mean=0, sd=1))
stan_code <- "
data{
real y[100];
}
parameters{
vector[2] alpha;
real<lower=0> sigma;
}
model{
real mu;
mu = alpha[1] + alpha[2];
y ~ normal(mu, sigma);
}
"
m_ <- stan(model_code=stan_code, data = dat, iter = 4000, chains = 2)
traceplot(m_, pars = c("alpha"), inc_warmup = F)
saveRDS(m_,file="unidentified1.rds")

dat <- list(y=rnorm(100, mean=0, sd=1))
stan_code <- "
data{
real y[100];
}
parameters{
vector[2] alpha;
real<lower=0> sigma;
}
model{
real mu;
alpha ~ normal(0, 10);
mu = alpha[1] + alpha[2];
y ~ normal(mu, sigma);
}
"
m_ <- stan(model_code=stan_code, data = dat, iter = 4000, chains = 2)
traceplot(m_, pars = c("alpha"), inc_warmup = F)
saveRDS(m_,file="unidentified2.rds")

###################################################
# GLMM
# fitting psychometric functions at the group level
# taking lapse probability into account
###################################################

# load dataset
# data are part of an upcoming study by Lisi, Solomon & Morgan
d <- read.table("bisection2.txt", sep="\t",header=T)

# bin data for plotting
# (stimulus intensities in the experiment were determined adaptively 
# using the Quest+ algorithm, see Watson, 2017, JVis for details)
d$dx_bin <- cut(d$dx, 4)
dag <- aggregate(cbind(rr,dx)~id+dx_bin+cond, d, mean)
dag$se <- aggregate(rr~id+dx_bin+cond, d, binomSEM)$rr
dag$n <- aggregate(rr~id+dx_bin+cond, d, length)$rr

ggplot(dag, aes(x=dx,y=rr,color=factor(cond)))+geom_errorbar(aes(ymin=rr-se,ymax=rr+se),width=0)+geom_point(pch=21)+facet_wrap(~id,ncol=7)+nice_theme+geom_smooth(data=d,aes(x=dx,y=rr),method="glm",method.args=list(family=binomial(probit)),size=0.8,fullrange=T,se=F)+scale_color_manual(values=c("black","blue"),guide=F)+labs(y="choice probability", x="distance offset [deg]")

# visualize priors and hyperpriors
# prior for lapse rate
x <- seq(0,1, length.out=500)
beta_x <- dbeta(x,1,5)
plot(x, beta_x, type="l",col="blue",ylab=expression(p(lambda)),xlab=expression(lambda),bty="n")
text(0.6,3,labels=expression(paste(lambda %~% "Beta(a=1,b=4)")),cex=1.5)

# hyperprior for beta distribution
x <- seq(0,10, length.out=1000)
shape_g1 <- dgamma(x, shape =4, rate=1)
plot(x, shape_g1, type="l", ylim=c(0,0.25),ylab="prior probability density", xlab="",col="blue",bty="n")
shape_g2 <- dgamma(x, shape =1, rate=0.2)
lines(x, shape_g2, col="red")
legend("topright",c("p(a)","p(b)"),lwd=1,col=c("red","blue"),bty="n",inset=0.2)

# prepare data for Stan fit
d_stan <- list(id=as.numeric(factor(d$id, labels=1:length(unique(d$id)))), ds=d$dx, rr = d$rr, N=nrow(d), J=length(unique(d$id)), cond=d$cond-1, ds_ppc=rep(seq(-4,4,length.out=20),2),N_ppc=40,cond_ppc=c(rep(0,20),rep(1,20)))
str(d_stan)

bisection_model <- stan(file = "bisection_stan_v1.stan", data = d_stan, iter = 4000, chains = 4)
#saveRDS(bisection_model,file="bisection_fit_ok.rds")
#bisection_model <- readRDS("bisection_fit_ppc.rds")

# inspect sampling outcomes
traceplot(bisection_model, pars = c("beta","sigma_u","beta_ab"), inc_warmup = F)
print(bisection_model, pars = c("beta","sigma_u","beta_ab"), probs = c(0.025, 0.975),digits=3)
print(bisection_model, pars = c("lambda"), probs = c(0.025, 0.975),digits=3)

# inference
beta_par <- extract(bisection_model, pars = "beta", inc_warmup = FALSE)$beta
jnd_1 <- 1/beta_par[,3]
jnd_2 <- 1/(beta_par[,3]+beta_par[,4])
jnd_diff <- jnd_2-jnd_1

# Highest Posterior Density intervals
# note the functions for the HPD intervals 
# are in the 'coda' package
jnd_ci_1 <- coda::HPDinterval(coda::as.mcmc(jnd_1), prob=0.95)
jnd_ci_2 <- coda::HPDinterval(coda::as.mcmc(jnd_2), prob=0.95)
jnd_ci_diff <- coda::HPDinterval(coda::as.mcmc(jnd_diff), prob=0.95)

library(bayesplot)
posterior2 <- extract(bisection_model, inc_warmup = F, permuted = FALSE)

mcmc_intervals(posterior2, pars = c("beta[1]","beta[2]","beta[3]","beta[4]","sigma_u[1]","sigma_u[2]","sigma_u[3]","sigma_u[4]"), prob=0.8, prob_outer = 0.95, point_est="mean")+labs(title="Posterior distributions",subtitle="with mean, 80% and 95% credible intervals")

mcmc_areas(posterior2, pars = c("beta[1]","beta[2]","beta[3]","beta[4]","sigma_u[1]","sigma_u[2]","sigma_u[3]","sigma_u[4]"), prob=0.80, point_est="mean")+labs(title="Posterior distributions",subtitle="with mean and 80% credible intervals")

traceplot(bisection_model, pars = c("beta","sigma_u","beta_ab"), inc_warmup = F)
pairs(bisection_model, pars = c("beta"))


#################################################################################
### hierarchical mixture model
#################################################################################

# load and prepare data
d <- read.table("cueingdd.txt",header=T,sep="\t")
d$alpha_rad <- d$alpha_deg/180*pi
d$rerr <- d$dir*angDiff(d$alpha_rad, d$resp)# transform in signed angular differences

# plot distributions
d$rerr_deg <- d$rerr/pi*180
d$cue_f <- ifelse(d$cue==1,"(1) pre-cue","(2) post-cue")
ggplot(d,aes(x=rerr_deg))+geom_histogram(binwidth=5,aes(y = ..density..),color="white",lwd=0.2)+geom_density(size=0.3,lty=1)+geom_vline(xintercept=0,lty=2,size=0.5)+facet_grid(id~cue_f)+nice_theme+labs(x=" response error [deg]")+scale_x_continuous(breaks=seq(-180,180,45))

# visualize hyperprior for lambda distribution
x <- seq(0,10, length.out=1000)
shape_g1 <- dgamma(x, shape =9, rate=1)
plot(x, shape_g1, type="l", ylim=c(0,0.25),ylab="prior probability density", xlab="",col="blue",bty="n")
shape_g2 <- dgamma(x, shape =1, rate=1)
lines(x, shape_g2, col="red")
legend("topright",c("p(a)","p(b)"),lwd=1,col=c("red","blue"),bty="n",inset=0.2)

# prepare data for stan
cue_stan <- list(id=as.numeric(d$id), J=length(unique(d$id)), rerr=d$rerr, cue=d$cue, N=nrow(d))
str(cue_stan)

library(rstan)
options(mc.cores = parallel::detectCores()) # indicate stan to use multiple cores if available
mix_model <- stan(file = "mixmodel_2.stan", data = cue_stan, iter = 1000, chains = 4)
# saveRDS(mix_model,file="mix_model_2.3.rds")
# mix_model <- readRDS("mixmodel_2.3.rds")

# check that sampling went fine
traceplot(mix_model, pars = c("mu","logk","logilambda"), inc_warmup = T)
traceplot(mix_model, pars = c("mu_j"), inc_warmup = F)
pairs(mix_model, pars = c("mu","logk","logilambda"), inc_warmup = F)
print(mix_model, pars = c("mu","logk","logilambda"), probs = c(0.025, 0.975),digits=3)

# compute 95% credible interval on bias differences
mu1 <- extract(mix_model, pars="mu[1]", inc_warmup = F, permuted = FALSE)
mu2 <- extract(mix_model, pars="mu[2]", inc_warmup = F, permuted = FALSE)
mu1 <- c(mu1[,1,1],mu1[,2,1])/pi*180
mu2 <- c(mu2[,1,1],mu2[,2,1])/pi*180
coda::HPDinterval(as.mcmc(mu2-mu1), prob=0.95)

# function to transform concentration parameter k into a standard deviation
k2SD <- function(k) sqrt(-2 *log(besselI(k, nu=1)/besselI(k, nu=0)))

# compute 95% credible interval on precision differences (quantified as standard deviations)
k1 <- extract(mix_model, pars="logk[1]", inc_warmup = F, permuted = FALSE)
k2 <- extract(mix_model, pars="logk[2]", inc_warmup = F, permuted = FALSE)
sd1 <- k2SD(exp(c(k1[,1,1],k1[,2,1])))/pi*180
sd2 <- k2SD(exp(c(k2[,1,1],k2[,2,1])))/pi*180
HPDinterval(as.mcmc((sd2-sd1), prob=0.95))
