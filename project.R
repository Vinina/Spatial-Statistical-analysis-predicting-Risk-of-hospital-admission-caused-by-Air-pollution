#load packages
library(sp)
library(ggplot2)
library(patchwork)
require(sf)
library(rgeos)
library(rgdal)
library(spdep)
library(nimble)
library(coda)
library(sf)
library(R2WinBUGS)
library(dplyr)
library(tidyverse)

#loading london data
load("london.RData")
#the data loads as london_3, checking the data first
head(london_3@data)

#creating a subset with the year 2007
london_3@data <- (subset(london_3@data, year%in%2007))
#converting to sf for ggplot
london_2007_sf <- st_as_sf(london_3)
#london_2007_sf <- subset(london_sf, year%in%2007)
head(london_2007_sf)
#calculating SIR
london_2007_sf$SIR <- london_2007_sf$Observed/london_2007_sf$Expected
#using ggplot 
ggplot(london_2007_sf) +
  geom_sf(aes(fill=SIR, colour=SIR)) +
  scale_fill_viridis_c(option="C") +
  scale_colour_viridis_c(option="C") +
  theme_void()
## Check distribution of SIR
hist(london_2007_sf$SIR, main = "Respiratory disease SIR across London", xlab = "Standardised incidence ratio")
#trying log of SIR
london_2007_sf$logSIR <- log(london_2007_sf$SIR)

ggplot(london_2007_sf) +
  geom_sf(aes(fill=logSIR, colour=logSIR)) + 
  scale_fill_viridis_c(option="C") +
  scale_colour_viridis_c(option="C")+
  theme_void() 
#checking scatter plots to understand effect of other variables on SIR
#plot(london_2007_sf$SIR~london_2007_sf$PM25, main="Scatterplot Example",
#     xlab="PM25 ", ylab="SIR", pch=19)

#plot(london_2007_sf$JSA, london_2007_sf$SIR, main="Scatterplot Example",
#     xlab="JSA ", ylab="SIR", pch=19)

#plot(london_2007_sf$Price, london_2007_sf$SIR, main="Scatterplot Example",
#     xlab="Price ", ylab="SIR", pch=19)

# Basic Scatterplot Matrix
#pairs(~SIR+PM25+JSA+Price,data=london_2007_sf,
#      main="Simple Scatterplot Matrix")

## Plot SIR against covariates, scatter plots:

ggplot(london_2007_sf) +
  geom_point(aes(PM25, SIR)) ->
  p1
ggplot(london_2007_sf) +
  geom_point(aes(JSA, SIR)) ->
  p2
ggplot(london_2007_sf) +
  geom_point(aes(Price, SIR)) ->
  p3

layout <-"
AABB
#CC#
"

p1+p2+p3 + plot_layout(design=layout)

#Question 2
## Fit the model #bugs

Data <- list(Y = london_2007_sf$Observed, E = london_2007_sf$Expected,
             N = nrow(london_2007_sf), JSA = london_2007_sf$JSA, PM25 = london_2007_sf$PM25, Price = london_2007_sf$Price)
Inits <- list(list(beta0 = 0, beta1 = 0, beta2 = 0, beta3 = 0))

preg.bugs <- bugs(data = Data,
                  model.file = "poisson-regression.txt",
                  parameters.to.save = c("beta0", "beta1", "beta2", "beta3"),
                  inits = Inits,
                  n.iter = 5000,
                  n.burnin = 2500,
                  n.thin = 1,
                  n.chains = 1,
                  bugs.directory="D:/Vinina Office/FinalLORs/StudentDocs/collegeInfo/MM915_stats/Assignment/winbugs14_full_patched/WinBUGS14")

sims <- as.data.frame(preg.bugs$sims.list)

ggplot(sims)+
  geom_line(aes(1:nrow(sims), beta0))+
  labs(x="Index")->
  p1

ggplot(sims)+
  geom_line(aes(1:nrow(sims), beta1))+
  labs(x="Index") ->
  p2

ggplot(sims)+
  geom_line(aes(1:nrow(sims), beta2))+
  labs(x="Index") ->
  p3

ggplot(sims)+
  geom_line(aes(1:nrow(sims), beta3))+
  labs(x="Index") ->
  p4

p1/p2

p3/p4

## Fit the model # nimble

model_code_poisson <- nimbleCode({
  for(i in 1:N){
    Y[i] ~dpois(mu[i])
    log(mu[i]) <- log(E[i]) + beta0 + beta1*JSA[i]+ beta2*PM25[i]+ beta3*Price[i]
  }
  
  beta0~dnorm(0,0.01)
  beta1~dnorm(0,0.01)
  beta2~dnorm(0,0.01)
  beta3~dnorm(0,0.01)
}) 

Data_poisson <- list(Y = london_2007_sf$Observed, 
             E = london_2007_sf$Expected, 
             JSA = london_2007_sf$JSA, 
             PM25 = london_2007_sf$PM25, 
             Price = london_2007_sf$Price)
Data_poisson <- list(Y = london_2007_sf$Observed, 
             E = london_2007_sf$Expected, 
             JSA = london_2007_sf$JSA - mean(london_2007_sf$JSA), 
             PM25 = london_2007_sf$PM25 - mean(london_2007_sf$PM25), 
             Price = london_2007_sf$Price - mean(london_2007_sf$Price))
Constants_poisson <- list(N = nrow(london_2007_sf))

Inits_poisson <- list(list(beta0 = 0, beta1 = 0, beta2 = 0, beta3 = 0),
              list(beta0 = 0, beta1 = 0, beta2 = 0, beta3 = 0))

preg_mod <- nimbleMCMC(data = Data_poisson,
                       constants=Constants_poisson,
                       code=model_code_poisson,
                       monitors = c("beta0", "beta1", "beta2", "beta3"),
                       inits = Inits_poisson,
                       niter = 5000,
                       nburnin = 2500,
                       thin = 1,
                       nchains = 2,
                       summary=TRUE,
                       samplesAsCodaMCMC=TRUE)

#with mu to calculate DIC
preg_mod_mean <- nimbleMCMC(data = Data_poisson,
                       constants=Constants_poisson,
                       code=model_code_poisson,
                       monitors = c("beta0", "beta1", "beta2", "beta3","mu"),
                       inits = Inits_poisson,
                       niter = 5000,
                       nburnin = 2500,
                       thin = 1,
                       nchains = 2,
                       summary=TRUE,
                       samplesAsCodaMCMC=TRUE)

#sims <- as.list(preg_mod$samples)
#sims <- Reduce("rbind", samples)
#sims <- as.data.frame(samples)
sims <- as.data.frame(Reduce("rbind", preg_mod$samples))
sims_mean <- as.data.frame(Reduce("rbind", preg_mod_mean$samples))


#sims <- as.data.frame(preg_mod$samples)
ggplot(sims)+
  geom_line(aes(1:nrow(sims), beta0))+
  labs(x="Iteration")->
  p1

ggplot(sims)+
  geom_line(aes(1:nrow(sims), beta1))+
  labs(x="Iteration") ->
  p2

ggplot(sims)+
  geom_line(aes(1:nrow(sims), beta2))+
  labs(x="Iteration") ->
  p3

ggplot(sims)+
  geom_line(aes(1:nrow(sims), beta3))+
  labs(x="Iteration") ->
  p4

p1+p2+p3+p4


# Gelman-Rubin diagnostic (Rhat) 
gelman.diag(preg_mod$samples)
# use coda package to get geweke values #Checking convergence
for_geweke <- coda::geweke.plot(preg_mod$samples) 

# This is a data frame with a row for each point on the geweke plot and a column for each parameter 
for_geweke <- as.data.frame(for_geweke) 

for_geweke %>% 
  # Transform to each row being for a particular point and a particular parameter (so we can use  facet_wrap) 
  pivot_longer(-start.iter, names_to=c("Variable","Chain"), values_to="z", names_prefix="z\\.", names_sep="\\.(?=c)")%>% 
  # Make the variable names nicer (Some beautiful regex there!) 
  mutate(Variable = str_replace(Variable, "\\.1$",""), 
         Variable=str_replace(Variable, "\\.", "["), 
         Variable=str_replace(Variable, "\\.", "]")) %>% 
  ggplot() + 
  # Plot the points 
  geom_point(aes(start.iter, z, color=Chain)) + 
  # Add dashed lines at +/-2 
  geom_hline(aes(yintercept=-2),linetype="dashed") + 
  geom_hline(aes(yintercept=2),linetype="dashed") + 
  # Separate plots for each parameter 
  facet_wrap(~Variable)

#check assumption-  Pearson residuals
# extract the beta samples 
# grepl() looks for the presence of "beta" in the data frame names
beta_samples <- apply(sims[,grepl("beta", names(sims))], 2, median)

#log_rate <- mean(sims$beta0) +  
#  mean(sims$beta1) * london_2007_sf$JSA+  
#  mean(sims$beta2) * london_2007_sf$PM25+  
#  mean(sims$beta3) * london_2007_sf$Price
# three covariates in this model
lp <- beta_samples[1] + beta_samples[2]*london_2007_sf$JSA + 
  beta_samples[3]*london_2007_sf$PM25 + beta_samples[4]*london_2007_sf$Price

mu <- london_2007_sf$Expected*exp(lp)
#mu <- london_2007_sf$Expected*exp(log_rate)

PREG_res <- (london_2007_sf$Observed - mu)/sqrt(mu)
#I want to check the mean=variance assumption first:

var(PREG_res)
ggplot() +
  geom_point(aes(mu, PREG_res))

#independence assumption
moran.test(PREG_res, nb2listw(poly2nb(london_3)), alternative="two.sided")


###Question 3 

# Produce the adjacency matrix, and various pieces of associated info BUGS needs
W <- nb2mat(poly2nb(london_3), style = "B")
inds <- lapply(1:nrow(W), function(i) which(W[i, ] == 1))
Adj <- Reduce("c", inds)
Num.Adj <- rowSums(W)
SumNumNeigh <- sum(Num.Adj)

# combine all of the data, and constants into a single list
pCAR_code <- nimbleCode({
  for(i in 1:N){
    Y[i] ~ dpois(mu[i])
    log(mu[i]) <- log(E[i]) + beta0 + beta1*JSA[i]+ beta2*PM25[i]+ beta3*Price[i] + phi[i]  
  }
  
  phi[1:N] ~ dcar_normal(Adj[1:L], weights[1:L], Num[1:N], tau, zero_mean=1)
  
  beta0~dnorm(0,0.01)
  beta1~dnorm(0,0.01)
  beta2~dnorm(0,0.01)
  beta3~dnorm(0,0.01)
  tau ~ dgamma(0.01,0.01)
})

Data_car <- list(Y = london_2007_sf$Observed, 
             E = london_2007_sf$Expected, 
             JSA = london_2007_sf$JSA - mean(london_2007_sf$JSA), 
             PM25 = london_2007_sf$PM25 - mean(london_2007_sf$PM25), 
             Price = london_2007_sf$Price - mean(london_2007_sf$Price))

Constants_car <- list(N = nrow(london_2007_sf), 
                  Adj = Adj, 
                  Num = Num.Adj,
                  L = SumNumNeigh, 
                  weights=rep(1, SumNumNeigh))
# Initial Values - it's a good idea to set initial values for this type of model
Inits_car <- list(list(beta0=0, beta1=0, beta2=0, beta3=0, tau=1, phi=rep(1, nrow(london_2007_sf))),
              list(beta0=0, beta1=0, beta2=0, beta3=0, tau=1, phi=rep(1, nrow(london_2007_sf))),
              list(beta0=0, beta1=0, beta2=0, beta3=0, tau=1, phi=rep(1, nrow(london_2007_sf))))

pCAR <- nimbleMCMC(data=Data_car,
                   constants = Constants_car,
                   code=pCAR_code,
                   monitors=c(paste0("beta", 0:3), "phi", "mu"),
                   inits=Inits_car,
                   nchains = 3,
                   niter=10000,
                   nburnin=5000,
                   summary=TRUE,
                   samplesAsCodaMCMC = TRUE)
car.sims <- as.data.frame(Reduce("rbind",pCAR$samples))
ggplot(car.sims) +
  geom_line(aes(1:nrow(car.sims), beta0)) +
  labs(x="Iteration") ->
  p1
ggplot(car.sims) +
  geom_line(aes(1:nrow(car.sims), beta1)) +
  labs(x="Iteration") ->
  p2
ggplot(car.sims) +
  geom_line(aes(1:nrow(car.sims), beta2)) +
  labs(x="Iteration") ->
  p3
ggplot(pCar_samples) +
  geom_line(aes(1:nrow(pCar_samples), beta3)) +
  labs(x="Iteration") ->
  p4
ggplot(car.sims) +
  geom_line(aes(1:nrow(car.sims), `phi[102]`)) +
  labs(x="Iteration") ->
  p5
ggplot(car.sims) +
  geom_line(aes(1:nrow(car.sims), `phi[230]`)) +
  labs(x="Iteration") ->
  p6
p1+p2+p3+p4+p5+p6
# Gelman-Rubin diagnostic (Rhat) 
for_diag <- lapply(pCAR$samples, function(x) x[,!grepl("mu", names(x[1,]))])
coda::gelman.diag(for_diag)


Neff <- function(samples, t=20){
  N <- nrow(samples)
  sum_cor <- apply(samples, 2,function(x){
    sum(sapply(1:t, function(y){
      abs(cor(x[-(1:y)], x[-((N-y+1):N)]))
    }))
  })
  return(N/(1+2*sum_cor))
}
Neff(as.data.frame(Reduce("rbind", pCAR$samples)))

# Run the model
pCAR <- nimbleMCMC(data=Data_car,
                   constants = Constants_car,
                   code=pCAR_code,
                   monitors=c(paste0("beta", 0:3), "phi", "mu"),
                   inits=Inits_car,
                   nchains = 3,
                   niter=170000,
                   nburnin=85000,
                   thin=17,
                   summary=TRUE,
                   samplesAsCodaMCMC = TRUE)

coda::gelman.diag(pCAR$samples, multivariate=FALSE)
pCar_samples <- as.data.frame(Reduce("rbind", pCAR$samples))
#lets check the trace plots
ggplot(pCar_samples) +
  geom_line(aes(1:nrow(pCar_samples), beta0)) +
  labs(x="Iteration") ->
  p1
ggplot(pCar_samples) +
  geom_line(aes(1:nrow(pCar_samples), beta1)) +
  labs(x="Iteration") ->
  p2
ggplot(pCar_samples) +
  geom_line(aes(1:nrow(pCar_samples), beta2)) +
  labs(x="Iteration") ->
  p3
ggplot(pCar_samples) +
  geom_line(aes(1:nrow(pCar_samples), beta3)) +
  labs(x="Iteration") ->
  p4
ggplot(pCar_samples) +
  geom_line(aes(1:nrow(pCar_samples), `phi[102]`)) +
  labs(x="Iteration") ->
  p5
ggplot(pCar_samples) +
  geom_line(aes(1:nrow(pCar_samples), `phi[230]`)) +
  labs(x="Iteration") ->
  p6
p1+p2+p3+p4+p5+p6

#let's check geweke plot for convergence
# use coda package to get geweke values #Checking convergence
for_diag <- lapply(pCAR$samples, function(x) x[,!grepl("mu", names(x[1,]))])
to_plot <- c("beta0","beta1","beta2","beta3", paste0("phi[", sample(1:nrow(london_2007_sf), 12, replace=FALSE), "]"))

for_geweke <- coda::geweke.plot(for_diag) 
# This is a data frame with a row for each point on the geweke plot and a column for each parameter 
for_geweke <- as.data.frame(for_geweke)
for_geweke %>% 
  # Transform to each row being for a particular point and a particular parameter (so we can use  facet_wrap) 
  pivot_longer(-start.iter, names_to=c("Variable","Chain"), values_to="z", names_prefix="z\\.", names_sep="\\.(?=c)")%>% 
  # Make the variable names nicer (Some beautiful regex there!) 
  mutate(Variable = str_replace(Variable, "\\.1$",""), 
         Variable=str_replace(Variable, "\\.", "["), 
         Variable=str_replace(Variable, "\\.", "]")) %>% 
  filter(Variable%in%to_plot) %>%
  ggplot() + 
  # Plot the points 
  geom_point(aes(start.iter, z, color=Chain)) + 
  # Add dashed lines at +/-2 
  geom_hline(aes(yintercept=-2),linetype="dashed") + 
  geom_hline(aes(yintercept=2),linetype="dashed") + 
  # Separate plots for each parameter 
  facet_wrap(~Variable)


#check assumption-  Pearson residuals
#I want to check the mean=variance assumption first:

mu_pCAR <- pCAR$summary$all.chains[grepl("mu", rownames(pCAR$summary$all.chains)),1]
res_pCAR <- (london_2007_sf$Observed - mu_pCAR)/sqrt(mu_pCAR)
var(res_pCAR)
ggplot() +
  geom_point(aes(mu_pCAR, res_pCAR))

#independence assumption
moran.test(res_pCAR, nb2listw(poly2nb(london_3)), alternative="two.sided")


##Question 4

DIC <- function(Observed, Fitted, name="mu"){
  Fitted <- Fitted[,grepl(name, names(Fitted))]
  deviance <- sum(log(dpois(Observed, lambda=apply(Fitted,2,median))))
  
  ave_dev <- mean(sapply(1:nrow(Fitted), function(x) sum(log(dpois(Observed, lambda=as.numeric(Fitted[x,]))))))
  
  pD <- 2*(deviance - ave_dev)
  DIC <- -2*deviance + 2*pD
  list(pD = pD, DIC = DIC)
}
preg_DIC <- DIC(london_2007_sf$Observed, sims_mean[grepl("mu", names(sims_mean))])
pCAR_DIC <- DIC(london_2007_sf$Observed, pCar_samples[grepl("mu", names(pCar_samples))])

exp(preg_mod$summary$all.chains[grepl("beta",rownames(preg_mod$summary$all.chains)),c(1,4,5)])
exp(pCAR$summary$all.chains[grepl("beta",rownames(pCAR$summary$all.chains)),c(1,4,5)])

