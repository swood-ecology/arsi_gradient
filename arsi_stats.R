## --------------------
#### Arsi analysis
## --------------------

## Load packages
library(tidyverse)  # For reading in data
library(arm)        # For standardize() function
library(rstan)      # For interfacing with Stan


## Read in data and generate new variables
arsi_full_data <- read_csv("~/Box Sync/Work/Writing/Manuscripts/Unsubmitted/arsi_negele/arsi_full_data.csv")
names(arsi_full_data)[4:24] <- c('Protein','Fe.crop','Zn.crop','Area','FoodYield','WheatYield','N.appl.amt','N.appl.p.ha','SIR','POM.C','MAOM.C','POM.N','MAOM.N','Sand','Silt','Clay','pH','Cu','Fe','Mn','Zn')

arsi_full_data$Protein <- arsi_full_data$Protein/100
arsi_full_data$MAOM.C.N <- arsi_full_data$MAOM.C/arsi_full_data$MAOM.N
arsi_full_data$POM.C.N <- arsi_full_data$POM.C/arsi_full_data$POM.N
arsi_full_data$Fines <- arsi_full_data$Silt + arsi_full_data$Clay
arsi_full_data$Total.N <- arsi_full_data$POM.N + arsi_full_data$MAOM.N


## Assess correlations
# Create function to determine correlations to plot on upper panels
panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  if(missing(cex.cor)) cex.cor <- 1/(strwidth(txt))
  cex.final = cex.cor * r
  if(cex.final < .5) cex.final <- .6
  text(0.5, 0.5, txt, cex = cex.final)
}
corr.vars <- c("WheatYield","N.appl.amt","POM.C","MAOM.C","POM.N","MAOM.N","Sand",
               "Silt","Clay","Cu","Fe","Mn","Zn","Dmun","Trd",
               "MAOM.C.N","POM.C.N","pH")
pairs(arsi_full_data[,corr.vars], lower.panel = panel.smooth, upper.panel = panel.cor)


## Yield model
lm.yield <- lm(log(WheatYield)~MAOM.C+POM.N+N.appl.amt+Trd,data=arsi_full_data)
summary(standardize(lm.yield))
# car::qqPlot(lm(log(WheatYield)~MAOM.C+POM.N+Trd+N.appl.amt,data=arsi_full_data))
# car::ncvTest(lm(log(WheatYield)~MAOM.C+POM.N+Trd+N.appl.amt,data=arsi_full_data))
# car::residualPlots(lm(log(WheatYield)~MAOM.C+POM.N+Trd+N.appl.amt,data=arsi_full_data))
# car::crPlots(lm(log(WheatYield)~MAOM.C+POM.N+Trd+N.appl.amt,data=arsi_full_data))
# car::avPlots(lm(log(WheatYield)~MAOM.C+POM.N+Trd+N.appl.amt,data=arsi_full_data))

# Drop NA values
yield.vars <- c("WheatYield","N.appl.amt","MAOM.C","POM.N","Trd")
yield.data <- arsi_full_data[complete.cases(arsi_full_data[,yield.vars]),c(yield.vars)]
yield.data <- as.data.frame(yield.data)

# Fit Stan model
yield.list <- list(N=nrow(yield.data),
                    K=4,
                    fert=yield.data$N.appl.amt,
                    maom=yield.data$MAOM.C,
                    pom=yield.data$POM.N,
                    trd=yield.data$Trd,
                    y=yield.data$WheatYield
                  )
yield.model <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/yield.stan", data = yield.list, control = list(adapt_delta=0.99,max_treedepth=15), chains = 4)

# Results
print(yield.model, probs=c(.025,.975),
      pars=c('beta_fert','beta_maom','beta_pom','beta_trd')) # Natural coefficients
plot(yield.model,pars=c('beta_std','sigma_std'),outer_level=0.975)
#traceplot(yield.model, inc_warmup = F,pars=c('beta','sigma','alpha'))

# Posterior predictive check
y_pred <- rstan::extract(yield.model,pars='y_tilde')
y_pred <- unlist(y_pred, use.names=FALSE)

yield.pp.data <- data.frame(c(y_pred,yield.list$y),c(rep("y_pred",length(y_pred)),rep("y_obs",length(yield.list$y))))
names(yield.pp.data) <- c("y","type")

ggplot(yield.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("Wheat yield") + ylab("Density") +
  scale_fill_manual(values=wesanderson::wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )


## Nutrient models
# Initial models
lm.zn <- lm(log(Zn.crop)~ N.appl.amt + pH + MAOM.N + POM.N + Silt + Trd,data=arsi_full_data)
summary(standardize(lm.zn))
lm.fe <- lm(log(Fe.crop)~ N.appl.amt + pH + MAOM.N + POM.N + Silt + Trd,data=arsi_full_data)
summary(standardize(lm.fe))

# Stan models
nutr.vars <- c("Zn.crop","Fe.crop","N.appl.amt","pH","Silt","MAOM.N","POM.N","Trd")
nutr.data <- arsi_full_data[complete.cases(arsi_full_data[,nutr.vars]),c(nutr.vars)]
nutr.data <- as.data.frame(nutr.data)

zn.list <- list(N=nrow(nutr.data),
                K=6,
                x=nutr.data[,c('N.appl.amt','pH','Silt','MAOM.N','POM.N','Trd')],
                y=nutr.data$Zn.crop
)

zn.model <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/crop_nutrients.stan", data = zn.list, control = list(adapt_delta=0.99,max_treedepth=17), chains=4)
print(zn.model, probs=c(.025,.975), pars=c('beta'))
plot(zn.model, pars=c("beta_std","sigma_std"),outer_level=0.975)
#traceplot(zn.model, inc_warmup = F, pars=c('alpha','beta','sigma'))

# Posterior predictive check
y_pred <- rstan::extract(zn.model,pars='y_tilde')
y_pred <- unlist(y_pred, use.names=FALSE)

zn.pp.data <- data.frame(c(y_pred,zn.list$y),c(rep("y_pred",length(y_pred)),rep("y_obs",length(zn.list$y))))
names(zn.pp.data) <- c("y","type")

ggplot(zn.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("Crop zinc concentration") + ylab("Density") +
  scale_fill_manual(values=wesanderson::wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )

fe.list <- list(N=nrow(nutr.data),
                K=6,
                x=nutr.data[,c('N.appl.amt','pH','Silt','MAOM.N','POM.N','Trd')],
                y=nutr.data$Fe.crop
)

fe.model <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/crop_nutrients.stan", data = fe.list, control = list(adapt_delta=0.999,max_treedepth=20), chains=4)

print(fe.model, probs=c(.025,.975), pars=c("beta"))
plot(fe.model, pars=c("beta_std","sigma_std"),outer_level=0.975)

# Posterior predictive check
y_pred <- rstan::extract(fe.model,pars='y_tilde')
y_pred <- unlist(y_pred, use.names=FALSE)

fe.pp.data <- data.frame(c(y_pred,fe.list$y),c(rep("y_pred",length(y_pred)),rep("y_obs",length(fe.list$y))))
names(fe.pp.data) <- c("y","type")

ggplot(fe.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("Crop iron concentration") + ylab("Density") +
  scale_fill_manual(values=wesanderson::wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )


## SOM models

# Set up data for Stan
# Define dummies
arsi_full_data <- cbind(arsi_full_data, dummies::dummy(arsi_full_data$Type)) %>%
  as.tibble()
names(arsi_full_data)[39:41] <- c('typeForest','typeHome','typeWheat')

# Compare models with different hierarchies
lm(log(POM.C) ~ rescale(typeHome) + rescale(typeWheat) + rescale(Dmun) + rescale(Fines) + rescale(pH),arsi_full_data) %>%
  summary()
lme4::lmer(log(POM.C) ~ rescale(Dmun) + rescale(Fines) + rescale(pH) + (1|Type),arsi_full_data) %>%
  summary()
lme4::lmer(log(POM.C) ~ rescale(Dmun) + rescale(Fines) + rescale(pH) + (1 + rescale(Dmun)|Type),arsi_full_data) %>%
  summary()

# Remove NAs
som.vars <- c('Location','POM.C','MAOM.C','SIR','typeHome','typeWheat','Dmun','Fines','pH')
som.data <- arsi_full_data[complete.cases(arsi_full_data[,som.vars]),c(som.vars)]

# Define data
pomList <- list(
  N = nrow(som.data),
  K = ncol(som.data[,c('typeHome','typeWheat','Dmun','Fines','pH')]),
  y = som.data$POM.C,
  x = som.data[,c('typeHome','typeWheat','Dmun','Fines','pH')]
)
maomList <- list(
  N = nrow(som.data),
  K = ncol(som.data[,c('typeHome','typeWheat','POM.C','Dmun','Fines','pH')]),
  y = som.data$MAOM.C,
  x = som.data[,c('typeHome','typeWheat','POM.C','Dmun','Fines','pH')]
)
sirList <- list(
  N = nrow(som.data),
  K = ncol(som.data[,c('typeHome','typeWheat','POM.C','MAOM.C','Dmun','Fines','pH')]),
  y = som.data$SIR,
  x = som.data[,c('typeHome','typeWheat','POM.C','MAOM.C','Dmun','Fines','pH')]
)

# Execute models
POM <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/som.stan",
     data = pomList,
     iter = 2000, chains = 4)
MAOM <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/som.stan",
            data = maomList,
            iter = 2000, chains = 4)
SIR <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/som.stan",
            data = sirList,
            iter = 2000, chains = 4)

# Print model
print(POM,pars=c('beta'))
plot(POM,pars=c('beta_std'),outer_level=0.975)
print(MAOM,pars=c('beta'))
plot(MAOM,pars=c('beta_std'),outer_level=0.975)
print(SIR,pars=c('beta'))
plot(SIR,pars=c('beta_std'),outer_level=0.975)


# Posterior predictive checks
# POM
y_pred_pom <- rstan::extract(POM,pars='y_tilde')
y_pred_pom <- unlist(y_pred_pom, use.names=FALSE)

pom.pp.data <- data.frame(c(y_pred_pom,pomList$y),c(rep("y_pred",length(y_pred_pom)),rep("y_obs",length(pomList$y))))
names(pom.pp.data) <- c("y","type")

ggplot(pom.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("POM") + ylab("Density") +
  scale_fill_manual(values=wesanderson::wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )

# MAOM
y_pred_maom <- rstan::extract(MAOM,pars='y_tilde')
y_pred_maom <- unlist(y_pred_maom, use.names=FALSE)

maom.pp.data <- data.frame(c(y_pred_maom,maomList$y),c(rep("y_pred",length(y_pred_maom)),rep("y_obs",length(maomList$y))))
names(maom.pp.data) <- c("y","type")

ggplot(maom.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("MAOM") + ylab("Density") +
  scale_fill_manual(values=wesanderson::wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )

# SIR
y_pred_sir <- rstan::extract(SIR,pars='y_tilde')
y_pred_sir <- unlist(y_pred_sir, use.names=FALSE)

sir.pp.data <- data.frame(c(y_pred_sir,sirList$y),c(rep("y_pred",length(y_pred_sir)),rep("y_obs",length(sirList$y))))
names(sir.pp.data) <- c("y","type")

ggplot(sir.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("Substrate-induced Respiration") + ylab("Density") +
  scale_fill_manual(values=wesanderson::wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )
