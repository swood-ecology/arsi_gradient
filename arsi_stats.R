## --------------------
#### Arsi analysis
## --------------------

## Load packages
library(tidyverse)  # For reading in data
library(arm)        # For standardize() function
library(rstan)      # For interfacing with Stan


## Relevant functions
# Determine correlations to plot on upper panels
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

## Read in data and generate new variables
arsi_full_data <- read_csv("~/Box Sync/Work/Writing/Manuscripts/Unsubmitted/arsi_negele/arsi_full_data.csv")
names(arsi_full_data)[4:24] <- c('Protein','Fe.crop','Zn.crop','Area','FoodYield','WheatYield','N.appl.amt','N.appl.p.ha','SIR','POM.C','MAOM.C','POM.N','MAOM.N','Sand','Silt','Clay','pH','Cu','Fe','Mn','Zn')

arsi_full_data$Protein <- arsi_full_data$Protein/100
arsi_full_data$MAOM.C.N <- arsi_full_data$MAOM.C/arsi_full_data$MAOM.N
arsi_full_data$POM.C.N <- arsi_full_data$POM.C/arsi_full_data$POM.N
arsi_full_data$Fines <- arsi_full_data$Silt + arsi_full_data$Clay
# Define dummies
arsi_full_data <- cbind(arsi_full_data, dummies::dummy(arsi_full_data$Type)) %>%
  as.tibble()
names(arsi_full_data)[39:41] <- c('typeForest','typeHome','typeWheat')


## Assess correlations
corr.vars <- c("WheatYield","N.appl.amt","POM.C","MAOM.C","POM.N","MAOM.N","Sand",
               "Silt","Clay","Cu","Fe","Mn","Zn","Dmun","Trd",
               "MAOM.C.N","POM.C.N","pH")
pairs(arsi_full_data[,corr.vars], lower.panel = panel.smooth, upper.panel = panel.cor)


# 1. Do yield and nutrients increase or decrease with distance to road?
# # 1.1. Generate vector of non-text location values
# loc.name <- as.vector(yield.data$Location)
# uniq <- unique(loc.name)
# J <- length(uniq)
# loc <- rep (NA,J)
# for (i in 1:J){
#   loc[loc.name==uniq[i]] <- i
# }

# 1.2. Define data
yield.vars <- c("WheatYield","Fe.crop","Zn.crop","Trd","Location","N.appl.amt","pH","POM.N","MAOM.C","MAOM.N")
yield.data <- arsi_full_data[complete.cases(arsi_full_data[,yield.vars]),c(yield.vars)]
yield.data <- as.data.frame(yield.data)

yield.list <- list(N=nrow(yield.data),
                   K=5,
                   fert=yield.data$N.appl.amt,
                   maom=yield.data$MAOM.C,
                   pom=yield.data$POM.N,
                   pH=yield.data$pH,
                   trd=yield.data$Trd,
                   y=yield.data$WheatYield
)
zn.list <- list(N=nrow(yield.data),
                K=5,
                fert=yield.data$N.appl.amt,
                maom=yield.data$MAOM.N,
                pom=yield.data$POM.N,
                pH=yield.data$pH,
                trd=yield.data$Trd,
                y=yield.data$Zn.crop)
fe.list <- list(N=nrow(yield.data),
                K=5,
                fert=yield.data$N.appl.amt,
                maom=yield.data$MAOM.N,
                pom=yield.data$POM.N,
                pH=yield.data$pH,
                trd=yield.data$Trd,
                y=yield.data$Fe.crop)

# 1.3. Call Stan models
yield.model <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/final_models/yield_and_nutrient.stan", 
                    data = yield.list, 
                    control = list(adapt_delta=0.99,max_treedepth=15), chains = 4)
zn.model <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/final_models/yield_and_nutrients.stan", 
                 data = zn.list, 
                 control = list(adapt_delta=0.99,max_treedepth=15), chains = 4)
fe.model <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/final_models/yield_and_nutrients.stan", 
                 data = fe.list, 
                 control = list(adapt_delta=0.99,max_treedepth=15), chains = 4)

# 1.4. Model results
print(yield.model,pars=c("beta_fert","beta_maom","beta_pom","beta_pH","beta_trd"),
      probs=c(0.05,0.95))
plot(yield.model,pars=c("beta_std"))

print(zn.model,pars=c("beta_fert","beta_maom","beta_pom","beta_pH","beta_trd"),
      probs=c(0.05,0.95))
plot(zn.model,pars=c("beta_std"))

print(fe.model,pars=c("beta_fert","beta_maom","beta_pom","beta_pH","beta_trd"),
      probs=c(0.05,0.95))
plot(fe.model,pars=c("beta_std"))

# 1.5. Posterior predictive checks
yield_pred <- rstan::extract(yield.model,pars='y_tilde')
yield_pred <- unlist(yield_pred, use.names=FALSE)
yield.pp.data <- data.frame(c(yield_pred,yield.list$y),c(rep("yield_pred",length(yield_pred)),rep("y_obs",length(yield.list$y))))
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


# 2. Does SOM vary along the same landscape gradient as yield? In other words, could it be a driver?
# 2.1. Set up data for Stan
som.vars <- c('Location','POM.C','MAOM.C','SIR','typeHome','typeWheat','Dmun','Fines','pH')
som.data <- arsi_full_data[complete.cases(arsi_full_data[,som.vars]),c(som.vars)]

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

# 2.2. Execute models
POM <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/final_models/som.stan",
            data = pomList,
            iter = 2000, chains = 4)
MAOM <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/final_models/som.stan",
             data = maomList,
             iter = 2000, chains = 4)
SIR <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/final_models/som.stan",
            data = sirList,
            iter = 2000, chains = 4)

# 2.3. Print model results
print(POM,pars=c('beta'),probs=c(0.05,0.95))
plot(POM,pars=c('beta_std'))
print(MAOM,pars=c('beta'),probs=c(0.05,0.95))
plot(MAOM,pars=c('beta_std'))
print(SIR,pars=c('beta'),probs=c(0.05,0.95))
plot(SIR,pars=c('beta_std'))


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

# OM models without other pools
maomList <- list(
  N = nrow(som.data),
  K = ncol(som.data[,c('typeHome','typeWheat','Dmun','Fines','pH')]),
  y = som.data$MAOM.C,
  x = som.data[,c('typeHome','typeWheat','Dmun','Fines','pH')]
)
sirList <- list(
  N = nrow(som.data),
  K = ncol(som.data[,c('typeHome','typeWheat','Dmun','Fines','pH')]),
  y = som.data$SIR,
  x = som.data[,c('typeHome','typeWheat','Dmun','Fines','pH')]
)

MAOM <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/final_models/som.stan",
             data = maomList,
             iter = 2000, chains = 4)
SIR <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/final_models/som.stan",
            data = sirList,
            iter = 2000, chains = 4)

print(MAOM,pars=c('beta'),probs=c(0.05,0.95))
print(SIR,pars=c('beta'),probs=c(0.05,0.95))
