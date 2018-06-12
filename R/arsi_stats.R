#####################################################################################
# Title:        Arsi soil and crop analysis                                         #
# Description:  Crop and soil data along a distance-to-forest gradient in Ethiopia  #
# Author:       Stephen Wood                                                        #
# Last updated: 5/31/18                                                             #
#####################################################################################

# LOAD PACKAGES
library(tidyverse)      # For reading in data
library(arm)            # For standardize() function
library(rstan)          # For interfacing with Stan
library(wesanderson)    # Color palette for plots
library(soiltexture)    # For plotting soil texture
source("R/panel_cor.R") # For plotting pairs plots with correlations

# READ DATA
arsi <- read_csv("data/arsi_full_data.csv")
names(arsi)[4:24] <- c('Protein','Fe.crop','Zn.crop','Area','FoodYield','WheatYield','N.appl.amt','N.appl.p.ha','SIR','POM.C','MAOM.C','POM.N','MAOM.N','Sand','Silt','Clay','pH','Cu','Fe','Mn','Zn')

# GENERATE NEW VARIABLES
arsi$pH <- 10^arsi$pH
arsi$Protein <- arsi$Protein/100
arsi$MAOM.C.N <- arsi$MAOM.C/arsi$MAOM.N
arsi$POM.C.N <- arsi$POM.C/arsi$POM.N
arsi$Fines <- arsi$Silt + arsi$Clay
arsi$Texture <- arsi$Fines/arsi$Sand
arsi$MAOM.C.perc <- arsi$MAOM.C*10
arsi$MAOM.N.perc <- arsi$MAOM.N*10
arsi$POM.N.perc <- arsi$POM.N*10

## Define dummies
arsi <- cbind(arsi, dummies::dummy(arsi$Type), dummies::dummy(arsi$Location)) %>% 
            as.tibble()
names(arsi)[43:49] <- c('typeForest','typeHome','typeWheat','locationForest',
                        'locationMiddle','locationNearForest','locationRoad')

# ASSESS CORRELATIONS
corr.vars <- c("WheatYield","N.appl.p.ha","POM.C","MAOM.C","POM.N","MAOM.N",
               "Texture","Cu","Fe","Mn","Zn","Dmun","MAOM.C.N","POM.C.N","pH")
pairs(arsi[,corr.vars], lower.panel = panel.smooth, upper.panel = panel_cor)
pairs(arsi[,c('WheatYield','Protein','Fe.crop','Zn.crop','N.appl.p.ha','POM.N','MAOM.C','MAOM.N')], lower.panel = panel.smooth, upper.panel = panel_cor)

# PLOT TEXTURE DATA
text.dat <- arsi[,c('Clay', 'Silt','Sand','SIR')]
names(text.dat) <- c("CLAY","SILT","SAND","SIR")
TT.plot(
  class.sys = "USDA.TT",
  tri.data = as.data.frame(text.dat),
  z.name="SIR",
  main = "Arsi Negele soil texture"
)
rm(text.dat)

# QUESTION 1
# Do yield and nutrients increase or decrease with distance to road?

## Define data
yield.vars <- c("WheatYield","Fe.crop","Zn.crop","Protein","N.appl.amt","N.appl.p.ha",
                "pH","Fe","Zn","POM.N","POM.C","MAOM.C","MAOM.N","MAOM.C.N","POM.C.N",
                "MAOM.C.perc","MAOM.N.perc","POM.N.perc","locationMiddle","locationNearForest",
                "locationRoad","typeWheat","typeHome")
yield.data <- arsi[complete.cases(arsi[,yield.vars]),c(yield.vars)]
yield.data <- as.data.frame(yield.data)

## List data for Stan model
### Use MAOM and POM N for yield and protein 
### because of expectation that aggrgate yield and amino acids depend on N
yield.list <- list(
                    N=nrow(yield.data),
                    K=ncol(yield.data[,c('locationNearForest','locationMiddle','N.appl.p.ha','MAOM.N','POM.N','pH')]),
                    y=yield.data$WheatYield,
                    x=yield.data[,c('locationNearForest','locationMiddle','N.appl.p.ha','MAOM.N','POM.N','pH')]
)
pro.list <- list(
  N=nrow(yield.data),
  K=ncol(yield.data[,c('locationNearForest','locationMiddle','N.appl.p.ha','MAOM.N','POM.N','pH')]),
  y=yield.data$Protein*100,
  x=yield.data[,c('locationNearForest','locationMiddle','N.appl.p.ha','MAOM.N','POM.N','pH')]
)
### Use MAOM C for Zn and Fe models because of expectation 
### that CEC is important for micronutrients
zn.list <- list(
                N=nrow(yield.data),
                K=ncol(yield.data[,c('locationNearForest','locationMiddle','N.appl.p.ha','MAOM.C','POM.N','pH')]),
                y=yield.data$Zn.crop,
                x=yield.data[,c('locationNearForest','locationMiddle','N.appl.p.ha','MAOM.C','POM.N','pH')]
)

fe.list <- list(
                N=nrow(yield.data),
                K=ncol(yield.data[,c('locationNearForest','locationMiddle','N.appl.p.ha','MAOM.C','POM.N','pH')]),
                y=yield.data$Fe.crop,
                x=yield.data[,c('locationNearForest','locationMiddle','N.appl.p.ha','MAOM.C','POM.N','pH')]
)

## Call Stan models
yield.model <- stan(file = "Stan/final_models/yield-nutrients.stan", 
                    data = yield.list, 
                    control = list(adapt_delta=0.99,max_treedepth=15), chains = 4)
pro.model <- stan(file = "Stan/final_models/yield-nutrients.stan", 
                  data = pro.list, 
                  control = list(adapt_delta=0.99,max_treedepth=15), chains = 4)
zn.model <- stan(file = "Stan/final_models/yield-nutrients.stan", 
                    data = zn.list, 
                    control = list(adapt_delta=0.99,max_treedepth=15), chains = 4)
fe.model <- stan(file = "Stan/final_models/yield-nutrients.stan", 
                    data = fe.list, 
                    control = list(adapt_delta=0.99,max_treedepth=15), chains = 4)

## Model results
### Wheat yield
print(yield.model,pars='beta',probs=c(0.05,0.95))
plot(yield.model,pars=c("beta_std"))

### Protein
print(pro.model,pars='beta',probs=c(0.05,0.95))
plot(pro.model,pars=c("beta_std"))

### Zinc
print(zn.model,pars='beta',probs=c(0.05,0.95))
plot(zn.model,pars=c("beta_std"))

### Iron
print(fe.model,pars='beta',probs=c(0.05,0.95))
plot(fe.model,pars=c("beta_std"))

## Posterior predictive checks
### Wheat yield
y_pred <- extract(yield.model,pars='y_tilde')
y_pred <- unlist(y_pred, use.names=FALSE)
yield.pp.data <- data.frame(
                  c(y_pred,yield.list$y),
                  c(rep("y_pred",length(y_pred)),
                    rep("y_obs",length(yield.list$y)))
                  )
names(yield.pp.data) <- c("y","type")
ggplot(yield.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("Wheat yield") + ylab("Density") +
  scale_fill_manual(values=wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )
rm(y_pred)

### Protein
y_pred <- extract(pro.model,pars='y_tilde')
y_pred <- unlist(y_pred, use.names=FALSE)
pro.pp.data <- data.frame(
                c(y_pred,pro.list$y),
                c(rep("y_pred",length(y_pred)),
                  rep("y_obs",length(pro.list$y)))
              )
names(pro.pp.data) <- c("y","type")
ggplot(pro.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("Grain protein") + ylab("Density") +
  scale_fill_manual(values=wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )
rm(y_pred)

### Zinc
y_pred <- extract(zn.model,pars='y_tilde')
y_pred <- unlist(y_pred, use.names=FALSE)
zn.pp.data <- data.frame(
                c(y_pred,zn.list$y),
                c(rep("y_pred",length(y_pred)),
                  rep("y_obs",length(zn.list$y)))
              )
names(zn.pp.data) <- c("y","type")
ggplot(zn.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("Crop zinc concentration") + ylab("Density") +
  scale_fill_manual(values=wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )
rm(y_pred)

### Iron
y_pred <- extract(fe.model,pars='y_tilde')
y_pred <- unlist(y_pred, use.names=FALSE)
fe.pp.data <- data.frame(
                c(y_pred,fe.list$y),
                c(rep("y_pred",length(y_pred)),
                  rep("y_obs",length(fe.list$y)))
              )
names(fe.pp.data) <- c("y","type")
ggplot(fe.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("Crop iron concentration") + ylab("Density") +
  scale_fill_manual(values=wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )
rm(y_pred)


# QUESTION 2
# Estimating nutrition impacts of nutrient gains from SOM
pro.list <- list(
  N=nrow(yield.data),
  K=6,
  nf=yield.data$locationNearForest,
  mid=yield.data$locationMiddle,
  maom=yield.data$MAOM.N,
  pom=yield.data$POM.N,
  pH=yield.data$pH,
  fert=yield.data$N.appl.p.ha,
  y=yield.data$Protein*100
)
pro.model <- stan(file = "Stan/final_models/nutrients.stan", 
                  data = pro.list, 
                  control = list(adapt_delta=0.99,max_treedepth=15), chains = 4)
print(pro,pars='beta',probs=c(0.05,0.95))
pro <- extract(pro.model,pars='pro_nourished') %>% unlist(use.names=FALSE)
hist(pro)



## Establish and convert parameters
protein.mean <- 0.20/10
zinc.mean <- 0.14/10
protein.5 <- -0.02/10
zinc.5 <- 0.07/10
protein.95 <- 0.43/10
zinc.95 <- 0.20/10

## Individual RDA for nutrients (Ethiopia specific estimate)
zn.rda <- 6.232   # mg / d
fe.rda <- 32.502  # mg / d
pro.rda <- 28.622 # g / d

## Total number of extra people nourished
ppl.nourished <- arsi[,c(1:3,46:49)]

ppl.nourished$pro.mean <- ((arsi$WheatYield * arsi$Area * 1000 * protein.mean)/365) / pro.rda
ppl.nourished$zn.mean <- ((arsi$WheatYield * arsi$Area * 1000 * zinc.mean)/365) / zn.rda

ppl.nourished$pro.min <- ((arsi$WheatYield * arsi$Area * 1000 * protein.5)/365) / pro.rda
ppl.nourished$zn.min <- ((arsi$WheatYield * arsi$Area * 1000 * zinc.5)/365) / zn.rda

ppl.nourished$pro.max <- ((arsi$WheatYield * 1000 * protein.95)/365) / pro.rda
ppl.nourished$zn.max <- ((arsi$WheatYield * 1000 * zinc.95)/365) / zn.rda

ppl.nourished$pro.mean <- ((arsi$WheatYield * 1000 * protein.mean)/365) / pro.rda
ppl.nourished$zn.mean <- ((arsi$WheatYield * 1000 * zinc.mean)/365) / zn.rda

ppl.nourished$pro.min <- ((arsi$WheatYield * 1000 * protein.5)/365) / pro.rda
ppl.nourished$zn.min <- ((arsi$WheatYield * 1000 * zinc.5)/365) / zn.rda

ppl.nourished$pro.max <- ((arsi$WheatYield * 1000 * protein.95)/365) / pro.rda
ppl.nourished$zn.max <- ((arsi$WheatYield * 1000 * zinc.95)/365) / zn.rda


## Plot effects
plot.data <- ppl.nourished %>% 
  group_by(Location) %>% 
  summarise(
    pro_mean = mean(pro.mean,na.rm=T),
    pro_sd = sd(pro.mean,na.rm=T),
    zn_mean = mean(zn.mean,na.rm=T),
    zn_sd = sd(zn.mean,na.rm=T),
    pro_min = mean(pro.min,na.rm=T),
    pro_min_sd = sd(pro.min,na.rm=T),
    zn_min = mean(zn.min,na.rm=T),
    zn_min_sd = sd(zn.min,na.rm=T),
    pro_max = mean(pro.max,na.rm=T),
    pro_max_sd = sd(pro.max,na.rm=T),
    zn_max = mean(zn.max,na.rm=T),
    zn_max_sd = sd(zn.max,na.rm=T)
    ) %>%
  filter(Location!='Forest')

plot.data$Location <- factor(plot.data$Location, levels=c("Road", "Middle", "Near Forest"))

ggplot(plot.data,aes(x=Location,y=pro_mean)) + 
  geom_bar(aes(x=Location,y=pro_max),
           stat="identity",
           color="red",
           fill="white",
           position="dodge",
           alpha=0.25) +
  geom_bar(aes(x=Location,y=pro_min),
           stat="identity",
           alpha=0.25,
           color="red",
           fill="white") +
  geom_bar(stat="identity",fill="grey30",color="grey30") +
  geom_errorbar(aes(ymin=pro_mean-pro_sd,ymax=pro_mean+pro_sd),
                size=0.25,
                width = 0.05,
                color="grey70",
                position = "dodge") +
  geom_errorbar(aes(ymin=pro_min-pro_sd,ymax=pro_min+pro_sd),
                size=0.25,
                width = 0.05,
                color="grey70",
                position = "dodge") +
  geom_errorbar(aes(ymin=pro_max-pro_sd,ymax=pro_max+pro_sd),
                size=0.25,
                width = 0.05,
                color="grey70",
                position = "dodge") +
  theme_bw() + ylab("Number of people\n") + xlab("") + 
  labs(title="Protein",subtitle="Extra people whose protein need could be met by increasing SOM by 1%") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )

ggplot(plot.data,aes(x=Location,y=zn_mean)) + 
  geom_bar(aes(x=Location,y=zn_max),
           stat="identity",
           color="red",
           fill="white",
           position="dodge",
           alpha=0.25) +
  geom_bar(aes(x=Location,y=zn_min),
           stat="identity",
           alpha=0.25,
           color="red",
           fill="white") +
  geom_bar(stat="identity",fill="grey30",color="grey30") +
  geom_errorbar(aes(ymin=zn_mean-zn_sd,ymax=zn_mean+zn_sd),
                size=0.25,
                width = 0.05,
                color="grey70",
                position = "dodge") +
  geom_errorbar(aes(ymin=zn_min-zn_sd,ymax=zn_min+zn_sd),
                size=0.25,
                width = 0.05,
                color="grey70",
                position = "dodge") +
  geom_errorbar(aes(ymin=zn_max-zn_sd,ymax=zn_max+zn_sd),
                size=0.25,
                width = 0.05,
                color="grey70",
                position = "dodge") +
  theme_bw() + ylab("Number of people\n") + xlab("") + 
  labs(title="Zinc",subtitle="Extra people whose zinc need could be met by increasing SOM by 1%") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )


# QUESTION 3
# Does SOM vary along the same landscape gradient as yield?

## Set up data for Stan
som.vars <- c('Location','POM.C','MAOM.C','SIR','typeHome','typeWheat','Fines','Texture',
              'locationForest','locationMiddle','locationNearForest','locationRoad','pH')
som.data <- arsi[complete.cases(arsi[,som.vars]),c(som.vars)]

## Convert data to list
pomList <- list(
                N = nrow(som.data),
                K = ncol(som.data[,c('typeHome','locationMiddle','locationNearForest','locationForest','Texture')]),
                y = som.data$POM.C,
                x = som.data[,c('typeHome','locationMiddle','locationNearForest','locationForest','Texture')]
            )
maomList <- list(
                N = nrow(som.data),
                K = ncol(som.data[,c('typeHome','locationMiddle','locationNearForest','locationForest','Texture')]),
                y = som.data$MAOM.C,
                x = som.data[,c('typeHome','locationMiddle','locationNearForest','locationForest','Texture')]
            )
sirList <- list(
                N = nrow(som.data),
                K = ncol(som.data[,c('typeHome','locationMiddle','locationNearForest','locationForest','Texture')]),
                y = som.data$SIR,
                x = som.data[,c('typeHome','locationMiddle','locationNearForest','locationForest','Texture')]
            )

## Execute models
POM <- stan(file = "Stan/final_models/som.stan",
                 data = pomList,
                 iter = 2000, chains = 4)
MAOM <- stan(file = "Stan/final_models/som.stan",
             data = maomList,
             iter = 2000, chains = 4)
SIR <- stan(file = "Stan/final_models/som.stan",
            data = sirList,
            iter = 2000, chains = 4)

## Print model results
print(POM,pars=c('beta'),probs=c(0.05,0.95))
plot(POM,pars=c('beta_std'))

print(MAOM,pars=c('beta'),probs=c(0.05,0.95))
plot(MAOM,pars=c('beta_std'))

print(SIR,pars=c('beta'),probs=c(0.05,0.95))
plot(SIR,pars=c('beta_std'))

## Posterior predictive checks
### POM
y_pred_pom <- extract(POM,pars='y_tilde')
y_pred_pom <- unlist(y_pred_pom, use.names=FALSE)

pom.pp.data <- data.frame(
                          c(y_pred_pom,pomList$y),
                          c(rep("y_pred",length(y_pred_pom)),
                            rep("y_obs",length(pomList$y)))
                          )
names(pom.pp.data) <- c("y","type")

ggplot(pom.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("POM") + ylab("Density") +
  scale_fill_manual(values=wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )

### MAOM
y_pred_maom <- extract(MAOM,pars='y_tilde')
y_pred_maom <- unlist(y_pred_maom, use.names=FALSE)

maom.pp.data <- data.frame(
                            c(y_pred_maom,maomList$y),
                            c(rep("y_pred",length(y_pred_maom)),
                              rep("y_obs",length(maomList$y)))
                            )
names(maom.pp.data) <- c("y","type")

ggplot(maom.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("MAOM") + ylab("Density") +
  scale_fill_manual(values=wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )

### SIR
y_pred_sir <- extract(SIR,pars='y_tilde')
y_pred_sir <- unlist(y_pred_sir, use.names=FALSE)

sir.pp.data <- data.frame(
                          c(y_pred_sir,sirList$y),
                          c(rep("y_pred",length(y_pred_sir)),
                            rep("y_obs",length(sirList$y)))
                          )
names(sir.pp.data) <- c("y","type")

ggplot(sir.pp.data, aes(x=y)) + 
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("Substrate-induced Respiration") + ylab("Density") +
  scale_fill_manual(values=wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )


# QUESTION 4
# How do results of Question 3 play out for wheat fields

## Filter data to wheat fields
som.data.wheat <- filter(som.data,typeHome==0 & locationForest == 0)

## Convert data to list
pomList <- list(
  N = nrow(som.data.wheat),
  K = ncol(som.data.wheat[,c('locationMiddle','locationNearForest','Texture')]),
  y = som.data.wheat$POM.C,
  x = som.data.wheat[,c('locationMiddle','locationNearForest','Texture')]
)
maomList <- list(
  N = nrow(som.data.wheat),
  K = ncol(som.data.wheat[,c('locationMiddle','locationNearForest','Texture')]),
  y = som.data.wheat$MAOM.C,
  x = som.data.wheat[,c('locationMiddle','locationNearForest','Texture')]
)
sirList <- list(
  N = nrow(som.data.wheat),
  K = ncol(som.data.wheat[,c('locationMiddle','locationNearForest','Texture')]),
  y = som.data.wheat$SIR,
  x = som.data.wheat[,c('locationMiddle','locationNearForest','Texture')]
)

## Execute models
POM <- stan(file = "Stan/final_models/som-wheat.stan",
            data = pomList,
            iter = 2000, chains = 4)
MAOM <- stan(file = "Stan/final_models/som-wheat.stan",
             data = maomList,
             iter = 2000, chains = 4)
SIR <- stan(file = "Stan/final_models/som-wheat.stan",
            data = sirList,
            iter = 2000, chains = 4)

## Print model results
print(POM,pars=c('beta'),probs=c(0.05,0.95))
plot(POM,pars=c('beta_std'))

print(MAOM,pars=c('beta'),probs=c(0.05,0.95))
plot(MAOM,pars=c('beta_std'))

print(SIR,pars=c('beta'),probs=c(0.05,0.95))
plot(SIR,pars=c('beta_std'))

## Posterior predictive checks
### POM
y_pred_pom <- extract(POM,pars='y_tilde')
y_pred_pom <- unlist(y_pred_pom, use.names=FALSE)

pom.pp.data <- data.frame(
  c(y_pred_pom,pomList$y),
  c(rep("y_pred",length(y_pred_pom)),
    rep("y_obs",length(pomList$y)))
)
names(pom.pp.data) <- c("y","type")

ggplot(pom.pp.data, aes(x=y)) +
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("POM") + ylab("Density") +
  scale_fill_manual(values=wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )

### MAOM
y_pred_maom <- extract(MAOM,pars='y_tilde')
y_pred_maom <- unlist(y_pred_maom, use.names=FALSE)

maom.pp.data <- data.frame(
  c(y_pred_maom,maomList$y),
  c(rep("y_pred",length(y_pred_maom)),
    rep("y_obs",length(maomList$y)))
)
names(maom.pp.data) <- c("y","type")

ggplot(maom.pp.data, aes(x=y)) +
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("MAOM") + ylab("Density") +
  scale_fill_manual(values=wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )

### SIR
y_pred_sir <- extract(SIR,pars='y_tilde')
y_pred_sir <- unlist(y_pred_sir, use.names=FALSE)

sir.pp.data <- data.frame(
  c(y_pred_sir,sirList$y),
  c(rep("y_pred",length(y_pred_sir)),
    rep("y_obs",length(sirList$y)))
)
names(sir.pp.data) <- c("y","type")

ggplot(sir.pp.data, aes(x=y)) +
  geom_density(aes(group=type, fill=type), alpha=0.75) + theme_bw() +
  xlab("Substrate-induced Respiration") + ylab("Density") +
  scale_fill_manual(values=wes_palette("Royal1",n=2)) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85,0.55),
    panel.grid = element_blank(),
    panel.background = element_rect(fill="white"),
    plot.background = element_rect(fill="white")
  )
