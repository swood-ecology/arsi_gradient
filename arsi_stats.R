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
# car::qqPlot(lm(log(WheatYield)~MAOM.C+POM.N+Trd))
# car::ncvTest(lm(log(WheatYield)~MAOM.C+POM.N+Trd))
# car::residualPlots(lm(log(WheatYield)~MAOM.C+POM.N+Trd))
# car::crPlots(lm(log(WheatYield)~MAOM.C+POM.N+Trd))
# car::avPlots(lm(log(WheatYield)~MAOM.C+POM.N+Trd))

# Drop NA values
yield.vars <- c("WheatYield","N.appl.amt","MAOM.C","POM.N","Trd")
yield.data <- arsi_full_data[complete.cases(arsi_full_data[,yield.vars]),c(yield.vars)]
yield.data <- as.data.frame(yield.data)

# Fit model
yield.list <- list(N=nrow(yield.data),
                    K=4,
                    fert=yield.data$N.appl.amt,
                    maom=yield.data$MAOM.C,
                    pom=yield.data$POM.N,
                    trd=yield.data$Trd,
                    y=yield.data$WheatYield
                  )
yield.model <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/yield.stan", data = yield.list, control = list(adapt_delta=0.99,max_treedepth=15), chains = 4)


# Look at results
print(yield.model, probs=c(.025,.975),
      pars=c('alpha_std','beta_std','sigma_std')) # Std coefficients
print(yield.model, probs=c(.025,.975),
      pars=c('beta_fert','beta_maom','beta_pom','beta_trd','sigma')) # Natural coefficients
plot(yield.model,pars=c('beta_std','sigma_std','alpha_std'))
traceplot(yield.model, inc_warmup = F,pars=c('beta','sigma','alpha'))

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

# Standardize variables
arsi_full_data$pH.z <- (arsi_full_data$pH - mean(arsi_full_data$pH,na.rm=T))/2*sd(arsi_full_data$pH,na.rm=T)
arsi_full_data$Silt.z <- (arsi_full_data$Silt - mean(arsi_full_data$Silt,na.rm=T))/2*sd(arsi_full_data$Silt,na.rm=T)
arsi_full_data$MAOM.N.z <- (arsi_full_data$MAOM.N - mean(arsi_full_data$MAOM.N,na.rm=T))/2*sd(arsi_full_data$MAOM.N,na.rm=T)

# Initial models
lm.zn <- lm(log(Zn.crop)~ N.appl.amt + pH + MAOM.N + POM.N + Silt,data=arsi_full_data)
summary(standardize(lm.zn))
lm.fe <- lm(log(Fe.crop)~ N.appl.amt + pH + MAOM.N + POM.N + Silt + Trd,data=arsi_full_data)
summary(standardize(lm.fe))

# Stan models
nutr.vars <- c("Zn.crop","Fe.crop","N.z","pH.z","Silt.z","MAOM.N.z","POM.N.z","Trd.z")
nutr.data <- arsi_full_data[complete.cases(arsi_full_data[,nutr.vars]),c(nutr.vars)]
nutr.data <- as.data.frame(nutr.data)

zn.list <- list(N=nrow(nutr.data),
                K=6,
                x=nutr.data[,c(3:8)],
                y=nutr.data$Zn.crop
)

zn.model <- stan(file = "~/Box Sync/Work/GitHub/arsi_gradient/Stan/crop_nutrients.stan", data = zn.list, control = list(adapt_delta=0.99,max_treedepth=17), chains=4)
print(zn.model, probs=c(.025,.975), pars=c('alpha','beta','sigma'))
plot(zn.model, pars=c('alpha','beta','sigma'))
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
                x=nutr.data[,c(3:8)],
                y=nutr.data$Fe.crop
)

fe.model <- stan(file = "~/Box Sync/Work/Writing/Manuscripts/Unsubmitted/arsi_negele/Stan/crop_nutrients.stan", data = fe.list, control = list(adapt_delta=0.999,max_treedepth=20), chains=4)

print(fe.model, probs=c(.025,.975), pars=c("alpha","beta","sigma"))
plot(fe.model, pars=c("alpha","beta","sigma"))

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

# Compare lm and lmm
lm(log(POM.C)~ MAOM.C,arsi_full_data) %>%
  standardize() %>%
    summary()
lme4::lmer(log(POM.C)~ MAOM.C + (1|Location),arsi_full_data) %>%
  standardize() %>%  
    summary()

# Set up data for Stan
# Remove NAs
som.data <- as.data.frame(arsi_full_data[complete.cases(arsi_full_data[,c(2,3,12:14,17,21:24,28,38)]),c(2,3,12:14,17,21:24,28,38)])

# Create numeric clustering variable
uniq <- unique(som.data$Location)
J <- length(uniq) 
loc <- rep(NA,J)

for (i in 1:J) {
  loc[som.data$Location == uniq[i]] <- i
}

# Define data
dataList <- list(
  N = length(som.data$POM.C),
  y = som.data$POM.C,
  x = som.data$MAOM.C,
  re = loc,
  J = max(loc)
)

# Execute model
POM <- stan(file = "~/Box Sync/Work/Writing/Manuscripts/Unsubmitted/arsi_negele/Stan/som_RE.stan",
     data = dataList,
     iter = 2000, chains = 2)

# Print model
print(POM)
plot(POM)
traceplot(POM,inc_warmup=F)






lm(log(POM.C)~ Type + Dmun + Zn + Cu + Mn + MAOM.C,arsi_full_data) %>%
  standardize() %>%
  summary()

lm.MAOM <- lm(log(MAOM.C)~ Type + Dmun + POM.C + Zn + Cu + Mn ,arsi_full_data )
lm.SIR <- lm(log(SIR)~ Type + Dmun + MAOM.C + POM.C + Sand + pH + Cu + Mn + Zn,arsi_full_data )

summary(standardize(lm.MAOM))
summary(standardize(lm.SIR))

stargazer::stargazer(lm.POM, lm.MAOM, lm.SIR, title="Organic matter model results", align=TRUE, omit.stat=c("LL","ser","f"), type="html")

POM <- stan_lm(log(POM.C)~ Type + Dmun + Cu + Mn + Zn + MAOM.C, prior=R2(location=0.5),data = arsi_full_data)
pp_check(POM)
