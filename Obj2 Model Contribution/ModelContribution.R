####------------------------------------------------------------------------------------#
#### Model Contribution by Covariates ####
####------------------------------------------------------------------------------------#

library(tidyverse)
library(lubridate)
library(scales)
library(R2jags)
library(MCMCvis)
library(loo)
library(viridis)
library(beepr)
library(GGally)
library(ggmcmc)



####-----------------------------------------#
#### Load and combine data ####
####-----------------------------------------#

# bootstrap-corrected GSI proportions
dat <- read_csv("GSI Analysis/By Section and Year/UpperSnake_GSI_SectionYear_BootstrappedProportions.csv") %>% 
  filter(mixture_collection != "snake_jldcattlemens") %>% 
  rename(mixture = mixture_collection, ppn_bscor = bs_corrected_repunit_ppn) %>%
  mutate(section = str_sub(mixture, end = -6), year = as.numeric(str_sub(mixture, start = -4))) %>%
  arrange(mixture, repunit)

# sample size per collection/mixture
ids <- read_csv("LabFieldIDs.csv")
data <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Methods/GSI sampling/Snake_GSI_field data_2020-2022_GenOnly_DropNoDrop_090823edit.csv")
metadata <- left_join(ids, data) %>% 
  select(indiv, CORRgenID, year, site, TL.mm, TL.in, FL.mm, weight.g, weight.lbs, collection.type, sizecat) %>% 
  filter(collection.type == "mixture", site != "snake_jldcattlemens")
metadata[duplicated(metadata$indiv),] # check for duplicates
metayr <- metadata %>% select(indiv, year)
mix_input <- read_csv("Baseline Data and Testing/UpperSnakeRiver_GTseq_InputData_NoSibs_clean_mixture.csv") %>%
  left_join(metayr) %>% mutate(collection = paste(collection, year, sep = "_")) %>% 
  group_by(collection) %>% summarize(total = n()) %>% ungroup() %>% rename(mixture = collection)

# groundwater metrics
gwmet <- read_csv("Landscape Covariates/Groundwater/GroundwaterMetrics_raw_RepUnits.csv") %>% 
  rename(repunit = site) %>% 
  mutate(logarea = log(areasqkm), loggwi = log(gwi_iew05km)) %>% 
  mutate(z_logarea = as.numeric(scale(logarea)), z_loggwi = as.numeric(scale(loggwi))) %>%
  arrange(repunit)
# gwmet$repunit <- str_replace_all(gwmet$repunit, "deadman _greys", "deadman_greys")

# flowline distances
dist <- read_csv("Landscape Covariates/Flowline Distance/SnakeRiverSections_RepUnits_FlowlineDistance.csv") %>%
  mutate(z_dist = as.numeric(scale(distkm)))
# dist$repunit <- str_replace_all(dist$repunit, "deadman _greys", "deadman_greys")

# scale tibble
scltib <- tibble(param = c("logarea", "loggwi", "dist"), 
                 mean = c(attributes(scale(gwmet$logarea))$`scaled:center`, attributes(scale(gwmet$loggwi))$`scaled:center`, attributes(scale(dist$distkm))$`scaled:center`), 
                 sd = c(attributes(scale(gwmet$logarea))$`scaled:scale`, attributes(scale(gwmet$loggwi))$`scaled:scale`, attributes(scale(dist$distkm))$`scaled:scale`))
write_csv(scltib, "Model Contribution/GSI_Covariates_MeansStDevs.csv")

# join tibbles
dat <- dat %>% left_join(mix_input) %>% left_join(gwmet) %>% left_join(dist)
dat <- dat %>% mutate(logarea = log(areasqkm), loggwi = log(gwi_iew05km))
view(dat)
write_csv(dat, "Model Contribution/GSI_ModelContribution_DataTable.csv")
dat <- read_csv("Model Contribution/GSI_ModelContribution_DataTable.csv")

# pairs plots - no strong correlations
ggpairs(dat %>% select(z_logarea, z_loggwi, z_dist))


####---------------------------------------------------#
#### Values for JAGS ####
####---------------------------------------------------#

# Bootstrap-corrected mixture "abundance"
# for multinomial-Dirichlet model
gsi_dat <- readRDS("GSI Analysis/By Section and Year/UpperSnake_GSI_SectionYear_output.RDS")
gsi_bscorr <- gsi_dat$bootstrapped_proportions %>%
  left_join(gsi_dat$indiv_posteriors %>% group_by(indiv, mixture_collection) %>% summarize(n = n()) %>% group_by(mixture_collection) %>% summarize(sampsize = n())) %>%
  mutate(bscorrnum = bs_corrected_repunit_ppn*sampsize, section = str_sub(mixture_collection, end = -6), year = str_sub(mixture_collection, start = -4)) %>%
  arrange(mixture_collection, repunit) %>% select(mixture_collection, section, year, repunit, bscorrnum) %>% spread(key = repunit, value = bscorrnum)
Y <- gsi_bscorr[,-c(1:3)]
Yint <- round(Y, digits = 0)
Nint <- rowSums(Yint) # sample size


# for Dirichlet only model
dat2 <- dat %>% arrange(mixture, repunit) %>% select(mixture, section, year, repunit, ppn_bscor) %>% spread(key = repunit, value = ppn_bscor)
Yppn <- (dat2[,-c(1:3)] + 0.0000001) / rowSums(dat2[,-c(1:3)] + 0.0000001)


# covariates
dist_cov <- dat %>% left_join(dist) %>% arrange(mixture, repunit) %>% select(mixture, repunit, z_dist) %>% spread(key = repunit, value = z_dist) %>% select(-mixture)
area_cov <- c(gwmet$z_logarea)
gw_cov <- c(gwmet$z_loggwi)

# misc values
nmix <- dim(Y)[1] # number of mixtures
nrgs <- dim(Y)[2] # number of reporting groups
nsect <- length(unique(gsi_bscorr$section)) # number of sections
sectid <- as.numeric(as.factor(gsi_bscorr$section)) # section IDs as numeric/factors
nyear <- length(unique(gsi_bscorr$year)) # number of years
yearid <- as.numeric(as.factor(gsi_bscorr$year))


####---------------------------------------------------#
#### Model in JAGS ####
####---------------------------------------------------#

#### 1. Multinomial-Dirichlet Regression using bootstrapped fish abundance ####
#### Dropped the overdispersion term (sensu Chong and Spencer (2018) and Ole Shelton (NOAA) suggestions) 
#### as it considerably affect model convergence and resulted in very inefficient MCMC sampling
#### and preliminary model runs indicated that data were not over or underdispersed (epsilon broadly overlapped 0 and transformation/exp broadly overlapped 1 for ~most mixtures)

# specify data 
jags.data <- list("nmix" = nmix, "nrgs" = nrgs, "nsect" = nsect, "sectid" = sectid, "nyear" = nyear, "yearid" = yearid, "Y" = Yint, "N" = Nint, "dist" = dist_cov, "area" = area_cov, "gw" = gw_cov)      # specify data
Sys.time()


# Model 0 - Dirichlet only to address Cliff's comment re: overparameterization
jags.data <- list("nmix" = nmix, "nrgs" = nrgs, "nsect" = nsect, "sectid" = sectid, "nyear" = nyear, "yearid" = yearid, "p" = Yppn, "dist" = dist_cov, "area" = area_cov, "gw" = gw_cov)      # specify data
jags.data <- list("nmix" = nmix, "nrgs" = nrgs, "nsect" = nsect, "sectid" = sectid, "nyear" = nyear, "yearid" = yearid, "Y" = Yint, "N" = rowSums(Yint), "dist" = dist_cov, "area" = area_cov, "gw" = gw_cov)      # specify data
params <- c("p", "a", "delta", "alpha", "beta1", "beta2", "beta3", "loglik", "res", "Y.pred", "theta")   # parameters to monitor
mod_d1 <- jags.parallel(jags.data, inits = NULL, parameters.to.save = params, model.file = "Model Contribution/Candidate models/JAGS models/Dirich_Simple_1.txt", n.chains = 3, n.thin = 10, n.burnin = 1000, n.iter = 16000, DIC = TRUE)
MCMCtrace(mod_d1, ind = TRUE, filename = "Model Contribution/Candidate models/Traceplots/Dirich_1_Traceplots.pdf")
# saveRDS(mod_d1, "Model Contribution/Candidate models/Model output/Dirich_1_modelout.RDS")
mod_d1$BUGSoutput$summary[,8][mod_d1$BUGSoutput$summary[,8] > 1.01]

Sys.time()


# Model 1
params <- c("p", "a", "delta", "alpha", "beta1", "beta2", "beta3", "loglik", "res", "Y.pred")   # parameters to monitor
mod_md1 <- jags.parallel(jags.data, inits = NULL, parameters.to.save = params, model.file = "Model Contribution/Candidate models/JAGS models/MultiDirich_Simple_1.txt", n.chains = 10, n.thin = 500, n.burnin = 10000, n.iter = 160000, DIC = TRUE)
MCMCtrace(mod_md1, ind = TRUE, filename = "Model Contribution/Candidate models/Traceplots/MultiDirich_1_Traceplots.pdf")
saveRDS(mod_md1, "Model Contribution/Candidate models/Model output/MultiDirich_1_modelout.RDS")
beep()
Sys.time()

# Model 2 
params <- c("p", "a", "delta", "alpha", "beta1", "beta2", "loglik", "res", "Y.pred")   # parameters to monitor
mod_md2 <- jags.parallel(jags.data, inits = NULL, parameters.to.save = params, model.file = "Model Contribution/Candidate models/JAGS models/MultiDirich_Simple_2.txt", n.chains = 10, n.thin = 500, n.burnin = 10000, n.iter = 160000, DIC = TRUE)
MCMCtrace(mod_md2, ind = TRUE, filename = "Model Contribution/Candidate models/Traceplots/MultiDirich_2_Traceplots.pdf")
saveRDS(mod_md2, "Model Contribution/Candidate models/Model output/MultiDirich_2_modelout.RDS")
beep()
Sys.time()

# Model 3
params <- c("p", "a", "delta", "alpha", "beta1", "beta3", "loglik", "res", "Y.pred")   # parameters to monitor
mod_md3 <- jags.parallel(jags.data, inits = NULL, parameters.to.save = params, model.file = "Model Contribution/Candidate models/JAGS models/MultiDirich_Simple_3.txt", n.chains = 10, n.thin = 500, n.burnin = 10000, n.iter = 160000, DIC = TRUE)
MCMCtrace(mod_md3, ind = TRUE, filename = "Model Contribution/Candidate models/Traceplots/MultiDirich_3_Traceplots.pdf")
saveRDS(mod_md3, "Model Contribution/Candidate models/Model output/MultiDirich_3_modelout.RDS")
beep()
Sys.time()

# Model 4
params <- c("p", "a", "delta", "alpha", "beta2", "beta3", "loglik", "res", "Y.pred")   # parameters to monitor
mod_md4 <- jags.parallel(jags.data, inits = NULL, parameters.to.save = params, model.file = "Model Contribution/Candidate models/JAGS models/MultiDirich_Simple_4.txt", n.chains = 10, n.thin = 500, n.burnin = 10000, n.iter = 160000, DIC = TRUE)
MCMCtrace(mod_md4, ind = TRUE, filename = "Model Contribution/Candidate models/Traceplots/MultiDirich_4_Traceplots.pdf")
saveRDS(mod_md4, "Model Contribution/Candidate models/Model output/MultiDirich_4_modelout.RDS")
beep()
Sys.time()

# Model 5
params <- c("p", "a", "delta", "alpha", "beta1", "loglik", "res", "Y.pred")   # parameters to monitor
mod_md5 <- jags.parallel(jags.data, inits = NULL, parameters.to.save = params, model.file = "Model Contribution/Candidate models/JAGS models/MultiDirich_Simple_5.txt", n.chains = 10, n.thin = 500, n.burnin = 10000, n.iter = 160000, DIC = TRUE)
MCMCtrace(mod_md5, ind = TRUE, filename = "Model Contribution/Candidate models/Traceplots/MultiDirich_5_Traceplots.pdf")
saveRDS(mod_md5, "Model Contribution/Candidate models/Model output/MultiDirich_5_modelout.RDS")
beep()
Sys.time()

# Model 6
params <- c("p", "a", "delta", "alpha", "beta2", "loglik", "res", "Y.pred")   # parameters to monitor
mod_md6 <- jags.parallel(jags.data, inits = NULL, parameters.to.save = params, model.file = "Model Contribution/Candidate models/JAGS models/MultiDirich_Simple_6.txt", n.chains = 10, n.thin = 500, n.burnin = 10000, n.iter = 160000, DIC = TRUE)
MCMCtrace(mod_md6, ind = TRUE, filename = "Model Contribution/Candidate models/Traceplots/MultiDirich_6_Traceplots.pdf")
saveRDS(mod_md6, "Model Contribution/Candidate models/Model output/MultiDirich_6_modelout.RDS")
beep()
Sys.time()

# Model 7
params <- c("p", "a", "delta", "alpha", "beta3", "loglik", "res", "Y.pred")   # parameters to monitor
mod_md7 <- jags.parallel(jags.data, inits = NULL, parameters.to.save = params, model.file = "Model Contribution/Candidate models/JAGS models/MultiDirich_Simple_7.txt", n.chains = 10, n.thin = 500, n.burnin = 10000, n.iter = 160000, DIC = TRUE)
MCMCtrace(mod_md7, ind = TRUE, filename = "Model Contribution/Candidate models/Traceplots/MultiDirich_7_Traceplots.pdf")
saveRDS(mod_md7, "Model Contribution/Candidate models/Model output/MultiDirich_7_modelout.RDS")
beep()
Sys.time()

# Model 8
params <- c("p", "a", "delta", "alpha", "loglik", "res", "Y.pred")   # parameters to monitor
mod_md8 <- jags.parallel(jags.data, inits = NULL, parameters.to.save = params, model.file = "Model Contribution/Candidate models/JAGS models/MultiDirich_Simple_8.txt", n.chains = 10, n.thin = 500, n.burnin = 10000, n.iter = 160000, DIC = TRUE)
MCMCtrace(mod_md8, ind = TRUE, filename = "Model Contribution/Candidate models/Traceplots/MultiDirich_8_Traceplots.pdf")
saveRDS(mod_md8, "Model Contribution/Candidate models/Model output/MultiDirich_8_modelout.RDS")
beep()
Sys.time()


# check for convergence, high R-hat, excluding residuals
mod_md1$BUGSoutput$summary[!str_detect(rownames(mod_md1$BUGSoutput$summary), "res"),8][mod_md1$BUGSoutput$summary[!str_detect(rownames(mod_md1$BUGSoutput$summary), "res"),8] > 1.01]
mod_md2$BUGSoutput$summary[!str_detect(rownames(mod_md2$BUGSoutput$summary), "res"),8][mod_md2$BUGSoutput$summary[!str_detect(rownames(mod_md2$BUGSoutput$summary), "res"),8] > 1.01]
mod_md3$BUGSoutput$summary[!str_detect(rownames(mod_md3$BUGSoutput$summary), "res"),8][mod_md3$BUGSoutput$summary[!str_detect(rownames(mod_md3$BUGSoutput$summary), "res"),8] > 1.01]
mod_md4$BUGSoutput$summary[!str_detect(rownames(mod_md4$BUGSoutput$summary), "res"),8][mod_md4$BUGSoutput$summary[!str_detect(rownames(mod_md4$BUGSoutput$summary), "res"),8] > 1.01]
mod_md5$BUGSoutput$summary[!str_detect(rownames(mod_md5$BUGSoutput$summary), "res"),8][mod_md5$BUGSoutput$summary[!str_detect(rownames(mod_md5$BUGSoutput$summary), "res"),8] > 1.01]
mod_md6$BUGSoutput$summary[!str_detect(rownames(mod_md6$BUGSoutput$summary), "res"),8][mod_md6$BUGSoutput$summary[!str_detect(rownames(mod_md6$BUGSoutput$summary), "res"),8] > 1.01]
mod_md7$BUGSoutput$summary[!str_detect(rownames(mod_md7$BUGSoutput$summary), "res"),8][mod_md7$BUGSoutput$summary[!str_detect(rownames(mod_md7$BUGSoutput$summary), "res"),8] > 1.01]
mod_md8$BUGSoutput$summary[!str_detect(rownames(mod_md8$BUGSoutput$summary), "res"),8][mod_md8$BUGSoutput$summary[!str_detect(rownames(mod_md8$BUGSoutput$summary), "res"),8] > 1.01]


# read in fitted models
mod_md1 <- readRDS("Model Contribution/Candidate models/Model output/MultiDirich_1_modelout.RDS")
mod_md2 <- readRDS("Model Contribution/Candidate models/Model output/MultiDirich_2_modelout.RDS")
mod_md3 <- readRDS("Model Contribution/Candidate models/Model output/MultiDirich_3_modelout.RDS")
mod_md4 <- readRDS("Model Contribution/Candidate models/Model output/MultiDirich_4_modelout.RDS")
mod_md5 <- readRDS("Model Contribution/Candidate models/Model output/MultiDirich_5_modelout.RDS")
mod_md6 <- readRDS("Model Contribution/Candidate models/Model output/MultiDirich_6_modelout.RDS")
mod_md7 <- readRDS("Model Contribution/Candidate models/Model output/MultiDirich_7_modelout.RDS")
mod_md8 <- readRDS("Model Contribution/Candidate models/Model output/MultiDirich_8_modelout.RDS")



psd1 <- mod_d1$BUGSoutput$summary[str_detect(rownames(mod_d1$BUGSoutput$summary), "a"),5]
psmd1 <- mod_md1$BUGSoutput$summary[str_detect(rownames(mod_md1$BUGSoutput$summary), "a"),5]

plot((psd1[1:884]) ~ log(psmd1[1:884]))

####---------------------------------------------------#
#### Model Selection Using LOO ####
####---------------------------------------------------#

fit_list <- list(mod_md1, mod_md2, mod_md3, mod_md4, mod_md5, mod_md6, mod_md7, mod_md8)
loo_list <- list()
nsamp <- 300 # number of samples retained per chain
for (i in 1:length(fit_list)) {
  m <- fit_list[[i]]
  loglik <- t(apply(m$BUGSoutput$sims.list$loglik, 1, c))
  reff <- relative_eff(exp(loglik), chain_id = c(rep(1,nsamp),rep(2,nsamp),rep(3,nsamp),rep(4,nsamp),rep(5,nsamp),rep(6,nsamp),rep(7,nsamp),rep(8,nsamp),rep(9,nsamp),rep(10,nsamp)))
  loo_list[[i]] <- loo(loglik, r_eff = reff)
  print(i)
}
lc <- loo_compare(loo_list)
print(lc, simplify = FALSE, digits = 2)
plot(loo_list[[5]])

# using DIC b/c can't figure out how to properly calculate log-likelihoods for LOO
dic_tib <- tibble(model = c(1:8), DIC = c(mod_md1$BUGSoutput$DIC,
                                          mod_md2$BUGSoutput$DIC,
                                          mod_md3$BUGSoutput$DIC,
                                          mod_md4$BUGSoutput$DIC,
                                          mod_md5$BUGSoutput$DIC,
                                          mod_md6$BUGSoutput$DIC,
                                          mod_md7$BUGSoutput$DIC,
                                          mod_md8$BUGSoutput$DIC)) %>% arrange(DIC)
# 
# 
# params <- c("p", "a", "delta", "alpha", "beta1", "beta2", "beta3", "loglik", "res", "Y.pred")  # parameters to monitor
# mod_md8 <- jags(jags.data, inits = NULL, parameters.to.save = params, 
#                          model.file = "Model Contribution/Candidate models/JAGS models/MultiDirich_OverdispSimple_Contributions_8.txt", 
#                          n.chains = 3, n.thin = 10, n.burnin = 100, n.iter = 200, DIC = TRUE)
# plot(as.numeric(mod_md8$BUGSoutput$mean$Y.pred) ~ as.numeric(as.matrix(Yint)))
# abline(a = 0, b = 1, col = "red")
# plot(as.numeric(mod_md8$BUGSoutput$mean$res) ~ as.numeric(as.matrix(Yint)))
# 
# mod_md8 <- rjags::jags.model(file = "Model Contribution/Candidate models/JAGS models/MultiDirich_OverdispSimple_Contributions_8.txt",
#                              n.chains = 3, n.adapt = 200, data = jags.data)
# update(mod_md8, 400)
# s <- rjags::coda.samples(mod_md8, params, n.iter = 400)
# plot(s[,"loglik[1,2]"])


####---------------------------------------------------#
#### Extract MCMC samples and parameter summary ####
####---------------------------------------------------#

top_mod <- mod_d1
full_mod <- readRDS("Model Contribution/Candidate models/Model output/MultiDirich_8_modelout.RDS")
top_mod <- readRDS("Model Contribution/Candidate models/Model output/MultiDirich_7_modelout.RDS")

# generate MCMC samples and store as an array
modelout <- top_mod$BUGSoutput
McmcList <- vector("list", length = dim(modelout$sims.array)[2])
for(i in 1:length(McmcList)) { McmcList[[i]] = as.mcmc(modelout$sims.array[,i,]) }
# rbind MCMC samples from 10 chains 
Mcmcdat <- rbind(McmcList[[1]], McmcList[[2]], McmcList[[3]], McmcList[[4]], McmcList[[5]], McmcList[[6]], McmcList[[7]], McmcList[[8]], McmcList[[9]], McmcList[[10]])
Mcmcdat <- rbind(McmcList[[1]], McmcList[[2]], McmcList[[3]])
param.summary <- modelout$summary


####---------------------------------------------------#
#### Model Diagnostic Plots ####
####---------------------------------------------------#

# subset expected and observed MCMC samples
ppdat_exp <- top_mod$BUGSoutput$mean$Y.pred
ppdat_obs <- as.matrix(Yint)

# Bayesian p-value
sum(ppdat_exp > ppdat_obs) / length(ppdat_obs)

jpeg("Model Contribution/Effect plots/MultiDirich_AdditiveEffects_Diagnostics.jpg", units = "in", width = 8, height = 4, res = 1000)
par(mfrow = c(1,2), mar = c(4,4,1,1), mgp = c(2.5,1,0))
# Posterior Predictive Check: plot median posterior expected length ~ observed length with Bayesian p-value
# jpeg("Growth by Covars/Figures/GrowthByCovars_PPCheck.jpg", units = "in", width = 5, height = 4.5, res = 1500)
# par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(ppdat_exp ~ ppdat_obs, xlab = "Observed", ylab = "Expected", bty = "l")
abline(a = 0, b = 1, lty = 2)
legend("topleft", "(a)", bty = "n")
# dev.off()

# histogram of residuals
# jpeg("Growth by Covars/Figures/GrowthByCovars_ResidHist.jpg", units = "in", width = 5, height = 4.5, res = 1500)
# par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
hist(top_mod$BUGSoutput$mean$res, main = "", xlab = "Model residuals")
box(bty = "l")
abline(v = mean(top_mod$BUGSoutput$mean$res), col = "red", lwd = 2)
legend("topleft", "(b)", bty = "n")
# abline(v = median(unlist(ppdat_obs - ppdat_exp)), col = "red", lwd = 2)
dev.off()

# RMSE
sqrt(sum((ppdat_exp - ppdat_obs)^2)/length(ppdat_exp))

####---------------------------------------------------#
#### Dot plots #### 
####---------------------------------------------------#

# Parameter (beta) estimates - dot plot
mod.gg <- ggs(as.mcmc(top_mod))
jpeg("Model Contribution/Effect plots/MultiDirich_AdditiveEffects_betas_dotplot.jpg", units = "in", width = 4, height = 2, res = 1000)
ggs_caterpillar(D = mod.gg, family = "^beta", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = c("Distance", "Area", "Groundwater"))
dev.off()


# full model
mod.gg <- ggs(as.mcmc(full_mod))
jpeg("Model Contribution/Effect plots/MultiDirich_Mod1_betas_dotplot.jpg", units = "in", width = 4, height = 3, res = 1000)
ggs_caterpillar(D = mod.gg, family = "^beta", thick_ci = c(0.25, 0.75), thin_ci = c(0.025, 0.975), sort = FALSE) + 
  theme_bw() + ylab("") + xlab("Posterior estimate") + geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = c("Distance", "Area", "Groundwater", "Distance*Area", "Distance*Groundwater", "Area*Groundwater"))
dev.off()

####---------------------------------------------------#
#### Marginal effects plots #### 
####---------------------------------------------------#

# marginal effects plots
scltib <- read_csv("Model Contribution/GSI_Covariates_MeansStDevs.csv")

# vectors of predictors for N = 52 reporting groups
pdist <- seq(from = min(dist_cov), to = max(dist_cov), length.out = nrgs)
parea <- seq(from = min(area_cov), to = max(area_cov), length.out = nrgs)
pgw <- seq(from = min(gw_cov), to = max(gw_cov), length.out = nrgs)

# color palette
mypal <- hcl.colors(2, "Purple-Green")


#### Flowline Distance

# calculate the linear predictor
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta1"]*pdist) }

# transform to proportions
pred_dist_tran <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_tran)) { pred_dist_tran[i,] <-  pred_dist_lin[i,] / rowSums(pred_dist_lin)[i] }
rowSums(pred_dist_tran) # check that rows sum to 1

# calculate mean and 95% confidence intervals by integrating over MCMC draws
pred_lower <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.5)

# Plot
jpeg("Model contribution/Effect plots/MultiDirich_Mod7_MarginalEffects_Distance.jpg", units = "in", width = 4.5, height = 4.25, res = 1000)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
x.axis <- seq(from = 0, to = 150, by = 50)
x.scaled <- (x.axis - scltib$mean[3]) / scltib$sd[3]
plot(pred_median ~ pdist, ylim = c(0,0.5), pch = NA, bty = "l", axes = F, xlab = "Flowline distance (km)", ylab = "Proportional contribution")
axis(1, at = x.scaled, labels = x.axis)
axis(2)
box(bty = "l")
points(as.matrix(dat2[,-c(1:3)]) ~ as.matrix(dist_cov), col = "grey80")
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2)
dev.off()



#### Marginal Effects on Proportional Contribution ####
jpeg("Model contribution/Effect plots/MultiDirich_AdditiveEffects_MarginalEffects.jpg", units = "in", width = 6, height = 2.25, res = 1000)

uplim <- 0.06
par(mar = c(4,1.5,0.5,0), mgp = c(2.5,1,0), mfrow = c(1,3), oma = c(0,2,0,1))

# distance
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta1"]*pdist) }
pred_dist_tran <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_tran)) { pred_dist_tran[i,] <-  pred_dist_lin[i,] / rowSums(pred_dist_lin)[i] }
pred_lower <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.5)
x.axis <- seq(from = 0, to = 150, by = 50)
x.scaled <- (x.axis - scltib$mean[3]) / scltib$sd[3]
plot(pred_median ~ pdist, ylim = c(0,uplim), pch = NA, bty = "l", axes = F, xlab = "Flowline distance (km)", ylab = "")
axis(1, at = x.scaled, labels = x.axis)
axis(2)
box(bty = "l")
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2)
legend("topright", "(a)", bty = "n")
mtext("Proportional contribution", 2, 2.5, cex = 0.7)

# area
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta2"]*parea) }
pred_dist_tran <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_tran)) { pred_dist_tran[i,] <-  pred_dist_lin[i,] / rowSums(pred_dist_lin)[i] }
pred_lower <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.5)
x.axis <- c(1,10,100,1000)
x.scaled <- (log(x.axis) - scltib$mean[1]) / scltib$sd[1]
plot(pred_median ~ parea, ylim = c(0,uplim), pch = NA, bty = "l", axes = F, xlab = expression(paste("Catchment area (km"^"2", ")")), ylab = "")
axis(1, at = x.scaled, labels = x.axis)
axis(2, labels = NA)
box(bty = "l")
polygon(c(parea, rev(parea)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.3), lty = 0)
lines(pred_median ~ parea, lwd = 2)
legend("topright", "(b)", bty = "n")

# groundwater
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta3"]*pgw) }
pred_dist_tran <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_tran)) { pred_dist_tran[i,] <-  pred_dist_lin[i,] / rowSums(pred_dist_lin)[i] }
pred_lower <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.5)
x.axis <- c(0.025,0.05,0.1,0.2,0.4)
x.scaled <- (log(x.axis) - scltib$mean[2]) / scltib$sd[2]
plot(pred_median ~ pgw, ylim = c(0,uplim), pch = NA, bty = "l", axes = F, xlab = "Groundwater index", ylab = "")
axis(1, at = x.scaled, labels = x.axis)
axis(2, labels = NA)
box(bty = "l")
polygon(c(pgw, rev(pgw)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.3), lty = 0)
lines(pred_median ~ pgw, lwd = 2)
legend("topright", "(c)", bty = "n")

dev.off()



#-----------------#

#### Additive Effects on Dirichlet a ####
jpeg("Model contribution/Effect plots/MultiDirich_Mod8_AdditiveEffects_DirichletA.jpg", units = "in", width = 4, height = 4, res = 1000)

par(mar = c(3.5,3.5,0.5,0.5), mgp = c(2,0.8,0))

# min Area
x.axis <- seq(from = 0, to = 150, by = 25)
x.scaled <- (x.axis - scltib$mean[3]) / scltib$sd[3]
plot(pred_median ~ pdist, xlim = c(pdist[1], 0.84), ylim = c(0,20), pch = NA, bty = "l", axes = F, xlab = "Flowline distance (km)", ylab = "Dirichlet a")
axis(1, at = x.scaled, labels = x.axis)
axis(2)
box(bty = "l")
# min gw
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta1"]*pdist + Mcmcdat[i,"beta2"]*min(parea) + Mcmcdat[i,"beta3"]*min(pgw)) }
pred_lower <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.5)
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mypal[1], 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2, col = mypal[1], lty = 3)
# max gw
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta1"]*pdist + Mcmcdat[i,"beta2"]*min(parea) + Mcmcdat[i,"beta3"]*max(pgw)) }
pred_lower <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.5)
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mypal[2], 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2, col = mypal[2], lty = 3)

# Max Area
# min gw
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta1"]*pdist + Mcmcdat[i,"beta2"]*max(parea) + Mcmcdat[i,"beta3"]*min(pgw)) }
pred_lower <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.5)
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mypal[1], 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2, col = mypal[1])
# max gw
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta1"]*pdist + Mcmcdat[i,"beta2"]*max(parea) + Mcmcdat[i,"beta3"]*max(pgw)) }
pred_lower <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.5)
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mypal[2], 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2, col = mypal[2])

# legend
legend(x = -0.35, y = 21, legend = c("Min. gw", "Max. gw"), fill = c(scales::alpha(mypal, 0.5)), bty = "n")
legend(x = -0.5, y = 18, legend = c("Min. area", "Max. area"), lty = c(3,1), lwd = 2, bty = "n")

dev.off()


#-----------------#
# boxplot approach: compare test RG with variable conditions to RG 
pdist <- seq(from = min(dist_cov), to = max(dist_cov), length.out = nrgs)
parea <- seq(from = min(area_cov), to = max(area_cov), length.out = nrgs)
pgw <- seq(from = min(gw_cov), to = max(gw_cov), length.out = nrgs)

mygw <- seq(from = min(gw_cov), to = max(gw_cov), length.out = 100)
mygw <- sort(c(mygw, -1, 1))
mycols <- rev(viridis(length(mygw)))
collogw <- mycols[which(mygw == -1)] 
colhigw <- mycols[which(mygw == 1)] 
mypal <- hcl.colors(2, "Purple-Green")

# scenarios differ by min/max of each variable
scenarios <- tibble(scenario = c(1,3,5,7,2,4,6,8),
                    dist = c(min(dist_cov), max(dist_cov), min(dist_cov), max(dist_cov), min(dist_cov), max(dist_cov), min(dist_cov), max(dist_cov)),
                    area = c(min(area_cov), min(area_cov), max(area_cov), max(area_cov), min(area_cov), min(area_cov), max(area_cov), max(area_cov)),
                    gw = c(min(gw_cov), min(gw_cov), min(gw_cov), min(gw_cov), max(gw_cov), max(gw_cov), max(gw_cov), max(gw_cov))) %>% arrange(scenario)
# scenarios differ by +/- 1 StDev of each variable
scenarios <- tibble(scenario = c(1,3,5,7,2,4,6,8),
                    dist = c(-1, 1, -1, 1, -1, 1, -1, 1),
                    area = c(-1, -1, 1, 1, -1, -1, 1, 1),
                    gw   = c(-1, -1, -1, -1, 1, 1, 1, 1)) %>% arrange(scenario)
# Distance and area differ by +/- 1 StDev, but Groundwater held at min/max
scenarios <- tibble(scenario = c(1,3,5,7,2,4,6,8),
                    dist = c(-1, 1, -1, 1, -1, 1, -1, 1),
                    area = c(-1, -1, 1, 1, -1, -1, 1, 1),
                    gw   = c(min(gw_cov), min(gw_cov), min(gw_cov), min(gw_cov), max(gw_cov), max(gw_cov), max(gw_cov), max(gw_cov))) %>% arrange(scenario)

pred_abs <- matrix(NA, nrow = nrow(Mcmcdat), ncol = dim(scenarios)[1])
pred_abs_ref <- matrix(NA, nrow = nrow(Mcmcdat), ncol = dim(scenarios)[1])
pred_diff <- matrix(NA, nrow = nrow(Mcmcdat), ncol = dim(scenarios)[1])
pred_diff_rel <- matrix(NA, nrow = nrow(Mcmcdat), ncol = dim(scenarios)[1])
for (j in 1:dim(scenarios)[1]) {
  pred_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
  for (i in 1:nrow(pred_lin)) { pred_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta1"]*c(scenarios$dist[j], rep(0, nrgs-1)) + Mcmcdat[i,"beta2"]*c(scenarios$area[j], rep(0, nrgs-1)) + Mcmcdat[i,"beta3"]*c(scenarios$gw[j], rep(0, nrgs-1))) }
  pred_tran <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
  for (i in 1:nrow(pred_tran)) { pred_tran[i,] <-  pred_lin[i,] / rowSums(pred_lin)[i] }
  pred_abs[,j] <- pred_tran[,1]
  pred_abs_ref[,j] <- pred_tran[,2]
  pred_diff[,j] <- pred_tran[,1] - pred_tran[,2]
  pred_diff_rel[,j] <- (pred_tran[,1] - pred_tran[,2]) / pred_tran[,2]
  print(j)
}

# paired t-tests to test for differences among groundwater treatments
t.test(x = pred_diff[,1], y = pred_diff[,2], alternative = "two.sided", paired = T) # p << 0.05, mean diff = -0.075
t.test(x = pred_diff[,3], y = pred_diff[,4], alternative = "two.sided", paired = T) # p << 0.05, mean diff = -0.00065
t.test(x = pred_diff[,5], y = pred_diff[,6], alternative = "two.sided", paired = T) # p << 0.05, mean diff = -0.151
t.test(x = pred_diff[,7], y = pred_diff[,8], alternative = "two.sided", paired = T) # p << 0.05, mean diff = -0.0023

# PLOT --- w/ broken axis
jpeg("Model contribution/Effect plots/MultiDirich_AdditiveEffects_Boxplots.jpg", units = "in", width = 4, height = 5, res = 1000)
par(mar = c(5, 4, 1, 1), mgp = c(3,2,0))
# postive differences
boxplot(pred_diff, col = rep(c(scales::alpha(mycols[c(1,102)], 0.7)), 4), axes = F, ylim = c(-0.02, 0.7))
axis(1, at = c(1.5, 3.5, 5.5, 7.5), labels = c("Near\nSmall", "Far\nSmall", "Near\nLarge", "Far\nLarge"))
par(mgp = c(2.5,1,0))
axis(2)
box()
# gridlines and legend
abline(v = c(2.5, 4.5, 6.5))
abline(h = 0, lty = 2)
legend("topright", legend = c("Min. GW", "Max. GW"), fill = scales::alpha(mycols[c(1,102)], 0.7), bty = "n", cex = 0.7)
# show t-test results
text(1, max(pred_diff[,1])+0.03, "*")
text(2, max(pred_diff[,2])+0.03, "**")
text(3, max(pred_diff[,3])+0.03, "*")
text(4, max(pred_diff[,4])+0.03, "**")
text(5, max(pred_diff[,5])+0.03, "*")
text(6, max(pred_diff[,6])+0.03, "**")
text(7, max(pred_diff[,7])+0.03, "*")
text(8, max(pred_diff[,8])+0.03, "**")
# axis labels
mtext("Tributary conditions", 1, 3.5)
mtext("Difference from average contribution", 2, 2.5)
dev.off()



# PLOT --- w/ broken axis
jpeg("Model contribution/Effect plots/MultiDirich_AdditiveEffects_Boxplots_BrokenAxis.jpg", units = "in", width = 4, height = 5, res = 1000)
layout(matrix(c(1, 1, 1, 1, 2, 2, 2), nrow = 7, ncol = 1))
par(mar = c(0.1, 6, 1, 1), oma = c(2,2,0,0), mgp = c(2.5,2,0))
# postive differences
boxplot(pred_diff, col = rep(c(scales::alpha(mycols[c(1,102)], 0.7)), 4), xaxt = "n", xaxs = "i", yaxs = "i", ylim = c(-0.02,0.75), las = 2)
abline(v = c(2.5, 4.5, 6.5))
legend("topright", legend = c("Low GW", "High GW"), fill = scales::alpha(mycols[c(1,102)], 0.7), bty = "n")
mtext("Difference from average contribution", 2, 5, at = 0.15)
# show t-test results
text(1, max(pred_diff[,1])+0.05, "*")
text(2, max(pred_diff[,2])+0.05, "**")
text(5, max(pred_diff[,5])+0.05, "*")
text(6, max(pred_diff[,6])+0.05, "**")
# negative differences
par(mar = c(4, 6, 0.1, 1))
boxplot(pred_diff, col = rep(c(scales::alpha(mycols[c(1,102)], 0.7)), 4), axes = F, xaxs = "i", yaxs = "i", ylim = c(-0.02,0), las = 2)
axis(2, at = seq(from = 0, to = -0.020, by = -0.005), labels = c(NA, "-0.005", "-0.010", "-0.015", "-0.020"), las = 2)
axis(1, at = c(1.5, 3.5, 5.5, 7.5), labels = c("Near\nSmall", "Far\nSmall", "Near\nLarge", "Far\nLarge"))
box()
abline(v = c(2.5, 4.5, 6.5))
# show t-test results
text(3, max(pred_diff[,3])+0.0025, "*")
text(4, max(pred_diff[,4])+0.0025, "**")
text(7, max(pred_diff[,7])+0.0025, "*")
text(8, max(pred_diff[,8])+0.0025, "**")
mtext("Tributary conditions", 1, 4)

dev.off()


# Some other plots...relative differences, absolute values, and absolute values of the reference groups
# just the negatives
boxplot(pred_diff[,c(3,4,7,8)], col = rep(c(collogw, colhigw), 4))
abline(h = 0, lty = 2)

# relative difference
boxplot(pred_diff_rel, col = rep(c(collogw, colhigw), 4))
abline(h = 0, lty = 2)

# absolute contributions
boxplot(pred_abs, col = rep(c(collogw, colhigw), 4))
abline(h = 0, lty = 2)

# absolute contributions of the references group(s)
boxplot(pred_abs_ref, col = rep(c(collogw, colhigw), 4))
abline(h = 0, lty = 2)

dev.off()



#-----------------------------#
# This is the additive effects plot, but definite not a clean/effective of a visual

# colors
mycols <- rev(viridis(2))

# number of reporting groups
nrgss <- 200

# common variables
pdist <- seq(from = min(dist_cov), to = max(dist_cov), length.out = nrgss/4)
parmin <- seq(from = min(area_cov), to = min(area_cov), length.out = nrgss/4)
parmax <- seq(from = max(area_cov), to = max(area_cov), length.out = nrgss/4)
pgwmin <- seq(from = min(gw_cov), to = min(gw_cov), length.out = nrgss/4)
pgwmax <- seq(from = max(gw_cov), to = max(gw_cov), length.out = nrgss/4)

parmin <- seq(from = -1, to = -1, length.out = nrgss/4)
parmax <- seq(from = 1, to = 1, length.out = nrgss/4)
pgwmin <- seq(from = min(gw_cov), to = min(gw_cov), length.out = nrgss/4)
pgwmax <- seq(from = max(gw_cov), to = max(gw_cov), length.out = nrgss/4)

# set up the plotting environment
x.axis <- seq(from = 0, to = 150, by = 50)
x.scaled <- (x.axis - scltib$mean[3]) / scltib$sd[3]
plot(seq(from = 0, to = 0.2, length.out = nrgss/4) ~ pdist, xlim = c(min(pdist), 0.84), ylim = c(0,0.18), pch = NA, bty = "l", xaxt = "n", xlab = "Flowline distance (km)", ylab = "Proportional contribution")
axis(1, at = x.scaled, labels = x.axis)

# linear predictor
pred_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgss)
for (i in 1:nrow(pred_lin)) { pred_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta1"]*c(pdist, pdist, pdist, pdist) + Mcmcdat[i,"beta2"]*c(parmin, parmin, parmax, parmax) + Mcmcdat[i,"beta3"]*c(pgwmin, pgwmax, pgwmin, pgwmax)) }
# transform to proportions
pred_tran <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgss)
for (i in 1:nrow(pred_tran)) { pred_tran[i,] <-  pred_lin[i,] / rowSums(pred_lin)[i] }

# SMALL BASINS
# low groundwater conditions
pred_tran_sub <- pred_tran[,c(1:50)]
pred_lower <- apply(pred_tran_sub, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_tran_sub, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_tran_sub, MARGIN = 2, quantile, prob = 0.5)
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mycols[1], 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2, col = mycols[1], lty = 2)
# high groundwater conditions
pred_tran_sub <- pred_tran[,c(51:100)]
pred_lower <- apply(pred_tran_sub, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_tran_sub, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_tran_sub, MARGIN = 2, quantile, prob = 0.5)
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mycols[2], 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2, col = mycols[2], lty = 2)

# LARGE BASINS
# low groundwater conditions
pred_tran_sub <- pred_tran[,c(101:150)]
pred_lower <- apply(pred_tran_sub, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_tran_sub, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_tran_sub, MARGIN = 2, quantile, prob = 0.5)
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mycols[1], 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2, col = mycols[1])
# high groundwater conditions
pred_tran_sub <- pred_tran[,c(151:200)]
pred_lower <- apply(pred_tran_sub, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_tran_sub, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_tran_sub, MARGIN = 2, quantile, prob = 0.5)
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mycols[2], 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2, col = mycols[2])




pred_lower <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.5)

axis(1, at = x.scaled, labels = x.axis)
axis(2)
box(bty = "l")
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha("black", 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2)

# min Area
x.axis <- seq(from = 0, to = 150, by = 25)
x.scaled <- (x.axis - scltib$mean[3]) / scltib$sd[3]
plot(pred_median ~ pdist, xlim = c(pdist[1], 0.84), ylim = c(0,20), pch = NA, bty = "l", axes = F, xlab = "Flowline distance (km)", ylab = "")
axis(1, at = x.scaled, labels = x.axis)
axis(2)
box(bty = "l")
# min gw
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta1"]*pdist + Mcmcdat[i,"beta2"]*min(parea) + Mcmcdat[i,"beta3"]*min(pgw)) }
pred_lower <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.5)
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mypal[1], 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2, col = mypal[1], lty = 3)
# max gw
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta1"]*pdist + Mcmcdat[i,"beta2"]*min(parea) + Mcmcdat[i,"beta3"]*max(pgw)) }
pred_lower <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_lin, MARGIN = 2, quantile, prob = 0.5)
polygon(c(pdist, rev(pdist)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mypal[2], 0.3), lty = 0)
lines(pred_median ~ pdist, lwd = 2, col = mypal[2], lty = 3)











#### Area x Groundwater interaction (Area focus)

jpeg("Model contribution/Effect plots/MultiDirich_Mod7_MarginalEffects_Area.jpg", units = "in", width = 4.5, height = 4.25, res = 1000)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
x.axis <- c(1,10,100,1000)
x.scaled <- (log(x.axis) - scltib$mean[1]) / scltib$sd[1]
plot(pred_median ~ parea, ylim = c(0,0.17), pch = NA, bty = "l", axes = F, xlab = expression(paste("Catchment area (km"^"2", ")")), ylab = "Proportional contribution")
axis(1, at = x.scaled, labels = x.axis)
axis(2)
box(bty = "l")

# min Groundwater
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta2"]*parea + Mcmcdat[i,"beta3"]*min(pgw) + Mcmcdat[i,"beta6"]*parea*min(pgw)) }
pred_dist_tran <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_tran)) { pred_dist_tran[i,] <-  pred_dist_lin[i,] / rowSums(pred_dist_lin)[i] }
pred_lower <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.5)
polygon(c(parea, rev(parea)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mypal[1], 0.3), lty = 0)
lines(pred_median ~ parea, lwd = 2, col = mypal[1])

# max Groundwater
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta2"]*parea + Mcmcdat[i,"beta3"]*max(pgw) + Mcmcdat[i,"beta6"]*parea*max(pgw)) }
pred_dist_tran <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_tran)) { pred_dist_tran[i,] <-  pred_dist_lin[i,] / rowSums(pred_dist_lin)[i] }
pred_lower <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.5)
polygon(c(parea, rev(parea)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mypal[10], 0.3), lty = 0)
lines(pred_median ~ parea, lwd = 2, col = mypal[10])

legend("topleft", legend = c("Groundwater", "Snowmelt"), fill = c(scales::alpha(mypal[10], 0.6), scales::alpha(mypal[1], 0.6)), bty = "n")

dev.off()



#### Area x Groundwater interaction (Groundwater focus)

jpeg("Model contribution/Effect plots/MultiDirich_Mod7_MarginalEffects_Groundwater.jpg", units = "in", width = 4.5, height = 4.25, res = 1000)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
x.axis <- c(0.025,0.05,0.1,0.2,0.4)
x.scaled <- (log(x.axis) - scltib$mean[2]) / scltib$sd[2]
plot(pred_median ~ pgw, ylim = c(0,0.17), pch = NA, bty = "l", axes = F, xlab = "Groundwater index", ylab = "Proportional contribution")
axis(1, at = x.scaled, labels = x.axis)
axis(2)
box(bty = "l")

# min Area
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta2"]*min(parea) + Mcmcdat[i,"beta3"]*pgw + Mcmcdat[i,"beta6"]*min(parea)*pgw) }
pred_dist_tran <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_tran)) { pred_dist_tran[i,] <-  pred_dist_lin[i,] / rowSums(pred_dist_lin)[i] }
pred_lower <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.5)
polygon(c(pgw, rev(pgw)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mypal[1], 0.3), lty = 0)
lines(pred_median ~ pgw, lwd = 2, col = mypal[1])

# max Groundwater
pred_dist_lin <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_lin)) { pred_dist_lin[i,] <- exp(Mcmcdat[i,"alpha"] + Mcmcdat[i,"beta2"]*max(parea) + Mcmcdat[i,"beta3"]*pgw + Mcmcdat[i,"beta6"]*max(parea)*pgw) }
pred_dist_tran <- matrix(NA, nrow = nrow(Mcmcdat), ncol = nrgs)
for (i in 1:nrow(pred_dist_tran)) { pred_dist_tran[i,] <-  pred_dist_lin[i,] / rowSums(pred_dist_lin)[i] }
pred_lower <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.025)
pred_upper <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.975)
pred_median <- apply(pred_dist_tran, MARGIN = 2, quantile, prob = 0.5)
polygon(c(pgw, rev(pgw)), c(pred_upper, rev(pred_lower)), col = scales::alpha(mypal[10], 0.3), lty = 0)
lines(pred_median ~ pgw, lwd = 2, col = mypal[10])

legend("topleft", legend = c("Max. area", "Min. area"), fill = c(scales::alpha(mypal[10], 0.6), scales::alpha(mypal[1], 0.6)), bty = "n")

dev.off()



# 
Y0 <- gsi_dat$bootstrapped_proportions %>%
  mutate(section = str_sub(mixture_collection, end = -6), year = str_sub(mixture_collection, start = -4)) %>%
  arrange(mixture_collection, repunit) %>% select(mixture_collection, section, year, repunit, bs_corrected_repunit_ppn) %>% 
  spread(key = repunit, value = bs_corrected_repunit_ppn)
Y1 <- Y0 %>% select(-c(mixture_collection, section, year))
Y1[Y1 < 0.05] <- 0

distinv <- 1/((dist_cov*scltib$sd[3])+scltib$mean[3])

rich <- specnumber(Y1) # richness
H <- diversity(Y1) # Shannon's H
evns <- H / log(rich)
gwwt <- c() # groundwater weighted mean
for (i in 1:17) { gwwt[i] <- weighted.mean(gwmet$gwi_iew05km, w = distinv[i,]) }

plot(rich ~ gwwt, type = "n")
text(gwwt, rich, labels = sectid)

plot(H ~ gwwt, type = "n")
text(gwwt, H, labels = sectid)

plot(evns ~ gwwt, type = "n")
text(gwwt, evns, labels = sectid)

plot(H ~ gwwt)
plot(evns ~ gwwt)










###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################

# simple diagnostic tools
rowSums(mod_multi$BUGSoutput$mean$p) # check to make sure estimated proportions sum to 1
plot(as.numeric(mod_multi$BUGSoutput$mean$p) ~ as.numeric(unlist(dat2[,-c(1:3)]))) # how do estimated proportions compare to bootstrap-corrected mixing proportions
sum(as.numeric(mod_multi$BUGSoutput$mean$p) > as.numeric(unlist(dat2[,-c(1:3)]))) / length(as.numeric(mod_multi$BUGSoutput$mean$p))
abline(a = 0, b = 1, col = "red")

plot(as.numeric(mod_multi$BUGSoutput$mean$res) ~ as.numeric(unlist(dist_cov)))




#### 2. Multinomial-Dirichlet Regression using bootstrapped fish abundance ####
#### this estimates an intensity paramater (functionally equivalent to overdispersion?) sensu Harrison et al (2020 Mol Ecol Res) and (???) Douma and Weedon (2019)

# specify data and parameters
jags.data <- list("nmix" = nmix, "nrgs" = nrgs, "nsect" = nsect, "sect" = sect, "Y" = Yint, "N" = Nint, "dist" = dist_cov, "area" = area_cov, "gw" = gw_cov)      # specify data
params <- c("p", "alpha", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "phi", "loglik")  # parameters to monitor
# Run the model...slow!
mod_dm <- jags(jags.data, inits = NULL, parameters.to.save = params, model.file = "Model Contribution/MultinomialDirichletContributions_covars_intensity.txt", 
                  n.chains = 3, n.thin = 10, n.burnin = 1000, n.iter = 3000, DIC = TRUE)
MCMCtrace(mod_dm, ind = TRUE, filename = "Model Contribution/MultinomialDirichletContribution_Intensity_Traceplots.pdf")
# simple diagnostic tools
rowSums(mod_dm$BUGSoutput$mean$p) # check to make sure estimated proportions sum to 1
plot(as.numeric(mod_dm$BUGSoutput$mean$p) ~ as.numeric(unlist(dat2[,-c(1:3)]))) # how do estimated proportions compare to bootstrap-corrected mixing proportions
abline(a = 0, b = 1, col = "red")




#### Dirichlet Regression using bootstrapped proportions ####

# specify data and parameters
jags.data <- list("nmix" = nmix, "nrgs" = nrgs, "Y" = Y, "dist" = dist_cov, "area" = area_cov, "gw" = gw_cov)      # specify data
params <- c("p", "alpha", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "sigma")  # parameters to monitor
# Run the model...
mod_dirich <- jags(jags.data, inits = NULL, parameters.to.save = params, model.file = "Model Contribution/DirichletContributions_covars.txt", 
                   n.chains = 3, n.thin = 10, n.burnin = 1000, n.iter = 4000, DIC = TRUE)
MCMCtrace(mod_dirich, ind = TRUE, filename = "Model Contribution/DirichletContribution_Traceplots.pdf")
# simple diagnostic tools
rowSums(mod_dirich$BUGSoutput$mean$p)
plot(as.numeric(mod_dirich$BUGSoutput$mean$p) ~ as.numeric(unlist(dat2[,-c(1:3)])), xlim = c(0, 0.5), ylim = c(0, 0.5)) 
abline(a = 0, b = 1, col = "red")
# this does a very poor job at estimating proportions



#### simple plots

mysum <- mod_multi$BUGSoutput$summary[c(1:7),]
area_new <- seq(from = min(dat$z_logarea), to = max(dat$z_logarea), length.out = 100)

par(mfrow = c(2,2))
# Min GW, Min Dist
pred <- exp(mysum[1,1] + mysum[2,1]*min(dat$z_dist) + mysum[3,1]*area_new + mysum[4,1]*min(dat$z_loggwi) + mysum[5,1]*min(dat$z_dist)*area_new + mysum[6,1]*min(dat$z_dist)*min(dat$z_loggwi) + mysum[7,1]*min(dat$z_logarea)*area_new )
plot(pred ~ area_new, type = "l", ylim = c(0,3))

# Max GW, Min Dist
pred <- exp(mysum[1,1] + mysum[2,1]*min(dat$z_dist) + mysum[3,1]*area_new + mysum[4,1]*max(dat$z_loggwi) + mysum[5,1]*min(dat$z_dist)*area_new + mysum[6,1]*min(dat$z_dist)*max(dat$z_loggwi) + mysum[7,1]*min(dat$z_logarea)*area_new )
plot(pred ~ area_new, type = "l", ylim = c(0,3))

# Min GW, Max Dist
pred <- exp(mysum[1,1] + mysum[2,1]*max(dat$z_dist) + mysum[3,1]*area_new + mysum[4,1]*min(dat$z_loggwi) + mysum[5,1]*max(dat$z_dist)*area_new + mysum[6,1]*max(dat$z_dist)*min(dat$z_loggwi) + mysum[7,1]*min(dat$z_logarea)*area_new )
plot(pred ~ area_new, type = "l", ylim = c(0,3))

# Max GW, Max Dist
pred <- exp(mysum[1,1] + mysum[2,1]*max(dat$z_dist) + mysum[3,1]*area_new + mysum[4,1]*max(dat$z_loggwi) + mysum[5,1]*max(dat$z_dist)*area_new + mysum[6,1]*max(dat$z_dist)*max(dat$z_loggwi) + mysum[7,1]*max(dat$z_logarea)*area_new )
plot(pred ~ area_new, type = "l", ylim = c(0,3))



dist_new <- seq(from = min(dat$z_dist), to = min(dat$z_dist), length.out = 100)
gw_new <- seq(from = mean(dat$z_loggwi), to = mean(dat$z_loggwi), length.out = 100)

pred1 <- 
pred2 <- 

plot(pred1 ~ newdat, type = "l", lty = 1, ylim = c(0,1))
lines(pred2 ~ newdat, type = "l", lty = 2)

plot(preds ~ c((newdat*scltib$sd[3]) + scltib$mean[3]), type = "l")





# 
V <- compositions::ilrBase(D = nrgs)
#generalized inverse of transpose of V, needed to do inverse of ilr in stan
#Argument: V, used in ilr transformation
#Value: generalized inverse of transpose of V
gettVinv <- function(V){
  return(MASS::ginv(t(V)))
}
tVinv <- gettVinv(V)

transect_dat <- list(nrgs = nrgs, 
                     nmix = nmix,
                     Yint = Yint,
                     area_cov = area_cov,
                     tVinv = tVinv)
niter <- 1e3 #number of stan iterations: used 1e4 in paper, but 1e3 to demonstrate code
fit <- stan(file = "./Model Contribution/multimod.stan", data = transect_dat, 
            iter = niter, chains = 1) 

plot(rstan::traceplot(fit, pars = c("beta0", "beta1", "beta2", "LSigma"), inc_warmup = TRUE))

print(summary(fit, pars = c("beta0", "beta1", "beta2", "LSigma"))$summary)




















####-----------------------------------------#
#### Create functions ####
####-----------------------------------------#

# logit function
logit <- function(x) { log(x / (1-x)) } 

# inverse logit function
invlogit <- function(x) { exp(x) / (1 + exp(x)) } 


p <- seq(0.001, 0.999, 0.001)
pLogit <- car::logit(p)
plot(p, pLogit, type='l', lwd=2, col='red', las=1, xlab='p', ylab='logit(p)')
plot(logit(p)~car::logit(p))

####-----------------------------------------#
#### Data visualization ####
####-----------------------------------------#

# Raw proportions ~ covariats
# area
plot(ppn_bscor ~ areasqkm, dat)
plot(ppn_bscor ~ logarea, dat)
# flowline distance
plot(ppn_bscor ~ distkm, dat)
# groundwater
plot(ppn_bscor ~ gwi_iew05km, dat)
plot(ppn_bscor ~ loggwi, dat)


# plots with logit transform
plot(logit(dat$ppn_bscor + 0.00001) ~ dat$logarea)
plot(logit(dat$ppn_bscor + 0.00001) ~ dat$distkm)
plot(logit(dat$ppn_bscor + 0.00001) ~ dat$loggwi)
# correlation
cor(logit(dat$ppn_bscor + 0.00001), dat$logarea)
cor(logit(dat$ppn_bscor + 0.00001), dat$distkm)
cor(logit(dat$ppn_bscor + 0.00001), dat$loggwi)


####-----------------------------------------#
#### Simple models ####
####-----------------------------------------#

# mod <- glm(ppn_bscor ~ logarea + distkm + loggwi + logarea*distkm + logarea*loggwi + distkm*loggwi , data = dat, family = "binomial", weight = total)
# plot(mod)
# summary(mod)

dat$logitppn <- car::logit(dat$ppn_bscor)
dat$ppndum <- dat$ppn_bscor + 0.00001
dat$num_bscor <- dat$ppn_bscor * dat$total

mod <- lm(logitppn ~ z_logarea + z_dist + z_loggwi, data = dat, weight = total)
plot(mod)
summary(mod)

mod <- glm(ppn_bscor ~ z_logarea + z_dist + z_loggwi, data = dat, weight = total, family = "quasibinomial")
plot(mod)
summary(mod)


bmod <- betareg(ppndum ~ z_logarea + z_dist + z_loggwi, data = dat, link = "logit", weight = total)
plot(bmod)
summary(bmod)

zbmod <- zoib(ppn_bscor ~ z_dist + z_loggwi + z_logarea*z_dist + z_logarea*z_loggwi + z_dist*z_loggwi, data = dat)

q = qt(0.975, df = df.residual(mod))

# marginal effect - logarea
newdata <- with(dat, data.frame(z_logarea = seq(min(z_logarea), max(z_logarea), len = 100),
                                z_loggwi = seq(mean(z_loggwi), mean(z_loggwi), len = 100),
                                z_dist = seq(mean(z_dist), mean(z_dist), len = 100)))
fit <- predict(mod, newdata = newdata, type = "response", interval = "confidence")
# newdata <- cbind(newdata, fit = binomial()$linkinv(fit$fit), lower = binomial()$linkinv(fit$fit - q * fit$se.fit), upper = binomial()$linkinv(fit$fit + q * fit$se.fit))
jpeg("Model contribution/Effect plots/Cont_logarea_marginal.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(invlogit(fit[,1]) ~ newdata$z_logarea, type = "n", ylim = c(0,0.12), bty = "l", xlab = "log(Basin aream, km^2)", ylab = "Proportional contribution")
polygon(c(newdata$z_logarea, rev(newdata$z_logarea)), c(invlogit(fit[,3]), rev(invlogit(fit[,2]))), col = scales::alpha("blue", 0.4), lty = 0)
lines(invlogit(fit[,1]) ~ newdata$z_logarea, lwd = 2)
dev.off()

# marginal effect - loggwi
newdata <- with(dat, data.frame(z_logarea = seq(min(dat$z_logarea), max(dat$z_logarea), len = 100),
                                z_loggwi = seq(0, 0, len = 100),
                                z_dist = seq(0, 0, len = 100)))
fit <- predict(mod, newdata = newdata, type = "link", se = TRUE)
newdata <- cbind(newdata, fit = binomial()$linkinv(fit$fit), lower = binomial()$linkinv(fit$fit - q * fit$se.fit), upper = binomial()$linkinv(fit$fit + q * fit$se.fit))
jpeg("Model contribution/Effect plots/Cont_loggwi_marginal.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(fit ~ z_logarea, newdata, type = "n", ylim = c(0,0.12), bty = "l", xlab = "log(Area)", ylab = "Proportional contribution")
polygon(c(newdata$z_logarea, rev(newdata$z_logarea)), c(newdata$upper, rev(newdata$lower)), col = scales::alpha("blue", 0.4), lty = 0)
lines(newdata$fit ~ newdata$z_logarea, lwd = 2)

newdata <- with(dat, data.frame(z_dist = seq(min(dat$z_dist), max(dat$z_dist), len = 100),
                                z_loggwi = seq(0, 0, len = 100),
                                z_logarea = seq(0, 0, len = 100)))
fit <- predict(mod, newdata = newdata, type = "link", se = TRUE)
newdata <- cbind(newdata, fit = binomial()$linkinv(fit$fit), lower = binomial()$linkinv(fit$fit - q * fit$se.fit), upper = binomial()$linkinv(fit$fit + q * fit$se.fit))
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(fit ~ z_dist, newdata, type = "n", ylim = c(0,0.12), bty = "l", xlab = "log(Area)", ylab = "Proportional contribution")
polygon(c(newdata$z_dist, rev(newdata$z_dist)), c(newdata$upper, rev(newdata$lower)), col = scales::alpha("blue", 0.4), lty = 0)
lines(newdata$fit ~ newdata$z_dist, lwd = 2)

dev.off()

newdata <- with(dat, data.frame(z_logarea = seq(mean(z_logarea), mean(z_logarea), len = 100),
                                z_loggwi = seq(min(z_loggwi), max(z_loggwi), len = 100),
                                z_dist = seq(mean(z_dist), mean(z_dist), len = 100)))
fit <- predict(mod, newdata = newdata, type = "response", interval = "confidence")
# newdata <- cbind(newdata, fit = binomial()$linkinv(fit$fit), lower = binomial()$linkinv(fit$fit - q * fit$se.fit), upper = binomial()$linkinv(fit$fit + q * fit$se.fit))
jpeg("Model contribution/Effect plots/Cont_logarea_marginal.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(invlogit(fit[,1]) ~ newdata$z_loggwi, type = "n", ylim = c(0,0.12), bty = "l", xlab = "log(Groundwater index)", ylab = "Proportional contribution")
polygon(c(newdata$z_loggwi, rev(newdata$z_loggwi)), c(invlogit(fit[,3]), rev(invlogit(fit[,2]))), col = scales::alpha("blue", 0.4), lty = 0)
lines(invlogit(fit[,1]) ~ newdata$z_loggwi, lwd = 2)
dev.off()



# marginal effect - distkm
newdata <- with(dat, data.frame(logarea = seq(mean(logarea), mean(logarea), len = 100),
                                loggwi = seq(mean(loggwi), mean(loggwi), len = 100),
                                distkm = seq(min(distkm), max(distkm), len = 100)))
fit <- predict(mod, newdata = newdata, type = "link", se = TRUE)
newdata <- cbind(newdata, fit = binomial()$linkinv(fit$fit), lower = binomial()$linkinv(fit$fit - q * fit$se.fit), upper = binomial()$linkinv(fit$fit + q * fit$se.fit))
jpeg("Model contribution/Effect plots/Cont_dist_marginal.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(fit ~ distkm, newdata, type = "n", ylim = c(0,0.12), bty = "l", xlab = "Flowline distance (km)", ylab = "Proportional contribution")
polygon(c(newdata$distkm, rev(newdata$distkm)), c(newdata$upper, rev(newdata$lower)), col = scales::alpha("blue", 0.4), lty = 0)
lines(newdata$fit ~ newdata$distkm, lwd = 2)
dev.off()

newdata <- with(dat, data.frame(z_logarea = seq(mean(z_logarea), mean(z_logarea), len = 100),
                                z_loggwi = seq(mean(z_loggwi), mean(z_loggwi), len = 100),
                                z_dist = seq(min(z_dist), max(z_dist), len = 100)))
fit <- predict(mod, newdata = newdata, type = "response", interval = "confidence")
# newdata <- cbind(newdata, fit = binomial()$linkinv(fit$fit), lower = binomial()$linkinv(fit$fit - q * fit$se.fit), upper = binomial()$linkinv(fit$fit + q * fit$se.fit))
jpeg("Model contribution/Effect plots/Cont_logarea_marginal.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
plot(invlogit(fit[,1]) ~ newdata$z_dist, type = "n", ylim = c(0,0.12), bty = "l", xlab = "Flowline distance (km)", ylab = "Proportional contribution")
polygon(c(newdata$z_dist, rev(newdata$z_dist)), c(invlogit(fit[,3]), rev(invlogit(fit[,2]))), col = scales::alpha("blue", 0.4), lty = 0)
lines(invlogit(fit[,1]) ~ newdata$z_dist, lwd = 2)
dev.off()


# 2-way interaction: gw * dist
jpeg("Model contribution/Effect plots/Cont_dist_gw_interaction.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
# minimum gw
# newdata <- with(dat, data.frame(logarea = seq(mean(logarea), mean(logarea), len = 100),
#                                 loggwi = seq(min(loggwi), min(loggwi), len = 100),
#                                 distkm = seq(min(distkm), max(distkm), len = 100)))
# fit <- predict(mod, newdata = newdata, type = "link", se = TRUE)
# newdata <- cbind(newdata, fit = binomial()$linkinv(fit$fit), lower = binomial()$linkinv(fit$fit - q * fit$se.fit), upper = binomial()$linkinv(fit$fit + q * fit$se.fit))
# plot(fit ~ distkm, newdata, type = "n", ylim = c(0,0.3), bty = "l", xlab = "Flowline distance (km)", ylab = "Proportional contribution")
# polygon(c(newdata$distkm, rev(newdata$distkm)), c(newdata$upper, rev(newdata$lower)), col = scales::alpha("blue", 0.4), lty = 0)
# lines(newdata$fit ~ newdata$distkm, lwd = 2, col = "blue")
newdata <- with(dat, data.frame(z_logarea = seq(mean(z_logarea), mean(z_logarea), len = 100),
                                z_loggwi = seq(min(z_loggwi), min(z_loggwi), len = 100),
                                z_dist = seq(min(z_dist), max(z_dist), len = 100)))
fit <- predict(bmod, newdata = newdata, type = "response")
plot(seq(0,1,length.out = 100) ~ newdata$z_dist, type = "n", ylim = c(0,0.12), bty = "l", xlab = "Flowline distance (km)", ylab = "Proportional contribution")
# polygon(c(newdata$z_dist, rev(newdata$z_dist)), c(invlogit(fit[,3]), rev(invlogit(fit[,2]))), col = scales::alpha("blue", 0.4), lty = 0)
lines(fit ~ newdata$z_dist, lwd = 2, col = "blue")
# maximum gw
# newdata <- with(dat, data.frame(logarea = seq(mean(logarea), mean(logarea), len = 100),
#                                 loggwi = seq(max(loggwi), max(loggwi), len = 100),
#                                 distkm = seq(min(distkm), max(distkm), len = 100)))
# fit <- predict(mod, newdata = newdata, type = "link", se = TRUE)
# newdata <- cbind(newdata, fit = binomial()$linkinv(fit$fit), lower = binomial()$linkinv(fit$fit - q * fit$se.fit), upper = binomial()$linkinv(fit$fit + q * fit$se.fit))
# polygon(c(newdata$distkm, rev(newdata$distkm)), c(newdata$upper, rev(newdata$lower)), col = scales::alpha("red", 0.4), lty = 0)
# lines(newdata$fit ~ newdata$distkm, lwd = 2, col = "red")
newdata <- with(dat, data.frame(z_logarea = seq(mean(z_logarea), mean(z_logarea), len = 100),
                                z_loggwi = seq(max(z_loggwi), max(z_loggwi), len = 100),
                                z_dist = seq(min(z_dist), max(z_dist), len = 100)))
fit <- predict(bmod, newdata = newdata, type = "response", interval = "confidence")
polygon(c(newdata$z_dist, rev(newdata$z_dist)), c(invlogit(fit[,3]), rev(invlogit(fit[,2]))), col = scales::alpha("red", 0.4), lty = 0)
lines(fit ~ newdata$z_dist, lwd = 2, col = "red")
# legend
legend("topright", legend = c("Min. GW", "Max. GW"), bty = "n", fill = c("blue", "red"))
dev.off()


# 2-way interaction: gw * dist
jpeg("Model contribution/Effect plots/Cont_area_gw_interaction.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
# minimum gw
newdata <- with(dat, data.frame(logarea = seq(min(logarea), max(logarea), len = 100),
                                loggwi = seq(min(loggwi), min(loggwi), len = 100),
                                distkm = seq(mean(distkm), mean(distkm), len = 100)))
fit <- predict(mod, newdata = newdata, type = "link", se = TRUE)
newdata <- cbind(newdata, fit = binomial()$linkinv(fit$fit), lower = binomial()$linkinv(fit$fit - q * fit$se.fit), upper = binomial()$linkinv(fit$fit + q * fit$se.fit))
plot(fit ~ logarea, newdata, type = "n", ylim = c(0,0.3), bty = "l", xlab = "log(Basin area, km^2)", ylab = "Proportional contribution")
polygon(c(newdata$logarea, rev(newdata$logarea)), c(newdata$upper, rev(newdata$lower)), col = scales::alpha("blue", 0.4), lty = 0)
lines(newdata$fit ~ newdata$logarea, lwd = 2, col = "blue")
# maximum gw
newdata <- with(dat, data.frame(logarea = seq(min(logarea), max(logarea), len = 100),
                                loggwi = seq(max(loggwi), max(loggwi), len = 100),
                                distkm = seq(mean(distkm), mean(distkm), len = 100)))
fit <- predict(mod, newdata = newdata, type = "link", se = TRUE)
newdata <- cbind(newdata, fit = binomial()$linkinv(fit$fit), lower = binomial()$linkinv(fit$fit - q * fit$se.fit), upper = binomial()$linkinv(fit$fit + q * fit$se.fit))
polygon(c(newdata$logarea, rev(newdata$logarea)), c(newdata$upper, rev(newdata$lower)), col = scales::alpha("red", 0.4), lty = 0)
lines(newdata$fit ~ newdata$logarea, lwd = 2, col = "red")
# legend
legend("topright", legend = c("Min. GW", "Max. GW"), bty = "n", fill = c("blue", "red"))
dev.off()


# 2-way interaction: area * dist
jpeg("Model contribution/Effect plots/Cont_area_dist_interaction.jpg", units = "in", width = 4.5, height = 4.25, res = 1500)
par(mar = c(4,4,1,1), mgp = c(2.5,1,0))
# minimum gw
newdata <- with(dat, data.frame(logarea = seq(min(logarea), max(logarea), len = 100),
                                loggwi = seq(mean(loggwi), mean(loggwi), len = 100),
                                distkm = seq(min(distkm), min(distkm), len = 100)))
fit <- predict(mod, newdata = newdata, type = "link", se = TRUE)
newdata <- cbind(newdata, fit = binomial()$linkinv(fit$fit), lower = binomial()$linkinv(fit$fit - q * fit$se.fit), upper = binomial()$linkinv(fit$fit + q * fit$se.fit))
plot(fit ~ logarea, newdata, type = "n", ylim = c(0,0.3), bty = "l", xlab = "log(Basin area, km^2)", ylab = "Proportional contribution")
polygon(c(newdata$logarea, rev(newdata$logarea)), c(newdata$upper, rev(newdata$lower)), col = scales::alpha("blue", 0.4), lty = 0)
lines(newdata$fit ~ newdata$logarea, lwd = 2, col = "blue")
# maximum gw
newdata <- with(dat, data.frame(logarea = seq(min(logarea), max(logarea), len = 100),
                                loggwi = seq(mean(loggwi), mean(loggwi), len = 100),
                                distkm = seq(max(distkm), max(distkm), len = 100)))
fit <- predict(mod, newdata = newdata, type = "link", se = TRUE)
newdata <- cbind(newdata, fit = binomial()$linkinv(fit$fit), lower = binomial()$linkinv(fit$fit - q * fit$se.fit), upper = binomial()$linkinv(fit$fit + q * fit$se.fit))
polygon(c(newdata$logarea, rev(newdata$logarea)), c(newdata$upper, rev(newdata$lower)), col = scales::alpha("red", 0.4), lty = 0)
lines(newdata$fit ~ newdata$logarea, lwd = 2, col = "red")
# legend
legend("topright", legend = c("Min. Distance", "Max. Distance"), bty = "n", fill = c("blue", "red"))
dev.off()

ggplot(data = newdata, aes(y = fit, x = distkm)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.3) + 
  geom_line() +
  scale_y_continuous("Proportional contribution") + theme_classic()

