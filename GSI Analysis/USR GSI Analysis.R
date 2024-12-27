####------------------------------------------------------------------------------#
#### Upper Snake Genetics Stock Identification Run GSI and examine z-scores ####
####------------------------------------------------------------------------------#

library(rubias)
library(tidyverse)


####---------------------------------------------------#
####  Read in formatted data #### 
####  See BaselineDataTesting.R
####---------------------------------------------------#

# formatted reference and mixture genotypes
ref_input <- read_csv("Baseline Data and Testing/UpperSnakeRiver_GTseq_InputData_NoSibs_clean_baseline.csv")
mix_input <- read_csv("Baseline Data and Testing/UpperSnakeRiver_GTseq_InputData_NoSibs_clean_mixture.csv") %>% filter(collection != "snake_jldcattlemens")
ref_input <- ref_input %>% mutate_if(is.double, as.integer) # make sure everything is an interger
mix_input <- mix_input %>% mutate_if(is.double, as.integer) # make sure everything is an interger

# get metadata
ids <- read_csv("LabFieldIDs.csv")
data <- read_csv("/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Methods/GSI sampling/Snake_GSI_field data_2020-2022_GenOnly_DropNoDrop_090823edit.csv")
metadata <- left_join(ids, data) %>% 
  select(indiv, CORRgenID, year, site, TL.mm, TL.in, FL.mm, weight.g, weight.lbs, collection.type, sizecat) %>% 
  filter(collection.type == "mixture", site != "snake_jldcattlemens")
metadata[duplicated(metadata$indiv),] # check for duplicates

# mixture global
mix_input_global <- mix_input %>% mutate(collection = "snake")
unique(mix_input_global$collection)

# mixtures by year
metayr <- metadata %>% select(indiv, year)
mix_input_yr <- mix_input %>% left_join(metayr) %>% mutate(collection = paste(collection, year, sep = "_")) %>% select(-year)
unique(mix_input_yr$collection)

# mixtures by size only
metasz <- metadata %>% select(indiv, sizecat)
mix_input_sz <- mix_input %>% left_join(metasz) %>% mutate(collection = paste("size", sizecat, sep = "_")) %>% select(-sizecat)
sort(unique(mix_input_sz$collection))
mix_input_sz %>% group_by(collection) %>% summarize(catch = n()) %>% ggplot() + geom_bar(aes(x = collection, y = catch), stat = "identity")


####---------------------------------------------------#
#### Run GSI - Global ####
####---------------------------------------------------#

### method BR does not work with more than 1000 markers!
### method PB provides bootstrapped corrected mixing proportions
gsi_global <- infer_mixture(reference = ref_input, mixture = mix_input_global, gen_start_col = 5, method = "PB")
saveRDS(gsi_global, "GSI Analysis/Global/UpperSnake_GSI_Global_output.RDS")
write_csv(gsi_global$bootstrapped_proportions, "GSI Analysis/Global/UpperSnake_GSI_Global_BootstrappedProportions.csv")


#### GSI Diagnostics: MAP, z-scores ####

# get the maximum-a-posteriori population for each individual
map_rows <- gsi_global$indiv_posteriors %>%
  group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()
ss <- map_rows %>% group_by(repunit) %>% summarise(n = n())

jpeg("GSI Analysis/Global/gsi_IndividualPosteriorAssignmentProbability.jpg", res = 1000, units = "in", width = 7, height = 5)
ggplot(map_rows, aes(x = repunit, y = PofZ)) +
  geom_boxplot(aes(colour = repunit)) +
  geom_text(data = ss, mapping = aes(y = 1.025, label = n), angle = 90, hjust = 0, vjust = 0.5, size = 3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9, vjust = 0.5)) +
  ylim(c(NA, 1.05)) +
  guides(colour = FALSE)
dev.off()

# plot density of observed z-scores (blue) against expected (black)...if all contributing pops are represented in the baseline, distributions should be similar
normo <- tibble(z_score = rnorm(1e06))
ks.test(normo, map_rows$z_score) # p << 0.05 indicates there are problematic individuals (distributions are not the same)

jpeg("GSI Analysis/Global/gsi_zscore_dist.jpg", res = 1000, units = "in", width = 5, height = 5)
ggplot(map_rows, aes(x = z_score)) +
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black")
dev.off()

# Following Bowersox, Hargrove, et al 2023 (NAJFM), individual likely originated from an unsampled pop if +/- 2 SDs from the mean z-score
meanz <- mean(map_rows$z_score)
zproblems <- map_rows %>% filter(z_score > meanz+2 | z_score < meanz-2)
dim(zproblems)[1] / dim(map_rows)[1] # fish originating from pops outside our baseline represent ~6.4% of all fish sampled in the mixture

jpeg("GSI Analysis/Global/zscore_comparison_boxplot.jpg", res = 1000, units = "in", width = 11, height = 11)
par(mfrow = c(6,7), mar = c(2,2.5,3,0.1))
sites <- sort(unique(map_rows$repunit))
for (i in 1:length(sites)) {
  mx <- map_rows %>% filter(repunit == sites[i])
  bl <- self_test_exam %>% filter(repunit == sites[i])
  ztib <- tibble(type = c(rep("mx", times = dim(mx)[1]), rep("bl", times = dim(bl)[1])), z = c(mx$z_score, bl$z_score))
  boxplot(z ~ type, ztib, col = c("blue", "red"), main = sites[i])
}
dev.off()

# Distribution of z-scores by population
jpeg("GSI Analysis/Global/gsi_zscore_dist_bypop.jpg", res = 1000, units = "in", width = 9, height = 5)
ggplot(map_rows, aes(x = repunit, y = z_score)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# view posterior density curves from raw individual data
# not super legit given that we will ultimate used bootstrapped proportions
# select the top 10 most abundant
# pp <- "snake_moosewilson"
top10 <- gsi_global$mixing_proportions %>% 
  group_by(mixture_collection, repunit) %>% 
  summarise(repprop = sum(pi)) %>% 
  # filter(mixture_collection == pp) %>%
  arrange(desc(repprop)) %>% slice(1:10)
# check how many MCMC sweeps were done
nsweeps <- max(mixture_analysis_init$mix_prop_traces$sweep)
# discard first 200 sweeps as burn-in and select traces from top 10
trace_subset <- gsi_global$mix_prop_traces %>%
  filter(sweep > 200) %>%
  group_by(sweep, repunit) %>%
  summarise(repprop = sum(pi)) %>% 
  filter(repunit %in% top10$repunit)
# plot density traces...fair amount of uncertainty
ggplot(trace_subset, aes(x = repprop, colour = repunit)) + geom_density()
# compute credible intervals
top10_cis <- trace_subset %>% group_by(repunit) %>%
  summarize(loci = quantile(repprop, probs = 0.025),
            medi = quantile(repprop, probs = 0.5),
            hici = quantile(repprop, probs = 0.975)) %>%
  arrange(desc(medi))
# view as bar plots
ggplot(top10_cis, aes(x = repunit, y = medi)) + geom_col() + 
  geom_errorbar(aes(ymin = loci, ymax = hici), width = 0.2, position=position_dodge(.9)) +
  scale_x_discrete(limits = c(top10_cis$repunit)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Reporting unit") + ylab(paste("Median contribution w/ 95% CI (", pp, ")", sep = ""))
  

#### Bootstrap-corrected mixing proportions ####
# because from simulations etc. we know that assignments will be biased

# summarize raw and bootstrap corrected mixing proportions by reporting group
tt <- gsi_global$mixing_proportions %>% 
  group_by(mixture_collection, repunit) %>% summarize(pi = sum(pi)) %>% ungroup() %>% 
  left_join(gsi_global$bootstrapped_proportions) %>%
  mutate(diff = bs_corrected_repunit_ppn - pi)
tt %>% group_by(mixture_collection) %>% summarize(mr = sum(bs_corrected_repunit_ppn - pi) / length(unique(repunit)))

# plot raw vs bootstraps
jpeg("GSI Analysis/Global/gsi_MixingProportions_RawVsBootstrapped.jpg", res = 1000, units = "in", width = 7, height = 5)
tt %>% ggplot() + 
  geom_point(aes(x = pi, y = bs_corrected_repunit_ppn)) + 
  geom_abline(intercept = 0, slope = 1, col = "red") + 
  # facet_wrap(~ factor(mixture_collection, levels = c("snake_pacificdeadmans", "snake_deadmansmoose", "snake_moosewilson", "snake_wilsonsouthpark", "snake_southparkastoria", "snake_astoriawesttable"))) + 
  xlab("Mixture Proportion") + ylab("Corrected Mixture Proportion") + theme_bw()
dev.off()

# view (bootstrap corrected) reporting group contributions by river section as bar plots
jpeg("GSI Analysis/Global/gsi_MixingProportions_bySection.jpg", res = 1000, units = "in", width = 8, height = 4)#9)
gsi_global$bootstrapped_proportions %>% 
  ggplot() +
  geom_bar(aes(x = repunit, y = bs_corrected_repunit_ppn), stat = "identity", fill = "grey70", color = "black") +
  # facet_wrap(~ factor(mixture_collection, levels = c("snake_pacificdeadmans", "snake_deadmansmoose", "snake_moosewilson", "snake_wilsonsouthpark", "snake_southparkastoria", "snake_astoriawesttable")), nrow = 6) +
  xlab("Reporting Group") + ylab("Corrected Mixture Proportion") + theme_bw() +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()


####---------------------------------------------------#
#### Run GSI - By Section ####
####---------------------------------------------------#

### method BR does not work with more than 1000 markers!
### method PB provides bootstrapped corrected mixing proportions
gsi_section <- infer_mixture(reference = ref_input, mixture = mix_input, gen_start_col = 5, method = "PB")
saveRDS(gsi_section, "GSI Analysis/By Section/UpperSnake_GSI_Section_output.RDS")
write_csv(gsi_section$bootstrapped_proportions, "GSI Analysis/By Section/UpperSnake_GSI_Section_BootstrappedProportions.csv")


#### GSI Diagnostics: MAP, z-scores ####

# get the maximum-a-posteriori population for each individual
map_rows <- gsi_section$indiv_posteriors %>%
  group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()
ss <- map_rows %>% group_by(repunit) %>% summarise(n = n())

jpeg("GSI Analysis/By Section/gsi_IndividualPosteriorAssignmentProbability.jpg", res = 1000, units = "in", width = 7, height = 5)
ggplot(map_rows, aes(x = repunit, y = PofZ)) +
  geom_boxplot(aes(colour = repunit)) +
  geom_text(data = ss, mapping = aes(y = 1.025, label = n), angle = 90, hjust = 0, vjust = 0.5, size = 3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9, vjust = 0.5)) +
  ylim(c(NA, 1.05)) +
  guides(colour = FALSE)
dev.off()

# plot density of observed z-scores (blue) against expected (black)...if all contributing pops are represented in the baseline, distributions should be similar
normo <- tibble(z_score = rnorm(1e06))
ks.test(normo, map_rows$z_score) # p << 0.05 indicates there are problematic individuals (distributions are not the same)

jpeg("GSI Analysis/By Section/gsi_zscore_dist.jpg", res = 1000, units = "in", width = 5, height = 5)
ggplot(map_rows, aes(x = z_score)) +
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black")
dev.off()

# Following Bowersox, Hargrove, et al 2023 (NAJFM), individual likely originated from an unsampled pop if +/- 2 SDs from the mean z-score
meanz <- mean(map_rows$z_score)
zproblems <- map_rows %>% filter(z_score > meanz+2 | z_score < meanz-2)
dim(zproblems)[1] / dim(map_rows)[1] # fish originating from pops outside our baseline represent ~6.6% of all fish sampled in the mixture

jpeg("GSI Analysis/By Section/zscore_comparison_boxplot.jpg", res = 1000, units = "in", width = 11, height = 11)
par(mfrow = c(6,7), mar = c(2,2.5,3,0.1))
sites <- sort(unique(map_rows$repunit))
for (i in 1:length(sites)) {
  mx <- map_rows %>% filter(repunit == sites[i])
  bl <- self_test_exam %>% filter(repunit == sites[i])
  ztib <- tibble(type = c(rep("mx", times = dim(mx)[1]), rep("bl", times = dim(bl)[1])), z = c(mx$z_score, bl$z_score))
  boxplot(z ~ type, ztib, col = c("blue", "red"), main = sites[i])
}
dev.off()

# Distribution of z-scores by population
jpeg("GSI Analysis/By Section/gsi_zscore_dist_bypop.jpg", res = 1000, units = "in", width = 9, height = 5)
ggplot(map_rows, aes(x = repunit, y = z_score)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# view posterior density curves from raw individual data
# not super legit given that we will ultimate used bootstrapped proportions
# select the top 10 most abundant
pp <- "snake_moosewilson"
top10 <- gsi_section$mixing_proportions %>% 
  group_by(mixture_collection, repunit) %>% 
  summarise(repprop = sum(pi)) %>% 
  filter(mixture_collection == pp) %>%
  arrange(desc(repprop)) %>% slice(1:10)
# check how many MCMC sweeps were done
nsweeps <- max(mixture_analysis_init$mix_prop_traces$sweep)
# discard first 200 sweeps as burn-in and select traces from top 10
trace_subset <- gsi_global$mix_prop_traces %>%
  filter(sweep > 200) %>%
  group_by(sweep, repunit) %>%
  summarise(repprop = sum(pi)) %>% 
  filter(repunit %in% top10$repunit)
# plot density traces...fair amount of uncertainty
ggplot(trace_subset, aes(x = repprop, colour = repunit)) + geom_density() #+ facet_wrap(~collection)
# compute credible intervals
top10_cis <- trace_subset %>% group_by(repunit) %>%
  summarize(loci = quantile(repprop, probs = 0.025),
            medi = quantile(repprop, probs = 0.5),
            hici = quantile(repprop, probs = 0.975)) %>%
  arrange(desc(medi))
# view as bar plots
ggplot(top10_cis, aes(x = repunit, y = medi)) + geom_col() + #facet_wrap(~collection) + 
  geom_errorbar(aes(ymin = loci, ymax = hici), width = 0.2, position=position_dodge(.9)) +
  scale_x_discrete(limits = c(top10_cis$repunit)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Reporting unit") + ylab(paste("Median contribution w/ 95% CI (", pp, ")", sep = ""))


#### Bootstrap-corrected mixing proportions ####
# because from simulations etc. we know that assignments will be biased

# summarize raw and bootstrap corrected mixing proportions by reporting group
tt <- gsi_section$mixing_proportions %>% 
  group_by(mixture_collection, repunit) %>% summarize(pi = sum(pi)) %>% ungroup() %>% 
  left_join(gsi_section$bootstrapped_proportions) %>%
  mutate(diff = bs_corrected_repunit_ppn - pi)
tt %>% group_by(mixture_collection) %>% summarize(mr = sum(bs_corrected_repunit_ppn - pi) / length(unique(repunit)))

# plot raw vs bootstraps
jpeg("GSI Analysis/By Section/gsi_MixingProportions_RawVsBootstrapped.jpg", res = 1000, units = "in", width = 7, height = 5)
tt %>% ggplot() + 
  geom_point(aes(x = pi, y = bs_corrected_repunit_ppn)) + 
  geom_abline(intercept = 0, slope = 1, col = "red") + 
  facet_wrap(~ factor(mixture_collection, levels = c("snake_pacificdeadmans", "snake_deadmansmoose", "snake_moosewilson", "snake_wilsonsouthpark", "snake_southparkastoria", "snake_astoriawesttable"))) + 
  xlab("Mixture Proportion") + ylab("Corrected Mixture Proportion") + theme_bw()
dev.off()

# view (bootstrap corrected) reporting group contributions by river section as bar plots
jpeg("GSI Analysis/By Section/gsi_MixingProportions_bySection.jpg", res = 1000, units = "in", width = 8, height = 9)
gsi_section$bootstrapped_proportions %>% 
  ggplot() +
  geom_bar(aes(x = repunit, y = bs_corrected_repunit_ppn), stat = "identity", fill = "grey70", color = "black") +
  facet_wrap(~ factor(mixture_collection, levels = c("snake_pacificdeadmans", "snake_deadmansmoose", "snake_moosewilson", "snake_wilsonsouthpark", "snake_southparkastoria", "snake_astoriawesttable")), nrow = 6) +
  xlab("Reporting Group") + ylab("Corrected Mixture Proportion") + theme_bw() +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()


####---------------------------------------------------#
#### Run GSI - By Section and Year ####
####---------------------------------------------------#

### method BR does not work with more than 1000 markers!
### method PB provides bootstrapped corrected mixing proportions
gsi_sectyear <- infer_mixture(reference = ref_input, mixture = mix_input_yr, gen_start_col = 5, method = "PB")
saveRDS(gsi_sectyear, "GSI Analysis/By Section and Year/UpperSnake_GSI_SectionYear_output.RDS")
write_csv(gsi_sectyear$bootstrapped_proportions, "GSI Analysis/By Section and Year/UpperSnake_GSI_SectionYear_BootstrappedProportions.csv")
gsi_sectyear <- readRDS("GSI Analysis/By Section and Year/UpperSnake_GSI_SectionYear_output.RDS")

#### GSI Diagnostics: MAP, z-scores ####

# get the maximum-a-posteriori population for each individual
map_rows <- gsi_sectyear$indiv_posteriors %>%
  group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()
ss <- map_rows %>% group_by(repunit) %>% summarise(n = n())

jpeg("GSI Analysis/By Section and Year/gsi_IndividualPosteriorAssignmentProbability.jpg", res = 1000, units = "in", width = 7, height = 5)
ggplot(map_rows, aes(x = repunit, y = PofZ)) +
  geom_boxplot(aes(colour = repunit)) +
  geom_text(data = ss, mapping = aes(y = 1.025, label = n), angle = 90, hjust = 0, vjust = 0.5, size = 3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9, vjust = 0.5)) +
  ylim(c(NA, 1.05)) +
  guides(colour = FALSE)
dev.off()

# plot density of observed z-scores (blue) against expected (black)...if all contributing pops are represented in the baseline, distributions should be similar
normo <- tibble(z_score = rnorm(1e06))
ks.test(normo, map_rows$z_score) # p << 0.05 indicates there are problematic individuals (distributions are not the same)

jpeg("GSI Analysis/By Section and Year/gsi_zscore_dist.jpg", res = 1000, units = "in", width = 5, height = 5)
ggplot(map_rows, aes(x = z_score)) +
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black")
dev.off()

# Following Bowersox, Hargrove, et al 2023 (NAJFM), individual likely originated from an unsampled pop if +/- 2 SDs from the mean z-score
meanz <- mean(map_rows$z_score)
zproblems <- map_rows %>% filter(z_score > meanz+2 | z_score < meanz-2)
dim(zproblems)[1] / dim(map_rows)[1] # fish originating from pops outside our baseline represent ~6.6% of all fish sampled in the mixture

jpeg("GSI Analysis/By Section and Year/zscore_comparison_boxplot.jpg", res = 1000, units = "in", width = 11, height = 11)
par(mfrow = c(6,7), mar = c(2,2.5,3,0.1))
sites <- sort(unique(map_rows$repunit))
for (i in 1:length(sites)) {
  mx <- map_rows %>% filter(repunit == sites[i])
  bl <- self_test_exam %>% filter(repunit == sites[i])
  ztib <- tibble(type = c(rep("mx", times = dim(mx)[1]), rep("bl", times = dim(bl)[1])), z = c(mx$z_score, bl$z_score))
  boxplot(z ~ type, ztib, col = c("blue", "red"), main = sites[i])
}
dev.off()

# Distribution of z-scores by population
jpeg("GSI Analysis/By Section and Year/gsi_zscore_dist_bypop.jpg", res = 1000, units = "in", width = 9, height = 5)
ggplot(map_rows, aes(x = repunit, y = z_score)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# view posterior density curves from raw individual data
# not super legit given that we will ultimate used bootstrapped proportions
# select the top 10 most abundant
pp <- "snake_moosewilson_2021"
top10 <- gsi_sectyear$mixing_proportions %>% 
  group_by(mixture_collection, repunit) %>% 
  summarise(repprop = sum(pi)) %>% 
  filter(mixture_collection == pp) %>%
  arrange(desc(repprop)) %>% slice(1:10)
# check how many MCMC sweeps were done
nsweeps <- max(gsi_sectyear$mix_prop_traces$sweep)
# discard first 200 sweeps as burn-in and select traces from top 10
trace_subset <- gsi_sectyear$mix_prop_traces %>%
  filter(sweep > 200) %>%
  group_by(sweep, repunit) %>%
  summarise(repprop = sum(pi)) %>% 
  filter(repunit %in% top10$repunit)
# plot density traces...fair amount of uncertainty
ggplot(trace_subset, aes(x = repprop, colour = repunit)) + geom_density() #+ facet_wrap(~collection)
# compute credible intervals
top10_cis <- trace_subset %>% group_by(repunit) %>%
  summarize(loci = quantile(repprop, probs = 0.025),
            medi = quantile(repprop, probs = 0.5),
            hici = quantile(repprop, probs = 0.975)) %>%
  arrange(desc(medi))
# view as bar plots
ggplot(top10_cis, aes(x = repunit, y = medi)) + geom_col() + #facet_wrap(~collection) + 
  geom_errorbar(aes(ymin = loci, ymax = hici), width = 0.2, position=position_dodge(.9)) +
  scale_x_discrete(limits = c(top10_cis$repunit)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Reporting unit") + ylab(paste("Median contribution w/ 95% CI (", pp, ")", sep = ""))


#### Bootstrap-corrected mixing proportions ####
# because from simulations etc. we know that assignments will be biased

# summarize raw and bootstrap corrected mixing proportions by reporting group
tt <- gsi_sectyear$mixing_proportions %>% 
  group_by(mixture_collection, repunit) %>% summarize(pi = sum(pi)) %>% ungroup() %>% 
  left_join(gsi_sectyear$bootstrapped_proportions) %>%
  mutate(diff = bs_corrected_repunit_ppn - pi,
         section = str_sub(mixture_collection, 7, -6),
         year = str_sub(mixture_collection, -4))
tt %>% group_by(mixture_collection) %>% summarize(mr = sum(bs_corrected_repunit_ppn - pi) / length(unique(repunit)))
tt <- tt %>% mutate(sectionID = recode(section, 
                                 "pacificdeadmans" = "A",
                                 "deadmansmoose" = "B",
                                 "moosewilson" = "C",
                                 "wilsonsouthpark" = "D",
                                 "southparkastoria" = "E",
                                 "astoriawesttable" = "F"),
                    year = as.numeric(year))
tt <- add_row(tt, year = 2020, sectionID = "A")
cor.test(tt$pi, tt$bs_corrected_repunit_ppn)

# plot raw vs bootstraps
jpeg("GSI Analysis/By Section and Year/gsi_MixingProportions_RawVsBootstrapped.jpg", res = 1000, units = "in", width = 6, height = 9)
tt %>% ggplot() + 
  geom_point(aes(x = pi, y = bs_corrected_repunit_ppn)) + 
  geom_abline(intercept = 0, slope = 1, col = "red") + 
  facet_grid(factor(sectionID) ~ year) + 
  xlab("Uncorrected mixing proportion") + ylab("Bootstrap-corrected mixing proportion") + theme_bw()
dev.off()

# view (bootstrap corrected) reporting group contributions by river section as bar plots
jpeg("GSI Analysis/By Section and Year/gsi_MixingProportions_bySectionYear.jpg", res = 1000, units = "in", width = 14, height = 9)
tt %>% 
  add_row(repunit = unique(tt$repunit), bs_corrected_repunit_ppn = rep(0, 52), section = rep("pacificdeadmans", 52), year = rep("2020", 52)) %>%
  ggplot() +
  geom_bar(aes(x = repunit, y = bs_corrected_repunit_ppn, fill = year), stat = "identity", position = "dodge", color = "black") +
  facet_wrap(~ factor(section, levels = c("pacificdeadmans", "deadmansmoose", "moosewilson", "wilsonsouthpark", "southparkastoria", "astoriawesttable")), nrow = 6) +
  xlab("Reporting Group") + ylab("Corrected Mixture Proportion") + theme_bw() +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()


####---------------------------------------------------#
#### Run GSI - By Size Class ####
####---------------------------------------------------#

### method BR does not work with more than 1000 markers!
### method PB provides bootstrapped corrected mixing proportions
gsi_size <- infer_mixture(reference = ref_input, mixture = mix_input_sz, gen_start_col = 5, method = "PB")
saveRDS(gsi_size, "GSI Analysis/By Size/UpperSnake_GSI_Global_output.RDS")
write_csv(gsi_size$bootstrapped_proportions, "GSI Analysis/By Size/UpperSnake_GSI_Size_BootstrappedProportions.csv")


#### GSI Diagnostics: MAP, z-scores ####

# get the maximum-a-posteriori population for each individual
map_rows <- gsi_size$indiv_posteriors %>%
  group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()
ss <- map_rows %>% group_by(repunit) %>% summarise(n = n())

jpeg("GSI Analysis/By Size/gsi_IndividualPosteriorAssignmentProbability.jpg", res = 1000, units = "in", width = 7, height = 5)
ggplot(map_rows, aes(x = repunit, y = PofZ)) +
  geom_boxplot(aes(colour = repunit)) +
  geom_text(data = ss, mapping = aes(y = 1.025, label = n), angle = 90, hjust = 0, vjust = 0.5, size = 3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9, vjust = 0.5)) +
  ylim(c(NA, 1.05)) +
  guides(colour = FALSE)
dev.off()

# plot density of observed z-scores (blue) against expected (black)...if all contributing pops are represented in the baseline, distributions should be similar
normo <- tibble(z_score = rnorm(1e06))
ks.test(normo, map_rows$z_score) # p << 0.05 indicates there are problematic individuals (distributions are not the same)

jpeg("GSI Analysis/By Size/gsi_zscore_dist.jpg", res = 1000, units = "in", width = 5, height = 5)
ggplot(map_rows, aes(x = z_score)) +
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black")
dev.off()

# Following Bowersox, Hargrove, et al 2023 (NAJFM), individual likely originated from an unsampled pop if +/- 2 SDs from the mean z-score
meanz <- mean(map_rows$z_score)
zproblems <- map_rows %>% filter(z_score > meanz+2 | z_score < meanz-2)
dim(zproblems)[1] / dim(map_rows)[1] # fish originating from pops outside our baseline represent ~6.6% of all fish sampled in the mixture

jpeg("GSI Analysis/By Size/zscore_comparison_boxplot.jpg", res = 1000, units = "in", width = 11, height = 11)
par(mfrow = c(6,7), mar = c(2,2.5,3,0.1))
sites <- sort(unique(map_rows$repunit))
for (i in 1:length(sites)) {
  mx <- map_rows %>% filter(repunit == sites[i])
  bl <- self_test_exam %>% filter(repunit == sites[i])
  ztib <- tibble(type = c(rep("mx", times = dim(mx)[1]), rep("bl", times = dim(bl)[1])), z = c(mx$z_score, bl$z_score))
  boxplot(z ~ type, ztib, col = c("blue", "red"), main = sites[i])
}
dev.off()

# Distribution of z-scores by population
jpeg("GSI Analysis/By Size/gsi_zscore_dist_bypop.jpg", res = 1000, units = "in", width = 9, height = 5)
ggplot(map_rows, aes(x = repunit, y = z_score)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# view posterior density curves from raw individual data
# not super legit given that we will ultimate used bootstrapped proportions
# select the top 10 most abundant
pp <- "size_6"
top10 <- gsi_size$mixing_proportions %>% 
  group_by(mixture_collection, repunit) %>% 
  summarise(repprop = sum(pi)) %>% 
  filter(mixture_collection == pp) %>%
  arrange(desc(repprop)) %>% slice(1:10)
# check how many MCMC sweeps were done
nsweeps <- max(gsi_size$mix_prop_traces$sweep)
# discard first 200 sweeps as burn-in and select traces from top 10
trace_subset <- gsi_size$mix_prop_traces %>%
  filter(sweep > 200) %>%
  group_by(sweep, repunit) %>%
  summarise(repprop = sum(pi)) %>% 
  filter(repunit %in% top10$repunit)
# plot density traces...fair amount of uncertainty
ggplot(trace_subset, aes(x = repprop, colour = repunit)) + geom_density() #+ facet_wrap(~collection)
# compute credible intervals
top10_cis <- trace_subset %>% group_by(repunit) %>%
  summarize(loci = quantile(repprop, probs = 0.025),
            medi = quantile(repprop, probs = 0.5),
            hici = quantile(repprop, probs = 0.975)) %>%
  arrange(desc(medi))
# view as bar plots
ggplot(top10_cis, aes(x = repunit, y = medi)) + geom_col() + #facet_wrap(~collection) + 
  geom_errorbar(aes(ymin = loci, ymax = hici), width = 0.2, position=position_dodge(.9)) +
  scale_x_discrete(limits = c(top10_cis$repunit)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Reporting unit") + ylab(paste("Median contribution w/ 95% CI (", pp, ")", sep = ""))


#### Bootstrap-corrected mixing proportions ####
# because from simulations etc. we know that assignments will be biased

# summarize raw and bootstrap corrected mixing proportions by reporting group
tt <- gsi_size$mixing_proportions %>% 
  group_by(mixture_collection, repunit) %>% summarize(pi = sum(pi)) %>% ungroup() %>% 
  left_join(gsi_size$bootstrapped_proportions) %>%
  mutate(diff = bs_corrected_repunit_ppn - pi)
tt %>% group_by(mixture_collection) %>% summarize(mr = sum(bs_corrected_repunit_ppn - pi) / length(unique(repunit)))

# plot raw vs bootstraps
jpeg("GSI Analysis/By Size/gsi_MixingProportions_RawVsBootstrapped.jpg", res = 1000, units = "in", width = 7, height = 7)
tt %>% ggplot() + 
  geom_point(aes(x = pi, y = bs_corrected_repunit_ppn)) + 
  geom_abline(intercept = 0, slope = 1, col = "red") + 
  facet_wrap(~ mixture_collection) + 
  xlab("Mixture Proportion") + ylab("Corrected Mixture Proportion") + theme_bw()
dev.off()

# view (bootstrap corrected) reporting group contributions by river section as bar plots
jpeg("GSI Analysis/By Size/gsi_MixingProportions_bySize.jpg", res = 1000, units = "in", width = 8, height = 9)
tt %>% 
  ggplot() +
  geom_bar(aes(x = repunit, y = bs_corrected_repunit_ppn), stat = "identity", fill = "grey70", color = "black") +
  facet_wrap(~ mixture_collection, nrow = 6) +
  xlab("Reporting Group") + ylab("Corrected Mixture Proportion") + theme_bw() +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()

