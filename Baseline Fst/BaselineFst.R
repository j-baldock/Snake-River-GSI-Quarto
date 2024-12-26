
#-----------------------------------------------------#
# Calculate pairwise Fst among baseline populations
#-----------------------------------------------------#

# see Willing et al. (2012) PLoS One for details on how #SNPs and sample size (# of individuals) affects the accuracy of Fst

library(tidyverse)
library(hierfstat)
library(plot.matrix)
library(viridis)


# load data
refdat <- read_csv("Baseline Data and Testing/UpperSnakeRiver_GTseq_InputData_NoSibs_clean_baseline.csv") %>% mutate_all(~replace_na(., 0))
ss <- refdat %>% group_by(repunit) %>% summarize(ss = n())
# rs1 <- refdat %>% group_by(repunit) %>% summarize(ss_full = n())
# drop_sibs <- read_csv("Baseline Relatedness/UpperSnakeRiver_GTseq_Baseline_DropFullSiblings.csv")

# # drop full siblings
# refdat <- refdat %>% filter(!indiv %in% drop_sibs$x)
# 
# # drop report groups with sample sizes <10 (goosewing, flagstaff, coburn)
# rs2 <- refdat %>% group_by(repunit) %>% summarize(ss_nosibs = n())
# refdat <- refdat %>% filter(repunit %in% ss$repunit)

# reformat genotype data for hierfstat
refdat_sub <- refdat[,c(5:ncol(refdat))] %>% unite(vnew, sep = "") %>% 
  separate(vnew, into = colnames(refdat)[c(seq(from = 5, to = ncol(refdat), by = 2))], sep = c(seq(from = 2, to = 532, by = 2))) 
refdat_sub <- as_tibble(data.frame(sapply(refdat_sub, as.numeric))) %>% mutate_all(~na_if(., 0))
refdat_sub$repunit <- refdat$repunit
refdat_sub <- refdat_sub %>% relocate(repunit)

# calculate basic stats
# pops <- unique(refdat_sub$repunit)
# bslist <- list()
# for (i in 1:length(pops)) {
#   d <- refdat_sub %>% filter(repunit == pops[i])
#   bs <- basic.stats(data.frame(d))
#   bslist[[i]] <- bs$overall
#   print(i)
#   }
basic.stats(data.frame(refdat_sub))

# calculate pairwise Fst following Weir and Cockerham (1984)...also see Harris et al (2022) Conservation Genetics
# NOTE: this takes ~1 hr to run
ref_fst <- pairwise.WCfst(data.frame(refdat_sub))
write.csv(data.frame(ref_fst), "GSI Analysis/Baseline Fst/UpperSnakeRiver_GTseq_BaselinePopsFst.csv", row.names = T)
ref_fst <- read.csv("GSI Analysis/Baseline Fst/UpperSnakeRiver_GTseq_BaselinePopsFst.csv")
gps <- read_csv("Landscape Covariates/RepUnit_LatLong.csv")
range(ref_fst[,-1], na.rm = T)

# plot histogram of pairwise Fst values
jpeg("GSI Analysis/Baseline Fst/BaselinePopsFst_hist.jpg", res = 1000, units = "in", width = 4.5, height = 4.25)
par(mar = c(4,4,1,1))
hist(ref_fst[lower.tri(ref_fst)], xlab = expression("F"["st"]), main = "", breaks = 50)
abline(v = 0.01, col = "red", lwd = 1.5)
dev.off()

# proportion of pairs for which Fst < 0.01
rrr <- ref_fst[lower.tri(ref_fst)]
length(rrr[rrr < 0.01]) / length(rrr)

# for which pairs is Fst < 0.01 
lapply(apply(ref_fst < 0.01, 1, which), names)
# UBBC, Cowboy Cabin, and Snake River Side Channel are all poorly differentiated from each other
# UBBC and Blacktail are poorly differentiated
# Lower and Upper Fall Ck are poorly differentiated

# plot color matrix of Fst
r2 <- ref_fst
r2[upper.tri(r2)] <- NA

jpeg("GSI Analysis/Baseline Fst/BaselinePopsFst_colormatrix.jpg", res = 1000, units = "in", width = 12, height = 12)
par(mar = c(12,12,4,4))
plot(r2, main = expression("F"["st"]), col = rev(magma(20)), las = 2, xlab = "", ylab = "", axes = F, border = "grey50")
axis(1, cex.axis = 0.2)
axis(2, cex.axis = 0.2)
dev.off()


# for Wy-ACT Fish-Flows MS, get range and median Fst among spring creeks with redd count data
unique(ref_fst$X)
ref_

ref_fst <- read.csv("GSI Analysis/Baseline Fst/UpperSnakeRiver_GTseq_BaselinePopsFst.csv", row.names = 1)
mypops <-  c("blacktail_NA", "cody_bluecrane", "cowboycabin_NA", "fish_NA", "flat_NA", "lowerbarbc_NA", "snakeriversidechannel_NA", "spring_tss", "threechannel_NA", "upperbarbc_NA")
ref_fst2 <- ref_fst[mypops, mypops]
ref_fst3 <- ref_fst2[lower.tri(ref_fst2)]

range(ref_fst3, na.rm = T)
median(ref_fst3, na.rm = T)

lower.tri()




