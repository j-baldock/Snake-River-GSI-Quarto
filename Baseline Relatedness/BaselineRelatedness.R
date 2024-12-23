#-----------------------------------------------------------------------------------------------------#
# Estimate relatedness with baseline populations
# this needs to run on an older version of R (e.g., 3.6.1)
# I could not install fts and Demerelate packages on my Mac, but was successful on the lab PC
#-----------------------------------------------------------------------------------------------------#

install.packages("sfsmisc")
install.packages("vegan")
install.packages("mlogit")

library(Demerelate)
library(tidyverse)


#-------------------------------------------------------------------#
# Data
#-------------------------------------------------------------------#

# Package test data
# data("demerelpop")
# demerelpop
# Demerelate(inputdata = demerelpop, value = "wang", object = TRUE)

# Upper Snake River Baseline/References Samples
data <- read_csv("G:/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Snake River GSI/GSI Analysis/UpperSnakeRiver_GTseq_InputData.csv") %>% 
  filter(sample_type == "reference", !repunit %in% c("bluecrane_ford", "schwabacher_NA", "southbuffalofork_NA", "bear_NA"))
unique(data$repunit)
data2 <- data %>% select(-c(1,3)) %>% relocate(indiv) #%>% filter(repunit %in% c("bailey_NA", "lowerbarbc_NA"))
colnames(data2)[1] <- "SampleID"
colnames(data2)[2] <- "Population"
data2$Population <- as.factor(data2$Population)
# does not seem to work with NAs, so replace with 0. Unclear how this affects relatedness estimates
data2 <- data2 %>% mutate_all(~replace_na(.,0))


#-------------------------------------------------------------------#
# Run Demerelate
#-------------------------------------------------------------------#

setwd("G:/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Analyses/Snake River GSI/Baseline Relatedness")
base_demerel <- Demerelate(inputdata = as.data.frame(data2), value = "wang", object = TRUE, file.output = TRUE)
write_rds(base_demerel, "UpperSnakeRiver_GTseq_Baseline_Relatedness.rds")

# histogram of relatedness statistics for a single population
hist(base_demerel[[1]]$Empirical_Relatedness$cliff_NA)

# convert wang relatedness estimates to csv
pops <- unique(data2$Population)
wanglist <- list()
for (i in 1:length(pops)) {
  wanglist[[i]] <- tibble(repunit = pops[i],
                          compar = names(unlist(base_demerel[[1]]$Empirical_Relatedness[i])),
                          wang = unlist(base_demerel[[1]]$Empirical_Relatedness[i]))
}
wangtib <- do.call(rbind, wanglist)
wangtib <- wangtib %>% 
  mutate(tmp = str_split(compar, "[.]")) %>% 
  mutate(compar2 = map_chr(tmp, 2)) %>% 
  select(-tmp) %>%
  mutate(tmp = str_split(compar2, "_")) %>%
  mutate(indiv1 = paste(map_chr(tmp, 1), map_chr(tmp, 2), sep = "_"),
         indiv2 = paste(map_chr(tmp, 3), map_chr(tmp, 4), sep = "_")) %>%
  select(-tmp) %>%
  select(repunit, compar2, indiv1, indiv2, wang) %>%
  rename(compar = compar2)
str(wangtib)
write_csv(wangtib, "UpperSnakeRiver_GTseq_Baseline_Relatedness_Wang.csv")

dim(wangtib %>% filter(wang >= 0.4))[1]
view(wangtib %>% filter(wang >= 0.4))


#-------------------------------------------------------------------#
# Filter and remove full siblings
#-------------------------------------------------------------------#

# data <- read_csv("GSI Analysis/UpperSnakeRiver_GTseq_InputData_clean_baseline.csv")
data <- read_csv("GSI Analysis/UpperSnakeRiver_GTseq_InputData.csv") %>% 
  filter(sample_type == "reference", !repunit %in% c("bluecrane_ford", "schwabacher_NA", "southbuffalofork_NA", "bear_NA"))
base_demerel <- read_rds("Baseline Relatedness/UpperSnakeRiver_GTseq_Baseline_Relatedness.rds")
wangtib <- read_csv("Baseline Relatedness/UpperSnakeRiver_GTseq_Baseline_Relatedness_Wang.csv")

# view full sib simulations
range(base_demerel[[1]]$Randomized_Populations_for_Relatedness_Statistics[[1]]$Randomized_Fullssibs)
boxplot(base_demerel[[1]]$Randomized_Populations_for_Relatedness_Statistics[[1]]$Randomized_Fullssibs)

dim(wangtib %>% filter(wang >= 0.4))[1]
view(wangtib %>% filter(wang >= 0.4))

# pull out only highly related individuals (full siblings)
wt2 <- wangtib %>% filter(wang >= 0.4)

# create full sibling famlies
group_pairs <- function(vector_1, vector_2) {
  groups_of_ids <- list()
  df <- cbind(vector_1, vector_2)
  ids_grouped <- c()
  for (id_group in unique(unlist(df))) {
    if (id_group %in% ids_grouped) {
      next
    }
    id <- id_group
    while (length(id) > 0) {
      newid <- unique(unlist(df[which(df[, 1] %in% id | df[, 2] %in% id), ]))
      id <- newid[!newid %in% id_group]
      id_group <- c(id_group, id)
    }
    ids_grouped <- c(ids_grouped, id_group)
    groups_of_ids[[length(groups_of_ids) + 1]] <- sort(unique(id_group))
  }
  return(groups_of_ids)
}
full_sibs_list <- group_pairs(wt2$indiv1, wt2$indiv2)
full_sibs_list
full_sibs <- unlist(full_sibs_list)

# these should match
length(unique(c(wt2$indiv1, wt2$indiv2))) # number of unique individuals from pairwise Wang estimates, where wang >= 0.4
length(unlist(full_sibs)) # number of unique individuals from list of full-sib families

# extract 2 sibling from each full-sib family 
# (see recommendations in Ostergren et al 2020 Mol Ecol Res)
fam_rep <- c(sapply(full_sibs_list, "[[", 1), sapply(full_sibs_list, "[[", 2))

# identify individuals to drop
drop_sibs <- full_sibs[!full_sibs %in% fam_rep]
write.csv(drop_sibs, "Baseline Relatedness/UpperSnakeRiver_GTseq_Baseline_DropFullSiblings.csv", row.names = FALSE)
