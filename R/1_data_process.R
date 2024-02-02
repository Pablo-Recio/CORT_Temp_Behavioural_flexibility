####################################
# 1_data_process
####################################
# Packages
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, forcats)
#
# Load data
data  <-  read.csv("./data/Learning.csv")
#
clean_df <- data %>%
  group_by(lizard_id) %>%
  filter(sum(is.na(FC_reversal)) <= 15) %>%
    mutate(temp = gsub("[AB]_", "", trt),
          cort = gsub("_[2][38]", "", trt))  %>%
    mutate(temp = factor(temp,
      levels = c("23", "28"),
      labels = c("23" = "Cold", "28" = "Hot"))) %>%
    mutate(cort = factor(cort,
      levels = c("B", "A"),
      labels = c("B" = "CORT", "A" = "Control"))) %>%
    mutate(trial_reversal=as.numeric(trial_reversal)) %>%
  data.frame()
write.csv(clean_df, "./output/databases_clean/clean_df.csv")


