####################################
# 1_data_process
####################################
# Packages
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, forcats)
#
# Load data
data  <-  read.csv("./data/Learning.csv")
###### Select individuals that pass the criterion of 80% of correct choices in the last 10 trials in the Associative
#
# Calculate the percentage of '1' values in FC_associative for the last 10 trials
percentage_per_lizard <- data %>%
  filter(trial_associative > 25 & trial_associative <= 35) %>%
  group_by(lizard_id) %>%
  summarise(percentage = mean(FC_associative == 1))
# Filter lizards where the percentage is at least 80%
filtered_lizards <- percentage_per_lizard %>%
  filter(percentage >= 0.8)
# Subset the original dataset based on the filtered lizard IDs
clean_df <- data %>%
  filter(lizard_id %in% filtered_lizards$lizard_id) %>%
  filter(sum(is.na(FC_reversal)) <= 15) ungroup()  %>%
    mutate(group = factor(group,
     levels = c("R_B", "B_R"),
     labels=c("R_B"="Red", "B_R"="Blue"))) %>%
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
write.csv(filtered_df, "./output/databases_clean/clean_database.csv")


