#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: setup
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, ggh4x, PupillometryR, cowplot, png, grid)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: fig-Methods
#| fig.cap: "Experimental design of early environment manipulation and learning tasks. Stages 1-3 indicate the different phases of the habituation phase. In the associative and reversal tasks, white lids show the ramps where the food reward was not accessible."

knitr::include_graphics("./Others/BEHFLEX_FIG_1.svg")

#
#
#
#
#
#
#
#
#
#| label: cleandata
# Obtain the main df using "./R/1_data_process.R"
source(here("R", "1_data_process.R"))
#
#
#
#| label: sampleSize
# List with the sample sizes from the database (here("output/databases_clean/data_asso.csv") as the sample size per species and group is the same on each task. We used function sample (see func.R) to estimate the sample size per treatment and species.
source(here("R", "func.R"))
#
specie <- c("delicata", "guichenoti")
hormone <- c("CORT", "Control")
temperature <- c("Cold", "Hot")
#
n_list <- list()
#
for(i in 1:length(specie)){
  for(k in 1:length(hormone)){
    for(l in 1:length(temperature)){
      n_list[[paste(specie[i], hormone[k], temperature[l], sep = "_")]] <- sample(clean_df, specie[i], hormone[k], temperature[l])
      }
    }
  }
#
#
#
#
#
#| label: models
# Fitting the model and extraction of posteriors for both types of task and species using fit_m function (see func.r in R folder). The result everytime the function is used is a df with the posteriors of the model. The functions saves the model automatically in output/models; and when the parameter refit = FALSE then the posteriors are extracted from the model previously written instead of fitting the model again each time.
source(here("R", "func.R"))
# Model formula: FC_reversal ~ trial_reversal*cort*temp + (1 + trial_reversal|lizard_id)
## A) L. delicata
deli <- fit_m(clean_df, "deli", "complete",  refit = FALSE)
write.csv(deli, file= "./output/Checking/deli.csv")
## B) L. guichenoti
guich <- fit_m(clean_df, "guich", "complete", refit = FALSE)
write.csv(guich, file= "./output/Checking/guich.csv")
#
#
#
# Rename some of the posteriors and make new estimates for the learning rate for the Reversal task doing the same thing we did in the chunk above.
## 1) L. delicata
deli_CORTCold <- deli$b_trial_reversal
deli_ControlCold <- (deli$'b_trial_reversal:cortControl' + deli$b_trial_reversal)
deli_CORTHot <- (deli$'b_trial_reversal:tempHot' + deli$b_trial_reversal)
deli_ControlHot <- (deli$'b_trial_reversal:cortControl:tempHot' + deli$b_trial_reversal + deli$'b_trial_reversal:cortControl' + deli$'b_trial_reversal:tempHot')
## 2) L. guichenoti
guich_CORTCold <- guich$b_trial_reversal
guich_ControlCold <- (guich$'b_trial_reversal:cortControl' + guich$b_trial_reversal)
guich_CORTHot <- (guich$'b_trial_reversal:tempHot' + guich$b_trial_reversal)
guich_ControlHot <- (guich$'b_trial_reversal:cortControl:tempHot' + guich$b_trial_reversal + guich$'b_trial_reversal:cortControl' + guich$'b_trial_reversal:tempHot')
#
#
#
#
#
#| label: fig-results
#| fig.cap: "Results for L. delicata (A,B) and L. guichenoti (C, D). A, C) Predicted probability of choosing the correct feeder first over trials. The lines represent the mean predicted probability of choosing the correct feeder first, and the shaded areas represent the standard deviation of the mean; both obtained by using the slope and intercept estimates from the posterior distributions. The different colours represent the different treatments. B, D) Distribution of the estimates of the slopes per each treatment. The x-axis represents the estimates of the slopes, and the y-axis represents the density of the estimates. The different colours represent the different treatments. Points and bars represent the mean and standard deviation of the mean of the estimates, respectively."
source(here("R", "func.R"))
# Create the dfs for the probability plot A and C) (df_prob)
# L. delicata
df_deliprob <- df_plotAC("deli")
# L. guichenoti
df_guichprob <- df_plotAC("guich")
#
# Create the df for the violin plots in plot B) (df_violin)
# L. delicata
df_deliV <- df_plotBD1("deli")
# L. guichenoti
df_guichV <- df_plotBD1("guich")
#
# Create the df for the points and bars in plot B) (df_points)
# L. delicata
df_deliP <- df_plotBD2("deli")
# L. guichenoti
df_guichP <- df_plotBD2("guich")
#
#
# Make plots deli and guich
source(here("R", "func.R"))
fig_deli <- plotting(sp = "deli", df_prob = df_deliprob, df_violin = df_deliV, df_points = df_deliP)
fig_guich <- plotting(sp = "guich", df_prob = df_guichprob, df_violin = df_guichV, df_points = df_guichP)
# Combine plots of both species into the final figure
fig_results <- plot_grid(fig_deli, fig_guich, ncol = 1)
ggsave("./output/figures/fig_results.png", plot=fig_results, width = 25, height = 18, units = "cm", dpi = 600)
knitr::include_graphics("./output/figures/fig_deli.png")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
cat("\\newpage")
#
#
#
#
# Chunk for calling all the models and making the residuals
## L. delicata
mod_deli <- readRDS(here("output/models/delicomplete.rds"))
resid_deli <- residuals(mod_deli)
## L. guichenoti
mod_guich <- readRDS(here("output/models/guichcomplete.rds"))
resid_guich <- residuals(mod_guich)
#
#
#
#
#
#
#
bayes_R2(mod_deli)
plot(mod_deli)
#
#
#
cat("\\newpage")
#
#
#
#
bayes_R2(mod_guich)
plot(mod_guich)
#
#
#
cat("\\newpage")
#
#
#
#
#| label: tbl-data
#| tbl-cap: "Estimates of Reversal learning slope for all the different treatments per each task, species, and group. Mean shows the arithmetic means of the estimates obtained from the posteriors of the model, and 95% CI indicates the 95% confidence interval of the mean. All pmcmc tested the hypothesis that the mean equals zero. In bold, those values that are significant (pmcmc <0.05)"
source(here("R", "func.R"))
#
############################## CREATING BIG DF FOR TABLE ##############################
# Building the vectors for titles of rows and columns
specie <- c("L. delicata", "L. guichenoti")
treatments <- c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")
values <- c("Mean", "95% CI", "pmcmc")
# Building the vectors for estimated means, co.intervals(95%), and pmcmc for the slopes obtained from posteriors, assuming a two-tailed test that testes the hypothesis that the value (slopes in this case) is 0.
#First get estimates for both tasks
estimates <- list(
  deli_CORTCold, deli_ControlCold, deli_CORTHot, deli_ControlHot, 
  guich_CORTCold, guich_ControlCold, guich_CORTHot, guich_ControlHot
)
#
# Then get the mean, co.intervals(95%), and pmcmc
mean <- format_dec(sapply(estimates, mean), 3)
interval_025 <- format_dec(sapply(estimates, function(x) quantile(x,0.025)), 3)
interval_975 <- format_dec(sapply(estimates, function(x) quantile(x,0.975)), 3)
intervals <- paste(interval_025, interval_975, sep = " , ")
Pmcmc <- format_p(sapply(estimates, pmcmc), 4)
#
# Building the df
table_df <- data.frame(
  Specie = rep(specie, each = length(treatments)),
  Treatment = rep(rep(treatments, each = 1), times = length(specie)),
  Mean = rep(mean, each = 1),
  CI = rep(intervals, each = 1),
  Pmcmc = rep(Pmcmc, each = 1)
)
#
write.csv(table_df, file= "./output/Checking/df.csv")
#
############################## ADDING SAMPLE SIZE TO DF FOR TABLE ##############################
# Make n_list into a df
n_df <- as.data.frame(do.call(rbind, n_list)) %>%
  rename("n" = V1) %>%
  rownames_to_column("model") %>%
  separate(model, into = c("Specie", "cort", "temp"), sep = "_") %>%
  unite("Treatment", c("cort", "temp"), sep = "-") %>%
  mutate(Specie = factor(Specie,
                  labels = c(delicata = "L. delicata", guichenoti = "L. guichenoti")),
        Treatment = factor(Treatment,
                   levels = c("CORT-Cold", "Control-Cold", "CORT-Hot","Control-Hot")))
# Estimate total sample size per species and number of observations (this will allow to include this data in the first column of the table)
n_deli <- n_list$delicata_CORT_Cold+n_list$delicata_Control_Cold+n_list$delicata_CORT_Hot+n_list$delicata_Control_Hot
ob_deli <- sum(clean_df$species == "delicata", na.rm = TRUE)
n_guich <- n_list$guichenoti_CORT_Cold+n_list$guichenoti_Control_Cold+n_list$guichenoti_CORT_Hot+n_list$guichenoti_Control_Hot
ob_guich <- sum(clean_df$species == "guichenoti", na.rm = TRUE)
# Merge both dfs, put sample size together with the treatment, and organize the new df to make it look like the table
new_table_df <- merge(table_df, n_df) %>%
  rename('pmcmc' = 'Pmcmc', '95% CI' = 'CI') %>% #Change the names of the columns for the table
  select(Specie, Treatment, Mean, `95% CI`, `pmcmc`, n) %>% #To order the columns in the way I want for the table
  mutate(Specie = factor(Specie,
                  levels = c("L. delicata", "L. guichenoti")),
        Treatment = factor(Treatment, 
                  levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")))%>%
  arrange(Specie, Treatment) %>% # To arrange the rows the way I want
  unite("Treatment", c("Treatment", "n"), sep = " (n = ") %>%
  mutate(Treatment = paste0(Treatment, ")"))
write.csv(new_table_df, file= "./output/Checking/new_table_data.csv")
#
############################## MAKING THE TABLE ##############################
## Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 fint.size = 10)
# Split the table_data df by task
real_table <- flextable(new_table_df) %>%
    bold(~ `pmcmc` < 0.05, ~ `pmcmc` + Mean + `95% CI`) %>%
    set_table_properties(width = 1) %>%
    align(align="center", part="all") %>% 
    set_header_labels(Mean = "Mean",
                      `95% CI` = "95% CI",
                      `pmcmc` = "pmcmc") %>%
    italic(i= c(1, 5), j = 1, italic = TRUE, part = "body") %>% # To have names of species in italics
    flextable::compose(i = c(2:4,6:8), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
    flextable::compose(i = 2, j = 1, value = as_paragraph(sprintf("N = %d", n_deli)), part = "body") %>% # To write the delicata sample size
    flextable::compose(i = 6, j = 1, value = as_paragraph(sprintf("N = %d", n_guich)), part = "body") %>% # To write the guichenoti sample size
    flextable::compose(i = 3, j = 1, value = as_paragraph(sprintf("Obs = %d", ob_deli)), part = "body") %>% # To write the delicata total observations
    flextable::compose(i = 7, j = 1, value = as_paragraph(sprintf("Obs = %d", ob_guich)), part = "body") %>% # To write the guichenoti total observations
    hline(i = 4, j = c(1:5), part = "body") %>% # To make some horizontal lines
    autofit() 
real_table
#
#
#
#
#
#| label: cleandata2
# Cleaning the data to select only those individuals that passed a learning criterion of 80% correct choices in the last 10 trials of the associative task.
data <- read.csv(here("data/Learning.csv"))
percentage_per_lizard <- data %>%
  filter(trial_associative > 25 & trial_associative <= 35) %>%
  group_by(lizard_id) %>%
  summarise(percentage = mean(FC_associative == 1))
# Filter lizards where the percentage is at least 80%
filtered_lizards <- percentage_per_lizard %>%
  filter(percentage >= 0.8)
# Subset the original dataset based on the filtered lizard IDs
clean_df_2 <- data %>%
  filter(lizard_id %in% filtered_lizards$lizard_id) %>%
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
#
#
#
#| label: sampleSize2
# List with the sample sizes from the database (here("output/databases_clean/data_asso.csv") as the sample size per species and group is the same on each task. We used function sample (see func.R) to estimate the sample size per treatment and species.
source(here("R", "func.R"))
#
specie <- c("delicata", "guichenoti")
hormone <- c("CORT", "Control")
temperature <- c("Cold", "Hot")
#
n_list_2 <- list()
for(i in 1:length(specie)){
  for(k in 1:length(hormone)){
    for(l in 1:length(temperature)){
      n_list_2[[paste(specie[i], hormone[k], temperature[l], sep = "_")]] <- sample(clean_df_2, specie[i], hormone[k], temperature[l])
      }
    }
  }
#
#
#
#| label: models2
# Fitting the model and extraction of posteriors for both types of task and species using fit_m function (see func.r in R folder). The result everytime the function is used is a df with the posteriors of the model. The functions saves the model automatically in output/models; and when the parameter refit = FALSE then the posteriors are extracted from the model previously written instead of fitting the model again each time.
source(here("R", "func.R"))
# Model formula: FC_reversal ~ trial_reversal*cort*temp + (1 + trial_reversal|lizard_id)
## A) L. delicata
deli_2 <- fit_m(clean_df_2, "deli", "suppl", refit = FALSE)
write.csv(deli_2, file= "./output/Checking/deli.csv")
## B) L. guichenoti
guich_2 <- fit_m(clean_df_2, "guich", "suppl", refit = FALSE)
write.csv(guich_2, file= "./output/Checking/guich.csv")
#
#
#
# Rename some of the posteriors and make new estimates for the learning rate for the Reversal task doing the same thing we did in the chunk above.
## 1) L. delicata
deli_CORTCold_2 <- deli_2$b_trial_reversal
deli_ControlCold_2 <- (deli_2$'b_trial_reversal:cortControl' + deli_2$b_trial_reversal)
deli_CORTHot_2 <- (deli_2$'b_trial_reversal:tempHot' + deli_2$b_trial_reversal)
deli_ControlHot_2 <- (deli_2$'b_trial_reversal:cortControl:tempHot' + deli_2$b_trial_reversal + deli_2$'b_trial_reversal:cortControl' + deli_2$'b_trial_reversal:tempHot')
## 2) L. guichenoti
guich_CORTCold_2 <- guich_2$b_trial_reversal
guich_ControlCold_2 <- (guich_2$'b_trial_reversal:cortControl' + guich_2$b_trial_reversal)
guich_CORTHot_2 <- (guich_2$'b_trial_reversal:tempHot' + guich_2$b_trial_reversal)
guich_ControlHot_2 <- (guich_2$'b_trial_reversal:cortControl:tempHot' + guich_2$b_trial_reversal + guich_2$'b_trial_reversal:cortControl' + guich_2$'b_trial_reversal:tempHot')
#
#
#
#| label: tbl-data_2
#| tbl-cap: "Estimates of Reversal learning slope for all the different treatments per each task, species, and group. Mean shows the arithmetic means of the estimates obtained from the posteriors of the model, and 95% CI indicates the 95% confidence interval of the mean. All pmcmc tested the hypothesis that the mean equals zero. In bold, those values that are significant (pmcmc <0.05)"
source(here("R", "func.R"))
#
############################## CREATING BIG DF FOR TABLE ##############################
# Building the vectors for titles of rows and columns
specie <- c("L. delicata", "L. guichenoti")
treatments <- c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")
values <- c("Mean", "95% CI", "pmcmc")
# Building the vectors for estimated means, co.intervals(95%), and pmcmcm for the slopes obtained from posteriors assuming a two-tailed test that testes the hypothesis that the value (slopes in this case) is 0.
#First get estimates for both tasks
estimates_2 <- list(
  deli_CORTCold_2, deli_ControlCold_2, deli_CORTHot_2, deli_ControlHot_2, 
  guich_CORTCold_2, guich_ControlCold_2, guich_CORTHot_2, guich_ControlHot_2
)
#
# Then get the mean, co.intervals(95%), and pmcmc
mean_2 <- format_dec(sapply(estimates_2, mean), 3)
interval_025_2 <- format_dec(sapply(estimates_2, function(x) quantile(x,0.025)), 3)
interval_975_2 <- format_dec(sapply(estimates_2, function(x) quantile(x,0.975)), 3)
intervals_2 <- paste(interval_025_2, interval_975_2, sep = " , ")
Pmcmc_2 <- format_p(sapply(estimates_2, pmcmc), 4)
#
# Building the df
table_df_2 <- data.frame(
  Specie = rep(specie, each = length(treatments)),
  Treatment = rep(rep(treatments, each = 1), times = length(specie)),
  Mean = rep(mean_2, each = 1),
  CI = rep(intervals_2, each = 1),
  Pmcmc = rep(Pmcmc_2, each = 1)
)
#
write.csv(table_df, file= "./output/Checking/df.csv")
#
############################## ADDING SAMPLE SIZE TO DF FOR TABLE ##############################
# Make n_list into a df
n_df_2 <- as.data.frame(do.call(rbind, n_list_2)) %>%
  rename("n" = V1) %>%
  rownames_to_column("model") %>%
  separate(model, into = c("Specie", "cort", "temp"), sep = "_") %>%
  unite("Treatment", c("cort", "temp"), sep = "-") %>%
  mutate(Specie = factor(Specie,
                  labels = c(delicata = "L. delicata", guichenoti = "L. guichenoti")),
        Treatment = factor(Treatment,
                   levels = c("CORT-Cold", "Control-Cold", "CORT-Hot","Control-Hot")))
# Estimate total sample size per species and number of observations (this will allow to include this data in the first column of the table)
n_deli_2 <- n_list_2$delicata_CORT_Cold+n_list_2$delicata_Control_Cold+n_list_2$delicata_CORT_Hot+n_list_2$delicata_Control_Hot
ob_deli_2 <- sum(clean_df_2$species == "delicata", na.rm = TRUE)
n_guich_2 <- n_list_2$guichenoti_CORT_Cold+n_list_2$guichenoti_Control_Cold+n_list_2$guichenoti_CORT_Hot+n_list_2$guichenoti_Control_Hot
ob_guich_2 <- sum(clean_df_2$species == "guichenoti", na.rm = TRUE)
# Merge both dfs, put sample size together with the treatment, and organize the new df to make it look like the table
new_table_df_2 <- merge(table_df_2, n_df_2) %>%
  rename('pmcmc' = 'Pmcmc', '95% CI' = 'CI') %>% #Change the names of the columns for the table
  select(Specie, Treatment, Mean, `95% CI`, `pmcmc`, n) %>% #To order the columns in the way I want for the table
  mutate(Specie = factor(Specie,
                  levels = c("L. delicata", "L. guichenoti")),
        Treatment = factor(Treatment, 
                  levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")))%>%
  arrange(Specie, Treatment) %>% # To arrange the rows the way I want
  unite("Treatment", c("Treatment", "n"), sep = " (n = ") %>%
  mutate(Treatment = paste0(Treatment, ")"))
write.csv(new_table_df, file= "./output/Checking/new_table_data.csv")
#
############################## MAKING THE TABLE ##############################
## Table format
set_flextable_defaults(
 font.family = "Times New Roman",
 fint.size = 10)
# Split the table_data df by task
real_table_2 <- flextable(new_table_df_2) %>%
    bold(~ `pmcmc` < 0.05, ~ `pmcmc` + Mean + `95% CI`) %>%
    set_table_properties(width = 1) %>%
    align(align="center", part="all") %>% 
    set_header_labels(Mean = "Mean",
                      `95% CI` = "95% CI",
                      `pmcmc` = "pmcmc") %>%
    italic(i = c(1,5), j = 1, italic = TRUE, part = "body") %>% # To have names od species in italics
    flextable::compose(i = c(2:4,6:8), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
    flextable::compose(i = 2, j = 1, value = as_paragraph(sprintf("N = %d", n_deli_2)), part = "body") %>% # To write the delicata sample size
    flextable::compose(i = 6, j = 1, value = as_paragraph(sprintf("N = %d", n_guich_2)), part = "body") %>% # To write the guichenoti sample size
    flextable::compose(i = 3, j = 1, value = as_paragraph(sprintf("Obs = %d", ob_deli_2)), part = "body") %>% # To write the delicata total observations
    flextable::compose(i = 7, j = 1, value = as_paragraph(sprintf("Obs = %d", ob_guich_2)), part = "body") %>% # To write the guichenoti total observations
    hline(i = 4, j = c(1:5), part = "body") %>% # To make some horizontal lines
    autofit() 
real_table_2
#
#
#
