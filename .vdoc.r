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
pacman::p_load(tidyverse, flextable, emmeans, DHARMa, brms, here, ggplot2, lme4, zoo, lmerTest, broom, tidybayes, ggh4x)
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
#| fig.cap: "Experimental design of early environment manipulation (**A**) and learning tasks (**B**). Stages 1-3 indicate the different phases of the habituation process. In the associative and reversal tasks, white lids show the ramps where the food reward was not attainable."

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
deli <- fit_m(clean_df, "deli", refit = FALSE)
write.csv(deli, file= "./output/Checking/deli.csv")
## B) L. guichenoti
guich <- fit_m(clean_df, "guich", refit = FALSE)
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
#| label: tbl-data
#| tbl-cap: "Estimates of Associative learning slope for all the different treatments per each task, species and group. Mean shows the aritmetic mean of the estimates obtained from the posteriors of the model, and 95% CI indicates the 95% confidence interval of the mean. All p-values were obtained using pmcmc and test the hypothesis that the mean is equal to zero. In bold, those values that are significant (p-value <0.05)"
source(here("R", "func.R"))
#
############################## CREATING BIG DF FOR TABLE ##############################
# Building the vectors for titles of rows and columns
specie <- c("L. delicata", "L. guichenoti")
treatments <- c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")
values <- c("Mean", "95% CI", "p-value")
# Building the vectors for estimated means, co.intervals(95%), and p-values for the slopes obtained from posteriors. p-values are obtained using pmcmc function (see func.R), assuming a two-tailed test that testes the hypothesis that the value (slopes in this case) is 0.
#First get estimates for both tasks
estimates <- list(
  deli_CORTCold, deli_ControlCold, deli_CORTHot, deli_ControlHot, 
  guich_CORTCold, guich_ControlCold, guich_CORTHot, guich_ControlHot
)
#
# Then get the mean, co.intervals(95%), and p-values
mean <- format_dec(sapply(estimates, mean), 3)
interval_025 <- format_dec(sapply(estimates, function(x) quantile(x,0.025)), 3)
interval_975 <- format_dec(sapply(estimates, function(x) quantile(x,0.975)), 3)
intervals <- paste(interval_025, interval_975, sep = " , ")
pvalue <- format_dec(sapply(estimates, pmcmc), 3)
#
# Building the df
table_df <- data.frame(
  Specie = rep(specie, each = length(treatments)),
  Treatment = rep(rep(treatments, each = 1), times = length(specie)),
  Mean = rep(mean, each = 1),
  CI = rep(intervals, each = 1),
  PValue = rep(pvalue, each = 1)
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
# Merge both dfs, put sample size together with the treatment, and organize the new df to make it look like the table
new_table_df <- merge(table_df, n_df) %>%
  rename('p-value' = 'PValue', '95% CI' = 'CI') %>% #Change the names of the columns for the table
  select(Specie, Treatment, Mean, `95% CI`, `p-value`, n) %>% #To order the columns in the way I want for the table
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
    bold(~ `p-value` < 0.05, ~ `p-value` + Mean + `95% CI`) %>%
    set_table_properties(width = 1) %>%
    align(align="center", part="all") %>% 
    set_header_labels(Mean = "Mean",
                      `95% CI` = "95% CI",
                      `p-value` = "p-value") %>%
    italic(j = 1, italic = TRUE, part = "body") %>% # To have names od species in italics
    flextable::compose(i = c(2:4,6:8), j = 1, value = as_paragraph(""), part = "body") %>% # To remove some of the values in the first column
    hline(i = 4, j = c(1:5), part = "body") %>% # To make some horizontal lines
    autofit() 
real_table
#
#
#
#
#
#| label: fig-results
#| fig.cap: "Predicted probability of choosing the correct feeder first over trials. The lines represent the mean predicted probability of choosing the correct feeder first, and the shaded areas represent the standard deviation of the mean both obtained by using the slope and intercept estimates from the posterior distributions. The different colours represent the different treatments. The different panels represent the different species and groups."
source(here("R", "func.R"))
# First step, create the dfs for all models
## 1) L. delicata
df_deli <- df_fig(as.data.frame(deli), "L. delicata")
write.csv(df_deli, file= "./output/Checking/df_deli.csv")
## 2) L. guichenoti
### Red
df_guich <- df_fig(as.data.frame(guich), "L. guichenoti")
#
# Second step, merge the df into a big Fig_df
Fig_df <- rbind(df_deli, df_guich) %>%
  mutate(Trial = gsub("X", "", Trial)) %>%
  mutate(Trial = as.numeric(Trial)) %>%
  group_by(Trial, Treatment, Species) %>%
  summarize(
    Mean_Predicted_prob = mean(Value),
    SE_Predicted_prob = sd(Value)
    ) %>%
  ungroup() %>%
  mutate(
    Treatment = factor(Treatment,
                        levels = c("CORT-Cold", "Control-Cold", "CORT-Hot", "Control-Hot")),
    ) %>%
data.frame()
write.csv(Fig_df, file= "./output/Checking/Fig_df.csv")
# Make the plot
figure_results <- plotting(Fig_df)
ggsave("./output/figures/figure_results.png", plot=figure_results, width = 18, height = 20, units = "cm", dpi = 3000)
knitr::include_graphics("./output/figures/figure_results.png")
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
mod_deli <- readRDS(here("output/models/deli.rds"))
resid_deli <- residuals(mod_deli)
## L. guichenoti
mod_guich <- readRDS(here("output/models/guich.rds"))
resid_guich <- residuals(mod_guich)
#
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
#
