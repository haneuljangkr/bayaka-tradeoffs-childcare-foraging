rm(list=ls())

## Load packages  
library(readr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(brms)
library(coda)
library(usdm)
library(rethinking)
library(dagitty)
library(truncnorm)

########################### 
## Data prep for models  ##
###########################
## Load data
d <- read.csv("data_tradeoffs_bayaka_mothers.csv", sep=",")
nrow(d) 

## Check data
length(unique(d$focalID)) # 23 focal women 
length(unique(d$dayID))   # 38 observation days 
length(unique(d$tripID))  # 523 foraging trips
length(unique(d$tripID[d$infant_presence==FALSE])) # 216 foraging trips without focal woman's infant in group
length(unique(d$tripID[d$infant_presence==TRUE]))  # 307 foraging trips with focal woman's infant in group
table(d$locID, useNA="always") # 3 locations (forest/garden/both) 
length(unique(d$groupID)) # 387 foraging groups
length(unique(d$food_items)) # 16 combinations of 7 unique items

# make factors
d$infant_presence <- ifelse(d$infant_presence==FALSE, 0, 1)
d$infant_presence <- as.factor(d$infant_presence)

# standardize continuous variables
d$n_adults <- as.numeric(scale(d$n_adult))
d$n_children <- as.numeric(scale(d$n_earlyC + d$n_middleC + d$n_adolescent))
d$n_females <- as.numeric(scale(d$n_females))
d$n_males <- as.numeric(scale(d$n_males))
d$n_kin <- as.numeric(scale(d$kin))
d$n_nonkin <- as.numeric(scale(d$nonkin))
d$n_girls_EC <- as.numeric(scale(d$earlyC_f))
d$n_girls_MC <- as.numeric(scale(d$middleC_f))
d$n_girls_AD <- as.numeric(scale(d$adolescent_f))
d$n_women <- as.numeric(scale(d$adult_f))

# Filter out missing outcomes
d1 <- d %>% 
  filter(!(is.na(kcal))) %>% 
  filter(!(is.na(trip_duration)))%>% 
  filter(!is.na(food_items)) 
nrow(d1) # 352 trips

# Compute group average of each predictor
d1$focalID = as.factor(d1$focalID)

d1 <- d1 %>%
  group_by(focalID) %>%
  
  # Get focalID-level average of predictors
  summarise(
    avg_prop_infant = mean(infant_presence == 1, na.rm = T),
    avg_num_adults = mean(n_adults, na.rm = T),
    avg_num_children = mean(n_children, na.rm = T),
    avg_num_females = mean(n_females, na.rm = T),
    avg_num_males = mean(n_males, na.rm = T),
    avg_num_kin = mean(n_kin, na.rm = T),
    avg_num_nonkin = mean(n_nonkin, na.rm = T),
    avg_num_girls_EC = mean(n_girls_EC, na.rm = T),
    avg_num_girls_MC = mean(n_girls_MC, na.rm = T),
    avg_num_girls_AD = mean(n_girls_AD, na.rm = T),
    avg_num_women = mean(n_women, na.rm = T)) %>%
  
  # merge with main dataframe
  right_join(d1) %>%
  
  mutate(
    n_adults_c = n_adults - avg_num_adults, # center on focalID average
    n_children_c = n_children - avg_num_children, # center on focalID average
    n_females_c = n_females - avg_num_females, # center on focalID average
    n_males_c = n_males - avg_num_males, # center on focalID average
    n_kin_c = n_kin - avg_num_kin, # center on focalID average
    n_nonkin_c = n_nonkin - avg_num_nonkin, # center on focalID average
    n_girls_EC_c = n_girls_EC - avg_num_girls_EC, # center on focalID average
    n_girls_MC_c = n_girls_MC - avg_num_girls_MC, # center on focalID average
    n_girls_AD_c = n_girls_AD - avg_num_girls_AD, # center on focalID average
    n_women_c = n_women - avg_num_women, # center on focalID average
    
    # scale responses by max value to kmeep in kcal of [0,1], makes it easier to set priors
    dur_s = trip_duration / max(trip_duration, na.rm=T),
    kcal_s = kcal / max(kcal, na.rm=T)
  )


########################### 
# Prior predictive checks #
###########################
# Weakly regularizing priors, understood in light of the data
hist(log(d1$kcal_s))

priors <- c(
  set_prior("normal(-1, 1)", class = "Intercept"),
  set_prior("exponential(5)", class = "sigma"),
  set_prior("normal(0,0.2)", class = "b"),
  set_prior("exponential(5)", class = "sd"),
  set_prior("exponential(5)", class = "sds")
)

m_prior_kcal_age <- brm(kcal_s ~
                           # focalID-level predictors (fixed effects)
                           avg_prop_infant + avg_num_adults + avg_num_children +
                           # trip-level predictors (fixed effects) and spline with an interaction
                           infant_presence + s(n_adults_c, by = infant_presence) + s(n_children_c, by = infant_presence) +
                           # trip-level predictors to control trip duration 
                           s(dur_s) +
                           # trip-level predictors (random effects)
                           (1 + infant_presence*(n_adults_c + n_children_c) | focalID) +
                           # additional random effect structure
                           (1 | dayID) + (1 | food_items),
                         
                         data = d1,
                         family = lognormal(),
                         chains = 4,
                         cores = 4,
                         prior = priors,
                         sample_prior = "only",
                         # backend = "cmdstanr", # optional, but faster
                         init = "0")


m_prior_kcal_kin <- brm(kcal_s ~
                           # focalID-level predictors (fixed effects)
                           avg_prop_infant + avg_num_kin + avg_num_nonkin +
                           # trip-level predictors (fixed effects) and spline with an interaction
                           infant_presence + s(n_kin_c, by = infant_presence) + s(n_nonkin_c, by = infant_presence) +
                           # trip-level predictors to control trip duration 
                           s(dur_s) +
                           # trip-level predictors (random effects)
                           (1 + infant_presence*(n_kin_c + n_nonkin_c) | focalID) +
                           # additional random effect structure
                           (1 | dayID) + (1 | food_items),
                         data = d1,
                         family = lognormal(),
                         chains = 4,
                         cores = 4,
                         prior = priors,
                         sample_prior = "only",
                         init = "0")

m_prior_kcal_sex <- brm(kcal_s ~
                           # focalID-level predictors (fixed effects)
                           avg_prop_infant + avg_num_females + avg_num_males +
                           # trip-level predictors (fixed effects) and spline with an interaction
                           infant_presence + s(n_females_c, by = infant_presence) + s(n_males_c, by = infant_presence) +
                           # trip-level predictors to control trip duration 
                           s(dur_s) +
                           # trip-level predictors (random effects)
                           (1 + infant_presence*(n_females_c + n_males_c) | focalID) +
                           # additional random effect structure
                           (1 | dayID) + (1 | food_items),
                         data = d1,
                         family = lognormal(),
                         chains = 4,
                         cores = 4,
                         prior = priors,
                         sample_prior = "only",
                         init = "0")

# check reasonableness of priors. Thin lines are the prior predictive densities, dark thicker line is the actual data
pp_check(m_prior_kcal_age, ndraws = 50) + 
  scale_x_continuous(limits = c(0,1.5)) +
  theme_bw(base_size = 18)

pp_check(m_prior_kcal_kin, ndraws = 50) + 
  scale_x_continuous(limits = c(0,1.5)) +
  theme_bw(base_size = 18)

pp_check(m_prior_kcal_sex, ndraws = 50) + 
  scale_x_continuous(limits = c(0,1.5)) +
  theme_bw(base_size = 18)

##########################
###  Fit the models  #####
##########################

# Baseline model for comparison
m_kcal_0 <- brm(kcal_s ~ s(dur_s) +
                   (1 | focalID) + (1 | dayID) + (1 | food_items),
                   data = d1,
                   family = lognormal(),
                   chains = 4,
                   cores = 4,
                   prior = c(
                   set_prior("normal(-1, 1)", class = "Intercept"),
                   set_prior("exponential(5)", class = "sigma"),
                   set_prior("exponential(5)", class = "sd")),
                   control = list(adapt_delta = 0.98),
                   init = "0")
summary(m_kcal_0)

# Model with predictors 
m_kcal_age <- brm(kcal_s ~
                     # focalID-level predictors (fixed effects)
                     avg_prop_infant + avg_num_adults + avg_num_children +
                     # trip-level predictors (fixed effects) and spline with an interaction
                     infant_presence + s(n_adults_c, by = infant_presence) + s(n_children_c, by = infant_presence) +
                     # trip-level predictors to control trip duration 
                     s(dur_s) +
                     # trip-level predictors (random effects)
                     (1 + infant_presence*(n_adults_c + n_children_c) | focalID) +
                     # additional random effect structure
                     (1 | dayID) + (1 | food_items),
                   data = d1,
                   family = lognormal(),
                   chains = 4,
                   cores = 4,
                   prior = priors,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.98),
                   init = "0")
summary(m_kcal_age)

m_kcal_kin <- brm(kcal_s ~
                     # focalID-level predictors (fixed effects)
                     avg_prop_infant + avg_num_kin + avg_num_nonkin +
                     # trip-level predictors (fixed effects) and spline with an interaction
                     infant_presence + s(n_kin_c, by = infant_presence) + s(n_nonkin_c, by = infant_presence) +
                     # trip-level predictors to control trip duration 
                     s(dur_s) +
                     # trip-level predictors (random effects)
                     (1 + infant_presence*(n_kin_c + n_nonkin_c) | focalID) +
                     # additional random effect structure
                     (1 | dayID) + (1 | food_items),
                   data = d1,
                   family = lognormal(),
                   chains = 4,
                   cores = 4,
                   prior = priors,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.98),
                   init = "0")
summary(m_kcal_kin)

m_kcal_sex <- brm(kcal_s ~
                     # focalID-level predictors (fixed effects)
                     avg_prop_infant + avg_num_females + avg_num_males +
                     # trip-level predictors (fixed effects) and spline with an interaction
                     infant_presence + s(n_females_c, by = infant_presence) + s(n_males_c, by = infant_presence) +
                     # trip-level predictors to control trip duration 
                     s(dur_s) +
                     # trip-level predictors (random effects)
                     (1 + infant_presence*(n_females_c + n_males_c) | focalID) +
                     # additional random effect structure
                     (1 | dayID) + (1 | food_items),      
                   data = d1,
                   family = lognormal(),
                   chains = 4,
                   cores = 4,
                   prior = priors,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.98),
                   init = "0")
summary(m_kcal_sex)

save.image("kcal_interaction.RData")

##########################
### Model comparison #####
##########################
m_kcal_0 <- add_criterion(m_kcal_0, criterion = "loo")
m_kcal_age <- add_criterion(m_kcal_age, criterion = "loo") # moment_match helps with influential observations, omitting that step for no
m_kcal_kin <- add_criterion(m_kcal_kin, criterion = "loo") # moment_match helps with influential observations, omitting that step for no
m_kcal_sex <- add_criterion(m_kcal_sex, criterion = "loo") # moment_match helps with influential observations, omitting that step for no

loo::loo_compare(loo(m_kcal_0), loo(m_kcal_age))
loo::loo_compare(loo(m_kcal_0), loo(m_kcal_kin))
loo::loo_compare(loo(m_kcal_0), loo(m_kcal_sex))

# Model weights interpreted as the probability that model will perform best on new data (compared to other models being compared)
loo::loo_model_weights(list(m0 = loo(m_kcal_0), m_age = loo(m_kcal_age)))
loo::loo_model_weights(list(m0 = loo(m_kcal_0), m_kin = loo(m_kcal_kin)))
loo::loo_model_weights(list(m0 = loo(m_kcal_0), m_sex = loo(m_kcal_sex)))

###########################################
###### Posterior predictive checks ########
###########################################
#The dark line is the observed distribution whereas the light lines are draws from the posterior.
pp_check(m_kcal_age, type = "dens_overlay") + scale_x_continuous(limits = c(0,1)) # I put a limit on the x-axis because we scaled the data [0,1].
pp_check(m_kcal_kin, type = "dens_overlay") + scale_x_continuous(limits = c(0,1)) # I put a limit on the x-axis because we scaled the data [0,1].
pp_check(m_kcal_sex, type = "dens_overlay") + scale_x_continuous(limits = c(0,1)) # I put a limit on the x-axis because we scaled the data [0,1].

# Empirical cumulative density function
pp_check(m_kcal_age, type = "ecdf_overlay") + scale_x_continuous(limits = c(0,1)) # I put a limit on the x-axis because we scaled the data [0,1].
pp_check(m_kcal_kin, type = "ecdf_overlay") + scale_x_continuous(limits = c(0,1)) # I put a limit on the x-axis because we scaled the data [0,1].
pp_check(m_kcal_sex, type = "ecdf_overlay") + scale_x_continuous(limits = c(0,1)) # I put a limit on the x-axis because we scaled the data [0,1].

# We can also check the model predictions by group witin the data, for example whether infants were present or not
pp_check(m_kcal_age, type = "dens_overlay_grouped", group = "infant_presence") + scale_x_continuous(limits = c(0,1))
pp_check(m_kcal_kin, type = "dens_overlay_grouped", group = "infant_presence") + scale_x_continuous(limits = c(0,1))
pp_check(m_kcal_sex, type = "dens_overlay_grouped", group = "infant_presence") + scale_x_continuous(limits = c(0,1))

# check for focal Id 
pp_check(m_kcal_age, type = "dens_overlay_grouped", group = "focalID") + scale_x_continuous(limits = c(0,1))

# ECDF plots for each focal woman
pp_check(m_kcal_age, type = "ecdf_overlay_grouped", group = "focalID") + scale_x_continuous(limits = c(0,1))

###########################
### Predictive plots #####
###########################

#################################################################
# (1) Predict effect of more adults, depending on infant presence
newdata_adults_kcal <- expand.grid(
  n_adults_c = seq(from = -2, to = 2, length.out = 30),
  infant_presence = unique(d1$infant_presence),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_children_c = 0
) %>% 
  mutate(condition = 1:n())

pred_adults <- posterior_epred(m_kcal_age, newdata = newdata_adults_kcal, re.form = NA)

pred_adults_long <- as.data.frame(pred_adults) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_adults_kcal) %>% 
  # put back onto original scale
  mutate(pred_kcal = pred * max(d1$kcal, na.rm=T))

pred_adults_long$infant_presence <- ifelse(pred_adults_long$infant_presence == 0, "infant absent", "infant present")

pred_adults_summary <- pred_adults_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_kcal), infant_presence = unique(infant_presence), n_adults_c = unique(n_adults_c))

# (2) Predict effect of more children, depending on infant presence
newdata_children_kcal <- expand.grid(
  n_children_c = seq(from = -2, to = 2, length.out = 30),
  infant_presence = unique(d1$infant_presence),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_adults_c = 0
) %>% 
  mutate(condition = 1:n())

pred_children <- posterior_epred(m_kcal_age, newdata = newdata_children_kcal, re.form = NA)

pred_children_long <- as.data.frame(pred_children) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_children_kcal) %>% 
  # put back onto original scale
  mutate(pred_kcal = pred * max(d1$kcal, na.rm=T))

pred_children_long$infant_presence <- ifelse(pred_children_long$infant_presence == 0, "infant absent", "infant present")

pred_children_summary <- pred_children_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_kcal), infant_presence = unique(infant_presence), n_children_c = unique(n_children_c))

###### OK, lets plot both #########
library(patchwork)

p_kcal_adults <- ggplot( data = filter(pred_adults_long, samp <= 100), aes(x = n_adults_c, y = pred_kcal, color = infant_presence, group = factor(samp)) ) +
  facet_wrap(~infant_presence) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_adults_summary, aes(x = n_adults_c, y = med, color = infant_presence, group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab("Food returns (kcal)") +
  xlab("Number of adults (z-score)") +
  theme(legend.position = "none", panel.spacing = unit(2, "lines"))

p_kcal_children <- ggplot( data = filter(pred_children_long, samp <= 100), aes(x = n_children_c, y = pred_kcal, color = infant_presence, group = factor(samp)) ) +
  facet_wrap(~infant_presence) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_children_summary, aes(x = n_children_c, y = med, color = infant_presence, group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab("Food returns (kcal)") +
  xlab("Number of children (z-score)") +
  theme(legend.position = "none", panel.spacing = unit(2, "lines"))

p_kcal_adults / p_kcal_children

# We can also plot the difference betwkcaln conditions (infant present/absent)
pred_adults_diff <- pred_adults_long %>% 
  group_by(samp, n_adults_c) %>% 
  summarise(diff = pred_kcal[infant_presence == "infant present"] - pred_kcal[infant_presence == "infant absent"]) 

pred_adults_diff_summary <- pred_adults_diff %>% group_by(n_adults_c) %>% 
  summarise(med = median(diff))

pred_children_diff <- pred_children_long %>% 
  group_by(samp, n_children_c) %>% 
  summarise(diff = pred_kcal[infant_presence == "infant present"] - pred_kcal[infant_presence == "infant absent"]) 

pred_children_diff_summary <- pred_children_diff %>% group_by(n_children_c) %>% 
  summarise(med = median(diff))

p_diff_adults <- ggplot( data = filter(pred_adults_diff, samp <= 100), aes(x = n_adults_c, y = diff, group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_adults_diff_summary, aes(x = n_adults_c, y = med, group = NULL), lwd = 2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", col = 'red', lwd = 1) +
  theme_minimal(base_size = 18) +
  ylab("Diff (kcal) [Infant present - Infant absent]") +
  xlab("Number of adults (z-score)") +
  theme(legend.position = "none", panel.spacing = unit(2, "lines"))

p_diff_children <- ggplot( data = filter(pred_children_diff, samp <= 100), aes(x = n_children_c, y = diff, group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_children_diff_summary, aes(x = n_children_c, y = med, group = NULL), lwd = 2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", col = 'red', lwd = 1) +
  theme_minimal(base_size = 18) +
  ylab("Diff (kcal) [Infant present - Infant absent]") +
  xlab("Number of children (z-score)") +
  theme(legend.position = "none", panel.spacing = unit(2, "lines"))

p_diff_adults / p_diff_children

#################################################################
# (2) Predict effect of more kin, depending on infant presence
newdata_kin_kcal <- expand.grid(
  n_kin_c = seq(from = -2, to = 2, length.out = 30),
  infant_presence = unique(d1$infant_presence),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_nonkin_c = 0
) %>% 
  mutate(condition = 1:n())

pred_kin <- posterior_epred(m_kcal_kin, newdata = newdata_kin_kcal, re.form = NA)

pred_kin_long <- as.data.frame(pred_kin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_kin_kcal) %>% 
  # put back onto original scale
  mutate(pred_kcal = pred * max(d1$kcal, na.rm=T))

pred_kin_long$infant_presence <- ifelse(pred_kin_long$infant_presence == 0, "infant absent", "infant present")

pred_kin_summary <- pred_kin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_kcal), infant_presence = unique(infant_presence), n_kin_c = unique(n_kin_c))

# (2) Predict effect of more nonkin, depending on infant presence
newdata_nonkin_kcal <- expand.grid(
  n_nonkin_c = seq(from = -2, to = 2, length.out = 30),
  infant_presence = unique(d1$infant_presence),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_kin_c = 0
) %>% 
  mutate(condition = 1:n())

pred_nonkin <- posterior_epred(m_kcal_kin, newdata = newdata_nonkin_kcal, re.form = NA)

pred_nonkin_long <- as.data.frame(pred_nonkin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_nonkin_kcal) %>% 
  # put back onto original scale
  mutate(pred_kcal = pred * max(d1$kcal, na.rm=T))

pred_nonkin_long$infant_presence <- ifelse(pred_nonkin_long$infant_presence == 0, "infant absent", "infant present")

pred_nonkin_summary <- pred_nonkin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_kcal), infant_presence = unique(infant_presence), n_nonkin_c = unique(n_nonkin_c))

###### OK, lets plot both #########
library(patchwork)

p_kcal_kin <- ggplot( data = filter(pred_kin_long, samp <= 100), aes(x = n_kin_c, y = pred_kcal, color = infant_presence, group = factor(samp)) ) +
  facet_wrap(~infant_presence) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_kin_summary, aes(x = n_kin_c, y = med, color = infant_presence, group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab("Food returns (kcal)") +
  xlab("Number of kin (z-score)") +
  theme(legend.position = "none", panel.spacing = unit(2, "lines"))

p_kcal_nonkin <- ggplot( data = filter(pred_nonkin_long, samp <= 100), aes(x = n_nonkin_c, y = pred_kcal, color = infant_presence, group = factor(samp)) ) +
  facet_wrap(~infant_presence) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_nonkin_summary, aes(x = n_nonkin_c, y = med, color = infant_presence, group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab("Food returns (kcal)") +
  xlab("Number of nonkin (z-score)") +
  theme(legend.position = "none", panel.spacing = unit(2, "lines"))

p_kcal_kin / p_kcal_nonkin

# We can also plot the difference betwkcaln conditions (infant present/absent)
pred_kin_diff <- pred_kin_long %>% 
  group_by(samp, n_kin_c) %>% 
  summarise(diff = pred_kcal[infant_presence == "infant present"] - pred_kcal[infant_presence == "infant absent"]) 

pred_kin_diff_summary <- pred_kin_diff %>% group_by(n_kin_c) %>% 
  summarise(med = median(diff))

pred_nonkin_diff <- pred_nonkin_long %>% 
  group_by(samp, n_nonkin_c) %>% 
  summarise(diff = pred_kcal[infant_presence == "infant present"] - pred_kcal[infant_presence == "infant absent"]) 

pred_nonkin_diff_summary <- pred_nonkin_diff %>% group_by(n_nonkin_c) %>% 
  summarise(med = median(diff))

p_diff_kin <- ggplot( data = filter(pred_kin_diff, samp <= 100), aes(x = n_kin_c, y = diff, group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_kin_diff_summary, aes(x = n_kin_c, y = med, group = NULL), lwd = 2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", col = 'red', lwd = 1) +
  theme_minimal(base_size = 18) +
  ylab("Diff (kcal) [Infant present - Infant absent]") +
  xlab("Number of kin (z-score)") +
  theme(legend.position = "none", panel.spacing = unit(2, "lines"))

p_diff_nonkin <- ggplot( data = filter(pred_nonkin_diff, samp <= 100), aes(x = n_nonkin_c, y = diff, group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_nonkin_diff_summary, aes(x = n_nonkin_c, y = med, group = NULL), lwd = 2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", col = 'red', lwd = 1) +
  theme_minimal(base_size = 18) +
  ylab("Diff (kcal) [Infant present - Infant absent]") +
  xlab("Number of nonkin (z-score)") +
  theme(legend.position = "none", panel.spacing = unit(2, "lines"))

p_diff_kin / p_diff_nonkin

#################################################################
# (3) Predict effect of more females, depending on infant presence
newdata_females_kcal <- expand.grid(
  n_females_c = seq(from = -2, to = 2, length.out = 30),
  infant_presence = unique(d1$infant_presence),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_males_c = 0
) %>% 
  mutate(condition = 1:n())

pred_females <- posterior_epred(m_kcal_sex, newdata = newdata_females_kcal, re.form = NA)

pred_females_long <- as.data.frame(pred_females) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_females_kcal) %>% 
  # put back onto original scale
  mutate(pred_kcal = pred * max(d1$kcal, na.rm=T))

pred_females_long$infant_presence <- ifelse(pred_females_long$infant_presence == 0, "infant absent", "infant present")

pred_females_summary <- pred_females_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_kcal), infant_presence = unique(infant_presence), n_females_c = unique(n_females_c))

# (2) Predict effect of more nonkin, depending on infant presence
newdata_males_kcal <- expand.grid(
  n_males_c = seq(from = -2, to = 2, length.out = 30),
  infant_presence = unique(d1$infant_presence),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_females_c = 0
) %>% 
  mutate(condition = 1:n())

pred_males <- posterior_epred(m_kcal_sex, newdata = newdata_males_kcal, re.form = NA)

pred_males_long <- as.data.frame(pred_males) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_males_kcal) %>% 
  # put back onto original scale
  mutate(pred_kcal = pred * max(d1$kcal, na.rm=T))

pred_males_long$infant_presence <- ifelse(pred_males_long$infant_presence == 0, "infant absent", "infant present")

pred_males_summary <- pred_males_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_kcal), infant_presence = unique(infant_presence), n_males_c = unique(n_males_c))

###### OK, lets plot both #########
library(patchwork)

p_kcal_females <- ggplot( data = filter(pred_females_long, samp <= 100), aes(x = n_females_c, y = pred_kcal, color = infant_presence, group = factor(samp)) ) +
  facet_wrap(~infant_presence) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_females_summary, aes(x = n_females_c, y = med, color = infant_presence, group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab("Food returns (kcal)") +
  xlab("Number of females (z-score)") +
  theme(legend.position = "none", panel.spacing = unit(2, "lines"))

p_kcal_males <- ggplot( data = filter(pred_males_long, samp <= 100), aes(x = n_males_c, y = pred_kcal, color = infant_presence, group = factor(samp)) ) +
  facet_wrap(~infant_presence) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_males_summary, aes(x = n_males_c, y = med, color = infant_presence, group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab("Food returns (kcal)") +
  xlab("Number of males (z-score)") +
  theme(legend.position = "none", panel.spacing = unit(2, "lines"))

p_kcal_females / p_kcal_males

# We can also plot the difference betwkcaln conditions (infant present/absent)
pred_females_diff <- pred_females_long %>% 
  group_by(samp, n_females_c) %>% 
  summarise(diff = pred_kcal[infant_presence == "infant present"] - pred_kcal[infant_presence == "infant absent"]) 

pred_females_diff_summary <- pred_females_diff %>% group_by(n_females_c) %>% 
  summarise(med = median(diff))

pred_males_diff <- pred_males_long %>% 
  group_by(samp, n_males_c) %>% 
  summarise(diff = pred_kcal[infant_presence == "infant present"] - pred_kcal[infant_presence == "infant absent"]) 

pred_males_diff_summary <- pred_males_diff %>% group_by(n_males_c) %>% 
  summarise(med = median(diff))

p_diff_females <- ggplot( data = filter(pred_females_diff, samp <= 100), aes(x = n_females_c, y = diff, group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_females_diff_summary, aes(x = n_females_c, y = med, group = NULL), lwd = 2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", col = 'red', lwd = 1) +
  theme_minimal(base_size = 18) +
  ylab("Diff (kcal) [Infant present - Infant absent]") +
  xlab("Number of females (z-score)") +
  theme(legend.position = "none", panel.spacing = unit(2, "lines"))

p_diff_males <- ggplot( data = filter(pred_males_diff, samp <= 100), aes(x = n_males_c, y = diff, group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_males_diff_summary, aes(x = n_males_c, y = med, group = NULL), lwd = 2) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", col = 'red', lwd = 1) +
  theme_minimal(base_size = 18) +
  ylab("Diff (kcal) [Infant present - Infant absent]") +
  xlab("Number of males (z-score)") +
  theme(legend.position = "none", panel.spacing = unit(2, "lines"))

p_diff_females / p_diff_males
