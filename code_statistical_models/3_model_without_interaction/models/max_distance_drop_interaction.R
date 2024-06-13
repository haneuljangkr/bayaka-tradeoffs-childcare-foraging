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

# Load data
d <- read.csv("data_tradeoffs_bayaka_mothers.csv", sep=",")
nrow(d) 

## Check data
length(unique(d$focalID)) # 22 focal women 
length(unique(d$dayID))   # 38 observation days 
length(unique(d$tripID))  # 503 foraging trips
length(unique(d$tripID[d$infant_presence==FALSE])) # 199 foraging trips without focal woman's infant in group
length(unique(d$tripID[d$infant_presence==TRUE]))  # 304 foraging trips with focal woman's infant in group
table(d$locID, useNA="always") # 3 locations (forest/garden/both) 
length(unique(d$groupID)) # 381 foraging groups
length(unique(d$food_items)) # 16 combinations of 7 unique items

# make factors
d$infant_presence <- ifelse(d$infant_presence==FALSE, 0, 1)
d$infant_presence <- as.factor(d$infant_presence)

# standardize continuous variables
d$n_adults_s <- as.numeric(scale(d$n_adult))
d$n_children <- d$n_earlyC + d$n_middleC + d$n_adolescent
d$n_children_s <- as.numeric(scale(d$n_children))
d$n_females_s <- as.numeric(scale(d$n_females))
d$n_males_s <- as.numeric(scale(d$n_males))
d$n_kin_s <- as.numeric(scale(d$kin))
d$n_nonkin_s <- as.numeric(scale(d$nonkin))
d$n_girls_EC_s <- as.numeric(scale(d$earlyC_f))
d$n_girls_MC_s <- as.numeric(scale(d$middleC_f))
d$n_girls_AD_s <- as.numeric(scale(d$adolescent_f))
d$n_women_s <- as.numeric(scale(d$adult_f))
d$n_nursing_s <- as.numeric(scale(d$n_nursing))
d$n_earlyC_s <- as.numeric(scale(d$n_earlyC))
d$n_middleC_s <- as.numeric(scale(d$n_middleC))
d$n_adolescent_s <- as.numeric(scale(d$n_adolescent))

# Filter out missing outcomes
d1 <- d %>% 
  filter(!(is.na(trip_duration)))

d1 <- d1 %>% 
  filter(!is.na(food_items)) 
nrow(d1) # 348 trips

# Compute group average of each predictor
d1$focalID = as.factor(d1$focalID)

d1 <- d1 %>%
  group_by(focalID) %>%
  
  # Get focalID-level average of predictors
  summarise(
    avg_prop_infant = mean(infant_presence == 1, na.rm = T),
    avg_num_adults = mean(n_adults_s, na.rm = T),
    avg_num_children = mean(n_children_s, na.rm = T),
    avg_num_females = mean(n_females_s, na.rm = T),
    avg_num_males = mean(n_males_s, na.rm = T),
    avg_num_kin = mean(n_kin_s, na.rm = T),
    avg_num_nonkin = mean(n_nonkin_s, na.rm = T),
    avg_num_girls_EC = mean(n_girls_EC_s, na.rm = T),
    avg_num_girls_MC = mean(n_girls_MC_s, na.rm = T),
    avg_num_girls_AD = mean(n_girls_AD_s, na.rm = T),
    avg_num_women = mean(n_women_s, na.rm = T),
    avg_num_nursing = mean(n_nursing_s, na.rm = T),
    avg_num_earlyC = mean(n_earlyC_s, na.rm = T),
    avg_num_middleC = mean(n_middleC_s, na.rm = T),
    avg_num_adolescent = mean(n_adolescent_s, na.rm = T))  %>%
  
  # merge with main dataframe
  right_join(d1) %>%
  
  mutate(
    n_adults_c = n_adults_s - avg_num_adults, # center on focalID average
    n_children_c = n_children_s - avg_num_children, # center on focalID average
    n_females_c = n_females_s - avg_num_females, # center on focalID average
    n_males_c = n_males_s - avg_num_males, # center on focalID average
    n_kin_c = n_kin_s - avg_num_kin, # center on focalID average
    n_nonkin_c = n_nonkin_s - avg_num_nonkin, # center on focalID average
    n_girls_EC_c = n_girls_EC_s - avg_num_girls_EC, # center on focalID average
    n_girls_MC_c = n_girls_MC_s - avg_num_girls_MC, # center on focalID average
    n_girls_AD_c = n_girls_AD_s - avg_num_girls_AD, # center on focalID average
    n_women_c = n_women_s - avg_num_women, # center on focalID average
    n_nursing_c = n_nursing_s - avg_num_nursing,
    n_earlyC_c = n_earlyC_s - avg_num_earlyC,
    n_middleC_c = n_middleC_s - avg_num_middleC,
    n_adolescent_c = n_adolescent_s - avg_num_adolescent,
    
    # scale responses by max value to keep in range of [0,1], makes it easier to set priors
    max_dist_s = maximum_distance / max(maximum_distance, na.rm=T),
  )

########################### 
# Prior predictive checks #
###########################
# Weakly regularizing priors, understood in light of the data
hist(log(d1$max_dist_s))

priors <- c(
  set_prior("normal(-1, 1)", class = "Intercept"),
  set_prior("exponential(5)", class = "sigma"),
  set_prior("normal(0,0.2)", class = "b"),
  set_prior("exponential(5)", class = "sd"))

m1_prior_max_age <- brm(max_dist_s ~
                        # focalID-level predictors (fixed effects)
                        avg_prop_infant + avg_num_adults + avg_num_children +
                        # trip-level predictors (fixed effects) 
                        infant_presence + n_adults_c + n_children_c +
                        # trip-level predictors (random effects)
                        (1 + infant_presence + n_adults_c + n_children_c | focalID) +
                        # additional random effect structure
                        (1 | dayID) + (1 | food_items),
                      data = d1,
                      family = lognormal(),
                      chains = 4,
                      cores = 4,
                      prior = priors,
                      sample_prior = "only",
                      init = "0")


m1_prior_max_kin <- brm(max_dist_s ~
                        # focalID-level predictors (fixed effects)
                        avg_prop_infant + avg_num_kin + avg_num_nonkin +
                        # trip-level predictors (fixed effects) 
                        infant_presence + n_kin_c + n_nonkin_c +
                        # trip-level predictors (random effects)
                        (1 + infant_presence + n_kin_c + n_nonkin_c | focalID) +
                        # additional random effect structure
                         (1 | dayID) + (1 | food_items),
                      data = d1,
                      family = lognormal(),
                      chains = 4,
                      cores = 4,
                      prior = priors,
                      sample_prior = "only",
                      init = "0")

m1_prior_max_sex <- brm(max_dist_s ~
                        # focalID-level predictors (fixed effects)
                        avg_prop_infant + avg_num_females + avg_num_males +
                        # trip-level predictors (fixed effects) 
                        infant_presence + n_females_c + n_males_c +
                        # trip-level predictors (random effects)
                        (1 + infant_presence + n_females_c + n_males_c | focalID) +
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
pp_check(m1_prior_max_age, ndraws = 50) + 
  scale_x_continuous(limits = c(0,1.5)) +
  theme_bw(base_size = 18)

pp_check(m1_prior_max_kin, ndraws = 50) + 
  scale_x_continuous(limits = c(0,1.5)) +
  theme_bw(base_size = 18)

pp_check(m1_prior_max_sex, ndraws = 50) + 
  scale_x_continuous(limits = c(0,1.5)) +
  theme_bw(base_size = 18)

##########################
###  Fit the models  #####
##########################

# Baseline model for comparison
m1_max_0 <- brm(max_dist_s ~
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
summary(m1_max_0)

# Model with predictors 
m1_max_age <- brm(max_dist_s ~
                  # focalID-level predictors (fixed effects)
                  avg_prop_infant + avg_num_adults + avg_num_children +
                  # trip-level predictors (fixed effects) 
                  infant_presence + n_adults_c + n_children_c +
                  # trip-level predictors (random effects)
                  (1 + infant_presence + n_adults_c + n_children_c | focalID) +
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
summary(m1_max_age)

m1_max_kin <- brm(max_dist_s ~
                  # focalID-level predictors (fixed effects)
                  avg_prop_infant + avg_num_kin + avg_num_nonkin +
                  # trip-level predictors (fixed effects) 
                  infant_presence + n_kin_c + n_nonkin_c +
                  # trip-level predictors (random effects)
                  (1 + infant_presence + n_kin_c + n_nonkin_c | focalID) +
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
summary(m1_max_kin)

m1_max_sex <- brm(max_dist_s ~
                  # focalID-level predictors (fixed effects)
                  avg_prop_infant + avg_num_females + avg_num_males +
                  # trip-level predictors (fixed effects) 
                  infant_presence + n_females_c + n_males_c +
                  # trip-level predictors (random effects)
                  (1 + infant_presence + n_females_c + n_males_c | focalID) +
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
summary(m1_max_sex)

m1_max_female_age <- brm(max_dist_s ~
                            # focalID-level predictors (fixed effects)
                            avg_prop_infant + avg_num_girls_EC + avg_num_girls_MC + avg_num_girls_AD + avg_num_women +
                            # trip-level predictors (fixed effects)
                            infant_presence + n_girls_EC_c + n_girls_MC_c + n_girls_AD_c + n_women_c +
                            # trip-level predictors (random effects)
                            (1 + infant_presence + n_girls_EC_c + n_girls_MC_c  + n_girls_AD_c + n_women_c | focalID) +
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
summary(m1_max_female_age)

m1_max_nursing <- brm(max_dist_s ~
                         # focalID-level predictors (fixed effects)
                         avg_prop_infant + avg_num_nursing + 
                         # trip-level predictors (fixed effects)
                         infant_presence + n_nursing_c + 
                         # trip-level predictors (random effects)
                         (1 + infant_presence + n_nursing_c | focalID) +
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
summary(m1_max_nursing)

save.image(file="max_distance_drop_interaction.RData")

##########################
### Model comparison #####
##########################
m1_max_0 <- add_criterion(m1_max_0, criterion = "loo")
m1_max_age <- add_criterion(m1_max_age, criterion = "loo") # moment_match helps with influential observations, omitting that step for no
m1_max_kin <- add_criterion(m1_max_kin, criterion = "loo") # moment_match helps with influential observations, omitting that step for no
m1_max_sex <- add_criterion(m1_max_sex, criterion = "loo") # moment_match helps with influential observations, omitting that step for no

loo::loo_compare(loo(m1_max_0), loo(m1_max_age))
loo::loo_compare(loo(m1_max_0), loo(m1_max_kin))
loo::loo_compare(loo(m1_max_0), loo(m1_max_sex))

# Model weights interpreted as the probability that model will perform best on new data (compared to other models being compared)
loo::loo_model_weights(list(m0 = loo(m1_max_0), m_age = loo(m1_max_age)))
loo::loo_model_weights(list(m0 = loo(m1_max_0), m_kin = loo(m1_max_kin)))
loo::loo_model_weights(list(m0 = loo(m1_max_0), m_sex = loo(m1_max_sex)))

###########################################
###### Posterior predictive checks ########
###########################################
#The dark line is the observed distribution whereas the light lines are draws from the posterior.
pp_check(m1_max_age, type = "dens_overlay") + scale_x_continuous(limits = c(0,1)) # I put a limit on the x-axis because we scaled the data [0,1].
pp_check(m1_max_kin, type = "dens_overlay") + scale_x_continuous(limits = c(0,1)) # I put a limit on the x-axis because we scaled the data [0,1].
pp_check(m1_max_sex, type = "dens_overlay") + scale_x_continuous(limits = c(0,1)) # I put a limit on the x-axis because we scaled the data [0,1].

# Empirical cumulative density function
pp_check(m1_max_age, type = "ecdf_overlay") + scale_x_continuous(limits = c(0,1)) # I put a limit on the x-axis because we scaled the data [0,1].
pp_check(m1_max_kin, type = "ecdf_overlay") + scale_x_continuous(limits = c(0,1)) # I put a limit on the x-axis because we scaled the data [0,1].
pp_check(m1_max_sex, type = "ecdf_overlay") + scale_x_continuous(limits = c(0,1)) # I put a limit on the x-axis because we scaled the data [0,1].

# per each focal woman
pp_check(m1_max_age, type = "dens_overlay_grouped", group = "focalID") + scale_x_continuous(limits = c(0,1))

# ECDF plots for each focal woman 
pp_check(m1_max_age, type = "ecdf_overlay_grouped", group = "focalID") + scale_x_continuous(limits = c(0,1))
