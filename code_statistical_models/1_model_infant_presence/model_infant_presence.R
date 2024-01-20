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

# NEW PART: other people presence in the group
d$n_others = d$n_earlyC + d$n_middleC + d$n_adolescent + d$n_adult -1
d$others_presence = ifelse(d$n_others==0, 0, 1)
table(d$others_presence)
d$others_presence <- as.factor(d$others_presence)

# Filter out missing outcomes
d1 <- d %>% 
  filter(!(is.na(trip_duration)))

d1 <- d1 %>% 
  filter(!is.na(food_items)) 
nrow(d1) # 359 trips

# Compute group average of each predictor
d1$focalID = as.factor(d1$focalID)

d1 <- d1 %>%
  group_by(focalID) %>%
  
  # Get focalID-level average of predictors
  summarise(
    avg_prop_infant = mean(infant_presence == 1, na.rm = T),
    avg_prop_others = mean(others_presence == 1, na.rm = T)) %>%
  
  # merge with main dataframe
  right_join(d1) %>%
  
  mutate(
    # scale responses by max value to keep in range of [0,1], makes it easier to set priors
    dur_s = trip_duration / max(trip_duration, na.rm=T),
    dist_s = trip_distance / max(trip_distance, na.rm = T), 
    max_dist_s = maximum_distance / max(maximum_distance, na.rm=T),
    range_s = trip_range_MCP / max(trip_range_MCP, na.rm=T),
    EE_s = Energy_expenditure / max(Energy_expenditure, na.rm = T),
    kcal_s = kcal / max(kcal, na.rm=T)
  )

##########################
###  Fit the models  #####
##########################
# Weakly regularizing priors, understood in light of the data
priors <- c(
  set_prior("normal(-1, 1)", class = "Intercept"),
  set_prior("exponential(5)", class = "sigma"),
  set_prior("normal(0,0.2)", class = "b"),
  set_prior("exponential(5)", class = "sd"))


### DURATION ###
m_dur_0 <- brm(dur_s ~ (1 | focalID) + (1 | dayID) + (1 | food_items),
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

m_dur <- brm(dur_s ~ # focalID-level predictors (fixed effects)
                     avg_prop_infant + avg_prop_others +
                     # trip-level predictors (fixed effects) 
                     infant_presence + others_presence + 
                     # trip-level predictors (random effects)
                     (1 + infant_presence + others_presence | focalID) +
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
summary(m_dur)
conditional_effects(m_dur)

m_dur_0 <- add_criterion(m_dur_0, criterion = "loo")
m_dur <- add_criterion(m_dur, criterion = "loo") 
loo::loo_compare(loo(m_dur_0), loo(m_dur))
loo::loo_model_weights(list(m0 = loo(m_dur_0), m_dur = loo(m_dur)))


### TOTAL DISTANCE ###
m_dist_0 <- brm(dist_s ~ (1 | focalID) + (1 | dayID) + (1 | food_items),
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

m_dist <- brm(dist_s ~ # focalID-level predictors (fixed effects)
                       avg_prop_infant + avg_prop_others +
                       # trip-level predictors (fixed effects) 
                       infant_presence + others_presence + 
                       # trip-level predictors (random effects)
                       (1 + infant_presence + others_presence | focalID) +
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
summary(m_dist)
conditional_effects(m_dist)

m_dist_0 <- add_criterion(m_dist_0, criterion = "loo")
m_dist <- add_criterion(m_dist, criterion = "loo") # moment_match helps with influential observations, omitting that step for no
loo::loo_compare(loo(m_dist_0), loo(m_dist_group_binary))
loo::loo_model_weights(list(m0 = loo(m_dist_0), m_dist = loo(m_dist)))

# when traveling with others, go farther trips, regardless of the baby presence

### MAX DISTANCE ###
m_max_0 <- brm(max_dist_s ~ (1 | focalID) + (1 | dayID) + (1 | food_items),
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

m_max <- brm(max_dist_s ~   # focalID-level predictors (fixed effects)
                            avg_prop_infant + avg_prop_others +
                            # trip-level predictors (fixed effects) 
                            infant_presence + others_presence + 
                            # trip-level predictors (random effects)
                            (1 + infant_presence + others_presence | focalID) +
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
summary(m_max)
conditional_effects(m_max)

m_max_0 <- add_criterion(m_max_0, criterion = "loo")
m_max <- add_criterion(m_max, criterion = "loo") # moment_match helps with influential observations, omitting that step for no
loo::loo_compare(loo(m_max_0), loo(m_max))
loo::loo_model_weights(list(m0 = loo(m_max_0), m_max = loo(m_max)))

# when traveling with others, go farther away, regardless of the baby presence


### range ###
m_range_0 <- brm(range_s ~ (1 | focalID) + (1 | dayID) + (1 | food_items),
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

m_range <- brm(range_s ~    # focalID-level predictors (fixed effects)
                            avg_prop_infant + avg_prop_others +
                            # trip-level predictors (fixed effects) 
                            infant_presence + others_presence + 
                            # trip-level predictors (random effects)
                            (1 + infant_presence + others_presence | focalID) +
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
summary(m_range)
conditional_effects(m_range)

m_range_0 <- add_criterion(m_range_0, criterion = "loo")
m_range_group_binary <- add_criterion(m_range, criterion = "loo") # moment_match helps with influential observations, omitting that step for no
loo::loo_compare(loo(m_range_0), loo(m_range))
loo::loo_model_weights(list(m0 = loo(m_range_0), m_range = loo(m_range)))

### EE ###
m_EE_0 <- brm(EE_s ~ (1 | focalID) + (1 | dayID) + (1 | food_items),
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

m_EE <- brm(EE_s ~   # focalID-level predictors (fixed effects)
                     avg_prop_infant + avg_prop_others +
                     # trip-level predictors (fixed effects) 
                     infant_presence + others_presence + 
                     # trip-level predictors (random effects)
                     (1 + infant_presence + others_presence | focalID) +
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
summary(m_EE)
conditional_effects(m_EE)

m_EE_0 <- add_criterion(m_EE_0, criterion = "loo")
m_EE <- add_criterion(m_EE, criterion = "loo") # moment_match helps with influential observations, omitting that step for no
loo::loo_compare(loo(m_EE_0), loo(m_EE))
loo::loo_model_weights(list(m0 = loo(m_EE_0), m_EE = loo(m_EE)))


### Kcal ###
m_kcal_0 <- brm(kcal_s ~  s(dur_s) + 
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

m_kcal <- brm(kcal_s ~       # trip-level predictors to control trip duration 
                             s(dur_s) + 
                             # focalID-level predictors (fixed effects)
                             avg_prop_infant + avg_prop_others +
                             # trip-level predictors (fixed effects) 
                             infant_presence + others_presence + 
                             # trip-level predictors (random effects)
                             (1 + infant_presence + others_presence | focalID) +
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
summary(m_kcal)
conditional_effects(m_kcal)

m_kcal_0 <- add_criterion(m_kcal_0, criterion = "loo")
m_kcal <- add_criterion(m_kcal, criterion = "loo") 
loo::loo_compare(loo(m_kcal_0), loo(m_kcal))
loo::loo_model_weights(list(m0 = loo(m_kcal_0), m_kcal = loo(m_kcal)))

# save.image("infant_presence.RData")

summary(m_dur)  # infant < other
summary(m_dist) # other
summary(m_max) # other
summary(m_range) # other
summary(m_EE)
summary(m_kcal) # nothing











