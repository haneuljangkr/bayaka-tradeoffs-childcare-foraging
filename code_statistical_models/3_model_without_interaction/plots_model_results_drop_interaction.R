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
library(patchwork)
library(ggpubr)
library(ggplot2)

####################################
### Posterior predictive plots #####
####################################

# Load the model results
load("duration_drop_interaction.RData")
load("total_distance_drop_interaction.RData")
load("max_distance_drop_interaction.RData")
load("range_drop_interaction.RData")
load("EE_drop_interaction.RData")
load("kcal_drop_interaction.RData")

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
nrow(d1) # 359 trips

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
    
    # scale responses by max value to kmeep in kcal of [0,1], makes it easier to set priors
    dur_s = trip_duration / max(trip_duration, na.rm=T),
    dist_s = trip_distance / max(trip_distance, na.rm=T),
    max_dist_s = maximum_distance / max(maximum_distance, na.rm=T),
    range_s = trip_range_MCP / max(trip_range_MCP, na.rm=T),
    EE_s = Energy_expenditure / max(Energy_expenditure, na.rm = T),
    kcal_s = kcal / max(kcal, na.rm=T)
    )


####################################
##### Travel Duration models #######
####################################

### (1) Predict effect of more adults

newdata_adults_dur <- expand.grid(
  n_adults_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_children_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_adults <- posterior_epred(m1_dur_age, newdata = newdata_adults_dur, re.form = NA)

pred_adults_long <- as.data.frame(pred_adults) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_adults_dur) %>% 
  # put back onto original scale
  mutate(pred_dur = pred * max(d1$trip_duration, na.rm=T)) 

pred_adults_summary <- pred_adults_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_dur), n_adults_c = unique(n_adults_c)) 


### (2) Predict effect of more children

newdata_children_dur <- expand.grid(
  n_children_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_adults_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_children <- posterior_epred(m1_dur_age, newdata = newdata_children_dur, re.form = NA)

pred_children_long <- as.data.frame(pred_children) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_children_dur) %>% 
  # put back onto original scale
  mutate(pred_dur = pred * max(d1$trip_duration, na.rm=T))

pred_children_summary <- pred_children_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_dur), n_children_c = unique(n_children_c))


### (3) Predict effect of more females 

newdata_females_dur <- expand.grid(
  n_females_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_males_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_females <- posterior_epred(m1_dur_sex, newdata = newdata_females_dur, re.form = NA)

pred_females_long <- as.data.frame(pred_females) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_females_dur) %>% 
  # put back onto original scale
  mutate(pred_dur = pred * max(d1$trip_duration, na.rm=T))

pred_females_summary <- pred_females_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_dur), n_females_c = unique(n_females_c))

### (4) Predict effect of more males 

newdata_males_dur <- expand.grid(
  n_males_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_females_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_males <- posterior_epred(m1_dur_sex, newdata = newdata_males_dur, re.form = NA)

pred_males_long <- as.data.frame(pred_males) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_males_dur) %>% 
  # put back onto original scale
  mutate(pred_dur = pred * max(d1$trip_duration, na.rm=T))

pred_males_summary <- pred_males_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_dur), n_males_c = unique(n_males_c))


### (5) Predict effect of more kin

newdata_kin_dur <- expand.grid(
  n_kin_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_nonkin_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_kin <- posterior_epred(m1_dur_kin, newdata = newdata_kin_dur, re.form = NA)

pred_kin_long <- as.data.frame(pred_kin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_kin_dur) %>% 
  # put back onto original scale
  mutate(pred_dur = pred * max(d1$trip_duration, na.rm=T))

pred_kin_summary <- pred_kin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_dur), n_kin_c = unique(n_kin_c))

### (6) Predict effect of more nonkin

newdata_nonkin_dur <- expand.grid(
  n_nonkin_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_kin_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_nonkin <- posterior_epred(m1_dur_kin, newdata = newdata_nonkin_dur, re.form = NA)

pred_nonkin_long <- as.data.frame(pred_nonkin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_nonkin_dur) %>% 
  # put back onto original scale
  mutate(pred_dur = pred * max(d1$trip_duration, na.rm=T))

pred_nonkin_summary <- pred_nonkin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_dur), n_nonkin_c = unique(n_nonkin_c))

#################################################################

pred_adults_long$orig = (pred_adults_long$n_adults_c)*sd(d1$n_adult) + mean(d1$n_adult)
pred_adults_summary$orig = (pred_adults_summary$n_adults_c)*sd(d1$n_adult) + mean(d1$n_adult)

p_dur_adults <- ggplot( data = filter(pred_adults_long, samp <= 100), 
                        aes(x = orig, y = pred_dur, group = factor(samp), color = "#fdae6b") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_adults_summary, aes(x = orig, y = med, group = NULL, color = "#fdae6b"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab("Trip duration (hr)") +
  xlab("Adults") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0.3), "lines")) +
  ylim(0.5, 8) +
  scale_colour_manual(values= "#fdae6b") + 
  scale_fill_manual(values= "#fdae6b")  
#p_dur_adults

pred_children_long$orig = (pred_children_long$n_children_c)*sd(d1$n_children) + mean(d1$n_children)
pred_children_summary$orig = (pred_children_summary$n_children_c)*sd(d1$n_children) + mean(d1$n_children)

p_dur_children <- ggplot(data = filter(pred_children_long, samp <= 100), 
                         aes(x = orig, y = pred_dur, group = factor(samp), color = "#2ca25f") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_children_summary, aes(x = orig, y = med, group = NULL, color = "#2ca25f"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Children") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 8) +
  scale_colour_manual(values= "#2ca25f") + 
  scale_fill_manual(values= "#2ca25f")  
#p_dur_children


pred_females_long$orig = (pred_females_long$n_females_c)*sd(d1$n_females) + mean(d1$n_females)
pred_females_summary$orig = (pred_females_summary$n_females_c)*sd(d1$n_females) + mean(d1$n_females)

p_dur_females <- ggplot( data = filter(pred_females_long, samp <= 100), 
                         aes(x = orig, y = pred_dur, group = factor(samp), color = "#fdae6b") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_females_summary, aes(x = orig, y = med, group = NULL, color = "#fdae6b"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Females") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  #  xlim(0, 5) + 
  ylim(0.5, 8) 
#p_dur_females 

pred_males_long$orig = (pred_males_long$n_males_c)*sd(d1$n_males) + mean(d1$n_males)
pred_males_summary$orig = (pred_males_summary$n_males_c)*sd(d1$n_males) + mean(d1$n_males)

p_dur_males <- ggplot( data = filter(pred_males_long, samp <= 100), 
                       aes(x = orig, y = pred_dur, group = factor(samp), color = "#473C8B") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_males_summary, aes(x = orig, y = med, group = NULL, color = "#473C8B"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Males") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 8) +
  scale_colour_manual(values= "#473C8B") + 
  scale_fill_manual(values= "#473C8B")
#p_dur_males 


pred_kin_long$orig = (pred_kin_long$n_kin_c)*sd(d1$kin) + mean(d1$kin)
pred_kin_summary$orig = (pred_kin_summary$n_kin_c)*sd(d1$kin) + mean(d1$kin)

p_dur_kin <- ggplot( data = filter(pred_kin_long, samp <= 100), 
                     aes(x = orig, y = pred_dur, color = "#8B5A2B", group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_kin_summary, aes(x = orig, y = med, color = "#8B5A2B", group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Kin") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 8) +
  scale_colour_manual(values= "#8B5A2B") + 
  scale_fill_manual(values= "#8B5A2B")
#p_dur_kin

pred_nonkin_long$orig = (pred_nonkin_long$n_nonkin_c)*sd(d1$nonkin) + mean(d1$nonkin)
pred_nonkin_summary$orig = (pred_nonkin_summary$n_nonkin_c)*sd(d1$nonkin) + mean(d1$nonkin)

p_dur_nonkin <- ggplot( data = filter(pred_nonkin_long, samp <= 100), 
                        aes(x = orig, y = pred_dur, color = "#36648B", group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_nonkin_summary, aes(x = orig, y = med, color = "#36648B", group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Nonkin") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 8) +
  scale_colour_manual(values= "#36648B") + 
  scale_fill_manual(values= "#36648B")
#p_dur_nonkin


figure1 = ggarrange(
  p_dur_adults, p_dur_children, p_dur_females, p_dur_males, p_dur_kin, p_dur_nonkin,
  ncol = 6,
  nrow = 1,
  label.x = 0,
  label.y = 1,
  hjust = -0.5,
  vjust = 1.5,
  font.label = list(size = 12, color = "black", face = "bold", family = NULL),
  align = c("none", "h", "v", "hv"),
  widths = 1.2,
  heights = 1,
  legend = NULL,
  common.legend = FALSE
)

figure1 = annotate_figure(figure1,
                          top = text_grob("(a) Travel duration ", color = "black", size = 15, face = "bold"),
                          left = text_grob("hour", color = "black", rot = 90, size = 15, face = "bold"))

figure1


####################################
##### Total travel distance  #######
####################################

### (1) Predict effect of more adults

newdata_adults_dist <- expand.grid(
  n_adults_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_children_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_adults <- posterior_epred(m1_dist_age, newdata = newdata_adults_dist, re.form = NA)

pred_adults_long <- as.data.frame(pred_adults) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_adults_dist) %>% 
  # put back onto original scale
  mutate(pred_km = (pred * max(d1$trip_distance, na.rm = T))/1000)

pred_adults_summary <- pred_adults_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_km), n_adults_c = unique(n_adults_c)) 


### (2) Predict effect of more children

newdata_children_dist <- expand.grid(
  n_children_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_adults_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_children <- posterior_epred(m1_dist_age, newdata = newdata_children_dist, re.form = NA)

pred_children_long <- as.data.frame(pred_children) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_children_dist) %>% 
  # put back onto original scale
  mutate(pred_km = (pred * max(d1$trip_distance, na.rm = T))/1000)

pred_children_summary <- pred_children_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_km), n_children_c = unique(n_children_c))


### (3) Predict effect of more females 

newdata_females_dist <- expand.grid(
  n_females_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_males_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_females <- posterior_epred(m1_dist_sex, newdata = newdata_females_dist, re.form = NA)

pred_females_long <- as.data.frame(pred_females) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_females_dist) %>% 
  # put back onto original scale
  mutate(pred_km = (pred * max(d1$trip_distance, na.rm = T))/1000)

pred_females_summary <- pred_females_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_km), n_females_c = unique(n_females_c))

### (4) Predict effect of more males 

newdata_males_dist <- expand.grid(
  n_males_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_females_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_males <- posterior_epred(m1_dist_sex, newdata = newdata_males_dist, re.form = NA)

pred_males_long <- as.data.frame(pred_males) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_males_dist) %>% 
  # put back onto original scale
  mutate(pred_km = (pred * max(d1$trip_distance, na.rm = T))/1000)

pred_males_summary <- pred_males_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_km), n_males_c = unique(n_males_c))


### (5) Predict effect of more kin

newdata_kin_dist <- expand.grid(
  n_kin_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_nonkin_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_kin <- posterior_epred(m1_dist_kin, newdata = newdata_kin_dist, re.form = NA)

pred_kin_long <- as.data.frame(pred_kin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_kin_dist) %>% 
  # put back onto original scale
  mutate(pred_km = (pred * max(d1$trip_distance, na.rm = T))/1000)

pred_kin_summary <- pred_kin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_km), n_kin_c = unique(n_kin_c))

### (6) Predict effect of more nonkin

newdata_nonkin_dist <- expand.grid(
  n_nonkin_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_kin_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_nonkin <- posterior_epred(m1_dist_kin, newdata = newdata_nonkin_dist, re.form = NA)

pred_nonkin_long <- as.data.frame(pred_nonkin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_nonkin_dist) %>% 
  # put back onto original scale
  mutate(pred_km = (pred * max(d1$trip_distance, na.rm = T))/1000)

pred_nonkin_summary <- pred_nonkin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_km), n_nonkin_c = unique(n_nonkin_c))

#################################################################

pred_adults_long$orig = (pred_adults_long$n_adults_c)*sd(d1$n_adult) + mean(d1$n_adult)
pred_adults_summary$orig = (pred_adults_summary$n_adults_c)*sd(d1$n_adult) + mean(d1$n_adult)

p_dist_adults <- ggplot( data = filter(pred_adults_long, samp <= 100), 
                         aes(x = orig, y = pred_km, group = factor(samp), color = "#fdae6b") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_adults_summary, aes(x = orig, y = med, group = NULL, color = "#fdae6b"), lwd = 2) +
  theme_minimal(base_size = 18) +
  xlab("Adults") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        #  axis.text.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0.3), "lines")) +
  ylim(0.5, 15) +
  scale_colour_manual(values= "#fdae6b") + 
  scale_fill_manual(values= "#fdae6b")  
#p_dist_adults

pred_children_long$orig = (pred_children_long$n_children_c)*sd(d1$n_children) + mean(d1$n_children)
pred_children_summary$orig = (pred_children_summary$n_children_c)*sd(d1$n_children) + mean(d1$n_children)

p_dist_children <- ggplot(data = filter(pred_children_long, samp <= 100), 
                          aes(x = orig, y = pred_km, group = factor(samp), color = "#2ca25f") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_children_summary, aes(x = orig, y = med, group = NULL, color = "#2ca25f"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Children") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 15) +
  scale_colour_manual(values= "#2ca25f") + 
  scale_fill_manual(values= "#2ca25f")  
#p_dist_children


pred_females_long$orig = (pred_females_long$n_females_c)*sd(d1$n_females) + mean(d1$n_females)
pred_females_summary$orig = (pred_females_summary$n_females_c)*sd(d1$n_females) + mean(d1$n_females)

p_dist_females <- ggplot( data = filter(pred_females_long, samp <= 100), 
                          aes(x = orig, y = pred_km, group = factor(samp), color = "#fdae6b") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_females_summary, aes(x = orig, y = med, group = NULL, color = "#fdae6b"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Females") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 15) 
#p_dist_females 

pred_males_long$orig = (pred_males_long$n_males_c)*sd(d1$n_males) + mean(d1$n_males)
pred_males_summary$orig = (pred_males_summary$n_males_c)*sd(d1$n_males) + mean(d1$n_males)

p_dist_males <- ggplot( data = filter(pred_males_long, samp <= 100), 
                        aes(x = orig, y = pred_km, group = factor(samp), color = "#473C8B") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_males_summary, aes(x = orig, y = med, group = NULL, color = "#473C8B"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Males") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 15) +
  scale_colour_manual(values= "#473C8B") + 
  scale_fill_manual(values= "#473C8B")
#p_dist_males 


pred_kin_long$orig = (pred_kin_long$n_kin_c)*sd(d1$kin) + mean(d1$kin)
pred_kin_summary$orig = (pred_kin_summary$n_kin_c)*sd(d1$kin) + mean(d1$kin)

p_dist_kin <- ggplot( data = filter(pred_kin_long, samp <= 100), 
                      aes(x = orig, y = pred_km, color = "#8B5A2B", group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_kin_summary, aes(x = orig, y = med, color = "#8B5A2B", group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Kin") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 15) +
  scale_colour_manual(values= "#8B5A2B") + 
  scale_fill_manual(values= "#8B5A2B")
#p_dist_kin

pred_nonkin_long$orig = (pred_nonkin_long$n_nonkin_c)*sd(d1$nonkin) + mean(d1$nonkin)
pred_nonkin_summary$orig = (pred_nonkin_summary$n_nonkin_c)*sd(d1$nonkin) + mean(d1$nonkin)

p_dist_nonkin <- ggplot( data = filter(pred_nonkin_long, samp <= 100), 
                         aes(x = orig, y = pred_km, color = "#36648B", group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_nonkin_summary, aes(x = orig, y = med, color = "#36648B", group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Nonkin") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 15) +
  scale_colour_manual(values= "#36648B") + 
  scale_fill_manual(values= "#36648B")
#p_dist_nonkin


figure2 = ggarrange(
  p_dist_adults, p_dist_children, p_dist_females, p_dist_males, p_dist_kin, p_dist_nonkin,
  ncol = 6,
  nrow = 1,
  label.x = 0,
  label.y = 1,
  hjust = -0.5,
  vjust = 1.5,
  font.label = list(size = 13, color = "black", face = "bold", family = NULL),
  align = c("none", "h", "v", "hv"),
  widths = 1.2,
  heights = 1,
  legend = NULL,
  common.legend = FALSE,
  legend.grob = NULL
)

figure2 = annotate_figure(figure2,
                          top = text_grob("(b) Total travel distance ", color = "black", size = 15, face = "bold"),
                          left = text_grob("km", color = "black", rot = 90, size = 15, face = "bold"))

figure2

####################################
#### Max distance from village #####
####################################

### (1) Predict effect of more adults

newdata_adults_max <- expand.grid(
  n_adults_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_children_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_adults <- posterior_epred(m1_max_age, newdata = newdata_adults_max, re.form = NA)

pred_adults_long <- as.data.frame(pred_adults) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_adults_max) %>% 
  # put back onto original scale
  mutate(pred_max = pred * max(d1$maximum_distance, na.rm=T)/1000)

pred_adults_summary <- pred_adults_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_max), n_adults_c = unique(n_adults_c)) 


### (2) Predict effect of more children

newdata_children_max <- expand.grid(
  n_children_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_adults_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_children <- posterior_epred(m1_max_age, newdata = newdata_children_max, re.form = NA)

pred_children_long <- as.data.frame(pred_children) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_children_max) %>% 
  # put back onto original scale
  mutate(pred_max = pred * max(d1$maximum_distance, na.rm=T)/1000)

pred_children_summary <- pred_children_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_max), n_children_c = unique(n_children_c))


### (3) Predict effect of more females 

newdata_females_max <- expand.grid(
  n_females_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_males_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_females <- posterior_epred(m1_max_sex, newdata = newdata_females_max, re.form = NA)

pred_females_long <- as.data.frame(pred_females) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_females_max) %>% 
  # put back onto original scale
  mutate(pred_max = pred * max(d1$maximum_distance, na.rm=T)/1000)

pred_females_summary <- pred_females_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_max), n_females_c = unique(n_females_c))

### (4) Predict effect of more males 

newdata_males_max <- expand.grid(
  n_males_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_females_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_males <- posterior_epred(m1_max_sex, newdata = newdata_males_max, re.form = NA)

pred_males_long <- as.data.frame(pred_males) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_males_max) %>% 
  # put back onto original scale
  mutate(pred_max = pred * max(d1$maximum_distance, na.rm=T)/1000)

pred_males_summary <- pred_males_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_max), n_males_c = unique(n_males_c))


### (5) Predict effect of more kin

newdata_kin_max <- expand.grid(
  n_kin_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_nonkin_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_kin <- posterior_epred(m1_max_kin, newdata = newdata_kin_max, re.form = NA)

pred_kin_long <- as.data.frame(pred_kin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_kin_max) %>% 
  # put back onto original scale
  mutate(pred_max = pred * max(d1$maximum_distance, na.rm=T)/1000)

pred_kin_summary <- pred_kin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_max), n_kin_c = unique(n_kin_c))

### (6) Predict effect of more nonkin

newdata_nonkin_max <- expand.grid(
  n_nonkin_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_kin_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_nonkin <- posterior_epred(m1_max_kin, newdata = newdata_nonkin_max, re.form = NA)

pred_nonkin_long <- as.data.frame(pred_nonkin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_nonkin_max) %>% 
  # put back onto original scale
  mutate(pred_max = pred * max(d1$maximum_distance, na.rm=T)/1000)

pred_nonkin_summary <- pred_nonkin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_max), n_nonkin_c = unique(n_nonkin_c))

#################################################################

pred_adults_long$orig = (pred_adults_long$n_adults_c)*sd(d1$n_adult) + mean(d1$n_adult)
pred_adults_summary$orig = (pred_adults_summary$n_adults_c)*sd(d1$n_adult) + mean(d1$n_adult)

p_max_adults <- ggplot( data = filter(pred_adults_long, samp <= 100), 
                        aes(x = orig, y = pred_max, group = factor(samp), color = "#fdae6b") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_adults_summary, aes(x = orig, y = med, group = NULL, color = "#fdae6b"), lwd = 2) +
  theme_minimal(base_size = 18) +
  # ylab(" ") +
  xlab("Adults") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        #  axis.text.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0.3), "lines")) +
  ylim(0.5, 3) +
  scale_colour_manual(values= "#fdae6b") + 
  scale_fill_manual(values= "#fdae6b")  
#p_max_adults

pred_children_long$orig = (pred_children_long$n_children_c)*sd(d1$n_children) + mean(d1$n_children)
pred_children_summary$orig = (pred_children_summary$n_children_c)*sd(d1$n_children) + mean(d1$n_children)

p_max_children <- ggplot(data = filter(pred_children_long, samp <= 100), 
                         aes(x = orig, y = pred_max, group = factor(samp), color = "#2ca25f") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_children_summary, aes(x = orig, y = med, group = NULL, color = "#2ca25f"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Children") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 3) +
  scale_colour_manual(values= "#2ca25f") + 
  scale_fill_manual(values= "#2ca25f")  
#p_max_children


pred_females_long$orig = (pred_females_long$n_females_c)*sd(d1$n_females) + mean(d1$n_females)
pred_females_summary$orig = (pred_females_summary$n_females_c)*sd(d1$n_females) + mean(d1$n_females)

p_max_females <- ggplot( data = filter(pred_females_long, samp <= 100), 
                         aes(x = orig, y = pred_max, group = factor(samp), color = "#fdae6b") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_females_summary, aes(x = orig, y = med, group = NULL, color = "#fdae6b"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Females") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 3) 
#p_max_females 

pred_males_long$orig = (pred_males_long$n_males_c)*sd(d1$n_males) + mean(d1$n_males)
pred_males_summary$orig = (pred_males_summary$n_males_c)*sd(d1$n_males) + mean(d1$n_males)

p_max_males <- ggplot( data = filter(pred_males_long, samp <= 100), 
                       aes(x = orig, y = pred_max, group = factor(samp), color = "#473C8B") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_males_summary, aes(x = orig, y = med, group = NULL, color = "#473C8B"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Males") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 3) +
  scale_colour_manual(values= "#473C8B") + 
  scale_fill_manual(values= "#473C8B")
#p_max_males 


pred_kin_long$orig = (pred_kin_long$n_kin_c)*sd(d1$kin) + mean(d1$kin)
pred_kin_summary$orig = (pred_kin_summary$n_kin_c)*sd(d1$kin) + mean(d1$kin)

p_max_kin <- ggplot( data = filter(pred_kin_long, samp <= 100), 
                     aes(x = orig, y = pred_max, color = "#8B5A2B", group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_kin_summary, aes(x = orig, y = med, color = "#8B5A2B", group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Kin") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 3) +
  scale_colour_manual(values= "#8B5A2B") + 
  scale_fill_manual(values= "#8B5A2B")
#p_max_kin

pred_nonkin_long$orig = (pred_nonkin_long$n_nonkin_c)*sd(d1$nonkin) + mean(d1$nonkin)
pred_nonkin_summary$orig = (pred_nonkin_summary$n_nonkin_c)*sd(d1$nonkin) + mean(d1$nonkin)

p_max_nonkin <- ggplot( data = filter(pred_nonkin_long, samp <= 100), 
                        aes(x = orig, y = pred_max, color = "#36648B", group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_nonkin_summary, aes(x = orig, y = med, color = "#36648B", group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Nonkin") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0.5, 3) +
  scale_colour_manual(values= "#36648B") + 
  scale_fill_manual(values= "#36648B")
#p_max_nonkin


figure3 = ggarrange(
  p_max_adults, p_max_children, p_max_females, p_max_males, p_max_kin, p_max_nonkin,
  ncol = 6,
  nrow = 1,
  label.x = 0,
  label.y = 1,
  hjust = -0.5,
  vjust = 1.5,
  font.label = list(size = 13, color = "black", face = "bold", family = NULL),
  align = c("none", "h", "v", "hv"),
  widths = 1.2,
  heights = 1,
  legend = NULL,
  common.legend = FALSE,
  legend.grob = NULL
)

figure3 = annotate_figure(figure3,
                          top = text_grob("(c) Maximum distance from the village ", color = "black", size = 15, face = "bold"),
                          left = text_grob("km", color = "black", rot = 90, size = 15, face = "bold"))

figure3

####################################
##### Exploration Range (km2) ######
####################################

### (1) Predict effect of more adults

newdata_adults_range <- expand.grid(
  n_adults_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_children_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_adults <- posterior_epred(m1_range_age, newdata = newdata_adults_range, re.form = NA)

pred_adults_long <- as.data.frame(pred_adults) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_adults_range) %>% 
  # put back onto original scale
  mutate(pred_range = pred * max(d1$trip_range_MCP, na.rm=T))

pred_adults_summary <- pred_adults_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_range), n_adults_c = unique(n_adults_c)) 


### (2) Predict effect of more children

newdata_children_range <- expand.grid(
  n_children_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_adults_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_children <- posterior_epred(m1_range_age, newdata = newdata_children_range, re.form = NA)

pred_children_long <- as.data.frame(pred_children) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_children_range) %>% 
  # put back onto original scale
  mutate(pred_range = pred * max(d1$trip_range_MCP, na.rm=T))

pred_children_summary <- pred_children_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_range), n_children_c = unique(n_children_c))


### (3) Predict effect of more females 

newdata_females_range <- expand.grid(
  n_females_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_males_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_females <- posterior_epred(m1_range_sex, newdata = newdata_females_range, re.form = NA)

pred_females_long <- as.data.frame(pred_females) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_females_range) %>% 
  # put back onto original scale
  mutate(pred_range = pred * max(d1$trip_range_MCP, na.rm=T))

pred_females_summary <- pred_females_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_range), n_females_c = unique(n_females_c))

### (4) Predict effect of more males 

newdata_males_range <- expand.grid(
  n_males_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_females_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_males <- posterior_epred(m1_range_sex, newdata = newdata_males_range, re.form = NA)

pred_males_long <- as.data.frame(pred_males) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_males_range) %>% 
  # put back onto original scale
  mutate(pred_range = pred * max(d1$trip_range_MCP, na.rm=T))

pred_males_summary <- pred_males_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_range), n_males_c = unique(n_males_c))


### (5) Predict effect of more kin

newdata_kin_range <- expand.grid(
  n_kin_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_nonkin_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_kin <- posterior_epred(m1_range_kin, newdata = newdata_kin_range, re.form = NA)

pred_kin_long <- as.data.frame(pred_kin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_kin_range) %>% 
  # put back onto original scale
  mutate(pred_range = pred * max(d1$trip_range_MCP, na.rm=T))

pred_kin_summary <- pred_kin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_range), n_kin_c = unique(n_kin_c))

### (6) Predict effect of more nonkin

newdata_nonkin_range <- expand.grid(
  n_nonkin_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_kin_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_nonkin <- posterior_epred(m1_range_kin, newdata = newdata_nonkin_range, re.form = NA)

pred_nonkin_long <- as.data.frame(pred_nonkin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_nonkin_range) %>% 
  # put back onto original scale
  mutate(pred_range = pred * max(d1$trip_range_MCP, na.rm=T))

pred_nonkin_summary <- pred_nonkin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_range), n_nonkin_c = unique(n_nonkin_c))

#################################################################

###### lets plot #########
library(patchwork)

pred_adults_long$orig = (pred_adults_long$n_adults_c)*sd(d1$n_adult) + mean(d1$n_adult)
pred_adults_summary$orig = (pred_adults_summary$n_adults_c)*sd(d1$n_adult) + mean(d1$n_adult)

p_range_adults <- ggplot( data = filter(pred_adults_long, samp <= 100), 
                          aes(x = orig, y = pred_range, group = factor(samp), color = "#fdae6b") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_adults_summary, aes(x = orig, y = med, group = NULL, color = "#fdae6b"), lwd = 2) +
  theme_minimal(base_size = 18) +
  # ylab(" ") +
  xlab("Adults") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        #  axis.text.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0.3), "lines")) +
  ylim(0, 2) +
  scale_colour_manual(values= "#fdae6b") + 
  scale_fill_manual(values= "#fdae6b")  
#p_range_adults

pred_children_long$orig = (pred_children_long$n_children_c)*sd(d1$n_children) + mean(d1$n_children)
pred_children_summary$orig = (pred_children_summary$n_children_c)*sd(d1$n_children) + mean(d1$n_children)

p_range_children <- ggplot(data = filter(pred_children_long, samp <= 100), 
                           aes(x = orig, y = pred_range, group = factor(samp), color = "#2ca25f") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_children_summary, aes(x = orig, y = med, group = NULL, color = "#2ca25f"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Children") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0, 2) +
  scale_colour_manual(values= "#2ca25f") + 
  scale_fill_manual(values= "#2ca25f")  
#p_range_children


pred_females_long$orig = (pred_females_long$n_females_c)*sd(d1$n_females) + mean(d1$n_females)
pred_females_summary$orig = (pred_females_summary$n_females_c)*sd(d1$n_females) + mean(d1$n_females)

p_range_females <- ggplot( data = filter(pred_females_long, samp <= 100), 
                           aes(x = orig, y = pred_range, group = factor(samp), color = "#fdae6b") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_females_summary, aes(x = orig, y = med, group = NULL, color = "#fdae6b"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Females") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0, 2) 
#p_range_females 

pred_males_long$orig = (pred_males_long$n_males_c)*sd(d1$n_males) + mean(d1$n_males)
pred_males_summary$orig = (pred_males_summary$n_males_c)*sd(d1$n_males) + mean(d1$n_males)

p_range_males <- ggplot( data = filter(pred_males_long, samp <= 100), 
                         aes(x = orig, y = pred_range, group = factor(samp), color = "#473C8B") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_males_summary, aes(x = orig, y = med, group = NULL, color = "#473C8B"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Males") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0, 2) +
  scale_colour_manual(values= "#473C8B") + 
  scale_fill_manual(values= "#473C8B")
#p_range_males 


pred_kin_long$orig = (pred_kin_long$n_kin_c)*sd(d1$kin) + mean(d1$kin)
pred_kin_summary$orig = (pred_kin_summary$n_kin_c)*sd(d1$kin) + mean(d1$kin)

p_range_kin <- ggplot( data = filter(pred_kin_long, samp <= 100), 
                       aes(x = orig, y = pred_range, color = "#8B5A2B", group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_kin_summary, aes(x = orig, y = med, color = "#8B5A2B", group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Kin") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0, 2) +
  scale_colour_manual(values= "#8B5A2B") + 
  scale_fill_manual(values= "#8B5A2B")
#p_range_kin

pred_nonkin_long$orig = (pred_nonkin_long$n_nonkin_c)*sd(d1$nonkin) + mean(d1$nonkin)
pred_nonkin_summary$orig = (pred_nonkin_summary$n_nonkin_c)*sd(d1$nonkin) + mean(d1$nonkin)

p_range_nonkin <- ggplot( data = filter(pred_nonkin_long, samp <= 100), 
                          aes(x = orig, y = pred_range, color = "#36648B", group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_nonkin_summary, aes(x = orig, y = med, color = "#36648B", group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Nonkin") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(0, 2) +
  scale_colour_manual(values= "#36648B") + 
  scale_fill_manual(values= "#36648B")
#p_range_nonkin


figure4 = ggarrange(
  p_range_adults, p_range_children, p_range_females, p_range_males, p_range_kin, p_range_nonkin,
  ncol = 6,
  nrow = 1,
  label.x = 0,
  label.y = 1,
  hjust = -0.5,
  vjust = 1.5,
  font.label = list(size = 13, color = "black", face = "bold", family = NULL),
  align = c("none", "h", "v", "hv"),
  widths = 1.2,
  heights = 1,
  legend = NULL,
  common.legend = FALSE,
  legend.grob = NULL
)

figure4 = annotate_figure(figure4,
                          top = text_grob("(d) Exploration range", color = "black", size = 15, face = "bold"),
                          left = text_grob("km2", color = "black", rot = 90, size = 15, face = "bold"))

figure4

####################################
### Energy Expenditure (kcal/min) ###
####################################

### (1) Predict effect of more adults

newdata_adults_EE <- expand.grid(
  n_adults_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_children_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_adults <- posterior_epred(m1_EE_age, newdata = newdata_adults_EE, re.form = NA)

pred_adults_long <- as.data.frame(pred_adults) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_adults_EE) %>% 
  # put back onto original scale
  mutate(pred_EE = pred * max(d1$Energy_expenditure, na.rm=T))

pred_adults_summary <- pred_adults_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_EE), n_adults_c = unique(n_adults_c)) 


### (2) Predict effect of more children

newdata_children_EE <- expand.grid(
  n_children_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_adults_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_children <- posterior_epred(m1_EE_age, newdata = newdata_children_EE, re.form = NA)

pred_children_long <- as.data.frame(pred_children) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_children_EE) %>% 
  # put back onto original scale
  mutate(pred_EE = pred * max(d1$Energy_expenditure, na.rm=T))

pred_children_summary <- pred_children_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_EE), n_children_c = unique(n_children_c))


### (3) Predict effect of more females 

newdata_females_EE <- expand.grid(
  n_females_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_males_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_females <- posterior_epred(m1_EE_sex, newdata = newdata_females_EE, re.form = NA)

pred_females_long <- as.data.frame(pred_females) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_females_EE) %>% 
  # put back onto original scale
  mutate(pred_EE = pred * max(d1$Energy_expenditure, na.rm=T))

pred_females_summary <- pred_females_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_EE), n_females_c = unique(n_females_c))

### (4) Predict effect of more males 

newdata_males_EE <- expand.grid(
  n_males_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_females_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_males <- posterior_epred(m1_EE_sex, newdata = newdata_males_EE, re.form = NA)

pred_males_long <- as.data.frame(pred_males) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_males_EE) %>% 
  # put back onto original scale
  mutate(pred_EE = pred * max(d1$Energy_expenditure, na.rm=T))

pred_males_summary <- pred_males_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_EE), n_males_c = unique(n_males_c))


### (5) Predict effect of more kin

newdata_kin_EE <- expand.grid(
  n_kin_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_nonkin_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_kin <- posterior_epred(m1_EE_kin, newdata = newdata_kin_EE, re.form = NA)

pred_kin_long <- as.data.frame(pred_kin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_kin_EE) %>% 
  # put back onto original scale
  mutate(pred_EE = pred * max(d1$Energy_expenditure, na.rm=T))

pred_kin_summary <- pred_kin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_EE), n_kin_c = unique(n_kin_c))

### (6) Predict effect of more nonkin

newdata_nonkin_EE <- expand.grid(
  n_nonkin_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_kin_c = 0,
  infant_presence = 0
) %>% 
  mutate(condition = 1:n())

pred_nonkin <- posterior_epred(m1_EE_kin, newdata = newdata_nonkin_EE, re.form = NA)

pred_nonkin_long <- as.data.frame(pred_nonkin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_nonkin_EE) %>% 
  # put back onto original scale
  mutate(pred_EE = pred * max(d1$Energy_expenditure, na.rm=T))

pred_nonkin_summary <- pred_nonkin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_EE), n_nonkin_c = unique(n_nonkin_c))

#################################################################

pred_adults_long$orig = (pred_adults_long$n_adults_c)*sd(d1$n_adult) + mean(d1$n_adult)
pred_adults_summary$orig = (pred_adults_summary$n_adults_c)*sd(d1$n_adult) + mean(d1$n_adult)

p_EE_adults <- ggplot( data = filter(pred_adults_long, samp <= 100), 
                       aes(x = orig, y = pred_EE, group = factor(samp), color = "#fdae6b") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_adults_summary, aes(x = orig, y = med, group = NULL, color = "#fdae6b"), lwd = 2) +
  theme_minimal(base_size = 18) +
  # ylab(" ") +
  xlab("Adults") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        #  axis.text.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0.3), "lines")) +
  ylim(1, 3) +
  scale_colour_manual(values= "#fdae6b") + 
  scale_fill_manual(values= "#fdae6b")  
#p_EE_adults

pred_children_long$orig = (pred_children_long$n_children_c)*sd(d1$n_children) + mean(d1$n_children)
pred_children_summary$orig = (pred_children_summary$n_children_c)*sd(d1$n_children) + mean(d1$n_children)

p_EE_children <- ggplot(data = filter(pred_children_long, samp <= 100), 
                        aes(x = orig, y = pred_EE, group = factor(samp), color = "#2ca25f") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_children_summary, aes(x = orig, y = med, group = NULL, color = "#2ca25f"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Children") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(1, 3) +
  scale_colour_manual(values= "#2ca25f") + 
  scale_fill_manual(values= "#2ca25f")  
#p_EE_children


pred_females_long$orig = (pred_females_long$n_females_c)*sd(d1$n_females) + mean(d1$n_females)
pred_females_summary$orig = (pred_females_summary$n_females_c)*sd(d1$n_females) + mean(d1$n_females)

p_EE_females <- ggplot( data = filter(pred_females_long, samp <= 100), 
                        aes(x = orig, y = pred_EE, group = factor(samp), color = "#fdae6b") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_females_summary, aes(x = orig, y = med, group = NULL, color = "#fdae6b"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Females") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(1, 3) 
#p_EE_females 

pred_males_long$orig = (pred_males_long$n_males_c)*sd(d1$n_males) + mean(d1$n_males)
pred_males_summary$orig = (pred_males_summary$n_males_c)*sd(d1$n_males) + mean(d1$n_males)

p_EE_males <- ggplot( data = filter(pred_males_long, samp <= 100), 
                      aes(x = orig, y = pred_EE, group = factor(samp), color = "#473C8B") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_males_summary, aes(x = orig, y = med, group = NULL, color = "#473C8B"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Males") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(1, 3) +
  scale_colour_manual(values= "#473C8B") + 
  scale_fill_manual(values= "#473C8B")
#p_EE_males 


pred_kin_long$orig = (pred_kin_long$n_kin_c)*sd(d1$kin) + mean(d1$kin)
pred_kin_summary$orig = (pred_kin_summary$n_kin_c)*sd(d1$kin) + mean(d1$kin)

p_EE_kin <- ggplot( data = filter(pred_kin_long, samp <= 100), 
                    aes(x = orig, y = pred_EE, color = "#8B5A2B", group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_kin_summary, aes(x = orig, y = med, color = "#8B5A2B", group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  #  ylab(" ") +
  xlab("Kin") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(1, 3) +
  scale_colour_manual(values= "#8B5A2B") + 
  scale_fill_manual(values= "#8B5A2B")
#p_EE_kin

pred_nonkin_long$orig = (pred_nonkin_long$n_nonkin_c)*sd(d1$nonkin) + mean(d1$nonkin)
pred_nonkin_summary$orig = (pred_nonkin_summary$n_nonkin_c)*sd(d1$nonkin) + mean(d1$nonkin)

p_EE_nonkin <- ggplot( data = filter(pred_nonkin_long, samp <= 100), 
                       aes(x = orig, y = pred_EE, color = "#36648B", group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_nonkin_summary, aes(x = orig, y = med, color = "#36648B", group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Nonkin") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(1, 3) +
  scale_colour_manual(values= "#36648B") + 
  scale_fill_manual(values= "#36648B")
#p_EE_nonkin


figure5 = ggarrange(
  p_EE_adults, p_EE_children, p_EE_females, p_EE_males, p_EE_kin, p_EE_nonkin,
  ncol = 6,
  nrow = 1,
  label.x = 0,
  label.y = 1,
  hjust = -0.5,
  vjust = 1.5,
  font.label = list(size = 13, color = "black", face = "bold", family = NULL),
  align = c("none", "h", "v", "hv"),
  widths = 1.2,
  heights = 1,
  legend = NULL,
  common.legend = FALSE,
  legend.grob = NULL
)

figure5 = annotate_figure(figure5,
                          top = text_grob("(e) Energy expenditure", color = "black", size = 15, face = "bold"),
                          left = text_grob("kcal/min", color = "black", rot = 90, size = 15, face = "bold"))

figure5

####################################
##### net food returns (kcal) ######
####################################

### (1) Predict effect of more adults

newdata_adults_kcal <- expand.grid(
  n_adults_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_children_c = 0,
  infant_presence = 0, 
  dur_s = 0
) %>% 
  mutate(condition = 1:n())

pred_adults <- posterior_epred(m1_kcal_age, newdata = newdata_adults_kcal, re.form = NA)

pred_adults_long <- as.data.frame(pred_adults) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_adults_kcal) %>% 
  # put back onto original scale
  mutate(pred_kcal = pred * max(d1$kcal, na.rm=T))

pred_adults_summary <- pred_adults_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_kcal), n_adults_c = unique(n_adults_c)) 


### (2) Predict effect of more children

newdata_children_kcal <- expand.grid(
  n_children_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_adults = 0,
  avg_num_children = 0,
  n_adults_c = 0,
  infant_presence = 0,
  dur_s = 0
) %>% 
  mutate(condition = 1:n())

pred_children <- posterior_epred(m1_kcal_age, newdata = newdata_children_kcal, re.form = NA)

pred_children_long <- as.data.frame(pred_children) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_children_kcal) %>% 
  # put back onto original scale
  mutate(pred_kcal = pred * max(d1$kcal, na.rm=T))

pred_children_summary <- pred_children_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_kcal), n_children_c = unique(n_children_c))


### (3) Predict effect of more females 

newdata_females_kcal <- expand.grid(
  n_females_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_males_c = 0,
  infant_presence = 0,
  dur_s = 0
) %>% 
  mutate(condition = 1:n())

pred_females <- posterior_epred(m1_kcal_sex, newdata = newdata_females_kcal, re.form = NA)

pred_females_long <- as.data.frame(pred_females) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_females_kcal) %>% 
  # put back onto original scale
  mutate(pred_kcal = pred * max(d1$kcal, na.rm=T))

pred_females_summary <- pred_females_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_kcal), n_females_c = unique(n_females_c))

### (4) Predict effect of more males 

newdata_males_kcal <- expand.grid(
  n_males_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_females = 0,
  avg_num_males = 0,
  n_females_c = 0,
  infant_presence = 0,
  dur_s = 0
) %>% 
  mutate(condition = 1:n())

pred_males <- posterior_epred(m1_kcal_sex, newdata = newdata_males_kcal, re.form = NA)

pred_males_long <- as.data.frame(pred_males) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_males_kcal) %>% 
  # put back onto original scale
  mutate(pred_kcal = pred * max(d1$kcal, na.rm=T))

pred_males_summary <- pred_males_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_kcal), n_males_c = unique(n_males_c))


### (5) Predict effect of more kin

newdata_kin_kcal <- expand.grid(
  n_kin_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_nonkin_c = 0,
  infant_presence = 0,
  dur_s = 0
) %>% 
  mutate(condition = 1:n())

pred_kin <- posterior_epred(m1_kcal_kin, newdata = newdata_kin_kcal, re.form = NA)

pred_kin_long <- as.data.frame(pred_kin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_kin_kcal) %>% 
  # put back onto original scale
  mutate(pred_kcal = pred * max(d1$kcal, na.rm=T))

pred_kin_summary <- pred_kin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_kcal), n_kin_c = unique(n_kin_c))

### (6) Predict effect of more nonkin

newdata_nonkin_kcal <- expand.grid(
  n_nonkin_c = seq(from = -2, to = 2, length.out = 30),
  # other predictors kept at their average values (0)
  avg_prop_infant = 0,
  avg_num_kin = 0,
  avg_num_nonkin = 0,
  n_kin_c = 0,
  infant_presence = 0,
  dur_s = 0
) %>% 
  mutate(condition = 1:n())

pred_nonkin <- posterior_epred(m1_kcal_kin, newdata = newdata_nonkin_kcal, re.form = NA)

pred_nonkin_long <- as.data.frame(pred_nonkin) %>% 
  mutate(samp = 1:n()) %>% 
  pivot_longer(-samp, names_to = "condition", values_to = "pred") %>% 
  mutate(condition = as.numeric(substr(condition, 2, nchar(condition)))) %>% 
  left_join(newdata_nonkin_kcal) %>% 
  # put back onto original scale
  mutate(pred_kcal = pred * max(d1$kcal, na.rm=T))

pred_nonkin_summary <- pred_nonkin_long %>% 
  group_by(condition) %>% 
  summarise(med = median(pred_kcal), n_nonkin_c = unique(n_nonkin_c))

#################################################################

pred_adults_long$orig = (pred_adults_long$n_adults_c)*sd(d1$n_adult) + mean(d1$n_adult)
pred_adults_summary$orig = (pred_adults_summary$n_adults_c)*sd(d1$n_adult) + mean(d1$n_adult)

p_kcal_adults <- ggplot( data = filter(pred_adults_long, samp <= 100), 
                       aes(x = orig, y = pred_kcal, group = factor(samp), color = "#fdae6b") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_adults_summary, aes(x = orig, y = med, group = NULL, color = "#fdae6b"), lwd = 2) +
  theme_minimal(base_size = 18) +
  # ylab(" ") +
  xlab("Adults") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0.3), "lines")) +
  ylim(2000, 10000) +
  scale_colour_manual(values= "#fdae6b") + 
  scale_fill_manual(values= "#fdae6b")  
#p_kcal_adults

pred_children_long$orig = (pred_children_long$n_children_c)*sd(d1$n_children) + mean(d1$n_children)
pred_children_summary$orig = (pred_children_summary$n_children_c)*sd(d1$n_children) + mean(d1$n_children)

p_kcal_children <- ggplot(data = filter(pred_children_long, samp <= 100), 
                        aes(x = orig, y = pred_kcal, group = factor(samp), color = "#2ca25f") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_children_summary, aes(x = orig, y = med, group = NULL, color = "#2ca25f"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Children") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(2000, 10000) +
  scale_colour_manual(values= "#2ca25f") + 
  scale_fill_manual(values= "#2ca25f")  
#p_kcal_children


pred_females_long$orig = (pred_females_long$n_females_c)*sd(d1$n_females) + mean(d1$n_females)
pred_females_summary$orig = (pred_females_summary$n_females_c)*sd(d1$n_females) + mean(d1$n_females)

p_kcal_females <- ggplot( data = filter(pred_females_long, samp <= 100), 
                        aes(x = orig, y = pred_kcal, group = factor(samp), color = "#fdae6b") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_females_summary, aes(x = orig, y = med, group = NULL, color = "#fdae6b"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Females") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
      ylim(2000, 10000) 
#p_kcal_females 

pred_males_long$orig = (pred_males_long$n_males_c)*sd(d1$n_males) + mean(d1$n_males)
pred_males_summary$orig = (pred_males_summary$n_males_c)*sd(d1$n_males) + mean(d1$n_males)

p_kcal_males <- ggplot( data = filter(pred_males_long, samp <= 100), 
                      aes(x = orig, y = pred_kcal, group = factor(samp), color = "#473C8B") ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_males_summary, aes(x = orig, y = med, group = NULL, color = "#473C8B"), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Males") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(2000, 10000) +
  scale_colour_manual(values= "#473C8B") + 
  scale_fill_manual(values= "#473C8B")
#p_kcal_males 


pred_kin_long$orig = (pred_kin_long$n_kin_c)*sd(d1$kin) + mean(d1$kin)
pred_kin_summary$orig = (pred_kin_summary$n_kin_c)*sd(d1$kin) + mean(d1$kin)

p_kcal_kin <- ggplot( data = filter(pred_kin_long, samp <= 100), 
                    aes(x = orig, y = pred_kcal, color = "#8B5A2B", group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_kin_summary, aes(x = orig, y = med, color = "#8B5A2B", group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  #  ylab(" ") +
  xlab("Kin") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(2000, 10000) +
  scale_colour_manual(values= "#8B5A2B") + 
  scale_fill_manual(values= "#8B5A2B")
#p_kcal_kin

pred_nonkin_long$orig = (pred_nonkin_long$n_nonkin_c)*sd(d1$nonkin) + mean(d1$nonkin)
pred_nonkin_summary$orig = (pred_nonkin_summary$n_nonkin_c)*sd(d1$nonkin) + mean(d1$nonkin)

p_kcal_nonkin <- ggplot( data = filter(pred_nonkin_long, samp <= 100), 
                       aes(x = orig, y = pred_kcal, color = "#36648B", group = factor(samp)) ) +
  geom_line(alpha = 0.3, lwd = 0.3) +
  geom_line(data = pred_nonkin_summary, aes(x = orig, y = med, color = "#36648B", group = NULL), lwd = 2) +
  theme_minimal(base_size = 18) +
  ylab(" ") +
  xlab("Nonkin") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 11, face= "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(1.5,0,1,0), "lines")) +
  ylim(2000, 10000) +
  scale_colour_manual(values= "#36648B") + 
  scale_fill_manual(values= "#36648B")
#p_kcal_nonkin


figure6 = ggarrange(
  p_kcal_adults, p_kcal_children, p_kcal_females, p_kcal_males, p_kcal_kin, p_kcal_nonkin,
  ncol = 6,
  nrow = 1,
  label.x = 0,
  label.y = 1,
  hjust = -0.5,
  vjust = 1.5,
  font.label = list(size = 13, color = "black", face = "bold", family = NULL),
  align = c("none", "h", "v", "hv"),
  widths = 1.2,
  heights = 1,
  legend = NULL,
  common.legend = FALSE,
  legend.grob = NULL
)

figure6 = annotate_figure(figure6,
                          top = text_grob("(f) Net food returns", color = "black", size = 15, face = "bold"),
                          left = text_grob("kcal", color = "black", rot = 90, size = 15, face = "bold"))

figure6

#################################################################
#################################################################

figure_all = ggarrange(
  figure1, figure2, figure3, figure4, figure5, figure6,
  ncol = 1,
  nrow = 6,
  label.x = 0,
  label.y = 1,
  hjust = -0.5,
  vjust = 1.5,
  font.label = list(size = 8, color = "black", face = "plain", family = NULL),
  align = c("none", "h", "v", "hv"),
  widths = 1.2,
  heights = 1,
  legend = NULL,
  common.legend = FALSE,
  legend.grob = NULL
)

figure_all  

ggsave("Rplot08.jpeg", plot = figure_all, 
       height = 15, width = 11.25, 
       dpi = 1200)




