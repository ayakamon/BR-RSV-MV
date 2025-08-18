### created by Ayaka Monoi
### initial input from Koltai, M., Moyes, J., Nyawanda, B. et al. Estimating the cost-effectiveness of maternal vaccination and monoclonal antibodies for respiratory syncytial virus in Kenya and South Africa. BMC Med 21, 120 (2023). https://doi.org/10.1186/s12916-023-02806-w



#### libraries ####
# library(MASS)
# library(dplyr)
# library(conflicted)
# library(here)
# conflicts_prefer(dplyr::filter)
# conflicts_prefer(dplyr::select)

set.seed(1234)
#################################Country-specific RSV-associated deaths ##############################
##Data
# to calculate death averted during the 1st year of life, sum up death_avert_i from 0 month to 11 month
# from RSV SARI hospitalized to RSV-associated deaths (point estimate & CIs) in Koltai et al 2023

# South Africa

## SARI
sari_sa <- read.csv(here("data","Koltai_figure_2B.csv")) %>%
  filter(disease_type_medic_status == "SARI hospitalised") %>%
  filter(age %in% 0:11) %>% 
  select(age, rate, rate_CI_lower, rate_CI_upper, popul_denom) #select age and rate

## death
burden_sa <- read.csv(here("data","Koltai_figure_2B.csv")) %>%
  filter(disease_type_medic_status == "SARI hospitalised") %>%
  filter(age %in% 0:11) %>% #under 1y
  mutate(in_hosp_cfr = c(.0121, .0121,.0121,.0121,.0121,.0121,
                         .0108,.0108,.0108,.0108,.0108,.0108)) %>%#Add in hospital CFR from Koltai et al Fig 3
  mutate(out_in = rep(.35, 12)) %>% #ratio of out vs in hospital CFR from Fig3 in Koltai et al 2023
  mutate(death_mean = rate  * in_hosp_cfr * (1+ out_in) /12)%>% #point estimate of deaths at each age in month # per 100000 person years
  mutate(death_ci_low = rate_CI_lower  * in_hosp_cfr * (1+ out_in) /12 ) %>% #lower CI
  mutate(death_ci_high = rate_CI_upper * in_hosp_cfr * (1+ out_in) /12) %>% #upper CI
  dplyr::select(age, death_mean, death_ci_high, death_ci_low)


# 
# 
##################
## Make posterior distribution
## generate 10000 samples of RSV-associated deaths at each age in month, assuming normal distribution

num_s = 10000

### severe ARI hospitalised
# list for age-specific result
sampled_sari_list <- list()

##sample by age in month
for (i in 1:nrow(sari_sa)){
  # get age information
  age <- sari_sa$age[i]
  
  # get mean and CIs from that age
  mean <- sari_sa$rate[i]
  lower <- sari_sa$rate_CI_lower[i]
  upper <- sari_sa$rate_CI_upper[i]
  
  # calculate SE from data
  initial_se <- (upper - lower) / (2 * qnorm(0.975))
  
  # draw samples assuming normal distribution
  sample <- rnorm(num_s, mean = mean, sd = initial_se)
  
  # return to the list
  sampled_sari_list[[i]] <- data.frame(age = age, sari = sample)
}

### death 
# list for age-specific result
sampled_death_list <- list()

##sample by age in month
for (i in 1:nrow(burden_sa)){
  # get age information
  age <- burden_sa$age[i]

  # get mean and CIs from that age
  mean <- burden_sa$death_mean[i]
  lower <- burden_sa$death_ci_low[i]
  upper <- burden_sa$death_ci_high[i]

  # calculate SE from data
  initial_se <- (upper - lower) / (2 * qnorm(0.975))

  # draw samples assuming normal distribution
  sample <- rnorm(num_s, mean = mean, sd = initial_se)

  # return to the list
  sampled_death_list[[i]] <- data.frame(age = age, death = sample)
}


# convert the list to dataframe
sampled_death <- do.call(rbind, sampled_death_list)
sampled_death%>%mutate(country = "ZAF") -> death_sa  ##posterior samples of RSV-associated deaths at age 0:11 in month


### underlying severe disease and death in South Africa ###
# under-1-year-old RSV SARI hospitalised
# samples for each age in month (0-11)
for (i in 1:12) {
  sampled_sari_list[[i]] %>% pull(sari) -> temp_sa
  if (i == 1) {
    total_sa <- temp_sa
  } else {
    total_sa <- total_sa + temp_sa
  }
}

# median
total_sa %>% median()
# lower 95% Crl
total_sa %>% quantile(.025)
# upper 95% Crl
total_sa %>% quantile(.975)



# under-1-year-old RSV deaths
# samples for each age in month (0-11)
for (i in 1:12) {
  sampled_death_list[[i]] %>% pull(death) -> temp
  if (i == 1) {
    total_death <- temp
  } else {
    total_death <- total_death + temp
  }
}

# median
total_death %>% median()
# lower 95% Crl
total_death %>% quantile(.025)
# upper 95% Crl
total_death %>% quantile(.975)
