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
## generate 1000 samples of RSV-associated deaths at each age in month, assuming normal distribution

num_s = 10000

# list for age-specific result
sampled_sari_list <- list()

##sample by age in month
#South Africa
for (i in 1:nrow(burden_sa)){
  # get age information
  age <- burden_sa$age[i]

  # get mean and CIs from that age
  mean <- burden_sa$death_mean[i]
  lower <- burden_sa$death_ci_low[i]
  upper <- burden_sa$death_ci_high[i]

  # calculate SE from data
  initial_se <- (upper - lower) / (2 * qnorm(0.975))

  # draw 1000 samples assuming normal distribution
  sample <- rnorm(num_s, mean = mean, sd = initial_se)

  # return to the list
  sampled_sari_list[[i]] <- data.frame(age = age, death = sample)
}


# convert the list to dataframe
sampled_sari <- do.call(rbind, sampled_sari_list)
sampled_sari%>%mutate(country = "ZAF") -> death_sa  ##posterior samples of RSV-associated deaths at age 0:11 in month

