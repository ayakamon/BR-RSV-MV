### created by Ayaka Monoi 



### Before this, ve_fit.R(waning protection estimates) and RSV_burden.R (estimates of country-specific RSV disease) need to be run.#####
set.seed(1234)
#### VE (from ve_fit.R)
##posterior samples of VE after birth
ve_post <-resp %>%
  rowwise() %>%
  mutate(VE_t = list(VE = erlang.decay(VE0 = VE0, T = T, t=1:365))) %>%
  unnest(VE_t) %>%
  mutate(t = rep(1:365, times=dim(resp)[1])) %>%
  group_by(group,t)

##select 1000 samples of VE vs severe RSV disease at age 0 to 11(month)
# prepare list
data_list <- list()

# from age 0 to 11(repeat from t = 15 for 12 times by 30 days)
for (i in seq(15, 365, by = 30)) {
  # sample VE vs severe at t = i
  ve <- ve_post %>%
#    filter(group == "severe") %>%
    filter(group == "Severe") %>%
    filter(t == i) %>%
    mutate(age = (i+15)/30 - 1) #age in month
  # return the sampled VE
  data_list[[(i+15)/30]] <- ve
}
#data_list provides 1000 samples of VE by age in month

#### RSV-associated deaths (from RSV_burden.R)
## 1000 samples of deaths by age in month
#SA
death_sa


#Combine SA
death_comb <- rbind(death_sa)

# countries (SA and KEN)
countries <- unique(death_comb$country)

################ main analysis if VE follows Erlang-2 distribution###########
### Calculate averted deaths by age in month (0to11 month)
dat_list <- list()
benef_by_country <- list()

### Generate 1000 samples for RSV-averted deaths
# calculate averted deaths by country
for (country in countries) {
  # filter by country
  death_comb_country <- death_comb %>% filter(country == !!country)

  # averted deaths at each age in month
  dat_list <- list()
  for (i in seq(1:12)) {
    ve_age <- data_list[[i]] %>% as_tibble() %>% dplyr::select(VE_t)
    death_age <- death_comb_country %>% filter(age == i - 1) %>% dplyr::select(death)
    death_avert_age <- ve_age$VE_t * death_age$death  #VE*RSV-deaths at that age in month
    dat_list[[i]] <- death_avert_age
  }

  #dat_list includes 10000 samples for averted deaths by age in month 
  
  #########

# 10000 samples for infant deaths averted 
  sum_samples <- numeric(num_s)
  for (j in 1:num_s) {
    sample_age <- numeric(12)
    for (i in 1:12) {
      samp_age <- sample(unlist(dat_list[[i]]), 1)
      sample_age[i] <- samp_age #averted deaths at each age in month
    }
    sum_samples[j] <- sum(sample_age) #sum up until 1 year old
  }

  benef_median <- quantile(sum_samples, .5)  #benefit result (median)
  benef_low <- quantile(sum_samples, .025) #lower 95% Crl
  benef_high <- quantile(sum_samples, .975) #upper 95% Crl

  benef_by_country[[country]] <- list(
    median = benef_median,
    low = benef_low,
    high = benef_high
  )
}

### Results ###
#Country-specific RSV-associated infant deaths averted through vaccination (median, 95% Crls)
print(benef_by_country)


############################################################
#### To caluclate benefit-risk ratio and for violin plot, re-generate posterior samples for benefits in ZAF in main scenario####
# Define the target country
target_country <- "ZAF"

# Filter data for the target country
death_comb_country <- death_comb %>% filter(country == target_country)

# Averted death by age in month
dat_list <- list()
for (i in seq(1:12)) {
  ve_age <- data_list[[i]] %>% as_tibble() %>% dplyr::select(VE_t)
  death_age <- death_comb_country %>% filter(age == i - 1) %>% dplyr::select(death)
  death_avert_age <- ve_age$VE_t * death_age$death
  dat_list[[i]] <- death_avert_age
}

# Under 1 year deaths averted by sample
post_benef <- numeric(num_s)
for (j in 1:num_s) {
  sample_age <- numeric(12)
  for (i in 1:12) {
    samp_age <- sample(unlist(dat_list[[i]]), 1)
    sample_age[i] <- samp_age
  }
  post_benef[j] <- sum(sample_age)
}
post_benef%>%median()
# this should be almost the same with benef_by_country
################### use this in benef_risk.R ##############################



###########Sensitivity analysis if VE = 80% for 1 year (Supplement) ###############
#Calculate deaths averted under 1 year (0to11 month)
da_list <- list()

benef_sens <- list()

# calculate averted death by country
for (country in countries) {
  # filter by country
  death_comb_country <- death_comb %>% filter(country == !!country)

  # averted death by month
  da_list <- list()
  for (i in seq(1:12)) {
    ve <- 0.8
    death_age <- death_comb_country %>% filter(age == i - 1) %>% dplyr::select(death)
    death_avert_age <- ve * death_age$death
    da_list[[i]] <- death_avert_age
  }

  #under 1 year deaths averted by sample
  sum_sample <- numeric(num_s)
  for (j in 1:num_s) {
    sample_age <- numeric(12)
    for (i in 1:12) {
      samp_age <- sample(unlist(da_list[[i]]), 1)
      sample_age[i] <- samp_age
    }
    sum_sample[j] <- sum(sample_age)
  }

  bene_median <- quantile(sum_sample, .5)
  bene_low <- quantile(sum_sample, .025)
  bene_high <- quantile(sum_sample, .975)

  benef_sens[[country]] <- list(
    median = bene_median,
    low = bene_low,
    high = bene_high
  )
}


### Results ###
#Country-specific RSV-associated deaths averted
print(benef_sens)


# main analysis
benef_sa <- data.frame(
  country = "ZAF",
  median = benef_by_country$ZAF$median[["50%"]],
  low_2.5 = benef_by_country$ZAF$low[["2.5%"]],
  high_97.5 = benef_by_country$ZAF$high[["97.5%"]]
)

write.csv(benef_sa, "output/benef_main.csv", row.names = FALSE)

# sensitivity analysis
benef_sa_sens <- data.frame(
  country = "ZAF",
  median = benef_sens$ZAF$median[["50%"]],
  low_2.5 = benef_sens$ZAF$low[["2.5%"]],
  high_97.5 = benef_sens$ZAF$high[["97.5%"]]
)

write.csv(benef_sa_sens, "output/benef_sens.csv", row.names = FALSE)




