### created by Ayaka Monoi 




set.seed(1234)
###############################risk-benefit ratio####################
### check the median of benefit and risk
# Posterior samples for RSV-deaths averted (benefit) from benef_calc.R
post_benef%>%median()
post_benef%>%quantile(.025)
post_benef%>%quantile(.975)

# Posterior samples for excess deaths (risk) from risk_calc.R
#### posterior distribution of excess deaths
# main analysis
post_risk <- calculate_sum(risk_nmr, riskM1) #change this from riskM1 to riskM5
post_risk%>%median()
post_risk%>%quantile(.025)
post_risk%>%quantile(.975)


#### benefit-risk ratio ####
#with resampled births in trial

# excess death per 1 live saved
rr_post <- post_risk/post_benef
# results
rr_post%>%median()
rr_post%>%quantile(.025)
rr_post%>%quantile(.975)



##################################
#### risk-benefit ratio (Table 1)####

# Calculate post_risk for analyses
post_risk_list <- lapply(paste0("riskM", 1:4), function(riskM) {
  calculate_sum(risk_nmr, get(riskM))
})

# Compute rr_post for each post_risk and calculate the median and 95% CrIs
rr_post_stats <- lapply(post_risk_list, function(post_risk) {
  rr_post <- post_risk / post_benef
  data.frame(
    median = median(rr_post),
    ci_low = quantile(rr_post, 0.025),
    ci_high = quantile(rr_post, 0.975)
  )
})


# Combine results into a data frame
rr_post_df <- do.call(rbind, rr_post_stats)
rr_post_df$analysis <- paste0("riskM", 1:4)
 
{cat("Risk-benefit ratio by scenario\n"); print(rr_post_df)}

#### excess deaths per one life saved for all analyses ####



# Compute risk-benfeit difference for each analysis and calculate the median and 95% CrIs
rd_post_stats <- lapply(post_risk_list, function(post_risk) {
  rd_post <- (post_risk - post_benef)
  data.frame(
    median = median(rd_post),
    ci_low = quantile(rd_post, 0.025),
    ci_high = quantile(rd_post, 0.975)
  )
})

#### excess deaths minus one life saved for all analyses ####


# Combine results into a data frame
rd_post_df <- do.call(rbind, rd_post_stats)
rd_post_df$analysis <- paste0("riskM", 1:4)
{cat("Risk-benefit difference by scenario\n"); print(rd_post_df)}
#### excess deaths per one life saved for all analyses ####



#### proportion of risk/benefit >1 ####
#risk per benefit
rr_post_calc <- lapply(post_risk_list, function(post_risk) {
  rr_post <- post_risk / post_benef
})

## count iterations in which ratio >1
count_over_one <- sum(rr_post_calc[[4]] > 1)
#
## percentage of iterations in which ratio >1
count_over_one/length(rr_post_calc[[1]])*100

#repeat this from 1 to 4
count <- numeric(4)
percentage <- numeric(4)

for (i in 1:4) {
  # count iterations in which risk/benefit >1
  count[i] <- sum(rr_post_calc[[i]] > 1)

  # percentage benefit/risk>1
  percentage[i] <- (1- (count[i] / length(rr_post_calc[[i]])) )* 100
}

# results
count
### percentage of iterations in which benefit exceeds risk
{cat("Pecentage of simulations where benefit exceeds risk by scenario\n main analysis, withou 1 earliest birth, without 5 earliest births, 27-36 weeks vaccination\n"); print(percentage)}

#### proportion of risk/benefit >5 ####
#risk per benefit
rr_post_calc_5 <- lapply(post_risk_list, function(post_risk) {
  rr_post_5 <- 5*post_risk / post_benef
})

## count iterations in which ratio >1
count_over_5 <- sum(rr_post_calc[[4]] > 1)
#
## percentage of iterations in which ratio >1
count_over_5/length(rr_post_calc_5[[1]])*100

#repeat this from 1 to 4
count_5 <- numeric(4)
percentage_5 <- numeric(4)

for (i in 1:4) {
  # count iterations in which risk/benefit >1
  count_5[i] <- sum(rr_post_calc_5[[i]] > 1)

  # percentage benefit/risk>1
  percentage_5[i] <- (1- (count_5[i] / length(rr_post_calc_5[[i]])) )* 100
}

# results
count_5
 ### percentage of iterations in which benefit exceeds risk
{cat("Pecentage of simulations where benefit exceeds 5 times risk by scenario\n main analysis, withou 1 earliest birth, without 5 earliest births, 27-36 weeks vaccination\n"); print(percentage_5)}


######### Table 1 ends ###################


#########################
#### Comparison between modelled risk and observed risk in South African component of trial ####
#### Applying the risk to the number of participants in SA in trial ####
#### from ptb_posterior.R
# posterior probabilities of GA-specific delivery
# for placebo
prp_p <- myM/sum(dt$ptb)
# for vaccine
prp_v <- myM_v/sum(dt_v$ptb)



# function to calculate posterior samples of excess neonatal deaths among vaccine group compared to unvaccinated group
calculate_death <- function(risk_nmr, matrix) {
  results <- numeric(length = num_sm)

  for (i in 1:num_sm) {
    nmr_risk <- as.vector(risk_nmr[[i]])
    ptb_riskdif <- as.vector(matrix[i,])
    result <- sum(nmr_risk *  ptb_riskdif)

    results[i] <- result
  }
  return(results)
}

####if apply the model to the number of mothers in the trial arms
#South Africa
total_v = 469
total_p = 471


## results
#vaccine arm
{calculate_death(risk_nmr, prp_v)*total_v}%>%median()
{calculate_death(risk_nmr, prp_v)*total_v}%>%quantile(.025)
{calculate_death(risk_nmr, prp_v)*total_v}%>%quantile(.975)

#placebo arm
{calculate_death(risk_nmr, prp_p)*total_p}%>%median()
{calculate_death(risk_nmr, prp_p)*total_p}%>%quantile(.025)
{calculate_death(risk_nmr, prp_p)*total_p}%>%quantile(.975)



#### Benefits
#### if benefit applied to SA in trials####


#Calculate deaths averted under 1 year (0to11 month)
da_list <- list()
bene_by_country <- list()

# calculate averted deaths by country
for (country in countries) {
  # RSV-associated deaths by age in month by country
  death_comb_country <- death_comb %>% filter(country == !!country)

  # expected deaths among vaccinees by age in month
  da_list <- list()
  for (i in seq(1:12)) {
    ve_age <- data_list[[i]] %>% as_tibble() %>% dplyr::select(VE_t) #VE at each age
    death_age <- death_comb_country %>% filter(age == i - 1) %>% dplyr::select(death) #death at each age
    rsv_death_age <- (1-ve_age$VE_t) * death_age$death #expected deaths after vaccination at each age
    da_list[[i]] <- rsv_death_age
  }

  #posterior samples for averted infant deaths
  su_samples <- numeric(1000)
  for (j in 1:1000) {
    sample_age <- numeric(12)
    for (i in 1:12) {
      samp_age <- sample(unlist(da_list[[i]]), 1)
      sample_age[i] <- samp_age #summarise by age in month
    }
    su_samples[j] <- sum(sample_age)
  }

  bene_median <- quantile(su_samples, .5)
  bene_low <- quantile(su_samples, .025)
  bene_high <- quantile(su_samples, .975)

  bene_by_country[[country]] <- list(
    median = bene_median,
    low = bene_low,
    high = bene_high
  )
}

#Expected RSV-associated deaths among vaccinees per 100000 by country
bene_by_country


####if applied the model to the number of participants in South African component of the trial
#Expected RSV deaths in vaccine arm
bene_by_country$ZAF$median/100000*total_v
bene_by_country$ZAF$low/100000*total_v
bene_by_country$ZAF$high/100000*total_v

#averted deaths per the participants in vaccine arm
benef_by_country$ZAF$median/100000*total_v
benef_by_country$ZAF$low/100000*total_v
benef_by_country$ZAF$high/100000*total_v



#########

## Expected deaths in placebo arm (without vaccination)
#Calculate deaths under 1 year (0to11 month)

ben_by_country <- list()

for (country in countries) {
  # filter by country
  death_comb_country <- death_comb %>% filter(country == !!country)
  # averted death by month
  d_list <- list()
  for (i in seq(1:12)) {
    ve_age <- data_list[[i]] %>% as_tibble() %>% dplyr::select(VE_t)
    death_age <- death_comb_country %>% filter(age == i - 1) %>% dplyr::select(death)
    rsv_death_age <-  death_age$death
    d_list[[i]] <- rsv_death_age
  }

  #under 1 year deaths averted by sample
  s_samples <- numeric(1000)
  for (j in 1:1000) {
    sampl_age <- numeric(12)
    for (i in 1:12) {
      sam_age <- sample(unlist(d_list[[i]]), 1)
      sampl_age[i] <- sam_age
    }
    s_samples[j] <- sum(sampl_age)
  }

  ben_median <- quantile(s_samples, .5)
  ben_low <- quantile(s_samples, .025)
  ben_high <- quantile(s_samples, .975)

  ben_by_country[[country]] <- list(
    median = ben_median,
    low = ben_low,
    high = ben_high
  )
}

#Expected RSV-associated deaths among un-vaccinees per 100000 by country
ben_by_country


####if applied the model to the number of mothers in the trial arms
#RSV-associated deaths in placebo arm
ben_by_country$ZAF$median/100000*total_p
ben_by_country$ZAF$low/100000*total_p
ben_by_country$ZAF$high/100000*total_p

