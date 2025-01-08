### created by Ayaka Monoi 



#### libraries ####
# library(ggplot2)
set.seed(1234)
### Before this, nmr.R(estimates of neonatal mortality(NMR) in South Africa) and ptb_posterior.R (GA distributions by trial arm) need to be run.#####

#### NMR (nmr.R)###########
#To check GA weeks included in data of South Africa's cohort study
dat_dch %>% mutate(GA=Tmin/7+1)

#### Generate posterior samples of NMR by GA
##from nmr.R after run MCMC
a <-nmr_bay %>%
  rowwise() %>%
  mutate(nmr_t = list(nmr = erlang.func(NMR0 = NMR0, T = T, t=min(dat_dch$Tmin):max(dat_dch$Tmax))))

#number of samples
b <- sample(a$nmr_t,num_sm) #number of samples defined in ptb_posterior.R

#posterior samples for estimated NMR at each GA week
risk_nmr <- list()
for (i in 1:num_sm) {
  risk_nmr[[i]] <- b[i]$nmr[c(seq(3, 126, by = 7))] #from 27 weeks to 44 weeks
}




#### Main and Sensitivity analyses (on early preterm infants) ####
## after ptb_posterior.R###


#risk difference of delivery from posterior samples (ptb_posterior.R)
riskM1 #main; all data in South Africa
riskM2 #wo 1 earliest birth in each arm in South Africa
riskM3 #wo 5 earliest births in each arm in South Africa

## combine posterior probabilities for delivery for main and sensitivity analyses
riskM_list <- list(riskM1, riskM2, riskM3)


#########################

#function
calculate_sum <- function(risk_nmr, matrix) {
  results <- numeric(length = num_sm)

  for (i in 1:num_sm) {
    nmr_risk <- as.vector(risk_nmr[[i]])
    ptb_riskdif <- as.vector(matrix[i,])
    result <- sum(nmr_risk *  ptb_riskdif) * 100000 #per 100000 live births

    results[i] <- result
  }
  return(results)
}

#results list
results_list <- list()

for (i in 1:length(riskM_list)) {
  matrix <- riskM_list[[i]]
  results <- calculate_sum(risk_nmr, matrix)

  risk_median <- quantile(results, probs = 0.5)
  risk_low <- quantile(results, probs = 0.025)
  risk_high <- quantile(results, probs = 0.975)

  results_list[[paste0("riskM", i)]] <- list(
    median = risk_median,
    low = risk_low,
    high = risk_high
  )
}


##estimated deaths potentially attributable to prematurity

{cat("Estimated neonatal deaths by scenario \n(riskM1:main analysis, riskM2: without 1 earliest birth, riskM3:without 5 earliest births)\n"); print(results_list)}


####to calculate benefit-risk ratio, go to "benef_risk.R"






########################Sensitivity analyses with fixed ptb (without resampling)#########
## after ptb_posterior.R###

# #risk difference of delivery from posterior samples (ptb.posterior.R)
# riskdif1 #main; all data in SA
# riskdif2 #wo 27w and 30w in SA
# riskdif3 #wo earliest 5 in each arm in SA

#
#
# ##############with fixed preterm birth risk###########
# ## combine risk difference of GA-delivery by sensitivity analysis
riskdif_list <- list(riskdif1, riskdif2, riskdif3)
#
# #calculate risk with fixed risk difference
calculat_sum <- function(risk_nmr, risk_ptb) {
  results <- numeric(length = num_sm)

  for (i in 1:num_sm) {
    nmr_risk <- as.vector(risk_nmr[[i]])
    result <- sum(nmr_risk *  risk_ptb) * 100000 #per 100000 live births

    results[i] <- result
  }
  return(results)
}

#all births in SA
calculat_sum(risk_nmr, riskdif1)%>%quantile(.5)
calculat_sum(risk_nmr, riskdif1)%>%quantile(.025)
calculat_sum(risk_nmr, riskdif1)%>%quantile(.975)


# ############posterior distribution of excess deaths########
 post_risk_fixed <- calculat_sum(risk_nmr, riskdif1)
#
#
# ###############

#results list
results_fix_list <- list()

#
for (i in 1:length(riskdif_list)) {
  riskdif <- riskdif_list[[i]]
  results <- calculat_sum(risk_nmr, riskdif)
  risk_median <- quantile(results, probs = 0.5)
  risk_low <- quantile(results, probs = 0.025)
  risk_high <- quantile(results, probs = 0.975)

  results_fix_list[[paste0("riskdif", i)]] <- list(
    median = risk_median,
    low = risk_low,
    high = risk_high
  )
}

# ##estimated deaths potentially attributable to prematurity with fixed ptb
results_fix_list %>%
{cat("Excess neonatal deaths(risk) without bootstrapping GA\n"); print(.)}


################Go to benefit_risk.R######





######### Plot Fig 2C##############

#### check which GA has large impact on risk estimates ####
#### with resampled ptb
ga_nmr_diff <- function(risk_nmr, riskM) {
  results <- list()
  for (i in 1:num_sm) {
    nmr_risk <- risk_nmr[[i]]
    ptb_riskdif <- riskM[i,]
    results[[i]] <- sapply(1:18, function(k) nmr_risk[k]*ptb_riskdif[k]*100000) #per 100000 live births
  }
  return(results)
}

#posterior samples of deaths born at each GA
sm <- ga_nmr_diff(risk_nmr,riskM1)

#one by one
i=1
elements <- sapply(sm, function(x) x[i])
stat_summary <- list(
  median = median(elements),
  quantile_2.5 = quantile(elements, probs = 0.025, na.rm=T),
  quantile_97.5 = quantile(elements, probs = 0.975, na.rm=T)
)
stat_summary

stat_summaries <- lapply(1:18, function(k) {
  elements <- sapply(sm, function(x) x[k])
  list(
    median = median(elements),
    quantile_2.5 = quantile(elements, probs = 0.025, na.rm=T),
    quantile_97.5 = quantile(elements, probs = 0.975, na.rm=T)
  )
})

#median and 95%CrI by GA week
stat_summaries

#check whether the sum is close to the risk estimate
total_sum <- 0
for (i in 1:18) {
  total_sum <- total_sum + stat_summaries[[i]]$median
}
total_sum
#### check ends #####

#### plot Fig 2C######

# dataframe
data <- data.frame(
  index = 1:18,
  median = c(
    stat_summaries[[1]]$median,
    stat_summaries[[2]]$median,
    stat_summaries[[3]]$median,
    stat_summaries[[4]]$median,
    stat_summaries[[5]]$median,
    stat_summaries[[6]]$median,
    stat_summaries[[7]]$median,
    stat_summaries[[8]]$median,
    stat_summaries[[9]]$median,
    stat_summaries[[10]]$median,
    stat_summaries[[11]]$median,
    stat_summaries[[12]]$median,
    stat_summaries[[13]]$median,
    stat_summaries[[14]]$median,
    stat_summaries[[15]]$median,
    stat_summaries[[16]]$median,
    stat_summaries[[17]]$median,
    stat_summaries[[18]]$median
  ),
  quantile_2.5 = c(
    stat_summaries[[1]]$quantile_2.5,
    stat_summaries[[2]]$quantile_2.5,
    stat_summaries[[3]]$quantile_2.5,
    stat_summaries[[4]]$quantile_2.5,
    stat_summaries[[5]]$quantile_2.5,
    stat_summaries[[6]]$quantile_2.5,
    stat_summaries[[7]]$quantile_2.5,
    stat_summaries[[8]]$quantile_2.5,
    stat_summaries[[9]]$quantile_2.5,
    stat_summaries[[10]]$quantile_2.5,
    stat_summaries[[11]]$quantile_2.5,
    stat_summaries[[12]]$quantile_2.5,
    stat_summaries[[13]]$quantile_2.5,
    stat_summaries[[14]]$quantile_2.5,
    stat_summaries[[15]]$quantile_2.5,
    stat_summaries[[16]]$quantile_2.5,
    stat_summaries[[17]]$quantile_2.5,
    stat_summaries[[18]]$quantile_2.5
  ),
  quantile_97.5 = c(
    stat_summaries[[1]]$quantile_97.5,
    stat_summaries[[2]]$quantile_97.5,
    stat_summaries[[3]]$quantile_97.5,
    stat_summaries[[4]]$quantile_97.5,
    stat_summaries[[5]]$quantile_97.5,
    stat_summaries[[6]]$quantile_97.5,
    stat_summaries[[7]]$quantile_97.5,
    stat_summaries[[8]]$quantile_97.5,
    stat_summaries[[9]]$quantile_97.5,
    stat_summaries[[10]]$quantile_97.5,
    stat_summaries[[11]]$quantile_97.5,
    stat_summaries[[12]]$quantile_97.5,
    stat_summaries[[13]]$quantile_97.5,
    stat_summaries[[14]]$quantile_97.5,
    stat_summaries[[15]]$quantile_97.5,
    stat_summaries[[16]]$quantile_97.5,
    stat_summaries[[17]]$quantile_97.5,
    stat_summaries[[18]]$quantile_97.5
  )
)

#add GA
data$label <- 27:44

#change color
ggplot(data, aes(x = factor(label), y = median )) +
  geom_bar(stat = "identity", aes(fill = ifelse(median < 0, "Negative", "Positive"))) +
  geom_errorbar(aes(ymin = quantile_2.5 , ymax = quantile_97.5 ), width = 0) +
  scale_y_continuous(breaks = seq(-100, 200, by = 20), labels = scales::comma) +
  labs(x = "Gestational age at birth (weeks)", y = "Excess neonatal deaths \n(per 100,000 live births)", fill = " ") +
  scale_fill_manual(values = c("Negative" = "#0C7BDC", "Positive" = "#FFC20A"),
                    labels = c("Negative" = "Excess deaths estimated among controls", "Positive" = "Excess deaths estimated among vaccinees"))+
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24),
    plot.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    legend.position = c(.6,0.7)
  )+
  annotate("text", x = 0, y = 210, 
           label = "(C)", size = 8, hjust = 0)

# ### save as pdf (Fig 2C)###
ggsave(
  filename = "output/fig2C.pdf",
  plot = last_plot(),
  device = "pdf",
  width = 11, height = 7,
  units = "in"
)


#proportion of deaths of infants born at 27wks out of overall excess death
stat_summaries[[1]]$median/results_list$riskM1$median*100




