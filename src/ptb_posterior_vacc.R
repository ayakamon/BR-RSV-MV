### created by Ayaka Monoi on 2024/9/4


### libraries ###
# library(MCMCpack)
# library(tidyr)
# library(dplyr)
# library(ggplot2)
# library(knitr)
# library(conflicted)
# conflicts_prefer(dplyr::select)
# conflicts_prefer(dplyr::filter)

### Data ###
### Number of births by trial arm by GA in Pfizer's Phase 3 trials (downloaded from the WHO website: https://terrance.who.int/mediacentre/data/sage/SAGE_Slidedeck_September-2024.pdf. on October 9, 2024 )########
set.seed(1234)
################ main analysis ###################################
## All births in South Africa
# vaccine arm
read.csv(here("data", "GA_birth_by_week_country_vacc_window.csv"))%>%
  filter(COUNTRY=="ZAF",
         ARM=='RSVpreF 120 ug') -> vacc_v
sum(vacc_v$ptb) #check the total number of births

# posterior distribution for vaccine
myM_v4 = matrix(NA,nrow=num_sm, ncol = length(27:44))
for (i in 1:num_sm)
  myM_v4[i,] = sample(vacc_v$GA,sum(vacc_v$ptb), replace=T, prob = vacc_v$ptb) %>%
  factor(levels =27:44) %>%
  table()


#placebo arm
read.csv(here("data", "GA_birth_by_week_country_vacc_window.csv")) %>%
  filter(COUNTRY=="ZAF",
         ARM=='Placebo')  -> vacc

#posterior distribution for placebo
myM4 = matrix(NA,nrow=num_sm, ncol = length(27:44))
for (i in 1:num_sm)
  myM4[i,] = sample(vacc$GA,sum(vacc$ptb), replace=T, prob = vacc$ptb) %>%
  factor(levels =27:44) %>%
  table()


#proportion probabilities for delivery
#placebo
prp_p4 <- myM4/sum(vacc$ptb)
#vaccine
prp_v4 <- myM_v4/sum(vacc_v$ptb)

# risk diff between arms
riskM4 <- matrix(NA, nrow=num_sm, ncol = length(27:44))
for(i in 1:num_sm){
  rd_ptb = prp_v4[i,] - prp_p4[i,]
  riskM4[i,] = rd_ptb
}

#10000 samples for risk diff in GA sets
riskM4
{cat("Estimated neonatal deaths with 27-36 vaccination")}

print(calculate_sum(risk_nmr, riskM4)%>%median())
print(calculate_sum(risk_nmr, riskM4)%>%quantile(.025))
print(calculate_sum(risk_nmr, riskM4)%>%quantile(.975))




#### Fig 3B ####

#posterior samples for placebo
vac_df_placebo <- data.frame(myM4)
colnames(vac_df_placebo) <- c(27:44)
#posterior samples for vaccine
vac_df_vaccine <- data.frame(myM_v4)
colnames(vac_df_vaccine) <-c(27:44)

# calculate median and 95%Crl
vac_stat_placebo <- data.frame(
  GA = c(27:44),
  Median = apply(vac_df_placebo, 2, median, na.rm = TRUE),
  Quantile_2.5 = apply(vac_df_placebo, 2, quantile, probs = 0.025, na.rm = TRUE),
  Quantile_97.5 = apply(vac_df_placebo, 2, quantile, probs = 0.975, na.rm = TRUE),
  ARM = "Placebo"
)

vac_stat_vaccine <- data.frame(
  GA = c(27:44),
  Median = apply(vac_df_vaccine, 2, median, na.rm = TRUE),
  Quantile_2.5 = apply(vac_df_vaccine, 2, quantile, probs = 0.025, na.rm = TRUE),
  Quantile_97.5 = apply(vac_df_vaccine, 2, quantile, probs = 0.975, na.rm = TRUE),
  ARM = "RSVpreF 120μg"
)

# combine data
combined_data_vac <- rbind(vac_stat_placebo, vac_stat_vaccine)


##################

### Plot Fig 3B
# change values in  placebo to negative
combined_data_vac <- combined_data_vac %>%
  mutate(Median = ifelse(ARM == "Placebo", -Median, Median))

# calculate diff of medians between placebo and vaccine
difference_data_vac <- combined_data_vac %>%
  group_by(GA) %>%
  summarise(Difference = sum(Median[ARM == "RSVpreF 120μg"]) + sum(Median[ARM == "Placebo"]))

#### plot
ggplot() +
  geom_bar(data = combined_data_vac, aes(x = GA, y = Median, fill = ARM), stat = "identity", position = "stack") +
  geom_errorbar(data = combined_data_vac, aes(x = GA, ymin = ifelse(ARM == "Placebo", -Quantile_97.5, Quantile_2.5),
                                          ymax = ifelse(ARM == "Placebo", -Quantile_2.5, Quantile_97.5)), width = 0) +
  geom_bar(data = difference_data_vac, aes(x = GA, y = Difference, fill = factor(sign(Difference), labels = c("negative", "zero", "positive"))), stat = "identity", alpha = 0.5) +
  xlab("Gestational age at birth (weeks)") +
  ylab("Number of births in the trial") +
  scale_y_continuous(breaks = seq(-150, 150, by = 30), limits = c(-160, 160),
                     labels = function(x) ifelse(x < 0, -x, x)) +
   scale_x_continuous(breaks = seq(27, 44, by = 1)) +
 scale_fill_manual(name = "", values = c("Placebo" = "#696969", "RSVpreF 120μg" = "gray", "positive" = "#FFDD44", "negative" = "#005AB5"), labels = c("Excess births in placebo arm",  "Excess births in intervention arm", "Placebo", "Intervention")) +
   guides(fill = guide_legend(reverse = TRUE)) +
  theme_minimal() +
  theme(
    legend.position = c(0.3, 0.8),
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    panel.grid.major.x = element_line(color = "grey", size = 0.1), 
    panel.grid.minor.x = element_blank() 
  ) +
  annotate("text", x = 26, y = 160, 
           label = "(B)", size = 8, hjust = 0)
  

### save as pdf (Fig 3B)###
ggsave(
  filename = "output/fig3B.pdf",
  plot = last_plot(),
  device = "pdf",
  width = 11, height = 7,
  units = "in"
)

######### Plot Fig 3C##############

#### check which GA has large impact on risk estimates ####

#### with resampled ptb
ga_nmr_diff <- function(risk_nmr, riskM) {
  results <- list()
  for (i in 1:num_sm) {
    nmr_risk <- risk_nmr[[i]]
    ptb_riskdif <- riskM[i,]
    #result <- (nmr_risk *  ptb_riskdif) * 1000 #per 1000 live births
    # results[[i]] <- result
    results[[i]] <- sapply(1:18, function(k) nmr_risk[k]*ptb_riskdif[k]*100000) #per 100000 live births
  }
  return(results)
}

#posterior samples of deaths born at each GA
smp <- ga_nmr_diff(risk_nmr,riskM4)

#one by one
i=11
elements <- sapply(smp, function(x) x[i])
stat_summary <- list(
  median = median(elements),
  quantile_2.5 = quantile(elements, probs = 0.025, na.rm=T),
  quantile_97.5 = quantile(elements, probs = 0.975, na.rm=T)
)
stat_summary

stat_summaries_vac <- lapply(1:18, function(k) {
  elements <- sapply(smp, function(x) x[k])
  list(
    median = median(elements),
    quantile_2.5 = quantile(elements, probs = 0.025, na.rm=T),
    quantile_97.5 = quantile(elements, probs = 0.975, na.rm=T)
  )
})

#median and 95%CrI by GA week
stat_summaries_vac

#check whether the sum is close to the risk estimate
total_sum <- 0
for (i in 1:18) {
  total_sum <- total_sum + stat_summaries_vac[[i]]$median
}


     

#### check ends #####

#### plot Fig 3C ######

# dataframe
data_vac <- data.frame(
  index = 1:18,
  median = c(
    stat_summaries_vac[[1]]$median,
    stat_summaries_vac[[2]]$median,
    stat_summaries_vac[[3]]$median,
    stat_summaries_vac[[4]]$median,
    stat_summaries_vac[[5]]$median,
    stat_summaries_vac[[6]]$median,
    stat_summaries_vac[[7]]$median,
    stat_summaries_vac[[8]]$median,
    stat_summaries_vac[[9]]$median,
    stat_summaries_vac[[10]]$median,
    stat_summaries_vac[[11]]$median,
    stat_summaries_vac[[12]]$median,
    stat_summaries_vac[[13]]$median,
    stat_summaries_vac[[14]]$median,
    stat_summaries_vac[[15]]$median,
    stat_summaries_vac[[16]]$median,
    stat_summaries_vac[[17]]$median,
    stat_summaries_vac[[18]]$median
  ),
  quantile_2.5 = c(
    stat_summaries_vac[[1]]$quantile_2.5,
    stat_summaries_vac[[2]]$quantile_2.5,
    stat_summaries_vac[[3]]$quantile_2.5,
    stat_summaries_vac[[4]]$quantile_2.5,
    stat_summaries_vac[[5]]$quantile_2.5,
    stat_summaries_vac[[6]]$quantile_2.5,
    stat_summaries_vac[[7]]$quantile_2.5,
    stat_summaries_vac[[8]]$quantile_2.5,
    stat_summaries_vac[[9]]$quantile_2.5,
    stat_summaries_vac[[10]]$quantile_2.5,
    stat_summaries_vac[[11]]$quantile_2.5,
    stat_summaries_vac[[12]]$quantile_2.5,
    stat_summaries_vac[[13]]$quantile_2.5,
    stat_summaries_vac[[14]]$quantile_2.5,
    stat_summaries_vac[[15]]$quantile_2.5,
    stat_summaries_vac[[16]]$quantile_2.5,
    stat_summaries_vac[[17]]$quantile_2.5,
    stat_summaries_vac[[18]]$quantile_2.5
  ),
  quantile_97.5 = c(
    stat_summaries_vac[[1]]$quantile_97.5,
    stat_summaries_vac[[2]]$quantile_97.5,
    stat_summaries_vac[[3]]$quantile_97.5,
    stat_summaries_vac[[4]]$quantile_97.5,
    stat_summaries_vac[[5]]$quantile_97.5,
    stat_summaries_vac[[6]]$quantile_97.5,
    stat_summaries_vac[[7]]$quantile_97.5,
    stat_summaries_vac[[8]]$quantile_97.5,
    stat_summaries_vac[[9]]$quantile_97.5,
    stat_summaries_vac[[10]]$quantile_97.5,
    stat_summaries_vac[[11]]$quantile_97.5,
    stat_summaries_vac[[12]]$quantile_97.5,
    stat_summaries_vac[[13]]$quantile_97.5,
    stat_summaries_vac[[14]]$quantile_97.5,
    stat_summaries_vac[[15]]$quantile_97.5,
    stat_summaries_vac[[16]]$quantile_97.5,
    stat_summaries_vac[[17]]$quantile_97.5,
    stat_summaries_vac[[18]]$quantile_97.5
  )
)

#add GA
data_vac$label <- 27:44

#change color
ggplot(data_vac, aes(x = factor(label), y = median )) +
  geom_bar(stat = "identity", aes(fill = ifelse(median < 0, "Negative", "Positive"))) +
  geom_errorbar(aes(ymin = quantile_2.5 , ymax = quantile_97.5 ), width = 0) +
  scale_y_continuous(breaks = seq(-100, 200, by = 20), labels = scales::comma) +
  labs(x = "Gestational age at birth (weeks)", y = "Excess neonatal deaths \n(per 100,000 live births)", fill = " ") +
  scale_fill_manual(values = c("Negative" = "#0072B2", "Positive" = "#FFC20A"),
                    labels = c("Negative" = "Excess deaths estimated among controls", "Positive" = "Excess deaths estimated among vaccinees"))+
   guides(fill = guide_legend(reverse = TRUE)) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24),
    plot.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 24),
    legend.title = element_text(size = 24),
    legend.position = c(.6,0.2),
    panel.grid.major.x = element_line(color = "grey", size = 0.1), 
    panel.grid.minor.x = element_blank() 
  )+
  annotate("text", x = 0, y = 45, 
           label = "(C)", size = 8, hjust = 0)


### save as pdf (Fig 3C)###
ggsave(
  filename = "output/fig3C.pdf",       
  plot = last_plot(),           
  device = "pdf",              
  width = 11, height = 7,       
  units = "in"                  
)
