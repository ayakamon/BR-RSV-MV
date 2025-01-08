### created by Ayaka Monoi



### libraries ###
# library(dplyr)
# library(ggplot2)
# library(ggtext)

### 
# analysis1 = "South Africa"
# analysis2 = "South Africa without 27-week in vaccine arm and 30-week in plaebo arm"
# analysis3 = "South Africa without 5 earliest births in each arm"
# analysis4 = "South Africa with mothers vaccinated at 27-36 weeks"

set.seed(1234)
### Fig S1 ###
# Calculate post_risk for analyses
post_risk_list <- lapply(paste0("riskM", 1:3), function(riskM) {
  calculate_sum(risk_nmr, get(riskM))
})


## data for plot
rr_post_calc <- lapply(post_risk_list, function(results_risk) {
  rr_post <- results_risk/ post_benef
})

### Violin plot ###
plo_data <- data.frame()

for (i in 1:3) {
  data <- rr_post_calc[[i]]
  temp_df <- data.frame(value = data, list_element = i)
  plo_data <- bind_rows(plo_data, temp_df)
}


### Fig S1 ###
ggplot(plo_data, aes(x = factor(list_element, levels = rev(unique(list_element))), y = value)) +
  geom_violin(trim = FALSE, fill = "#4895D1") +
  geom_hline(yintercept = 1) +
  labs(title = "", x = "", y = "Risk-benefit ratio") +
  scale_x_discrete(labels = rev(c("All births in South Africa",
                                  "South Africa without \n1 earliest birth in each arm",
                                  "South Africa without \n5 earliest births in each arm"))) +
  scale_y_continuous(
    breaks = seq(floor(min(plo_data$value)), ceiling(max(plo_data$value)), by = 2), 
    limits = c(floor(min(plo_data$value)), ceiling(max(plo_data$value))) 
  ) +
  theme_minimal() +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 20, angle = 0),
    strip.text = element_text(size = 20),
    axis.text.y = element_text(size = 20),
      axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  )

##save
ggsave(
  filename = "output/figS1.pdf",
  plot = last_plot(),
  device = "pdf",
  width = 11, height = 7,
  units = "in"
)



### Fig S2 ###
# Calculate post_risk for analyses
post_risk_list_fix <- lapply(paste0("riskdif", 1:3), function(riskdif) { 
  calculat_sum(risk_nmr, get(riskdif))
})


## data for plot
rr_post_calc_fix <- lapply(post_risk_list_fix, function(results_risk_fix) {
  rr_post_fix <- results_risk_fix / post_benef
})

### Violin plot ###
plo_data_fix <- data.frame()

for (i in 1:3) {
  data_fix <- rr_post_calc_fix[[i]]
  temp_df <- data.frame(value = data_fix, list_element = i)
  plo_data_fix <- bind_rows(plo_data_fix, temp_df)
}


### Fig S2 ###
ggplot(plo_data_fix, aes(x = factor(list_element, levels = rev(unique(list_element))), y = value)) +
  geom_violin(trim = FALSE, fill = "#4895D1") +
  geom_hline(yintercept = 1) +
  labs(title = "", x = "", y = "Risk-benefit ratio") +
  scale_x_discrete(labels = rev(c("All births in South Africa",
                                  "South Africa without \n1 earliest birth in each arm",
                                  "South Africa without \n5 earliest births in each arm"))) +
  scale_y_continuous(breaks = c(0, 1, 2, 3,4)) +
  theme_minimal() +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 20, angle = 0),
    strip.text = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # strip.text.x = element_text(size = 30),
    # strip.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  )

##save
ggsave(
  filename = "output/figS2.pdf",
  plot = last_plot(),
  device = "pdf",
  width = 11, height = 7,
  units = "in"
)


### Percentage of simulations where benefit>risk ###

#risk per benefit
rr_post_calc_fix 
# rr_post_calc_fix <- lapply(post_risk_list_fix, function(post_risk_fix) {
#   rr_post <- post_risk / post_benef
# })

## count iterations in which Risk-benefit ratio >1
count_over_one_fix <- sum(rr_post_calc_fix[[1]] > 1)
#
## percentage of iterations in which ratio >1
count_over_one_fix/length(rr_post_calc_fix[[1]])*100

#repeat this from 1 to 3
count_fix <- numeric(3)
percentage_fix <- numeric(3)

for (i in 1:3) {
  # count iterations in which risk/benefit >1
  count_fix[i] <- sum(rr_post_calc_fix[[i]] > 1)
  
  # percentage benefit/risk>1
  percentage_fix[i] <- (1- (count_fix[i] / length(rr_post_calc_fix[[i]])) )* 100
}

# results
count_fix
### percentage of iterations in which benefit exceeds risk
{cat("Pecentage of simulations where benefit exceeds risk by scenario without bootstrapping trial birth data\n main analysis, without 1 earliest birth, without 5 earliest births\n"); print(percentage_fix)}





### Fig 4 ###
excess_vacc <- calculate_sum(risk_nmr, riskM4)

rr_vacc = excess_vacc/post_benef
table(rr_vacc>1)
table(rr_vacc>0.2)

##plot
data_vacc_rb = as.data.frame(cbind(risk = excess_vacc, benefit = post_benef))

##South Africa with 27+ vacc window
ggplot(data = data_vacc_rb, aes(x=risk, y=benefit)) + geom_point() +
  geom_abline(slope=1, intercept=0, color = 'darkgrey', lty='dashed', size = 1) +
  geom_abline(slope=5, intercept=0, color = 'darkgrey', lty='dotdash', size = 1) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlim(c(-275, 90)) + ylim(c(0, 50)) +
  labs( x='Risk\n(Estimated excess neonatal deaths \npotentially attributable to vaccine-associated preterm birth)', y = 'Benefit \n(Estimated vaccine-preventable RSV-associated infant deaths)') +
  theme_minimal()  +
  theme(axis.line = element_line(colour = "grey50", linewidth = 1))+
  theme(
    axis.title.x = element_text(size = 19, lineheight = 1.2),
    axis.title.y = element_text(size = 19, lineheight = 1.2) ,
    axis.text.x=element_text(size=18),
    axis.text.y=element_text(size=18),
  ) +
  coord_flip()

### save
ggsave(
  filename = "output/fig4.pdf",
  plot = last_plot(),
  device = "pdf",
  width = 11, height = 8,
  units = "in"
)

