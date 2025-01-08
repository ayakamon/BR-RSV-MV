### created by Ayaka Monoi 


### libraries ###
# library(MCMCpack)
# library(tidyr)
# library(dplyr)
# library(ggplot2)
# library(knitr)
# library(conflicted)
# library(here)
# 
# conflicts_prefer(dplyr::select)
# conflicts_prefer(dplyr::filter)

set.seed(1234)
### Data ###

################ main analysis ###################################
## All births in South Africa
#vaccine arm
dt_v <- read.csv(here("data", 'GA_birth_by_week_country.csv')) %>%
  filter(COUNTRY=="ZAF",ARM=='RSVpreF 120 ug')
sum(dt_v$ptb) #check the total number of births

#placebo arm
dt <- read.csv(here("data", 'GA_birth_by_week_country.csv')) %>%
  filter(COUNTRY=="ZAF",
         ARM=='Placebo')
sum(dt$ptb)


#### resampling ####
num_sm = 10000 

## placebo arm
myM = matrix(NA,nrow=num_sm, ncol = length(27:44))
for (i in 1:num_sm)
  myM[i,] = sample(dt$GA,sum(dt$ptb), replace=T, prob = dt$ptb) %>%
  factor(levels =27:44) %>%
  table()

## vaccine arm
myM_v = matrix(NA,nrow=num_sm, ncol = length(27:44))
for (i in 1:num_sm)
  myM_v[i,] = sample(dt_v$GA,sum(dt_v$ptb), replace=T, prob = dt_v$ptb) %>%
  factor(levels =27:44) %>%
  table()

#### Visualisation of resampled births in trial (Fig 2B) ####

#posterior samples for placebo arm
df_placebo <- data.frame(myM)
colnames(df_placebo) <- c(27:44)

#posterior samples for vaccine arm
df_vaccine <- data.frame(myM_v)
colnames(df_vaccine) <- c(27:44)

# median and 95%Crl
#placebo
stat_placebo <- data.frame(
  GA = c(27:44),
  Median = apply(df_placebo, 2, median, na.rm = TRUE),
  Quantile_2.5 = apply(df_placebo, 2, quantile, probs = 0.025, na.rm = TRUE),
  Quantile_97.5 = apply(df_placebo, 2, quantile, probs = 0.975, na.rm = TRUE),
  ARM = "Placebo"
)

#vaccine
stat_vaccine <- data.frame(
  GA = c(27:44),
  Median = apply(df_vaccine, 2, median, na.rm = TRUE),
  Quantile_2.5 = apply(df_vaccine, 2, quantile, probs = 0.025, na.rm = TRUE),
  Quantile_97.5 = apply(df_vaccine, 2, quantile, probs = 0.975, na.rm = TRUE),
  ARM = "RSVpreF 120μg"
)

# combine data
combined_data <- rbind(stat_placebo, stat_vaccine)


## Plot
# change values in  placebo to negative
combined_data <- combined_data %>%
  mutate(Median = ifelse(ARM == "Placebo", -Median, Median))

# calculate diff of medians between placebo and vaccine
difference_data <- combined_data %>%
  group_by(GA) %>%
  summarise(Difference = sum(Median[ARM == "RSVpreF 120μg"]) + sum(Median[ARM == "Placebo"]))

#### plot Fig 2B ###
ggplot() +
  geom_bar(data = combined_data, aes(x = GA, y = Median, fill = ARM), stat = "identity", position = "stack") +
  geom_errorbar(data = combined_data, aes(x = GA, ymin = ifelse(ARM == "Placebo", -Quantile_97.5, Quantile_2.5),
                                          ymax = ifelse(ARM == "Placebo", -Quantile_2.5, Quantile_97.5)), width = 0) +
  geom_bar(data = difference_data, aes(x = GA, y = Difference, fill = factor(sign(Difference), labels = c("negative", "zero", "positive"))), stat = "identity", alpha = 0.5) +
  xlab("Gestational age at birth (weeks)") +
  ylab("Number of births in the trial") +
  scale_y_continuous(breaks = seq(-150, 150, by = 30), limits = c(-160, 160),
                     labels = function(x) ifelse(x < 0, -x, x)) +
  scale_x_continuous(breaks = seq(27, 44, by = 1)) +
  scale_fill_manual(name = "", values = c("Placebo" = "#696969", "RSVpreF 120μg" = "gray", "positive" = "#FFDD44", "negative" = "#005AB5"), labels = c("Excess births in placebo arm",  "Excess births in intervention arm", "Placebo", "Vaccine")) +
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




### save as pdf (Fig 2B)###
ggsave(
  filename = "output/fig2B.pdf",
  plot = last_plot(),
  device = "pdf",
  width = 11, height = 7,
  units = "in"
)



#### Generate posterior distributions for GA distributions in the trial####

#proportion risk
prp_p <- myM/sum(dt$ptb)
prp_v <- myM_v/sum(dt_v$ptb)


riskM1 <- matrix(NA, nrow=num_sm, ncol = length(27:44))
for(i in 1:num_sm){
  rd_ptb = prp_v[i,] - prp_p[i,]
  riskM1[i,] = rd_ptb
}

#10000 samples for risk diff in GA sets
riskM1




#### Sensitivity analysis ####

#### without birth at 27 week in vaccine arm and 30 week in placebo arm in South African trial component ####

#### data prep ####
# vaccine
dt_v2 <- dt_v
dt_v2$ptb[1] <-0 #27 wks

# placebo
dt2 <- dt
dt2$ptb[4]<-0 #30 wks

#posterior distribution for vaccine
myM_v2 = matrix(NA,nrow=num_sm, ncol = length(27:44))
for (i in 1:num_sm)
  myM_v2[i,] = sample(dt_v2$GA,sum(dt_v2$ptb), replace=T, prob = dt_v2$ptb) %>%
  factor(levels =27:44) %>%
  table()

#posterior distribution for placebo
myM2 = matrix(NA,nrow=num_sm, ncol = length(27:44))
for (i in 1:num_sm)
  myM2[i,] = sample(dt2$GA,sum(dt2$ptb), replace=T, prob = dt2$ptb) %>%
  factor(levels =27:44) %>%
  table()

#proportion risk
prp_p2 <- myM2/sum(dt2$ptb)
prp_v2 <- myM_v2/sum(dt_v2$ptb)

riskM2 <- matrix(NA, nrow=num_sm, ncol = length(27:44))
for(i in 1:num_sm){
  rd_ptb2 = prp_v2[i,] - prp_p2[i,]
  riskM2[i,] = rd_ptb2
}


#### without 5 earliest births in each arm in South Africa ####

#### data prep ####
# vaccine
dt_v3 <- dt_v
dt_v3$ptb[1:7] <-0

# placebo
dt3 <- dt
dt3$ptb[1:7]<-0
dt3$ptb[8]<-4

#posterior distribution for vaccine
myM_v3 = matrix(NA,nrow=num_sm, ncol = length(27:44))
for (i in 1:num_sm)
  myM_v3[i,] = sample(dt_v3$GA,sum(dt_v3$ptb), replace=T, prob = dt_v3$ptb) %>%
  factor(levels =27:44) %>%
  table()

#posterior distribution for placebo
myM3 = matrix(NA,nrow=num_sm, ncol = length(27:44))
for (i in 1:num_sm)
  myM3[i,] = sample(dt3$GA,sum(dt3$ptb), replace=T, prob = dt3$ptb) %>%
  factor(levels =27:44) %>%
  table()

#proportion risk
prp_p3 <- myM3/sum(dt3$ptb)
prp_v3 <- myM_v3/sum(dt_v3$ptb)

riskM3 <- matrix(NA, nrow=num_sm, ncol = length(27:44))
for(i in 1:num_sm){
  rd_ptb3 = prp_v3[i,] - prp_p3[i,]
  riskM3[i,] = rd_ptb3
}



#### without birth at 27 week in vaccine arm in South African trial component ####

#### data prep ####
# vaccine
dt_v5 <- dt_v
dt_v5$ptb[1] <-0 #27 wks

# placebo
dt5 <- dt

#posterior distribution for vaccine

myM_v5 = matrix(NA,nrow=num_sm, ncol = length(27:44))
for (i in 1:num_sm)
  myM_v5[i,] = sample(dt_v5$GA,sum(dt_v5$ptb), replace=T, prob = dt_v5$ptb) %>%
  factor(levels =27:44) %>%
  table()

#posterior distribution for placebo
myM5 = matrix(NA,nrow=num_sm, ncol = length(27:44))
for (i in 1:num_sm)
  myM5[i,] = sample(dt5$GA,sum(dt5$ptb), replace=T, prob = dt5$ptb) %>%
  factor(levels =27:44) %>%
  table()

#proportion risk
prp_p5 <- myM5/sum(dt5$ptb)
prp_v5 <- myM_v5/sum(dt_v5$ptb)


riskM5 <- matrix(NA, nrow=num_sm, ncol = length(27:44))
for(i in 1:num_sm){
  rd_ptb5 = prp_v5[i,] - prp_p5[i,]
  riskM5[i,] = rd_ptb5
}




############fixed delivery risk without resampling ##################
####main analysis#######
##all births in SA
riskdif1 =  (dt_v$ptb/sum(dt_v$ptb) - dt$ptb/sum(dt$ptb))

# ##without 27 week in vaccine arm and 30 week in placebo arm in South Africa

riskdif2 = (dt_v2$ptb/sum(dt_v2$ptb) - dt2$ptb/sum(dt2$ptb))

##without 5 earliest births in each arm in South Africa
riskdif3 = (dt_v3$ptb/sum(dt_v3$ptb) - dt3$ptb/sum(dt3$ptb))



####################### To caluclate excess neonatal deaths, go to risk_calc.R ################


