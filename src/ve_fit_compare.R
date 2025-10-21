### created by Ayaka Monoi 
### initial input from Stefan Flasche


#### libraries ####
# library(tidyverse)
# library(epitools)
# library(BayesianTools)
# library(knitr)
# library(ggplot2)
 # library(fitR)
 # library(lattice)
 # library(gridExtra)
# library(MCMCpack)


set.seed(1234)

start_time <- Sys.time()
#### Data preparation ####
#### Phase 3 trial of Pfizer's RSVpreF (Munjai et al. 2024/2. RSVVW'24)
#### VE vs RSV-positive severe MA-LRTI
N_placebo = 3563
N_intervention = 3585

tibble(Tmin = c(0,30,60,90,120,150),
       Tmax = c(30,60,90,120,150,180),
       Intervention = diff(c(0,1,4,6,13,18,21)),
       Placebo = diff(c(0,10,28,34,49,61,70))) %>%
  rowwise() %>%
  mutate(Tmid = (Tmin+Tmax)/2,
         VEmid = 1-riskratio.wald(rbind(c(N_placebo-Placebo,Placebo), c(N_intervention-Intervention,Intervention)))$measure[2,1],
         VEhi = 1-riskratio.wald(rbind(c(N_placebo-Placebo,Placebo), c(N_intervention-Intervention,Intervention)))$measure[2,2],
         VElo = 1-riskratio.wald(rbind(c(N_placebo-Placebo,Placebo), c(N_intervention-Intervention,Intervention)))$measure[2,3],
         group = 'Severe') -> dat_mat_sev



##### VE vs RSV-positive less severe MA-LRTI
tibble(Tmin = c(0,30,60,90,120,150),
       Tmax = c(30,60,90,120,150,180),
       Intervention = diff(c(0,2,14,25,40,55,67)),
       Placebo = diff(c(0,15,38,59,88,110,132)) ) %>%
  rowwise() %>%
  mutate(Tmid = (Tmin+Tmax)/2,
         VEmid = 1-riskratio.wald(rbind(c(N_placebo-Placebo,Placebo), c(N_intervention-Intervention,Intervention)))$measure[2,1],
         VEhi = 1-riskratio.wald(rbind(c(N_placebo-Placebo,Placebo), c(N_intervention-Intervention,Intervention)))$measure[2,2],
         VElo = 1-riskratio.wald(rbind(c(N_placebo-Placebo,Placebo), c(N_intervention-Intervention,Intervention)))$measure[2,3],
         group = 'Less severe') -> dat_mat_LRTI


#### get corresponding VE with Erlang-2 decay ####
erlang.decay = function(VE0, T, k=2, t=1:150) {
  res = 0
  for(n in 0:(k-1)){
    res = res + 1/factorial(n)*exp(-T*t)*(T*t)^n }
  return( VE0 * res)
}




#### define LL ####
LL = function(x){
  VE0_s = x[1]
  VE0_l = x[2]
  T = x[3]
  dat_mat_sev %>%
    rowwise() %>%
    mutate(VE_est = erlang.decay(VE0 = VE0_s, T = T, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    mutate(LL = log(dbinom(Intervention,N_intervention,(1-VE_est)*Placebo/N_placebo))) -> datatmp_s
  dat_mat_LRTI %>%
    rowwise() %>%
    mutate(VE_est = erlang.decay(VE0 = VE0_l, T = T, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    mutate(LL = log(dbinom(Intervention,N_intervention,(1-VE_est)*Placebo/N_placebo))) -> datatmp_l
  return(sum(datatmp_s$LL) + sum(datatmp_l$LL))
}

# To calculate WAIC add pont-wise likelihood
LL = function(x, sum = TRUE){
  VE0_s = x[1]
  VE0_l = x[2]
  T = x[3]
  dat_mat_sev %>%
    rowwise() %>%
    mutate(VE_est = erlang.decay(VE0 = VE0_s, T = T, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    mutate(LL = log(dbinom(Intervention,N_intervention,(1-VE_est)*Placebo/N_placebo))) -> datatmp_s
  dat_mat_LRTI %>%
    rowwise() %>%
    mutate(VE_est = erlang.decay(VE0 = VE0_l, T = T, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    mutate(LL = log(dbinom(Intervention,N_intervention,(1-VE_est)*Placebo/N_placebo))) -> datatmp_l
  logLik <- c(datatmp_s$LL, datatmp_l$LL)
  
  if(sum)
    return(sum(logLik))
   else 
    return(logLik)
  
  }


LL_freq = function(x) -LL(x)


# Bayesian
bayesianSetup = createBayesianSetup(likelihood = LL, lower = c(0.001,0.001,0.001), upper = c(.99,.99,.99))
iter = 500000
settings = list(iterations = iter, message = FALSE, nrChains = 2)
out_ba <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)


summary(out_ba) %>% print()
plot(out_ba) %>% print()

waic_result <- WAIC(out_ba)
cat("=== Erlang-2 model ===\n")
print(waic_result)

#### get model predictions ####
out_ba[[1]]$chain[seq(1,iter,by=10),1:3] %>%  ##can try until this line
  as_tibble() %>%
  pivot_longer(-'3', names_to = 'group', values_to = 'VE0') %>%
  rename('T'='3') %>%
  mutate(group = case_match(group, '1'~'Severe', '2'~'Less severe')) -> resp

##### To calculate the benefits, go to benef_calc.R #####




#### get model predictions ####
out_ba[[1]]$chain[seq(1,iter,by=10),1:3] %>%  ##can try until this line
  as_tibble() %>%
  pivot_longer(-'3', names_to = 'group', values_to = 'VE0') %>%
  rename('T'='3') %>%
  mutate(group = case_match(group, '1'~'Severe', '2'~'Less severe')) -> resp


#########

#### summary of predicted VE ####
resp %>%
  group_by(group) %>%
  summarise(VE0mid = median(VE0),
            VE0lo = quantile(VE0, probs=.025),
            VE0hi = quantile(VE0, probs=.975),
            Tmid = 2/median(T),
            Tlo = 2/quantile(T, probs=.975),
            Thi = 2/quantile(T, probs=.025)) %>%
  kable(format = "markdown", digits = 2)

#### VE vs severe RSV disease at day t after birth
model_ve <-resp %>%
  rowwise() %>%
  mutate(VE_t = list(VE = erlang.decay(VE0 = VE0, T = T, t=1:365))) %>%
  unnest(VE_t) %>%
  mutate(t = rep(1:365, times=dim(resp)[1])) %>%
  group_by(group,t)%>%
  summarise(VE_t_mid = median(VE_t),
            VE_t_lo = quantile(VE_t, probs = .025),
            VE_t_hi = quantile(VE_t, probs = .975))%>%
  filter(group == "Severe")

#### VE vs RSV-positive severe MA-LRTI at each age after birth
model_ve%>%
  filter(t == 150) %>% #1 day
  {cat("VE (Erlang-2) vs severe disease at 150 days of life\n"); print(.)}


# model_ve%>%
#   filter(t == 365) %>% #end of 1st year
#   {cat("VE vs severe disease at 1 year of life\n"); print(.)}


#### VE vs less severe

model_ve_less <-resp %>%
  rowwise() %>%
  mutate(VE_t = list(VE = erlang.decay(VE0 = VE0, T = T, t=1:365))) %>%
  unnest(VE_t) %>%
  mutate(t = rep(1:365, times=dim(resp)[1])) %>%
  group_by(group,t)%>%
  summarise(VE_t_mid = median(VE_t),
            VE_t_lo = quantile(VE_t, probs = .025),
            VE_t_hi = quantile(VE_t, probs = .975))%>%
  filter(group == "Less severe")


model_ve_less%>%
  filter(t == 150) %>% #5 month
  {cat("VE vs less severe disease at 150 days of life\n"); print(.)}


# model_ve_less%>%
#   filter(t == 1) %>% #1 day
#   {cat("VE vs less severe disease at 1 day of life\n"); print(.)}

#### plot ####

#### vs RSV-positive severe MA-LRTI 
ggplot(data = dat_mat_sev, aes(x=Tmid, y=VEmid, ymin=VElo, ymax=VEhi)) +
  geom_pointrange(color=grey(.5), size = 1) +
  geom_ribbon(data = model_ve, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi),
              lwd = 1.5, alpha = .2) +
  geom_line(data = model_ve , aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi), lwd = 1.5, alpha = .5) +
  geom_hline(yintercept = 0, size = .5) +  # Adding the line at y = 0% and making it thicker
  geom_hline(yintercept = 1, size = .5) +  # Adding the line at y = 100% and making it thicker
  xlim(0,365) +
  xlab('Time in days') + ylab('Vaccine efficacy') +
  scale_y_continuous(breaks = seq(0, 1, by = .2), labels = scales::percent) +
  coord_cartesian(ylim=c(0,1)) +
  theme_minimal()+
  theme(axis.text.x=element_text(size=36,angle=0),strip.text=element_text(size=24),
        axis.text.y=element_text(size=36),
        axis.title.x = element_text(size = 36),
        axis.title.y = element_text(size = 36),
        legend.text = element_text(size = 36))

#### vs RSV-positive less severe MA-LRTI
ggplot(data = dat_mat_LRTI, aes(x=Tmid, y=VEmid, ymin=VElo, ymax=VEhi)) +
  geom_pointrange(color=grey(.5), size = 1) +
  geom_ribbon(data = model_ve_less, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi),
              lwd = 1.5, alpha = .2) +
  geom_line(data = model_ve_less , aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi), lwd = 1.5, alpha = .5) +
  geom_hline(yintercept = 0, size = .5) +  # Adding the line at y = 0% and making it thicker
  geom_hline(yintercept = 1, size = .5) +  # Adding the line at y = 100% and making it thicker
  xlim(0,365) +
  xlab('Time in days') + ylab('Vaccine efficacy') +
  scale_y_continuous(breaks = seq(0, 1, by = .2), labels = scales::percent) +
  coord_cartesian(ylim=c(0,1)) +
  theme_minimal()+
  theme(axis.text.x=element_text(size=36,angle=0),strip.text=element_text(size=24),
        axis.text.y=element_text(size=36),
        axis.title.x = element_text(size = 36),
        axis.title.y = element_text(size = 36),
        legend.text = element_text(size = 36))


### Plot joint fit
# Combined data
model_comb_erlang <-resp %>%
  rowwise() %>%
  mutate(VE_t = list(VE = erlang.decay(VE0 = VE0, T = T, t=1:365))) %>%
  unnest(VE_t) %>%
  mutate(t = rep(1:365, times=dim(resp)[1])) %>%
  group_by(group,t) %>%
  summarise(VE_t_mid = median(VE_t),
            VE_t_lo = quantile(VE_t, probs = .025),
            VE_t_hi = quantile(VE_t, probs = .975)) %>% 
  mutate(model = "Erlang-2 distribution")

# plot 
custom_labeller <- labeller(group = c("Less severe" = "(A) Less severe", 
                                      "Severe" = "(B) Severe"))
rbind(dat_mat_sev,dat_mat_LRTI) %>%
  ggplot(aes(x=Tmid, y=VEmid, ymin=VElo, ymax=VEhi)) +
  geom_pointrange(color=grey(.5)) + 
  geom_ribbon(data = model_comb_erlang, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi), 
              lwd = 1.5, alpha = .2) +
  geom_line(data = model_comb_erlang, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi), 
            lwd = 1.5, alpha = .5) +
  facet_grid(~group, labeller = custom_labeller) +
  xlim(0, 365) +
  xlab('Time in days') + ylab('Vaccine efficacy') +
  scale_y_continuous(breaks = seq(0, 1, by = .2), labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 0),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    legend.text = element_text(size = 24),
    strip.text = element_text(size = 24),
    axis.line.x = element_blank()
  ) +
  annotate(
    "segment",
    x = 0, xend = 365,
    y = 0, yend = 0,
    arrow = arrow(length = unit(0.5, "cm"), ends = "last", type = "open"),
    color = "black", size = .5
  ) +
  annotate(
    "segment",
    x = 0, xend = 0,
    y = -.1, yend = 1,
    arrow = arrow(length = unit(0.5, "cm"), ends = "last", type = "open"),
    color = "black", size = .5
  )






####################### exponential decay ###########################

#### get corresponding VE with exponential decay ####
exp.decay = function(VE0, T,  t=1:150) {
  res =  exp(-T*t)
  return( VE0 * res)
}


LL = function(x, sum = TRUE){
  VE0_s = x[1]
  VE0_l = x[2]
  T = x[3]
  dat_mat_sev %>%
    rowwise() %>%
    mutate(VE_est = exp.decay(VE0 = VE0_s, T = T, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    mutate(LL = log(dbinom(Intervention,N_intervention,(1-VE_est)*Placebo/N_placebo))) -> datatmp_s
  dat_mat_LRTI %>%
    rowwise() %>%
    mutate(VE_est = exp.decay(VE0 = VE0_l, T = T, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    mutate(LL = log(dbinom(Intervention,N_intervention,(1-VE_est)*Placebo/N_placebo))) -> datatmp_l
  logLik <- c(datatmp_s$LL, datatmp_l$LL)
  
  if(sum)
    return(sum(logLik))
  else 
    return(logLik)
  
}


#### run fitting ####
# frequentist check
out_freq <- optim(c(.5,.5,1/365), LL_freq, lower = c(0.001,0.001,0.001), upper = c(1,1,1), method= 'L-BFGS-B')
# par 1: VE0_s, ...initial value of VE for severe RSV disease,
# par 2: VE0_l,...initial value of VE for less severe RSV disease
# par 3: T... inverse of scale parameter of Weibull distribution

# Bayesian
bayesianSetup = createBayesianSetup(likelihood = LL, lower = c(0.001,0.001,0.001), upper = c(.99,.99,.99))
iter = 500000

settings = list(iterations = iter, message = FALSE, nrChains = 1, thin = 3)
out_ba2 <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out_ba2)  %>% print()
plot(out_ba2)

waic_result_2 <- WAIC(out_ba2)
cat("=== Exponential model ===\n")
print(waic_result_2)

#########

getSample(out_ba2, coda = TRUE)[3] -> bay 
BurnBay <- burnAndThin(bay, burn = 110) #
plot(BurnBay)


### trace plot and density plot for the result without initial XX values ###
## plot all chains together
bay1 <- getSample(out_ba2, coda = TRUE)[1]
bay2 <- getSample(out_ba2, coda = TRUE)[2]
bay3 <- getSample(out_ba2, coda = TRUE)[3]

BurnBay1 <- burnAndThin(bay1, burn = 585)

BurnBay2 <- burnAndThin(bay2, burn = 585)
BurnBay3 <- burnAndThin(bay3, burn = 585)

trace1 <- mcmc(BurnBay1[[1]])
trace2 <- mcmc(BurnBay2[[1]])
trace3 <- mcmc(BurnBay3[[1]])

colnames(trace1) <- c("VE0_s", "VE0_l", "T")
colnames(trace2) <- c("VE0_s", "VE0_l", "T")
colnames(trace3) <- c("VE0_s", "VE0_l", "T")

trace <- mcmc.list(list(trace1, trace2, trace3))

head(trace, 3)

1 - rejectionRate(trace)

## for the 'xyplot' command
xyplot(trace)
densityplot(trace)

plot1 <- xyplot(trace)
plot2 <- densityplot(trace)

p_mc <- grid.arrange(plot1, plot2, ncol = 2)
print(p_mc)


#### get model predictions ####
#bay <- getSample(out_ba1, coda = TRUE)[[1]]
bay <- getSample(out_ba2, coda = TRUE)[[1]]
num = 1000
as_tibble(bay)[(nrow(bay) - (num -1)):nrow(bay),]-> bay
bay %>%
  pivot_longer(-'par 3', names_to = 'group', values_to = 'VE0') %>%
  rename('T'='par 3') %>%
  mutate(group = case_match(group, 'par 1'~'Severe', 'par 2'~'Less severe')) -> resp_ex


#### summary of predicted VE ####
resp_ex %>%
  group_by(group) %>%
  summarise(VE0mid = median(VE0),
            VE0lo = quantile(VE0, probs=.025),
            VE0hi = quantile(VE0, probs=.975),
            Tmid = 2/median(T),
            Tlo = 2/quantile(T, probs=.975),
            Thi = 2/quantile(T, probs=.025)) %>%
  kable(format = "markdown", digits = 2)

#### VE vs severe RSV disease at day t after birth
n_d = 365
model_ve_ex <-resp_ex %>%
  rowwise() %>%
  mutate(VE_t = list(VE = exp.decay(VE0 = VE0, T = T, t=1:n_d))) %>%
  unnest(VE_t) %>%
  mutate(t = rep(1:n_d, times=dim(resp_ex)[1])) %>%
  group_by(group,t)%>%
  summarise(VE_t_mid = median(VE_t),
            VE_t_lo = quantile(VE_t, probs = .025),
            VE_t_hi = quantile(VE_t, probs = .975))%>%
  filter(group == "Severe")

#### VE vs RSV-positive severe MA-LRTI at each age after birth
 model_ve_ex%>%
  filter(t == 150) %>% #5 month
  {cat("VE (exponential) vs severe disease at 150 days of life\n"); print(.)}


#### vs RSV-positive severe MA-LRTI 
ggplot(data = dat_mat_sev, aes(x=Tmid, y=VEmid, ymin=VElo, ymax=VEhi)) +
  geom_pointrange(color=grey(.5), size = 1) +
  geom_ribbon(data = model_ve_ex, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi),
              lwd = 1.5, alpha = .2) +
  geom_line(data = model_ve_ex , aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi), lwd = 1.5, alpha = .5) +
  geom_hline(yintercept = 0, size = .5) +  # Adding the line at y = 0% and making it thicker
  geom_hline(yintercept = 1, size = .5) +  # Adding the line at y = 100% and making it thicker
#  xlim(0,365) +
  xlim(0,n_d) +
  xlab('Time in days') + ylab('Vaccine efficacy') +
  scale_y_continuous(breaks = seq(0, 1, by = .2), labels = scales::percent) +
  coord_cartesian(ylim=c(0,1)) +
  theme_minimal()+
  theme(axis.text.x=element_text(size=36,angle=0),strip.text=element_text(size=24),
        axis.text.y=element_text(size=36),
        axis.title.x = element_text(size = 36),
        axis.title.y = element_text(size = 36),
        legend.text = element_text(size = 36))



### Plot joint fit
# Combined data
model_comb_exp <-resp_ex %>%
  rowwise() %>%
#  mutate(VE_t = list(VE = erlang.decay(VE0 = VE0, T = T, t=1:365))) %>%
  mutate(VE_t = list(VE = exp.decay(VE0 = VE0, T = T, t=1:n_d))) %>%
  unnest(VE_t) %>%
 # mutate(t = rep(1:365, times=dim(resp)[1])) %>%
  mutate(t = rep(1:n_d, times=dim(resp_ex)[1])) %>%
  group_by(group,t) %>%
  summarise(VE_t_mid = median(VE_t),
            VE_t_lo = quantile(VE_t, probs = .025),
            VE_t_hi = quantile(VE_t, probs = .975)) %>% 
  mutate(model = "Exponential distribution")

# plot 
custom_labeller <- labeller(group = c("Less severe" = "(A) Less severe", 
                                      "Severe" = "(B) Severe"))
ggplot() +
  geom_pointrange(data = rbind(dat_mat_sev,dat_mat_LRTI), aes(x=Tmid, y=VEmid, ymin=VElo, ymax=VEhi), color=grey(.5)) +
  geom_ribbon(data = model_comb_exp, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi), 
              lwd = 1.5, alpha = .2) +
  geom_line(data = model_comb_exp, aes(x=t, y=VE_t_mid),
            lwd = 1.5, alpha = .5) +
  facet_grid(~group, labeller = custom_labeller) +
  xlim(0, 365) +
  xlab('Time in days') + ylab('Vaccine efficacy') +
  scale_y_continuous(breaks = seq(0, 1, by = .2), labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 0),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    legend.text = element_text(size = 24),
    strip.text = element_text(size = 24),
    axis.line.x = element_blank()
  ) +
  annotate(
    "segment",
    x = 0, xend = 365,
    y = 0, yend = 0,
    arrow = arrow(length = unit(0.5, "cm"), ends = "last", type = "open"),
    color = "black", size = .5
  ) +
  annotate(
    "segment",
    x = 0, xend = 0,
    y = -.1, yend = 1,
    arrow = arrow(length = unit(0.5, "cm"), ends = "last", type = "open"),
    color = "black", size = .5
  )


### Comparison plot ###
combined_models <- rbind(model_comb_erlang, model_comb_exp)

# plot
custom_labeller <- labeller(group = c("Less severe" = "(A) Less severe",
                                      "Severe" = "(B) Severe"))
ggplot() +
  geom_pointrange(data = rbind(dat_mat_sev,dat_mat_LRTI), aes(x=Tmid, y=VEmid, ymin=VElo, ymax=VEhi), color=grey(.5)) +
  geom_ribbon(data = combined_models, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi, fill = model),
              lwd = 1.5, alpha = .2) +
  geom_line(data = combined_models, aes(x=t, y=VE_t_mid, color = model),
            lwd = 1.5, alpha = .5) +
  facet_grid(~group, labeller = custom_labeller) +
  xlim(0, 365) +
  xlab('Time in days') + ylab('Vaccine efficacy') +
  scale_y_continuous(breaks = seq(0, 1, by = .2), labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(values = c("Erlang-2 distribution" = "#0072B2", "Exponential distribution" = "#D55E00")) +
  scale_fill_manual(values = c("Erlang-2 distribution" = "#0072B2", "Exponential distribution" = "#D55E00")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 0),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    legend.position = "bottom", # Add legend
    legend.title = element_blank(), # No legend title
    legend.text = element_text(size = 24),
    strip.text = element_text(size = 24),
    axis.line.x = element_blank()
  ) +
  annotate(
    "segment",
    x = 0, xend = 365,
    y = 0, yend = 0,
    arrow = arrow(length = unit(0.5, "cm"), ends = "last", type = "open"),
    color = "black", size = .5
  ) +
  annotate(
    "segment",
    x = 0, xend = 0,
    y = -.1, yend = 1,
    arrow = arrow(length = unit(0.5, "cm"), ends = "last", type = "open"),
    color = "black", size = .5
  )

#### Exponential distribution model###

#### Weibull distribution model ###


###########
weibull.decay =  function(VE0, T, k=2, t=1:150) {
    res =  exp(-(T*t)^k)
    return( VE0 * res)
  }


LL = function(x, sum = TRUE){
  VE0_s = x[1]
  VE0_l = x[2]
  T = x[3]
  dat_mat_sev %>%
    rowwise() %>%
    mutate(VE_est = weibull.decay(VE0 = VE0_s, T = T, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    mutate(LL = log(dbinom(Intervention,N_intervention,(1-VE_est)*Placebo/N_placebo))) -> datatmp_s
  dat_mat_LRTI %>%
    rowwise() %>%
    mutate(VE_est = weibull.decay(VE0 = VE0_l, T = T, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    mutate(LL = log(dbinom(Intervention,N_intervention,(1-VE_est)*Placebo/N_placebo))) -> datatmp_l
  logLik <- c(datatmp_s$LL, datatmp_l$LL)
  
  if(sum)
    return(sum(logLik))
  else 
    return(logLik)
  
}


LL_freq = function(x) -LL(x)


#### run fitting ####
# frequentist check
out_freq <- optim(c(.5,.5,1/365), LL_freq, lower = c(0.001,0.001,0.001), upper = c(1,1,1), method= 'L-BFGS-B')
# par 1: VE0_s, ...initial value of VE for severe RSV disease,
# par 2: VE0_l,...initial value of VE for less severe RSV disease
# par 3: T... inverse of scale parameter of Weibull distribution


weibull.decay(out_freq$par[1],  out_freq$par[3], t=1:150)%>%
  plot
 
# Bayesian
bayesianSetup = createBayesianSetup(likelihood = LL, lower = c(0.001,0.001,0.001), upper = c(.99,.99,.99))
iter = 500000

settings = list(iterations = iter, message = FALSE, nrChains = 1, thin = 3)
out_ba3 <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
summary(out_ba3) %>% print()
plot(out_ba3)

waic_result_3 <- WAIC(out_ba3)
cat("=== Weibull model ===\n")
print(waic_result_3)

#getSample(out_ba1, coda = TRUE)[2] -> bay #change 1 to 3
getSample(out_ba3, coda = TRUE)[2] -> bay #change 1 to 3
BurnBay <- burnAndThin(bay, burn = 110) 
plot(BurnBay)


### trace plot and density plot for the result without initial XX values ###
## plot all chains together
bay1 <- getSample(out_ba3, coda = TRUE)[1]
bay2 <- getSample(out_ba3, coda = TRUE)[2]
bay3 <- getSample(out_ba3, coda = TRUE)[3]

BurnBay1 <- burnAndThin(bay1, burn = 585)

BurnBay2 <- burnAndThin(bay2, burn = 585)
BurnBay3 <- burnAndThin(bay3, burn = 585)

trace1 <- mcmc(BurnBay1[[1]])
trace2 <- mcmc(BurnBay2[[1]])
trace3 <- mcmc(BurnBay3[[1]])

colnames(trace1) <- c("VE0_s", "VE0_l", "T")
colnames(trace2) <- c("VE0_s", "VE0_l", "T")
colnames(trace3) <- c("VE0_s", "VE0_l", "T")

trace <- mcmc.list(list(trace1, trace2, trace3))

head(trace, 3)

1 - rejectionRate(trace)

## for the 'xyplot' command
xyplot(trace)
densityplot(trace)

plot1 <- xyplot(trace)
plot2 <- densityplot(trace)

p_mc <- grid.arrange(plot1, plot2, ncol = 2)
print(p_mc)


#### get model predictions ####
bay <- getSample(out_ba3, coda = TRUE)[[1]]
num = 1000
as_tibble(bay)[(nrow(bay) - (num -1)):nrow(bay),]-> bay
bay %>%
  pivot_longer(-'par 3', names_to = 'group', values_to = 'VE0') %>%
  rename('T'='par 3') %>%
  mutate(group = case_match(group, 'par 1'~'Severe', 'par 2'~'Less severe')) -> resp_wei


#### summary of predicted VE ####
resp_wei %>%
  group_by(group) %>%
  summarise(VE0mid = median(VE0),
            VE0lo = quantile(VE0, probs=.025),
            VE0hi = quantile(VE0, probs=.975),
            Tmid = 2/median(T),
            Tlo = 2/quantile(T, probs=.975),
            Thi = 2/quantile(T, probs=.025)) %>%
  kable(format = "markdown", digits = 2)

#### VE vs severe RSV disease at day t after birth
n_d = 365
model_ve_wei <-resp_wei %>%
  rowwise() %>%
  #  mutate(VE_t = list(VE = erlang.decay(VE0 = VE0, T = T, t=1:365))) %>%
  mutate(VE_t = list(VE = weibull.decay(VE0 = VE0, T = T, t=1:n_d))) %>%
  unnest(VE_t) %>%
  #  mutate(t = rep(1:365, times=dim(resp)[1])) %>%
  mutate(t = rep(1:n_d, times=dim(resp_wei)[1])) %>%
  group_by(group,t)%>%
  summarise(VE_t_mid = median(VE_t),
            VE_t_lo = quantile(VE_t, probs = .025),
            VE_t_hi = quantile(VE_t, probs = .975))%>%
  filter(group == "Severe")

model_ve_less_wei <-resp_wei %>%
  rowwise() %>%
  #  mutate(VE_t = list(VE = erlang.decay(VE0 = VE0, T = T, t=1:365))) %>%
  mutate(VE_t = list(VE = weibull.decay(VE0 = VE0, T = T, t=1:n_d))) %>%
  unnest(VE_t) %>%
  #  mutate(t = rep(1:365, times=dim(resp)[1])) %>%
  mutate(t = rep(1:n_d, times=dim(resp_wei)[1])) %>%
  group_by(group,t)%>%
  summarise(VE_t_mid = median(VE_t),
            VE_t_lo = quantile(VE_t, probs = .025),
            VE_t_hi = quantile(VE_t, probs = .975))%>%
  filter(group == "Less severe")


model_ve_wei%>%
  filter(group == "Severe") %>% 
  filter(t == 150) %>% #1 day
  {cat("VE (Weibull) vs severe disease at 150 days of life\n"); print(.)}


#### vs RSV-positive severe MA-LRTI 
ggplot(data = dat_mat_sev, aes(x=Tmid, y=VEmid, ymin=VElo, ymax=VEhi)) +
  geom_pointrange(color=grey(.5), size = 1) +
  geom_ribbon(data = model_ve_wei, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi),
              lwd = 1.5, alpha = .2) +
  geom_line(data = model_ve_wei , aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi), lwd = 1.5, alpha = .5) +
  geom_hline(yintercept = 0, size = .5) +  # Adding the line at y = 0% and making it thicker
  geom_hline(yintercept = 1, size = .5) +  # Adding the line at y = 100% and making it thicker
  xlim(0,365) +
  xlab('Time in days') + ylab('Vaccine efficacy') +
  scale_y_continuous(breaks = seq(0, 1, by = .2), labels = scales::percent) +
  coord_cartesian(ylim=c(0,1)) +
  theme_minimal()+
  theme(axis.text.x=element_text(size=36,angle=0),strip.text=element_text(size=24),
        axis.text.y=element_text(size=36),
        axis.title.x = element_text(size = 36),
        axis.title.y = element_text(size = 36),
        legend.text = element_text(size = 36))

#### vs RSV-positive less severe MA-LRTI
ggplot(data = dat_mat_LRTI, aes(x=Tmid, y=VEmid, ymin=VElo, ymax=VEhi)) +
  geom_pointrange(color=grey(.5), size = 1) +
  geom_ribbon(data = model_ve_less_wei, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi),
              lwd = 1.5, alpha = .2) +
  geom_line(data = model_ve_less_wei, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi), lwd = 1.5, alpha = .5) +
  geom_hline(yintercept = 0, size = .5) +  # Adding the line at y = 0% and making it thicker
  geom_hline(yintercept = 1, size = .5) +  # Adding the line at y = 100% and making it thicker
  xlim(0,365) +
  xlab('Time in days') + ylab('Vaccine efficacy') +
  scale_y_continuous(breaks = seq(0, 1, by = .2), labels = scales::percent) +
  coord_cartesian(ylim=c(0,1)) +
  theme_minimal()+
  theme(axis.text.x=element_text(size=36,angle=0),strip.text=element_text(size=24),
        axis.text.y=element_text(size=36),
        axis.title.x = element_text(size = 36),
        axis.title.y = element_text(size = 36),
        legend.text = element_text(size = 36))


model_comb_wei <-resp_wei %>%
  rowwise() %>%
  mutate(VE_t = list(VE = weibull.decay(VE0 = VE0, T = T, t=1:365))) %>%
  unnest(VE_t) %>%
  mutate(t = rep(1:365, times=dim(resp_wei)[1])) %>%
  group_by(group,t) %>%
  summarise(VE_t_mid = median(VE_t),
            VE_t_lo = quantile(VE_t, probs = .025),
            VE_t_hi = quantile(VE_t, probs = .975)) %>% 
  mutate(model = "Weibull distribution")


##########################

############## plot every distiburion together
### Comparison plot ###
combined_models <- rbind(model_comb_erlang, model_comb_exp, model_comb_wei)

# plot
custom_labeller <- labeller(group = c("Less severe" = "(A) Less severe",
                                      "Severe" = "(B) Severe"))
p_combined <- ggplot() +
  geom_pointrange(data = rbind(dat_mat_sev,dat_mat_LRTI), aes(x=Tmid, y=VEmid, ymin=VElo, ymax=VEhi), color=grey(.5)) +
  geom_ribbon(data = combined_models, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi, fill = model),
              lwd = 1.5, alpha = .2) +
  geom_line(data = combined_models, aes(x=t, y=VE_t_mid, color = model),
            lwd = 1.5, alpha = .5) +
  facet_grid(~group, labeller = custom_labeller) +
  xlim(0, 365) +
  xlab('Time in days') + ylab('Vaccine efficacy') +
  scale_y_continuous(breaks = seq(0, 1, by = .2), labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(values = c("Erlang-2 distribution" = "#0072B2", "Exponential distribution" = "#D55E00","Weibull distribution" = "green" )) +
  scale_fill_manual(values = c("Erlang-2 distribution" = "#0072B2", "Exponential distribution" = "#D55E00","Weibull distribution" = "green" )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, angle = 0),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    legend.position = "bottom", # Add legend
    legend.title = element_blank(), # No legend title
    legend.text = element_text(size = 24),
    strip.text = element_text(size = 24),
    axis.line.x = element_blank()
  ) +
  annotate(
    "segment",
    x = 0, xend = 365,
    y = 0, yend = 0,
    arrow = arrow(length = unit(0.5, "cm"), ends = "last", type = "open"),
    color = "black", size = .5
  ) +
  annotate(
    "segment",
    x = 0, xend = 0,
    y = -.1, yend = 1,
    arrow = arrow(length = unit(0.5, "cm"), ends = "last", type = "open"),
    color = "black", size = .5
  )

print(p_combined)

### save as pdf ###
ggsave(
  filename = "output/suppl_ve.pdf",
  plot = last_plot(),
  device = "pdf",
  width = 14, height = 8,
  units = "in"
)

end_time <- Sys.time()

durtion <- end_time - start_time
print(durtion)
