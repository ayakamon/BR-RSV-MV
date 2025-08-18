### created by Ayaka Monoi 
### initial input from Stefan Flasche


#### libraries ####
# library(tidyverse)
# library(epitools)
# library(BayesianTools)
# library(knitr)
# library(ggplot2)


set.seed(1234)
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


#### get corresponding VE with Erlang decay ####
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
LL_freq = function(x) -LL(x)

#### run fitting ####
# frequentist check
out_freq <- optim(c(.5,.5,1/365), LL_freq, lower = c(0.001,0.001,0.001), upper = c(1,1,1), method= 'L-BFGS-B')

# Bayesian
bayesianSetup = createBayesianSetup(likelihood = LL, lower = c(0.001,0.001,0.001), upper = c(.99,.99,.99))
iter = 10000
settings = list(iterations = iter, message = FALSE, nrChains = 2)
out_ba <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
summary(out_ba)
plot(out_ba)


#### get model predictions ####
out_ba[[1]]$chain[seq(1,iter,by=10),1:3] %>%  ##can try until this line
  as_tibble() %>%
  pivot_longer(-'3', names_to = 'group', values_to = 'VE0') %>%
  rename('T'='3') %>%
  mutate(group = case_match(group, '1'~'Severe', '2'~'Less severe')) -> resp

##### To calculate the benefits, go to benef_calc.R #####


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
sev_ve_1 <- model_ve%>%
  filter(t == 1) #1 day

sev_ve_365 <- model_ve%>%
  filter(t == 365) #end of 1st year

sev_ve_1 %>% bind_rows(sev_ve_365)


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


less_ve_365 <- model_ve_less%>%
  filter(t == 365)  #12 month


less_ve_1 <- model_ve_less%>%
  filter(t == 1) #1 day

result_ve <- less_ve_1 %>% bind_rows(less_ve_365) %>% bind_rows(sev_ve_1) %>% bind_rows(sev_ve_365)
write.csv(result_ve, "output/ve_summary.csv", row.names = FALSE)


#### plot Fig 1 ####

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
model_comb <-resp %>%
  rowwise() %>%
  mutate(VE_t = list(VE = erlang.decay(VE0 = VE0, T = T, t=1:365))) %>%
  unnest(VE_t) %>%
  mutate(t = rep(1:365, times=dim(resp)[1])) %>%
  group_by(group,t) %>%
  summarise(VE_t_mid = median(VE_t),
            VE_t_lo = quantile(VE_t, probs = .025),
            VE_t_hi = quantile(VE_t, probs = .975))

# plot 
custom_labeller <- labeller(group = c("Less severe" = "(A) Less severe", 
                                      "Severe" = "(B) Severe"))
rbind(dat_mat_sev,dat_mat_LRTI) %>%
  ggplot(aes(x=Tmid, y=VEmid, ymin=VElo, ymax=VEhi)) +
  geom_pointrange(color=grey(.5)) + 
  geom_ribbon(data = model_comb, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi), 
              lwd = 1.5, alpha = .2) +
  geom_line(data = model_comb, aes(x=t, y=VE_t_mid, ymin=VE_t_lo, ymax=VE_t_hi), 
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




### save as pdf (Fig 1)###
ggsave(
  filename = "output/fig1.pdf",
  plot = last_plot(),
  device = "pdf",
  width = 14, height = 8,
  units = "in"
)





## check convergence ######
out_ba[[1]]

bay1 <- getSample(out_ba, coda = TRUE)[1]

num = iter/settings$thin

trace_ve <- mcmc(bay1[[1]][seq(1,iter, by=10),])
colnames(trace_ve) <- c("VE0_s", "VE0_l", "T_v")

## for the 'xyplot' command
ess_ve <- effectiveSize(trace_ve)
write.csv(ess_ve, file = "output/ess_ve.csv")


xyplot(trace_ve)
densityplot(trace_ve)

plot_ve1 <- xyplot(trace_ve)
plot_ve2 <- densityplot(trace_ve)

p_mc_ve <- grid.arrange(plot_ve1, plot_ve2, ncol = 2)
print(p_mc_ve)

# ### save as pdf ###
ggsave(
  filename = "output/suppl_mcmc_ve.pdf",
  plot = p_mc_ve,
  device = "pdf",
  width = 10, height = 10,
  units = "in"
)
