### created by Ayaka Monoi 


### libraries ###
# #install.packages(c("dplyr", "BayesianTools", "ggplot2", "knitr", "tidyverse", "conflicted"))
# library(dplyr)
# library(BayesianTools)
# library(ggplot2)
# library(knitr)
# library(tidyverse)
# library(conflicted)
# conflicts_prefer(dplyr::filter)

set.seed(1234)
### Data ###
### Number of births and deaths at each gestational age

### Dummy data 
tibble(Tmin = seq(26*7,43*7, by=7), #from 27 ga weeks (days)
       Tmax = seq(27*7,44*7, by=7), #until 44 ga weeks (days)
       observed_deaths = c(1,1,1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1,1,1), #number of deaths at each GA
       total_population = c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,10,10,10,10,10,1000))%>% #births
  rowwise() %>%
  mutate(Tmid = (Tmin+Tmax)/2, #days after birth
         nmr = observed_deaths/total_population)-> dat_dch

confidence_intervals <- apply(dat_dch[, c("observed_deaths", "total_population")], #calculate CI using exact binomial test
                              1,
                              function(x) binom.test(x[1], x[2])$conf.int)

dat_dch <- cbind(dat_dch, t(confidence_intervals)) #point estimate and CIs
colnames(dat_dch)[7:8] <- c("CI95_low", "CI95_high")
dat_dch%>%
  mutate(GA = (Tmin/7+1))%>%
  kable(format = "markdown", digits = 4)

sum(dat_dch$total_population)
sum(dat_dch$observed_deaths)

dat_dch


###########Estimate GA-specific NMR following Erlang-2 distribution######
##NMR function for 28-36 GA weeks and constant VE at <28 or 37+
erlang.func <- function(NMR0, T, k = 2, t = min(dat_dch$Tmin):max(dat_dch$Tmax)) {
  re <- numeric(length(t))
  for (i in seq_along(t)) {
    re[i] <- ifelse(t[i] < 36 * 7,
                     sum(sapply(0:(k - 1), function(n) 1 / factorial(n) * exp(-T * (t[i] - 181)) * (T * (t[i] - 181))^n)), #for <37weeks
                     sum(sapply(0:(k - 1), function(n) 1 / factorial(n) * exp(-T * (36 * 7 - 181)) * (T * (36 * 7 - 181))^n))) #for 37+ weeks
  }
  return(NMR0 * re)
}


#define likelihood function
LL <- function(x) {
  a <- x[1] #par1
  b <- x[2] #par2
  dat_dch %>%
    rowwise() %>%
    dplyr::mutate(nmr_est = erlang.func(a, b, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    dplyr::mutate(LL = log(dbinom(observed_deaths,total_population,nmr_est))) -> datatmp
  return(sum(datatmp$LL))
}

#frequentist check
LL_freq = function(x) -LL(x)
out_freq <- optim(c(.5,.01), LL_freq)
#check NMR by plotting
erlang.func(out_freq$par[1], out_freq$par[2], t=min(dat_dch$Tmin):max(dat_dch$Tmax))%>%plot


##Bayesian
bayesianSetup = createBayesianSetup(likelihood = LL, lower = c(.1, .01), upper = c(.5,.1))
iter = 50000
burnin = 1000
settings = list(iterations = iter,  burnin=burnin,message = FALSE, nrChains = 4, thin =4)
out_bay <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
#check
summary(out_bay)
plot(out_bay)


num = (iter - burnin)/settings$nrChains

out_bay[[1]]$chain[seq(1,num),1:2]%>%
  as_tibble() %>%
  rename('NMR0'='1', 'T'='2') -> nmr_bay


#######To calculate the risk, go to risk_cal.R ####################

#Plot estimated NMR with observed NMR
model <-nmr_bay %>%
  rowwise() %>%
  mutate(nmr_t = list(nmr = erlang.func(NMR0 = NMR0, T = T, t=min(dat_dch$Tmin):max(dat_dch$Tmax)))) %>%
  unnest(nmr_t) %>%
  mutate(t = rep(min(dat_dch$Tmin):max(dat_dch$Tmax), times=dim(nmr_bay)[1])) %>%
  group_by(t) %>%
  summarise(nmr_t_mid = median(nmr_t),
            nmr_t_lo = quantile(nmr_t, probs = .025),
            nmr_t_hi = quantile(nmr_t, probs = .975))

#NMR at mean day of each GA weeks
model %>% filter(t%in%seq(25*7-4, 44*7-4, by=7))%>% mutate(GA=(t-3)/7+1)%>%kable()%>% 
  {cat("Estiamted neonatal mortality\n"); print(.)}

#plot (Fig 2A)
ggplot(data = dat_dch, aes(x=Tmid, y=nmr)) +
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = CI95_low, ymax = CI95_high), width = 0.2)+
  geom_ribbon(data = model, aes(x=t, y=nmr_t_mid, ymin=nmr_t_lo, ymax=nmr_t_hi),
              lwd = 1.5, alpha = .2) +
  geom_line(data = model, aes(x=t, y=nmr_t_mid),
            lwd = 1.5, alpha = .8) +
  xlab('Gestational age at births (weeks)') +
  ylab('Neonatal mortality (per 100,000 births)') + 
  scale_y_continuous(labels = function(x) x * 100000)+
  scale_x_continuous(breaks = dat_dch$Tmid, labels=c(27:44))+
  coord_cartesian(ylim=c(0,.9)) +
  theme_minimal()+
 # theme_classic() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold") ,
    panel.grid.major.x = element_line(color = "grey", size = 0.1), 
    panel.grid.minor.x = element_blank() 
  )+
  annotate("text", x = 178, y = 0.8, 
                 label = "(A)", size = 7, hjust = 0)

# ### save as pdf (Fig 2A)###
ggsave(
  filename = "output/fig2A.pdf",
  plot = last_plot(),
  device = "pdf",
  width = 10, height = 7,
  units = "in"
)
