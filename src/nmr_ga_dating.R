### created by Ayaka Monoi

set.seed(1234)

### Dummy data of births and deaths at each gestational age from the Drakenstein Child Health Study (Zar HJ, et al. Maternal health and birth outcomes in a South African birth cohort study. Hill B, editor. PLOS ONE. 2019 Nov 21;14(11):e0222399. )
tibble(Tmin = seq(26*7,36*7, by=7),
       Tmax = c(seq(27*7,37*7, by=7)), #tentatively consider 37 weeks to
       deaths_obs = c(1,1,1,1,1,1,1,1,1,1,1),
       births_obs = c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 1000))%>%
  rowwise() %>%
  mutate(Tmid = (Tmin+Tmax)/2,
         nmr = deaths_obs/births_obs)-> dat_dch

confidence_intervals <- apply(dat_dch[, c("deaths_obs", "births_obs")],
                              1,
                              function(x) binom.test(x[1], x[2])$conf.int)

dat_dch <- cbind(dat_dch, t(confidence_intervals))
colnames(dat_dch)[7:8] <- c("CI95_low", "CI95_high")
dat_dch

dat_dch %>%
  rowwise %>%
  mutate(GA =Tmin/7+1) ->dat_dch



##########################

##resample 5 times
birth_rsmpl = matrix(0, nrow = 5, ncol = length(25:39))
death_rsmpl = matrix(0, nrow = 5, ncol = length(25:39))
data_list <- list()

for(k in 27:37){
  dat_dch %>%
    filter(GA == k) %>%
    {tibble(
      GA_obs = rep(k, .$births_obs),
      die = c(rep(1, .$deaths_obs), rep(0, .$births_obs - .$deaths_obs)),
      GA = sample(c((k - 2):(k + 2)), .$births_obs, replace = TRUE, rep(1, 5))
    )} -> ga_rsmple

  data_list[[k - 26]] <- ga_rsmple
}
ga_rsmpl<- bind_rows(data_list)



##resample
num_rsmpl =5
birth_rsmpl = matrix(0, nrow = num_rsmpl, ncol = length(25:39))
death_rsmpl = matrix(0, nrow = num_rsmpl, ncol = length(25:39))
data_list <- list()

for(i in 1:num_rsmpl){

  for(k in 27:37){
    dat_dch %>%
      filter(GA == k) %>%
      {tibble(
        GA_obs = rep(k, .$births_obs),
        die = c(rep(1, .$deaths_obs), rep(0, .$births_obs - .$deaths_obs)),
        GA = sample(c((k - 2):(k + 2)), .$births_obs, replace = TRUE, rep(1, 5))
      )} -> ga_rsmple

    data_list[[k - 26]] <- ga_rsmple
  }
  ga_rsmpl<- bind_rows(data_list)

  ga_rsmpl$GA %>% factor(levels = 25:39)%>%table()%>%as.vector() -> birth_rsmpl[i,]
  ga_rsmpl%>%filter(die ==1) %>% .$GA %>%factor(levels = 25:39)%>%table()%>%as.vector() -> death_rsmpl[i,]

}


# number of live births at GA 25:39 after resampling
birth_rsmpl[1,]%>%sum()
# number of deaths born at GA 25:39 after resampling
death_rsmpl%>%sum()




tibble(GA = c(25:39),
       Tmin = seq(24*7,38*7, by=7),
       Tmax = seq(25*7,39*7, by=7), #need to be changed to 37+ constant
       death = death_rsmpl[i,],
       birth = birth_rsmpl[i,])-> dat_rsm

#NMR following Erlang-2 (without constant NMR

i = 5 #change to resample from 1 to 5
tibble(GA = c(25:39),
       Tmin = seq(24*7,38*7, by=7),
       Tmax = seq(25*7,39*7, by=7), #need to be changed to 37+ constant
       death = death_rsmpl[i,],
       birth = birth_rsmpl[i,])-> dat_rsm



#Smooth function for NMR
erlang.func_dat <- function(NMR0, T, k = 2, t = 168:308) {
  re <- numeric(length(t))
  for (i in seq_along(t)) {
    re[i] <- ifelse(t[i] < 36 * 7,
                     sum(sapply(0:(k - 1), function(n) 1 / factorial(n) * exp(-T * (t[i] - (24*7-1))) * (T * (t[i] - (24*7-1)))^n)), #181 = 26*7-1
                     sum(sapply(0:(k - 1), function(n) 1 / factorial(n) * exp(-T * (36 * 7 - (24*7-1))) * (T * (36 * 7 - (24*7-1)))^n)))
  }
  return(NMR0 * re)
}


LL <- function(x) {
  # Extract parameters
  a <- x[1]
  b <- x[2]
  dat_rsm%>%
    rowwise() %>%
    dplyr::mutate(nmr_est = erlang.func_dat(a, b, t=Tmin:Tmax) %>% mean()) %>%
    rowwise() %>%
    dplyr::mutate(LL = log(dbinom(death,birth,nmr_est))) -> datatmp
  return(sum(datatmp$LL))
}


#####MLE estimate###############
LL_freq = function(x) -LL(x)
out <- optim(c(.5,.01), LL_freq)

#plot
erlang.func_dat(out$par[1], out$par[2], t=min(dat_rsm$Tmin):max(dat_rsm$Tmax))%>%plot




###Bayesian ################
bayesianSetup = createBayesianSetup(likelihood = LL, lower = c(.1, .01), upper = c(.7,.1))
iter = 10000
burnin = 1000
settings = list(iterations = iter,  burnin=burnin,message = FALSE, nrChains = 4, thin =4)

out_bay_dat <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
summary(out_bay_dat)
plot(out_bay_dat)


out_bay_dat[[1]]$chain[seq(1,2250),1:2]%>%
  as_tibble() %>%
  rename('NMR0'='1', 'T'='2') -> nmr_bay_dat

##To calculate the risk, go to risk_cal_ci.R



##Make posterior samples of NMR
#from nmr.R after run MCMC
a <-nmr_bay_dat %>%
  rowwise() %>%
  mutate(nmr_t = list(nmr = erlang.func_dat(NMR0 = NMR0, T = T, t=182:308)))

#number of samples
num_smp <- 2000

b <- sample(a$nmr_t,num_smp)

#posterior samples for estimated NMR at each GA week
risk_nmr <- list()
for (i in 1:num_smp) {
  risk_nmr[[i]] <- b[i]$nmr[c(seq(3, 126, by = 7))] #from 27 weeks to 44 weeks
}


calculate_sum <- function(risk_nmr, riskM) {
  results <- numeric(length = num_smp)

  for (i in 1:num_smp) {
    nmr_risk <- as.vector(risk_nmr[[i]])
    ptb_riskdif <- as.vector(riskM[i,])
    result <- sum(nmr_risk *  ptb_riskdif) * 100000 #per 100000 live births

    results[i] <- result
  }
  return(results)
}


#main analysis
calculate_sum(risk_nmr, riskM1)%>%mean() %>%
{cat("Excess deaths incorporating uncertainty of GA dating\n"); print(.)}
calculate_sum(risk_nmr, riskM1)%>%quantile(.025) %>%
  {cat("Excess deaths incorporating uncertainty of GA dating\n"); print(.)}
calculate_sum(risk_nmr, riskM1)%>%quantile(.975) %>%
  {cat("Excess deaths incorporating uncertainty of GA dating\n"); print(.)}

post_risk <- calculate_sum(risk_nmr, riskM1)

###benefit-risk ratio#####

#risk-benefit ratio (excess death per 1 live saved)
rr_post <- post_risk/post_benef
rr_post%>%median()%>%
  {cat("Risk per benefit incorporating uncertainty of GA dating\n"); print(.)}
rr_post%>%quantile(.025)%>%
  {cat("Risk per benefit incorporating uncertainty of GA dating\n"); print(.)}
rr_post%>%quantile(.975)%>%
  {cat("Risk per benefit incorporating uncertainty of GA dating\n"); print(.)}

post_benef%>%median()
post_risk%>%median()



############################
#summary
model <-nmr_bay_dat %>%
  rowwise() %>%
  mutate(nmr_t = list(nmr = erlang.func(NMR0 = NMR0, T = T, t=182:308))) %>%
  unnest(nmr_t) %>%
  mutate(t = rep(182:308, times=dim(nmr_bay_dat)[1])) %>%
  group_by(t) %>%
  summarise(nmr_t_mid = median(nmr_t),
            nmr_t_lo = quantile(nmr_t, probs = .025),
            nmr_t_hi = quantile(nmr_t, probs = .975))


#plot
ggplot(data = dat_dch, aes(x=Tmid, y=nmr)) +
  geom_point()+
  geom_errorbar(aes(ymin = CI95_low, ymax = CI95_high), width = 0.2)+
  geom_ribbon(data = model, aes(x=t, y=nmr_t_mid, ymin=nmr_t_lo, ymax=nmr_t_hi),
              lwd = 1.5, alpha = .2) +
  geom_line(data = model, aes(x=t, y=nmr_t_mid),
            lwd = 1.5, alpha = .5) +
  xlab('Gestational age at births (weeks)') +
  ylab('Neonatal mortality (per 1000 births)') +
  scale_y_continuous(labels = function(x) x * 1000)+
  scale_x_continuous(breaks = dat_dch$Tmid, labels=c(c(27:37)))+
  coord_cartesian(ylim=c(0,.9)) +
  theme_minimal()+
  theme(legend.text = element_text(size = 12))





