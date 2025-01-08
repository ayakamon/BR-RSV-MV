### created by Ayaka Monoi



### libraries
#library(readxl)
set.seed(1234)
### data prep
##From Hazel EA, et al. Neonatal mortality risk of vulnerable newborns by fine stratum of gestational age and birth weight for 230 679 live births in nine low- and middle-income countries, 2000-2017. BJOG. 2024 Jan 16. doi: 10.1111/1471-0528.17743. Epub ahead of print. PMID: 38228570.
d <- read.csv(here("data","nmr_hazel.csv")) %>%
  mutate(GA = c(27:44))%>%
  filter(GA >= 37 & GA <= 44)

### reference RR in dataset in Hazel et al
sum(d$Died)/sum(d$Live.births)*1000-> hazel_ref 


### calculate relative risk of mortality for infants born at 37+
#d%>% filter (GA==37)%>% pull(NMR)/hazel_ref -> rr_37

rr_list <- lapply(37:44, function(i) {
  d %>%
    filter(GA == i) %>%
    pull(NMR) / hazel_ref
})

names(rr_list) <- paste0("rr_", 37:44)


### calculate NMR for infants born at 37+ in the South African cohort study using relative risk from the pooled LMIC dataset
### data from South African cohort study
read.csv(here("data","nmr.csv"))%>%filter(GA == "37+") -> d_dchs

### reference NMR in the South African study
death_count <- d_dchs %>%
  filter(Outcome == "death") %>%
  pull(Count)

live_birth_count <- d_dchs %>%
  filter(Outcome == "live birth") %>%
  pull(Count)

dchs_ref <- death_count / live_birth_count

# apply relative mortality risk among 37+ infants from pooled LMIC study to the South African study
nmr_list <- lapply(37:44, function(ga) {
  rr_list[[paste0("rr_", ga)]] * dchs_ref
})

names(nmr_list) <- paste0("nmr_", 37:44)
###estimated neonatal mortality for infants born at 37-44 weeks in South African cohort study (from nmr.R)
model %>% filter(t%in%seq(25*7-4, 44*7-4, by=7))%>% mutate(GA=(t-3)/7+1)%>%kable()

model_updated <- model %>%
  filter(t %in% seq(25 * 7 - 4, 44 * 7 - 4, by = 7)) %>%
  mutate(GA = (t - 3) / 7 + 1) %>%
  select(GA, nmr_t_mid)

# NMR from pooled LMIC study
nmr_df <- data.frame(
  GA = 37:44,
  nmr_t_mid_replacement = unlist(nmr_list)
)

# combine
model_updated <- model_updated %>%
  left_join(nmr_df, by = "GA") %>%
  mutate(nmr_t_mid = ifelse(!is.na(nmr_t_mid_replacement), nmr_t_mid_replacement, nmr_t_mid)) %>%
  select(-nmr_t_mid_replacement)  

nmr_t_mid_vector <- model_updated$nmr_t_mid



### excess deaths using estimated neonatal mortality from South African cohort study together with relative mortality risk among 37+ infants from pooled LMIC study
median(sapply(1:10000, function(i) {
  sum(nmr_t_mid_vector * riskM1[i, ]) * 100000
})) %>%
{cat("Excess deaths using estimated neonatal mortality from South African cohort study together with relative mortality risk among 37+ infants from pooled LMIC study\n"); print(.)}
