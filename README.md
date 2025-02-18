# The benefits and risks of maternal RSV vaccination on mortality in South Africa: a modelling study
Ayaka Monoi, Akira Endo,  Simon R Procter, Sequoia I. Leuba, Stefan Flasche,  Mark Jit, Maternal RSV vaccine benefit-risk advisory group

## Overview of the code files in src folder

## 0) main.R
* main.R is the script to run.

### 1) ve_fit.R
* Codes to esimate waning protection of RSVpreF during infancy using the phase 3 trial data (Mujai et al. RSVVW'24).
* Data used are embedded in the codes.

### 2) RSV_burden.R
*  Codes to estimate country-specific RSV-associated deaths with initial input from Koltai M, et al. Estimating the cost-effectiveness of maternal vaccination and monoclonal antibodies for respiratory syncytial virus in Kenya and South Africa. BMC Med. 2023 Mar 31;21(1):120. doi: 10.1186/s12916-023-02806-w. PMID: 37004062; PMCID: PMC10064962.
*  Data used are Koltai_figure_2B.csv

### 3) nmr.R
* Codes to estimate gestaional age-specific neonatal mortality in South Africa.
* Dammy data are embeded in the codes.

### 4) ptb_posterior.R
* Codes to obtain posterior samples of gestational ages in South African compornent of the phase 3 trial by trial arm.
* Data used are data/GA_birth_by_week_country.csv

### 4) ptb_posterior_vacc.R
* Codes to obtain posterior samples of gestational ages in South African compornent of the phase 3 trial by trial arm using data for infants born to mothers vaccinated at 27-36 weeks.
* Data used are data/GA_birth_by_week_country_vacc_window.csv

## Overview of data

### 1) GA_birth_by_week_country.csv
* Numbers of live births by gestaional age in South African compornent of the phase 3 trial by trial arm.

### 2) GA_birth_by_week_country_vacc_window.csv
* Numbers of live births by gestaional age in South African compornent of the phase 3trial by trial arm.

### 3) Koltai_figure_2B.csv
* RSV hospitalisations in South Africa from Koltai M, et al. Estimating the cost-effectiveness of maternal vaccination and monoclonal antibodies for respiratory syncytial virus in Kenya and South Africa. BMC Med. 2023 Mar 31;21(1):120. doi: 10.1186/s12916-023-02806-w. PMID: 37004062; PMCID: PMC10064962.

### 4) nmr.csv
* Dammy data for neonatal mortality in South Africa.

### 5) nmr_hazel.csv
* Neonatal mortality data from Hazel EA, et al. Neonatal mortality risk of vulnerable newborns by fine stratum of gestational age and birthweight for 230 679 live births in nine low- and middle-income countries, 2000-2017. BJOG. 2024 Jan 16. doi: 10.1111/1471-0528.17743. Epub ahead of print. PMID: 38228570.
