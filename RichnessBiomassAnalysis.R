# loading in packages
library(ggplot2)
library(MASS)
library(reshape2)
library(tidyverse)
library(lme4)
library(glmmTMB)
library(AICcmodavg)
library(emmeans)
library(dplyr)
library(AER)

# reading in the data
emergebiomass <- read.csv ("EmergeBiomassTotaled.csv")

#  totaling biomass and density per sample
biomasspersite <- emergebiomass %>%
  group_by(Wetland, DateCollected, Type, HLength, CoV, Sample) %>% 
  summarize(Biomass = sum (tot_biomass, na.rm = TRUE), 
            Density = sum (Density, na.rm = TRUE),
            .groups = "drop") 

# creating flux values 
biomasspersite <- biomasspersite %>%
  mutate(DateCollected = as.Date(DateCollected, format = "%m/%d/%Y")) %>% # Ensure Date format
  group_by(DateCollected) %>%
  mutate(Flux = case_when(
    DateCollected == as.Date("2023-02-21") ~ Biomass / 7,
    DateCollected == as.Date("2024-02-21") ~ Biomass / 7,
    DateCollected == as.Date("2023-03-06") ~ Biomass / 13,
    DateCollected == as.Date("2023-03-05") ~ Biomass / 13,
    DateCollected == as.Date("2023-03-21") ~ Biomass / 15,
    DateCollected == as.Date("2023-04-05") ~ Biomass / 15,
    DateCollected == as.Date("2023-04-20") ~ Biomass / 15,
    DateCollected == as.Date("2023-05-03") ~ Biomass / 13,
    DateCollected == as.Date("2024-05-03") ~ Biomass / 13,
    DateCollected == as.Date("2023-05-12") ~ Biomass / 9,
    DateCollected == as.Date("2023-05-17") ~ Biomass / 14,
    DateCollected == as.Date("2023-05-31") ~ Biomass / 14,
    DateCollected == as.Date("2023-06-14") ~ Biomass / 14,
    DateCollected == as.Date("2023-06-15") ~ Biomass / 15,
    DateCollected == as.Date("2023-06-23") ~ Biomass / 23,
    DateCollected == as.Date("2023-06-22") ~ Biomass / 22,
    DateCollected == as.Date("2023-07-05") ~ Biomass / 16,
    DateCollected == as.Date("2023-07-06") ~ Biomass / 16,
    DateCollected == as.Date("2023-07-19") ~ Biomass / 14,
    DateCollected == as.Date("2024-07-19") ~ Biomass / 14,
    DateCollected == as.Date("2023-08-01") ~ Biomass / 13,
    TRUE ~ Biomass
  ))

# renaming for simpler names
across <- biomasspersite %>%
  rename(Hydrology= HLength)%>%
  rename(Vegetation= Type)

# testing for normality, looking at the relationship between biomass flux and hydrology
# and vegetation
qqnorm(across$Flux)
qqline(across$Flux)

shapiro.test(across$Flux)

# testing for overdispersion
dispmodel <- glm(Flux~Hydrology, data = across, family = poisson)
dispersiontest(dispmodel)

dispmodel <- glm(Flux~Vegetation, data = across, family = poisson)
dispersiontest(dispmodel)



## now averaging flux and emergence across the replicates to get one value per site per
# collection date, fixing dates
replicates<- across %>%
  group_by(Wetland, DateCollected, Vegetation, Hydrology) %>% 
  summarize(Flux = mean(Flux, na.rm = TRUE), 
            Emergence = mean(Emergence, na.rm = TRUE),
            .groups = "drop") %>%
  rename(Hydroperiod= Hydrology)%>%
  mutate(DateCollected = as.Date(DateCollected, format = "%m/%d/%Y")) %>%
  mutate(DateCollected = case_when(
    DateCollected == as.Date("2023-02-21") ~ "Feb",
    DateCollected == as.Date("2024-02-21") ~ "Feb",
    DateCollected == as.Date("2023-03-06") ~ "Mar",
    DateCollected == as.Date("2023-03-05") ~ "Mar",
    DateCollected == as.Date("2023-03-21") ~ "Mar",
    DateCollected == as.Date("2023-04-05") ~ "Apr",
    DateCollected == as.Date("2023-04-20") ~ "Apr",
    DateCollected == as.Date("2023-05-03") ~ "May",
    DateCollected == as.Date("2024-05-03") ~ "May",
    DateCollected == as.Date("2023-05-12") ~ "May",
    DateCollected == as.Date("2023-05-17") ~ "May",
    DateCollected == as.Date("2023-05-31") ~ "May",
    DateCollected == as.Date("2023-06-14") ~ "Jun",
    DateCollected == as.Date("2023-06-15") ~ "Jun",
    DateCollected == as.Date("2023-06-23") ~ "Jun",
    DateCollected == as.Date("2023-06-22") ~ "Jun",
    DateCollected == as.Date("2023-07-05") ~ "Jul",
    DateCollected == as.Date("2023-07-06") ~ "Jul",
    DateCollected == as.Date("2023-07-19") ~ "Jul",
    DateCollected == as.Date("2024-07-19") ~ "Jul",
    DateCollected == as.Date("2023-08-01") ~ "Aug",
    TRUE ~ as.character(DateCollected)
  )) %>%
  relocate(DateCollected, .before = 1)

### anova of RAW data, replicates/across time for vegetation and hydrology as predictors
# vegetation as a predictor of flux across time
fluxveg<- aov(Flux ~ Vegetation, data = replicates)
summary(fluxveg)
# p = 0.869
emmeans(fluxveg, pairwise~Vegetation, type="response")
# Vegetation emmean   SE df lower.CL upper.CL
# Marsh        2.19 0.336 134     1.53     2.85
# Swamp        2.11 0.367 134     1.38     2.83
# p = 0.8686

# hydrology as a predictor of flux across time
fluxhydro<- aov(Flux ~ Hydroperiod, data = replicates)
summary(fluxhydro)
# p = 0.0052
emmeans(fluxhydro, pairwise~Hydroperiod, type="response")
# Hydroperiod  emmean   SE df lower.CL upper.CL
# Intermediate   1.38 0.406 133    0.577     2.18
# Long           3.14 0.386 133    2.377     3.91
# Short          1.75 0.458 133    0.842     2.65
# intermediate/long: p = 0.0059
# intermediate/short: p = 0.8211
# long/short: p = 0.0555


# looking at date collected as a predictor of flux across time
fluxdate<- aov(Flux ~ DateCollected, data = replicates)
summary(fluxdate)
# p = 0.454
emmeans(fluxdate, pairwise~DateCollected, type="response")
# DateCollected emmean   SE df lower.CL upper.CL
# Apr             2.34 0.614 129   1.1211     3.55
# Aug             2.57 0.868 129   0.8559     4.29
# Feb             1.35 0.831 129  -0.2953     2.99
# Jul             1.88 0.588 129   0.7140     3.04
# Jun             2.82 0.831 129   1.1741     4.46
# Mar             1.27 0.600 129   0.0864     2.46
# May             2.78 0.509 129   1.7706     3.78
# no significance between dates  
# flux highest in june and lowest in march


# looking at total flux, summing replicates to get total flux per site
totals<- replicates %>%
  group_by(Wetland, Vegetation, Hydroperiod) %>% 
  summarize(Biomass = sum(Flux, na.rm = TRUE), 
            .groups = "drop") 

# looking at hydroperiod as a predictor of total flux
fluxhydroT<- aov(Biomass ~ Hydroperiod, data = totals)
summary(fluxhydroT)
emmeans(fluxhydroT, pairwise~Hydroperiod, type="response")
# Hydroperiod  emmean   SE df lower.CL upper.CL
# Intermediate   16.2 6.11  8     2.13     30.3
# Long           40.8 6.11  8    26.74     54.9
# Short          21.5 7.06  8     5.28     37.8

# looking at vegetation as a predictor of total flux
fluxvegT<- aov(Biomass ~ Vegetation, data = totals)
summary(fluxvegT)
emmeans(fluxvegT, pairwise~Vegetation, type="response")
# Vegetation emmean   SE df lower.CL upper.CL
# Marsh        27.0 6.82  9    11.60     42.5
# Swamp        26.1 7.47  9     9.25     43.0

# creating a null model for comparison of other models
null_model_fr <- glm(Flux ~ 1, family = Gamma (link = "log"), data = replicates)
summary (null_model_fr)

## HYDROLOGY
# now using a Gamma distribution to correct for skewness in the data

# replicates: sampling date data, avg values from replicates so one value per
# sampling date per wetland

# looking at the relationship between hydroperiod length and flux per collection per site
hydrology_fr <- glm (Flux ~ Hydroperiod, family = Gamma (link = "log"), data = replicates)
summary (hydrology_fr)

# doing the anova comparing our model to the null model
anova(null_model_fr, hydrology_fr, test = "Chisq")
# the significant p-value (0.00511) indicates that including hydroperiod significantly improves model fit
emmeans(hydrology_fr, pairwise~Hydroperiod, type="response")
# predicted average daily flux for each hydroperiod category:
# Intermediate: 1.38 (95% CI:  0.947 to 2.01)
# Long: 3.14 (95% CI: 2.196 to 4.49)
# Short: 1.75 (95% CI: 1.143 to 2.67)
# Intermediate vs Long: Significant difference (p =  0.0061); flux is lower for intermediate hydroperiods than for long ones.
# Intermediate vs Short: No significant difference (p = 0.6912).
# Long vs Short: No significant difference (p = 0.0961); flux is higher for long hydroperiods than for short ones.


# looking at date collected as well
date_fr <- glm (Flux ~ DateCollected, family = Gamma (link = "log"), data = replicates)
summary (date_fr)
# anova comparing our model to the null
anova(null_model_fr, date_fr, test = "Chisq")



## VEGETATION
# vegetation and average flux per sampling date
# looking at the relationship between vegetation and richness
vegetation_fr <- glm (Flux ~ Vegetation, family = Gamma (link = "log"), data = replicates)
summary (vegetation_fr)

# anova comparing our model to the null
anova(null_model_fr, vegetation_fr, test = "Chisq")
# the insignificant p-value (0.8677) indicates that including vegetation does not improve model fit.
emmeans(vegetation_fr, pairwise~Vegetation, type="response")
# Predicted average daily flux for each vegetation category:
# Marsh: 2.19 (95% CI: 1.61 to 2.98)
# Swamp: 2.11 (95% CI: 1.51 to 2.95)
# Marsh vs swamp: insignificant difference (p = 0.8679)



## Now looking at total flux 
# creating a null model for comparison of other models
null_model_bt <- glm(Biomass ~ 1, family = Gamma (link = "log"), data = totals)
summary (null_model_bt)

# Hydrology and total biomass per site
hydrology_bt <- glm (Biomass ~ Hydroperiod, family = Gamma (link = "log"), data = totals)
summary (hydrology_bt)

# anova
anova(null_model_bt, hydrology_bt, test = "Chisq")
# the significant p-value (0.008107) indicates that including hydroperiod significantly improves model fit.
emmeans(hydrology_bt, pairwise~Hydroperiod, type="response")


# vegetation and total biomass per sampling dates
vegetation_bt <- glm (Biomass ~ Vegetation, family = Gamma (link = "log"), data = totals)
summary (vegetation_bt)

# anova
anova(null_model_bt, vegetation_bt, test = "Chisq")
# insignificant p-value (0.9325) 
emmeans(vegetation_bt, pairwise~Vegetation, type="response")

### Richness analysis
library(ggplot2)
library(MASS)
library(reshape2)
library(tidyverse)
library(lme4)
library(glmmTMB)
library(AICcmodavg)
library(emmeans)
library(dplyr)
library(AER)

across<- read.csv("AnalysisDataPrepped.csv")
across <- across %>%
  rename(Hydrology= HLength)%>%
  rename(Vegetation= Type)

# testing for normality
qqnorm(across$Richness)
qqline(across$Richness)

shapiro.test(across$Richness)

# measuring dispersion in the data
dispmodel <- glm(Richness~Hydrology, data = across, family = poisson)
dispersiontest(dispmodel)

dispmodel <- glm(Richness~Vegetation, data = across, family = poisson)
dispersiontest(dispmodel)

## analysis of replicates (replicates averaged for each sampling date)
replicates<- across %>%
  group_by(Wetland, DateCollected, Vegetation, Hydrology, CoV) %>% 
  summarize(Biomass = mean(Biomass, na.rm = TRUE), 
            Density = mean(Density, na.rm = TRUE),
            Richness = mean(Richness, na.rm = TRUE),
            .groups = "drop") %>%
  rename(Hydroperiod= Hydrology) %>%
  mutate(DateCollected = as.Date(DateCollected, format = "%m/%d/%Y")) %>%
  mutate(DateCollected = case_when(
    DateCollected == as.Date("2023-02-21") ~ "Feb",
    DateCollected == as.Date("2024-02-21") ~ "Feb",
    DateCollected == as.Date("2023-03-06") ~ "Mar",
    DateCollected == as.Date("2023-03-05") ~ "Mar",
    DateCollected == as.Date("2023-03-21") ~ "Mar",
    DateCollected == as.Date("2023-04-05") ~ "Apr",
    DateCollected == as.Date("2023-04-20") ~ "Apr",
    DateCollected == as.Date("2023-05-03") ~ "May",
    DateCollected == as.Date("2024-05-03") ~ "May",
    DateCollected == as.Date("2023-05-12") ~ "May",
    DateCollected == as.Date("2023-05-17") ~ "May",
    DateCollected == as.Date("2023-05-31") ~ "May",
    DateCollected == as.Date("2023-06-14") ~ "Jun",
    DateCollected == as.Date("2023-06-15") ~ "Jun",
    DateCollected == as.Date("2023-06-23") ~ "Jun",
    DateCollected == as.Date("2023-06-22") ~ "Jun",
    DateCollected == as.Date("2023-07-05") ~ "Jul",
    DateCollected == as.Date("2023-07-06") ~ "Jul",
    DateCollected == as.Date("2023-07-19") ~ "Jul",
    DateCollected == as.Date("2024-07-19") ~ "Jul",
    DateCollected == as.Date("2023-08-01") ~ "Aug",
    TRUE ~ as.character(DateCollected)
  )) %>%
  relocate(DateCollected, .before = 1)

## anova of raw data, replicates/across time
richnessveg<- aov(Richness ~ Vegetation, data = replicates)
summary(richnessveg)
# p = 0.392
emmeans(richnessveg, pairwise~Vegetation, type="response")
# Vegetation emmean   SE df lower.CL upper.CL
# Marsh        3.04 0.223 134     2.60     3.48
# Swamp        3.32 0.243 134     2.84     3.81
# p = 0.3919

richnesshydro<- aov(Richness ~ Hydroperiod, data = replicates)
summary(richnesshydro)
# p = 0.000323
emmeans(richnesshydro, pairwise~Hydroperiod, type="response")
# Hydroperiod  emmean   SE df lower.CL upper.CL
# Intermediate   2.64 0.265 133     2.11     3.16
# Long           3.99 0.252 133     3.49     4.49
# Short          2.70 0.299 133     2.11     3.29
# intermediate/long: p = 0.0009
# intermediate/short: p = 0.9869
# long/short: p = 0.0035

replicatesS <- replicates %>%
  filter(Hydroperiod == "Short")
replicatesI <- replicates %>%
  filter(Hydroperiod == "Intermediate")
replicatesL <- replicates %>%
  filter(Hydroperiod == "Long")


richnessdate<- aov(Richness ~ DateCollected, data = replicates)
summary(richnessdate)
# p = 0.145
emmeans(richnessdate, pairwise~DateCollected, type="response")
# DateCollected emmean   SE df lower.CL upper.CL
# Apr             3.64 0.403 129    2.844     4.44
# Aug             3.47 0.569 129    2.343     4.60
# Feb             1.93 0.545 129    0.852     3.01
# Jul             3.21 0.385 129    2.446     3.97
# Jun             3.88 0.545 129    2.796     4.95
# Mar             2.73 0.394 129    1.953     3.51
# May             3.23 0.334 129    2.569     3.89
# no significance between dates 

richnessvegT<- aov(Richness ~ Vegetation, data = totals)
summary(richnessvegT)
emmeans(richnessvegT, pairwise~Vegetation, type="response")

richnesshydroT<- aov(Richness ~ Hydroperiod, data = totals)
summary(richnesshydroT)
emmeans(richnesshydroT, pairwise~Hydroperiod, type="response")



# creating a null model for comparison of other models
null_model_rr <- glm.nb(Richness ~ 1, data = replicates, link = log)
summary (null_model_rr)


## HYDROLOGY

# RICHNESS: sampling date data, avg values from replicates so one value per
# sampling date per wetland

# looking at the relationship between hydroperiod length and richness
hydrology_rr <- glm.nb(Richness ~ Hydroperiod, data = replicates, link = log)
summary (hydrology_rr)

# anova
anova(null_model_rr, hydrology_rr, test = "Chisq")
# significant p-value (0.000308779) indicates that including hydroperiod significantly improves model fit.
emmeans(hydrology_rr, pairwise~Hydroperiod, type="response")
# Predicted richness for each hydroperiod category:
# Intermediate: 2.64 (95% CI: 2.21 to 3.15)
# Long: 3.99 (95% CI: 3.48 to 4.57)
# Short: 2.70 (95% CI: 2.21 to 3.29)
# Intermediate vs Long: Significant difference (p = 0.0009); richness is lower for intermediate hydroperiods than for long ones.
# Intermediate vs Short: No significant difference (p = 0.984).
# Long vs Short: Significant difference (p = 0.0042); richness is higher for long hydroperiods than for short ones.

date_rr <- glm.nb(Richness ~ DateCollected, data = replicates, link = log)
summary (date_rr)

# anova
anova(null_model_rr, date_rr, test = "Chisq")


## analysis of totals
totals<- replicates %>%
  group_by(Wetland, Vegetation, Hydroperiod, CoV) %>% 
  summarize(Biomass = sum(Biomass, na.rm = TRUE), 
            Density = sum(Density, na.rm = TRUE),
            .groups = "drop") %>%
  # bringing in richness values from the richness calculations in TotalRichnessBiomass.R
  mutate(Richness = case_when(
    Wetland == "W15" ~ "37",
    Wetland == "W21" ~ "32",
    Wetland == "W53" ~ "29",
    Wetland == "W46" ~ "31",
    Wetland == "W37" ~ "22",
    Wetland == "W42" ~ "31",
    Wetland == "W11" ~ "35",
    Wetland == "W68" ~ "23",
    Wetland == "W58" ~ "34",
    Wetland == "W32" ~ "32",
    Wetland == "W52" ~ "11",
  )) 

# creating a null model for comparison of other models
null_model_rt <- glm.nb(Richness ~ 1, data = totals, link = log)
summary (null_model_rt)

# Hydrology and total richness per site
# looking at the relationship between hydroperiod length and richness
hydrology_rt <- glm.nb(Richness ~ Hydroperiod, data = totals, link = log)
summary (hydrology_rt)

# anova
anova(null_model_rt, hydrology_rt, test = "Chisq")
# insignificant p-value (0.1611078) indicates that including hydroperiod does not significantly improve model fit.
emmeans(hydrology_rt, pairwise~Hydroperiod, type="response")
# Predicted richness for each hydroperiod category:
# Intermediate: 23.8 (95% CI: 18.9 to 29.8)
# Long: 32.8 (95% CI: 26.8 to 40.1)
# Short: 30.3 (95% CI: 23.9 to 38.6)
# Intermediate vs Long: insignificant difference (p = 0.0964)
# Intermediate vs Short: No significant difference (p = 0.3153).
# Long vs Short: Significant difference (p = 0.8811)


## VEGETATION

# vegetation and average richness per sampling date
# looking at the relationship between vegetation and richness
vegetation_rr <- glm.nb(Richness ~ Vegetation, data = replicates, link = log)
summary (vegetation_rr)

# anova
anova(null_model_rr, vegetation_rr, test = "Chisq")
# The insignificant p-value (0.38444) indicates that including vegetation does not improve model fit.
emmeans(vegetation_rr, pairwise~Vegetation, type="response")
# Predicted richness for each vegetation category:
# Marsh: 3.04 (95% CI: 2.65 to 3.49)
# Swamp: 3.32 (95% CI: 2.87 to 3.84)
# Marsh vs swamp: insignificant difference (p = 0.3835)

## analysis of totals

# vegetation and average richness per sampling date
# looking at the relationship between vegetation and richness
vegetation_rt <- glm.nb(Richness ~ Vegetation, data = totals, link = log)
summary (vegetation_rt)

# anova
anova(null_model_rt, vegetation_rt, test = "Chisq")
# insignificant p-value (0.4663191) indicates that including hydroperiod significantly improves model fit.
emmeans(vegetation_rt, pairwise~Vegetation, type="response")
# Predicted richness for each vegetation category:
# Marsh: 30.3 (95% CI: 24.7 to 37.2)
# Swamp: 27.0 (95% CI: 21.4 to 34.0)
# Marsh vs swamp: insignificant difference (p = 0.4599)

# vegetation is NOT a significant predictor of total richness OR richness across time.
# HOWEVER we know from nmds that COMPOSITION differs


### Benthic Analysis
benthicrichness <- read.csv ("BenthicBiomassTotaled.csv")

across <- benthicrichness %>%
  group_by(Wetland, Date_collected, Sample, HLength, Type, CoV) %>%
  summarise(Richness = n_distinct(Taxa)) %>%
  ungroup()

across <- across %>%
  rename(Hydrology= HLength)%>%
  rename(Vegetation= Type)%>%
  rename(DateCollected = Date_collected)

# testing for normality 
qqnorm(across$Richness)
qqline(across$Richness)

shapiro.test(across$Richness)

# looking at dispersion in the data
dispmodel <- glm(Richness~Hydrology, data = across, family = poisson)
dispersiontest(dispmodel)

dispmodel <- glm(Richness~Vegetation, data = across, family = poisson)
dispersiontest(dispmodel)


## analysis of replicates (replicates averaged for each sampling date)
replicates<- across %>%
  group_by(Wetland, DateCollected, Vegetation, Hydrology, CoV) %>% 
  summarize(Richness = mean(Richness, na.rm = TRUE),
            .groups = "drop") %>%
  rename(Hydroperiod= Hydrology)

## anova of raw data, replicates/across time
richnessveg<- aov(Richness ~ Vegetation, data = replicates)
summary(richnessveg)
# p = 0.00101
emmeans(richnessveg, pairwise~Vegetation, type="response")
# Vegetation emmean   SE df lower.CL upper.CL
# Marsh        9.75 0.757 44     8.22    11.27
# Swamp        6.05 0.725 44     4.59     7.52
# p = 0.0010


richnesshydro<- aov(Richness ~ Hydroperiod, data = replicates)
summary(richnesshydro)
# p = 0.484
emmeans(richnesshydro, pairwise~Hydroperiod, type="response")
# Hydroperiod  emmean   SE df lower.CL upper.CL
# Intermediate   7.33 0.970 43     5.37     9.28
# Long           8.62 0.894 43     6.82    10.43
# Short          6.96 1.333 43     4.27     9.65
# intermediate/long: p = 0.5925
# intermediate/short: p = 0.9733
# long/short: p = 0.5597


richnessdate<- aov(Richness ~ DateCollected, data = replicates)
summary(richnessdate)
# p = 0.902
emmeans(richnessdate, pairwise~DateCollected, type="response")
# DateCollected emmean   SE df lower.CL upper.CL
# 2/9/2023        7.76 1.24 41     5.25    10.26
# 3/7/2023        6.83 1.30 41     4.21     9.46
# 4/12/2023       8.68 1.37 41     5.91    11.45
# 5/12/2023       8.23 1.56 41     5.09    11.37
# 6/22/2023       7.81 1.37 41     5.05    10.58
# no significance between dates 



# creating a null model for comparison of other models
null_model_rr <- glm.nb(Richness ~ 1, data = replicates, link = log)
summary (null_model_rr)

## HYDROLOGY

# RICHNESS: sampling date data, avg values from replicates so one value per
# sampling date per wetland

# looking at the relationship between hydroperiod length and richness
hydrology_rr <- glm.nb(Richness ~ Hydroperiod, data = replicates, link = log)
summary (hydrology_rr)

# anova
anova(null_model_rr, hydrology_rr, test = "Chisq")
# The insignificant p-value (0.4374649) indicates that including hydroperiod significantly improves model fit.
emmeans(hydrology_rr, pairwise~Hydroperiod, type="response")
# Predicted richness for each hydroperiod category:
# Intermediate     7.33 0.866 Inf      5.81      9.24
# Long             8.62 0.898 Inf      7.03     10.58
# Short            6.96 1.147 Inf      5.04      9.62
# Intermediate vs Long: no Significant difference (p = 0.5557); richness is lower for intermediate hydroperiods than for long ones.
# Intermediate vs Short: No significant difference (p = 0.9655).
# Long vs Short: no Significant difference (p = 0.5157); richness is higher for long hydroperiods than for short ones.


## analysis of totals
totals<- benthicrichness %>%
  group_by(Wetland) %>%
  summarise(Richness = n_distinct(Taxa)) %>%
  ungroup()

# adding in metadata
totals <- totals %>%
  mutate(Type = case_when(
    Wetland == "W15" ~ "Marsh",
    Wetland == "W21" ~ "Marsh",
    Wetland == "W53" ~ "Marsh",
    Wetland == "W46" ~ "Marsh",
    Wetland == "W37" ~ "Marsh",
    Wetland == "W42" ~ "Marsh",
    Wetland == "W11" ~ "Swamp",
    Wetland == "W68" ~ "Swamp",
    Wetland == "W58" ~ "Swamp",
    Wetland == "W32" ~ "Swamp",
    Wetland == "W52" ~ "Swamp",
  )) %>%
  relocate(Type, .before = 1) %>%
  mutate(HLength = case_when(
    Wetland == "W15" ~ "Short",
    Wetland == "W21" ~ "Intermediate",
    Wetland == "W53" ~ "Intermediate",
    Wetland == "W46" ~ "Long",
    Wetland == "W37" ~ "Short",
    Wetland == "W42" ~ "Long",
    Wetland == "W11" ~ "Long",
    Wetland == "W68" ~ "Intermediate",
    Wetland == "W58" ~ "Long",
    Wetland == "W32" ~ "Short",
    Wetland == "W52" ~ "Intermediate",
  )) %>%
  relocate(Type, .before = 1) %>%
  mutate(DaysInund = case_when(
    Wetland == "W15" ~ "113",
    Wetland == "W21" ~ "165",
    Wetland == "W53" ~ "189",
    Wetland == "W46" ~ "229",
    Wetland == "W37" ~ "142",
    Wetland == "W42" ~ "327",
    Wetland == "W11" ~ "345",
    Wetland == "W68" ~ "199",
    Wetland == "W58" ~ "314",
    Wetland == "W32" ~ "139",
    Wetland == "W52" ~ "286",
  )) %>%
  relocate(HLength, .before = 1)%>%
  mutate(CoV = case_when(
    Wetland == "W15" ~ "0.5944965",
    Wetland == "W21" ~ "0.3968312",
    Wetland == "W53" ~ "0.3195769",
    Wetland == "W46" ~ "0.3835618",
    Wetland == "W37" ~ "0.5491605",
    Wetland == "W42" ~ "0.2360643",
    Wetland == "W11" ~ "0.3916141",
    Wetland == "W68" ~ "0.3681222",
    Wetland == "W58" ~ "0.3163063",
    Wetland == "W32" ~ "0.4723248",
    Wetland == "W52" ~ "0.4596531",
  )) %>%
  relocate(CoV, .before = 1)

totals <- totals %>%
  rename(Hydroperiod= HLength)%>%
  rename(Vegetation= Type)

## anova of raw data, totals
richnessveg<- aov(Richness ~ Vegetation, data = totals)
summary(richnessveg)
# p =  0.615
emmeans(richnessveg, pairwise~Vegetation, type="response")
# Vegetation emmean   SE df lower.CL upper.CL
# Marsh        40.8 6.51  9     26.1     55.6
# Swamp        35.8 7.14  9     19.7     51.9
# p = 0.6150

richnesshydro<- aov(Richness ~ Hydroperiod, data = totals)
summary(richnesshydro)
# p = 0.234
emmeans(richnesshydro, pairwise~Hydroperiod, type="response")
# Hydroperiod  emmean   SE df lower.CL upper.CL
# Intermediate   37.8 7.16  8    21.23     54.3
# Long           47.8 7.16  8    31.23     64.3
# Short          27.3 8.27  8     8.26     46.4
# intermediate/long: p = 0.6047
# intermediate/short: p = 0.6251
# long/short: p = 0.2101


# creating a null model for comparison of other models
null_model_rt <- glm.nb(Richness ~ 1, data = totals, link = log)
summary (null_model_rt)

# Hydrology and total richness per site
# looking at the relationship between hydroperiod length and richness
hydrology_rt <- glm.nb(Richness ~ Hydroperiod, data = totals, link = log)
summary (hydrology_rt)

# anova
anova(null_model_rt, hydrology_rt, test = "Chisq")
# The insignificant p-value (0.1428474) indicates that including hydroperiod does not significantly improve model fit.
emmeans(hydrology_rt, pairwise~Hydroperiod, type="response")
# Predicted richness for each hydroperiod category:
# Intermediate: 37.8 (95% CI: 27.3 to 52.2)
# Long: 47.8 (95% CI: 34.8 to 65.5)
# Short: 27.3 (95% CI: 18.5 to 40.4)
# Intermediate vs Long: insignificant difference (p = 0.5656)
# Intermediate vs Short: No significant difference (p = 0.4262).
# Long vs Short: Significant difference (p = 0.0755)


## VEGETATION

# vegetation and average richness per sampling date
# looking at the relationship between vegetation and richness
vegetation_rr <- glm.nb(Richness ~ Vegetation, data = replicates, link = log)
summary (vegetation_rr)

# anova
anova(null_model_rr, vegetation_rr, test = "Chisq")
# The significant p-value (0.0004540495) indicates that including vegetation does improve model fit.
emmeans(vegetation_rr, pairwise~Vegetation, type="response")
# Predicted richness for each vegetation category:
# Marsh: 9.75 (95% CI: 8.25 to 11.51)
# Swamp: 6.05 (95% CI: 5.02 to 7.31)
# Marsh vs swamp: significant difference (p = 0.0002)

## analysis of totals

# vegetation and average richness per sampling date
# looking at the relationship between vegetation and richness
vegetation_rt <- glm.nb(Richness ~ Vegetation, data = totals, link = log)
summary (vegetation_rt)

# anova
anova(null_model_rt, vegetation_rt, test = "Chisq")
# insignificant p-value (0.5809764) indicates that including hydroperiod insignificantly improves model fit.
emmeans(vegetation_rt, pairwise~Vegetation, type="response")
# Predicted richness for each vegetation category:
# Marsh: 40.8 (95% CI: 29.9 to 55.7)
# Swamp: 35.8 (95% CI: 25.4 to 50.5)
# Marsh vs swamp: insignificant difference (p = 0.5778)


## testing additive and interactive effects to determine whether veg and hydrology are related
# Additive model
additiveR <- glm.nb(Richness ~ Vegetation + Hydroperiod, data = replicates, link = "log")
additiveR
library(car)
vif(additiveR)

# Interactive model
interactiveR <- glm.nb(Richness ~ Vegetation * Hydroperiod, data = replicates, link = "log")
interactiveR

emmeans(interactiveR, pairwise~Vegetation * Hydroperiod, type="response")

# Compare models
anova(additiveR, interactiveR, test = "Chisq")

# there is no significant interaction between vegetation and hydrology, (0.1770409),
# they both influence richness independently

## collection date
# looking at the relationship between collection date and density
date_rr <- glm.nb(Richness ~ DateCollected, data = replicates, link = log)
summary (date_rr)
# (Intercept)             2.048670   0.145118  14.117   <2e-16 ***
# DateCollected3/7/2023  -0.126857   0.214401  -0.592    0.554    
# DateCollected4/12/2023  0.111882   0.212796   0.526    0.599    
# DateCollected5/12/2023  0.059521   0.230408   0.258    0.796    
# DateCollected6/22/2023  0.007351   0.216086   0.034    0.973 

# anova
anova(null_model_rr, date_rr, test = "Chisq")
# p-value (0.8699036) 
emmeans(date_rr, pairwise~DateCollected, type="response")
# no significant differences

replicatesM <- replicates %>%
  filter(Vegetation == "Marsh")
replicatesS <- replicates %>%
  filter(Vegetation == "Swamp")
null_model_rrm <- glm.nb (Richness ~ 1, link = "log", data = replicatesM)
null_model_rrs <- glm.nb (Richness ~ 1, link = "log", data = replicatesS)


date_rrm <- glm.nb (Richness ~ DateCollected, link = "log", data = replicatesM)
summary (date_rrm)

# anova
anova(null_model_rrm, date_rrm, test = "Chisq")
# The significant p-value (0.01219) indicates that including hydroperiod does not significantly improves model fit.
emmeans(date_rrm, pairwise~DateCollected, type="response")


date_rrs <- glm.nb (Richness ~ DateCollected, link = "log", data = replicatesM)
summary (date_rrs)

# anova
anova(null_model_rrs, date_rrs, test = "Chisq")
# The significant p-value (0.01219) indicates that including hydroperiod does not significantly improves model fit.
emmeans(date_rrs, pairwise~DateCollected, type="response")


### Benthic biomass analysis
# importing the data and totaling biomass per sample
across<- read.csv ("BenthicBiomassTotaled.csv")
across <- across %>%
  group_by(Wetland, Date_collected, Type, HLength, CoV, Sample) %>% 
  summarize(Biomass = sum (tot_biomass, na.rm = TRUE), 
            Density = sum(Count, na.rm = TRUE),
            .groups = "drop") 

# renaming some things for ease
across <- across %>%
  mutate(Density = round(Density))%>%
  rename(DateCollected= Date_collected)%>%
  rename(Hydrology= HLength)%>%
  rename(Vegetation= Type)

# testing for normality 
qqnorm(across$Biomass)
qqline(across$Biomass)

shapiro.test(across$Biomass)

# checking for overdispersion in the data
dispmodel <- glm(Biomass~Hydrology, data = across, family = poisson)
dispersiontest(dispmodel)

dispmodel <- glm(Biomass~Vegetation, data = across, family = poisson)
dispersiontest(dispmodel)

## analysis of replicates (averaging biomass per collection date)
replicates<- across %>%
  group_by(Wetland, DateCollected, Vegetation, Hydrology, CoV) %>% 
  summarize(Biomass = mean(Biomass, na.rm = TRUE), 
            Density = mean(Density, na.rm = TRUE),
            .groups = "drop") %>%
  rename(Hydroperiod= Hydrology)

## anova of raw data, replicates/across time
biomassveg<- aov(Biomass ~ Vegetation, data = replicates)
summary(biomassveg)
# p = 0.0022
emmeans(biomassveg, pairwise~Vegetation, type="response")
# Vegetation emmean   SE df lower.CL upper.CL
# Marsh        99.1 14.5 44    69.99    128.3
# Swamp        34.0 13.8 44     6.12     61.9
# p = 0.0022

biomasshydro<- aov(Biomass ~ Hydroperiod, data = replicates)
summary(biomasshydro)
# p = 0.00409
emmeans(biomasshydro, pairwise~Hydroperiod, type="response")
# Hydroperiod  emmean   SE df lower.CL upper.CL
#  Intermediate   38.2 16.3 43     5.31     71.1
#  Long          105.0 15.0 43    74.63    135.3
#  Short          27.7 22.4 43   -17.53     72.9
# intermediate/long: p = 0.0119
# intermediate/short: p = 0.9237
# long/short: p = 0.0174

biomassdate<- aov(Biomass ~ DateCollected, data = replicates)
summary(biomassdate)
# p = 0.0415
emmeans(biomassdate, pairwise~DateCollected, type="response")
# DateCollected emmean   SE df lower.CL upper.CL
# 2/9/2023        54.1 21.0 41    11.77     96.4
# 3/7/2023        77.5 22.0 41    33.06    121.9
# 4/12/2023       47.7 23.2 41     0.92     94.5
# 5/12/2023      135.6 26.3 41    82.52    188.7
# 6/22/2023       27.7 23.2 41   -19.10     74.5

# no significance between dates besides may and june, 
# biomass highest in may and lowest in june
#  (5/12/2023) - (6/22/2023), p = 0.0287



## creating models and model analysis
# creating a null model for comparison of other models
null_model_br <- glm(Biomass ~ 1, family = Gamma (link = "log"), data = replicates)
summary (null_model_br)

## HYDROLOGY

# replicates: sampling date data, avg values from replicates so one value per
# sampling date per wetland

# looking at the relationship between hydroperiod length and biomass in per collection per site
hydrology_br <- glm (Biomass ~ Hydroperiod, family = Gamma (link = "log"), data = replicates)
summary (hydrology_br)

# anova
anova(null_model_br, hydrology_br, test = "Chisq")
# The significant p-value (0.0002171) indicates that including hydroperiod does significantly improve model fit.
emmeans(hydrology_br, pairwise~Hydroperiod, type="response")
# Predicted biomass for each hydroperiod category:
# Intermediate: 38.2 (95% CI: 24.3     60.1)
# Long: 105.0 (95% CI: 69.1    159.4)
# Short: 27.7 (95% CI: 14.8     51.6)
# Intermediate vs Long: Insignificant difference (p =  0.0054); biomass is lower for intermediate hydroperiods than for long ones.
# Intermediate vs Short: Significant difference (p = 0.6783).
# Long vs Short: Insignificant difference (p = 0.0024); biomass is higher for long hydroperiods than for short ones.


## analysis of averages
across <- across %>%
  rename(Hydroperiod= Hydrology)
totals<- across %>%
  group_by(Wetland, Vegetation, Hydroperiod, CoV) %>% 
  summarize(Biomass = mean(Biomass, na.rm = TRUE), 
            Density = mean(Density, na.rm = TRUE),
            .groups = "drop") 

## anova of raw data
biomassveg<- aov(Biomass ~ Vegetation, data = totals)
summary(biomassveg)
# p =  0.0812
emmeans(biomassveg, pairwise~Vegetation, type="response")


biomasshydro<- aov(Biomass ~ Hydroperiod, data = totals)
summary(biomasshydro)
# p = 0.153
emmeans(biomasshydro, pairwise~Hydroperiod, type="response")
# Hydroperiod  emmean   SE df lower.CL upper.CL
#  Intermediate   38.2 16.3 43     5.31     71.1
#  Long          105.0 15.0 43    74.63    135.3
#  Short          27.7 22.4 43   -17.53     72.9
# intermediate/long: p = 0.0119
# intermediate/short: p = 0.9237
# long/short: p = 0.0174



# creating a null model for comparison of other models
null_model_bt <- glm(Biomass ~ 1, family = Gamma (link = "log"), data = totals)
summary (null_model_bt)

# Hydrology and total richness per site
# looking at the relationship between hydroperiod length and richness
hydrology_bt <- glm (Biomass ~ Hydroperiod, family = Gamma (link = "log"), data = totals)
summary (hydrology_bt)

# anova
anova(null_model_bt, hydrology_bt, test = "Chisq")
# The significant p-value (0.06132) indicates that including hydroperiod significantly model fit.
emmeans(hydrology_bt, pairwise~Hydroperiod, type="response")
# Hydroperiod  response    SE df lower.CL upper.CL
# Intermediate      162  44.5  8     86.3      305
# Long              525 143.7  8    279.0      987
# Short              83  26.3  8     40.0      172

#contrast             ratio    SE df null t.ratio p.value
# Intermediate / Long  0.309 0.120  8    1  -3.028  0.0389
# Intermediate / Short 1.956 0.818  8    1   1.603  0.2989
# Long / Short         6.321 2.645  8    1   4.407  0.0057

## VEGETATION

# vegetation and average biomass per sampling date
# looking at the relationship between vegetation and biomass
vegetation_br <- glm (Biomass ~ Vegetation, family = Gamma (link = "log"), data = replicates)
summary (vegetation_br)

# anova
anova(null_model_br, vegetation_br, test = "Chisq")
# The significant p-value (0.0002533) indicates that including vegetation improves model fit.
emmeans(vegetation_br, pairwise~Vegetation, type="response")
# Predicted biomass for each vegetation category:
# Marsh: 99.1 (95% CI: 65.2 to 150.7)
# Swamp: 34.0 (95% CI: 22.8 to 50.8)
# Marsh vs swamp: significant difference (p = 0.0006)

## analysis of totals

# vegetation and average biomass per sampling date
# looking at the relationship between vegetation and biomass
vegetation_bt <- glm (Biomass ~ Vegetation, family = Gamma (link = "log"), data = totals)
summary (vegetation_bt)

# Perform the likelihood ratio test
anova(null_model_bt, vegetation_bt, test = "Chisq")
# Compares the null model (without hydroperiod as a predictor) to the model with hydroperiod as a predictor.
# insignificant p-value (0.1231) indicates that including vegetation does not improve model fit.
emmeans(vegetation_bt, pairwise~Vegetation, type="response")
# Predicted biomass for each vegetation category:
# Marsh: 364 (95% CI: 167.9 to 787)
# Swamp: 163 (95% CI: 70.1 to 381)
# Marsh vs swamp: insignificant difference (p =  0.1486)

## testing additive and interactive effects to determine whether veg and hydrology are related
# Additive model
additiveB <- glm.nb(Biomass ~ Vegetation + Hydroperiod, data = replicates, link = "log")
additiveB
library(car)
vif(additiveB)

# Interactive model
interactiveB <- glm.nb(Biomass ~ Vegetation * Hydroperiod, data = replicates, link = "log")
interactiveB

emmeans(interactiveB, pairwise~Vegetation * Hydroperiod, type="response")

# Compare models
anova(additiveB, interactiveB, test = "Chisq")

## collection date
# looking at the relationship between collection date and biomass
replicatesM <- replicates %>%
  filter(Vegetation == "Marsh")
replicatesS <- replicates %>%
  filter(Vegetation == "Swamp")
null_model_brm <- glm(Biomass ~ 1, family = Gamma (link = "log"), data = replicatesM)
null_model_brs <- glm(Biomass ~ 1, family = Gamma (link = "log"), data = replicatesS)

date_br <- glm (Biomass ~ DateCollected, family = Gamma (link = "log"), data = replicates)
summary (date_br)

# anova
anova(null_model_br, date_br, test = "Chisq")

date_brm <- glm (Biomass ~ DateCollected, family = Gamma (link = "log"), data = replicatesM)
summary (date_brm)

# anova
anova(null_model_brm, date_brm, test = "Chisq")
# The significant p-value (0.01219) indicates that including hydroperiod does not significantly improves model fit.
emmeans(date_brm, pairwise~DateCollected, type="response")


date_brs <- glm (Biomass ~ DateCollected, family = Gamma (link = "log"), data = replicatesS)
summary (date_brs)

# anova
anova(null_model_brs, date_brs, test = "Chisq")
# The significant p-value (0.01219) indicates that including hydroperiod does not significantly improves model fit.
emmeans(date_brs, pairwise~DateCollected, type="response")
