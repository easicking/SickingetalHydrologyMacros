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


