### Proportion of insects to non-insects across time, Figure 6

# loading in packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(emmeans)

# reading in the data
benthicbiomass <- read.csv ("BenthicBiomassTotaled.csv")

# total biomass per sample per collection date
insects <- benthicbiomass %>%
  group_by(Wetland, Date_collected, Type, HLength, Sample, Insect) %>% 
  summarize(Biomass = sum (tot_biomass, na.rm = TRUE), 
            Density = sum(Count, na.rm = TRUE),
            .groups = "drop") 

# average biomass per collection date
insects <- insects %>%
  group_by(Wetland, Date_collected, Type, HLength, Insect) %>% 
  summarize(Biomass = mean (Biomass, na.rm = TRUE), 
            Density = mean (Density, na.rm = TRUE),
            .groups = "drop") 
# fixing dates
insects<- insects %>%
  mutate(Date_collected = as.Date(Date_collected, format = "%m/%d/%Y")) %>%
  mutate(DateCollected = case_when(
    Date_collected == as.Date("2023-02-09") ~ "Feb",
    Date_collected == as.Date("2023-03-07") ~ "Mar",
    Date_collected == as.Date("2023-04-12") ~ "Apr",
    Date_collected == as.Date("2023-05-12") ~ "May",
    Date_collected == as.Date("2023-06-22") ~ "Jun",
    TRUE ~ as.character(Date_collected)
  ))%>%
  relocate(Date_collected, .before = 1)


# Making figure 6
# making sure date and other factors are in the right order before plotting
insects$DateCollected = factor(insects$DateCollected, levels = c("Feb", "Mar", "Apr", "May", "Jun"))
insects$Insect = factor(insects$Insect, levels = c("NonInsect", "Insect"))

# plotting
Figure6 <- ggplot(insects, aes(x = DateCollected, y = Density, fill = Insect)) +
  geom_bar(position="fill", stat="identity")+
  labs(x = "Collection Period", y = "Proportion of Insects/Non-Insects") +
  scale_fill_manual(values = c("Insect" = "#6A88B4", "NonInsect" = "#CDD6DF"), 
                    name = "Type") +
  facet_wrap(~ Type, scales = "free_x") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "grey", fill = NA),
        axis.line = element_line(color = "grey"),
        axis.text = element_text(size = 14), # Increase axis text size
        axis.title = element_text(size = 16), # Increase axis title text size
        strip.text = element_text(size = 14), # Increase facet label text size
        legend.title = element_text(size = 14), # Increase legend title text size
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))
Figure6

## Analysis for Figure 6

# creating proportion dataset for analysis
proportion_data <- insects %>%
  group_by(Wetland, Date_collected, Type, HLength) %>%
  summarise(
    Insect_Biomass = sum(Biomass[Insect == "Insect"]),
    NonInsect_Biomass = sum(Biomass[Insect == "NonInsect"]),
    Insect_Density = sum(Density[Insect == "Insect"]),
    NonInsect_Density = sum(Density[Insect == "NonInsect"]),
    Biomass_Proportion = Insect_Biomass / NonInsect_Biomass * 100,
    Density_Proportion = Insect_Density / NonInsect_Density * 100,
  )
proportion_data <- proportion_data %>%
  filter(!is.infinite(Biomass_Proportion) & !is.infinite(Density_Proportion))
proportion_data<- proportion_data %>%
  mutate(Date_collected = as.Date(Date_collected, format = "%m/%d/%Y")) %>%
  mutate(DateCollected = case_when(
    Date_collected == as.Date("2023-02-09") ~ "Feb",
    Date_collected == as.Date("2023-03-07") ~ "Mar",
    Date_collected == as.Date("2023-04-12") ~ "Apr",
    Date_collected == as.Date("2023-05-12") ~ "May",
    Date_collected == as.Date("2023-06-22") ~ "Jun",
    TRUE ~ as.character(Date_collected)
  ))%>%
  relocate(Date_collected, .before = 1)

# separating into swamp and marsh datasets
proportion_dataM <- proportion_data %>%
  filter(Type == "Marsh")
proportion_dataS <- proportion_data %>%
  filter(Type == "Swamp")

# creating a null model for comparison of other models
nullS <- glm (Density_Proportion ~ 1, family = Gamma (link = "log"), data = proportion_dataS)
nullM <- glm (Density_Proportion ~ 1, family = Gamma (link = "log"), data = proportion_dataM)

# creating the models with date collected as a predictor
propmodelS <- glm (Density_Proportion ~ DateCollected, family = Gamma (link = "log"), data = proportion_dataS)
summary (propmodelS)
propmodelM <- glm (Density_Proportion ~ DateCollected, family = Gamma (link = "log"), data = proportion_dataM)
summary (propmodelM)

# comparing the null model (without date as a predictor) to the model with date as a predictor
anova(nullS, propmodelS, test = "Chisq")
# The significant p-value (2.057e-09) indicates that including collection date significantly improves model fit.

# looking at differences in proportion between each pair of months
emmeans(propmodelS, pairwise~DateCollected, type="response")

# comparing the null model (without date as a predictor) to the model with date as a predictor
anova(nullM, propmodelM, test = "Chisq")
# The significant p-value (0.007027) indicates that including collection date significantly improves model fit.

# looking at differences in proportion between each pair of months
emmeans(propmodelM, pairwise~ DateCollected, type="response")


## looking at everything combined, not separating by marsh and swamp
# creating null model
null <- glm (Density_Proportion ~ 1, family = Gamma (link = "log"), data = proportion_data)
summary (null)
# creating model with date collected as a predictor
propmodel <- glm (Density_Proportion ~ DateCollected, family = Gamma (link = "log"), data = proportion_data)
summary (propmodel)

# comparing the null model (without date as a predictor) to the model with date as a predictor
anova(null, propmodel, test = "Chisq")
# The significant p-value (0.007027) indicates that including collection date significantly improves model fit.
# (aka collection data is predictive of the proportion of insects to non-insects regardless of veg)

### looking at how density and biomass proportions of insects vs non-insects change over time
## this information is referenced in the results section "Proportion of insects to non-insects across time"
# loading in packages
library(tidyr)
library(dplyr)
library(ggplot2)

# reading in data
benthicbiomass <- read.csv ("BenthicBiomassTotaled.csv")

# total biomass per sample per collection date
insects <- benthicbiomass %>%
  group_by(Wetland, Date_collected, Type, HLength, Sample, Insect) %>% 
  summarize(Biomass = sum (tot_biomass, na.rm = TRUE), 
            Density = sum(Count, na.rm = TRUE),
            .groups = "drop") 

# average biomass per collection date
insects <- insects %>%
  group_by(Wetland, Date_collected, Type, HLength, Insect) %>% 
  summarize(Biomass = mean (Biomass, na.rm = TRUE), 
            Density = mean (Density, na.rm = TRUE),
            .groups = "drop") 

proportion_data <- insects %>%
  group_by(Date_collected, Type) %>%
  summarise(
    Insect_Biomass = sum(Biomass[Insect == "Insect"]),
    NonInsect_Biomass = sum(Biomass[Insect == "NonInsect"]),
    Insect_Density = sum(Density[Insect == "Insect"]),
    NonInsect_Density = sum(Density[Insect == "NonInsect"]),
    Biomass_Proportion = Insect_Biomass / NonInsect_Biomass * 100,
    Density_Proportion = Insect_Density / NonInsect_Density * 100,
    Biomass_T = Insect_Biomass + NonInsect_Biomass,
    Density_T = Insect_Density + NonInsect_Density,
    Biomass_Percent = NonInsect_Biomass/Biomass_T * 100,
    Density_Percent = NonInsect_Density/Density_T * 100,
    Biomass_PercentI = Insect_Biomass/Biomass_T * 100,
    Density_PercentI = Insect_Density/Density_T * 100,
  )