## Changes in emergence biomass and richness over time, figures 10 and 8
library(tidyr)
library(dplyr)
library(ggplot2)

# reading in the data
emergebiomass <- read.csv ("EmergeBiomassTotaled.csv")

biomasspersite <- emergebiomass %>%
  group_by(Wetland, DateCollected, Type, HLength, CoV, Sample) %>% 
  summarize(Biomass = sum (tot_biomass, na.rm = TRUE), 
            Density = sum (Density, na.rm = TRUE),
            .groups = "drop") 

# dividing the biomass collected on a given date by the number of days the trap
# had been collecting, therefore measuring the rate of flux (biomass/m^2/day)
biomasspersite <- biomasspersite %>%
  mutate(DateCollected = as.Date(DateCollected, format = "%m/%d/%Y")) %>% 
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

# dividing the density of insects collected on a given date by the number of days the trap
# had been collecting, therefore measuring the rate of emergence (density/m^2/day)
biomasspersite <- biomasspersite %>%
  mutate(DateCollected = as.Date(DateCollected, format = "%m/%d/%Y")) %>% 
  group_by(DateCollected) %>%
  mutate(Emergence = case_when(
    DateCollected == as.Date("2023-02-21") ~ Density / 7,
    DateCollected == as.Date("2024-02-21") ~ Density / 7,
    DateCollected == as.Date("2023-03-06") ~ Density / 13,
    DateCollected == as.Date("2023-03-05") ~ Density / 13,
    DateCollected == as.Date("2023-03-21") ~ Density / 15,
    DateCollected == as.Date("2023-04-05") ~ Density / 15,
    DateCollected == as.Date("2023-04-20") ~ Density / 15,
    DateCollected == as.Date("2023-05-03") ~ Density / 13,
    DateCollected == as.Date("2024-05-03") ~ Density / 13,
    DateCollected == as.Date("2023-05-12") ~ Density / 9,
    DateCollected == as.Date("2023-05-17") ~ Density / 14,
    DateCollected == as.Date("2023-05-31") ~ Density / 14,
    DateCollected == as.Date("2023-06-14") ~ Density / 14,
    DateCollected == as.Date("2023-06-15") ~ Density / 15,
    DateCollected == as.Date("2023-06-23") ~ Density / 23,
    DateCollected == as.Date("2023-06-22") ~ Density / 22,
    DateCollected == as.Date("2023-07-05") ~ Density / 16,
    DateCollected == as.Date("2023-07-06") ~ Density / 16,
    DateCollected == as.Date("2023-07-19") ~ Density / 14,
    DateCollected == as.Date("2024-07-19") ~ Density / 14,
    DateCollected == as.Date("2023-08-01") ~ Density / 13,
    TRUE ~ Density
  ))

# simplifying dates into month vs day of month
biomasspersite <- biomasspersite %>%
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

# averaging biomass, density, flux, and emergence values per collection date and wetland
biomassacross<- biomasspersite %>%
  group_by(Wetland, DateCollected, Type, HLength) %>% 
  summarize(Biomass = mean (Biomass, na.rm = TRUE), 
            Density = mean (Density, na.rm = TRUE),
            Flux = mean (Flux, na.rm = TRUE),
            Emergence = mean (Emergence, na.rm = TRUE),
            .groups = "drop") 

# ensuring months are in the correct order
biomassacross$DateCollected = factor(biomassacross$DateCollected, levels = c("Feb", "Mar", "Apr",
                                                                             "May", "Jun", "Jul",
                                                                             "Aug"))
# ensuring hydroperiod is in the correct order
biomassacross$HLength = factor(biomassacross$HLength, levels=c("Short", "Intermediate",
                                                               "Long"))

# plotting figure 10
Figure10 <- ggplot(biomassacross, aes(x = DateCollected, y = Flux, fill = HLength)) +
  geom_boxplot() +
  labs(x = "Collection Period", y = bquote('Flux'~(mg/m^2/day))) +
  scale_fill_manual(values = c("Short" = "#DBDEE6", "Intermediate" = "#6A88B4",
                               "Long" = "#354C82"), name = "Hydrology") +
  facet_wrap(~ HLength, scales = "free_x") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "grey", fill = NA),
        axis.line = element_line(color = "grey"),
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        strip.text = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
)
Figure10

## RICHNESS ACROSS TIME VISUALIZATIONS (Figure 8)
# reading in richness data
richness<- read.csv("AnalysisDataPrepped.csv")

# simplifying dates for visualizations
richnessacross <- richness %>%
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

# making sure dates and hydroperiods are ordered correctly
richnessacross$DateCollected = factor (richnessacross$DateCollected, levels = c("Feb", "Mar", "Apr",
                                                                                "May", "Jun", "Jul",
                                                                                "Aug"))
richnessacross$HLength = factor (richnessacross$HLength, levels=c("Short", "Intermediate",
                                                                  "Long"))

# visualizing richness across the hydroperiod, grouped by hydrology
Figure8 <- ggplot(richnessacross, aes(x = DateCollected, y = Richness, fill = HLength)) +
  geom_boxplot() +
  labs(x = "Collection Period", y = "Taxa Richness") +
  scale_fill_manual(values = c("Short" = "#DBDEE6", "Intermediate" = "#6A88B4",
                               "Long" = "#354C82"), name = "Hydrology") +
  facet_wrap(~ HLength, scales = "free_x") +
  coord_cartesian(ylim = c(1, 13)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "grey", fill = NA),
        axis.line = element_line(color = "grey"),
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        strip.text = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))
Figure8

### Benthic biomass and richness across time, figures 7 and 9
library(tidyr)
library(dplyr)
library(ggplot2)

benthicbiomass <- read.csv ("BenthicBiomassTotaled.csv")

# total biomass per sample per collection date
biomasspersite <- benthicbiomass %>%
  group_by(Wetland, Date_collected, Type, HLength, CoV, Sample) %>% 
  summarize(Biomass = sum (tot_biomass, na.rm = TRUE), 
            Density = sum(Count, na.rm = TRUE),
            .groups = "drop") 

# average biomass per collection date
biomasspersite <- biomasspersite %>%
  group_by(Wetland, Date_collected, Type, HLength, CoV) %>% 
  summarize(Biomass = mean (Biomass, na.rm = TRUE), 
            Density = mean (Density, na.rm = TRUE),
            .groups = "drop") 
# making dates into
biomassacross <- biomasspersite %>%
  mutate(Date_collected = as.Date(Date_collected, format = "%m/%d/%Y")) %>%
  mutate(Date_collected = case_when(
    Date_collected == as.Date("2023-02-09") ~ "Feb",
    Date_collected == as.Date("2023-03-07") ~ "Mar",
    Date_collected == as.Date("2023-04-12") ~ "Apr",
    Date_collected == as.Date("2023-05-12") ~ "May",
    Date_collected == as.Date("2023-06-22") ~ "Jun",
    TRUE ~ as.character(Date_collected)
  ))%>%
  relocate(Date_collected, .before = 1)

# making sure collection month and hydroperiod are ordered correctly
biomassacross$DateCollected = factor(biomassacross$Date_collected, levels = c("Feb", "Mar", "Apr", "May", "Jun"))
biomassacross$HLength = factor(biomassacross$HLength, levels=c("Short", "Intermediate",
                                                               "Long"))

# visualizing biomass across the hydroperiod, grouped by veg (Figure 9)
Figure9 <- ggplot(biomassacross, aes(x = DateCollected, y = Biomass, fill = Type)) +
  geom_boxplot() +
  labs(x = "Collection Period", y = bquote('Biomass'~(mg/m^2))) +
  scale_fill_manual(values = c("Marsh" = "#79B4A9",
                               "Swamp" = "#22577A"), name = "Vegetation") +
  facet_wrap(~ Type, scales = "free_x") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "grey", fill = NA),
        axis.line = element_line(color = "grey"),
        axis.text = element_text(size = 14), # Increase axis text size
        axis.title = element_text(size = 16), # Increase axis title text size
        strip.text = element_text(size = 14), # Increase facet label text size
        legend.title = element_text(size = 14), # Increase legend title text size
        legend.text = element_text(size = 12))
Figure9


## RICHNESS ACROSS THE HYDROPERIOD
#reading in data
benthicrichness <- read.csv ("BenthicBiomassTotaled.csv")

# calculating richness based on distinct number of taxa per wetland/collection date/sample
richness <- benthicrichness %>%
  group_by(Wetland, Date_collected, Sample) %>%
  summarise(richness = n_distinct(Taxa)) %>%
  ungroup()

# adding in metadata
richnesspersite <- richness %>%
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
  relocate(HLength, .before = 1)

# simplifying date into month collected
richnessacross<- richnesspersite %>%
  mutate(Date_collected = as.Date(Date_collected, format = "%m/%d/%Y")) %>%
  mutate(Date_collected = case_when(
    Date_collected == as.Date("2023-02-09") ~ "Feb",
    Date_collected == as.Date("2023-03-07") ~ "Mar",
    Date_collected == as.Date("2023-04-12") ~ "Apr",
    Date_collected == as.Date("2023-05-12") ~ "May",
    Date_collected == as.Date("2023-06-22") ~ "Jun",
    TRUE ~ as.character(Date_collected)
  ))%>%
  relocate(Date_collected, .before = 1)

# ensuring dates and hydroperiods are ordered correctly
richnessacross$DateCollected = factor(richnessacross$Date_collected, levels = c("Feb", "Mar", "Apr", "May", "Jun"))
richnessacross$HLength = factor(richnessacross$HLength, levels=c("Short", "Intermediate",
                                                                 "Long"))

# visualizing richness across the hydroperiod, grouped by veg (Figure 7)
Figure7 <- ggplot(richnessacross, aes(x = DateCollected, y = richness, fill = Type)) +
  geom_boxplot() +
  labs(x = "Collection Period", y = "Taxa Richness") +
  scale_fill_manual(values = c("Marsh" = "#79B4A9",
                               "Swamp" = "#22577A"), name = "Hydrology") +
  facet_wrap(~ Type, scales = "free_x") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "grey", fill = NA),
        axis.line = element_line(color = "grey"),
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        strip.text = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))
Figure7

