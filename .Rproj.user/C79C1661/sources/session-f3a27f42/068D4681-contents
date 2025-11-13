### Total emergence biomass and richness per site, Figure 11
# loading packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpmisc)

## starting with biomass first
# reading in the data
emergebiomass <- read.csv ("EmergeBiomassTotaled.csv")

# visualizing biomass by wetland
# need to start by averaging replicates then totaling replicates
# totaling biomass per sample
totalbiomasspersite <- emergebiomass %>%
  group_by(Wetland, DateCollected, Type, HLength, Sample) %>% 
  summarize(Biomass = sum (tot_biomass, na.rm = TRUE), 
            Density = sum (Density, na.rm = TRUE),
            .groups = "drop") 
# averaging biomass per collection date (taking an average across replicates)
totalbiomasspersite <- totalbiomasspersite %>%
  group_by(Wetland, DateCollected, Type, HLength) %>% 
  summarize(avg_biomass = mean(Biomass, na.rm = TRUE), 
            Density = mean(Density, na.rm = TRUE),
            .groups = "drop") 
# totaling biomass
totalbiomasspersite <- totalbiomasspersite %>%
  group_by(Wetland) %>% 
  summarize(Biomass = sum (avg_biomass, na.rm = TRUE), 
            Density = sum (Density, na.rm = TRUE),
            .groups = "drop") 

# getting an annual estimation of biomass/flux 
totalbiomasspersite <- totalbiomasspersite %>%
  group_by(Wetland) %>%
  mutate(YearBiomass = case_when (
    Wetland == "W15" ~ (Biomass / 168) *364,
    Wetland == "W21" ~ (Biomass / 168) *364,
    Wetland == "W53" ~ (Biomass / 168) *364,
    Wetland == "W46" ~ (Biomass / 168) *364,
    Wetland == "W37" ~ (Biomass / 168) *364,
    Wetland == "W42" ~ (Biomass / 168) *364,
    Wetland == "W11" ~ (Biomass / 168) *364,
    Wetland == "W68" ~ (Biomass / 168) *364,
    Wetland == "W58" ~ (Biomass / 168) *364,
    Wetland == "W32" ~ (Biomass / 168) *364,
    Wetland == "W52" ~ (Biomass / 168) *364
  ))

# getting an annual estimation of emergence density
totalbiomasspersite <- totalbiomasspersite %>%
  group_by(Wetland) %>%
  mutate(YearDensity = case_when (
    Wetland == "W15" ~ (Density / 168) *364,
    Wetland == "W21" ~ (Density / 168) *364,
    Wetland == "W53" ~ (Density / 168) *364,
    Wetland == "W46" ~ (Density / 168) *364,
    Wetland == "W37" ~ (Density / 168) *364,
    Wetland == "W42" ~ (Density / 168) *364,
    Wetland == "W11" ~ (Density / 168) *364,
    Wetland == "W68" ~ (Density / 168) *364,
    Wetland == "W58" ~ (Density / 168) *364,
    Wetland == "W32" ~ (Density / 168) *364,
    Wetland == "W52" ~ (Density / 168) *364
  ))

# adding in metadata
totalbiomasspersite <- totalbiomasspersite %>%
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
    Wetland == "W15" ~ "106.0556",
    Wetland == "W21" ~ "163.9706",
    Wetland == "W53" ~ "188.7879",
    Wetland == "W46" ~ "224.4815",
    Wetland == "W37" ~ "142.8125",
    Wetland == "W42" ~ "325.7391",
    Wetland == "W11" ~ "345.5263",
    Wetland == "W68" ~ "201.2759",
    Wetland == "W58" ~ "329.1905",
    Wetland == "W32" ~ "136.3939",
    Wetland == "W52" ~ "284.5000",
  )) %>%
  relocate(HLength, .before = 1)

# convert DaysInund to numeric during creation
totalbiomasspersite <- totalbiomasspersite %>%
  mutate(DaysInund = case_when(
      Wetland == "W15" ~ 106.0556,
      Wetland == "W21" ~ 163.9706,
      Wetland == "W53" ~ 188.7879,
      Wetland == "W46" ~ 224.4815,
      Wetland == "W37" ~ 142.8125,
      Wetland == "W42" ~ 325.7391,
      Wetland == "W11" ~ 345.5263,
      Wetland == "W68" ~ 201.2759,
      Wetland == "W58" ~ 329.1905,
      Wetland == "W32" ~ 136.3939,
      Wetland == "W52" ~ 284.5000
    ))

# plotting the flux side of figure 11
Figure11Flux <- ggplot(totalbiomasspersite, aes(x = DaysInund, y = YearBiomass, 
                                      color = Type, shape = Type)) +
  geom_point(size = 4) +
  # creating regressions with confidence bands
  geom_smooth(method = "lm", aes(fill = Type), alpha = 0.15, color = NA) + 
  geom_smooth(method = "lm", se = FALSE, size = 1) + 
  # regression equations to be displayed on the figure
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    size = 5,
    label.x.npc = "center",
    label.y.npc = c(0.9, 0.8)
  ) +
  labs(x = "Average Inundation Length (Days)", 
       y = bquote('Total Flux'~(mg/m^2/year))) +
  scale_color_manual(values = c("Marsh" = "#79B4A9", "Swamp" = "#354C82")) +
  scale_fill_manual(values = c("Marsh" = "#79B4A9", "Swamp" = "#354C82")) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "grey", fill = NA),
        axis.line = element_line(color = "grey"),
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        strip.text = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))
Figure11Flux

# making an alternative figure without the regressions and confidence bands
Fluxalternate <- ggplot(totalbiomasspersite, aes(x = DaysInund, y = YearBiomass)) +
  geom_point(size =4, aes(color = Type, shape = Type)) +
  labs(x = "Average Inundation Length (Days)", 
       y = bquote('Total Flux'~(mg/m^2/year))) +
  scale_color_manual(values = c("Marsh" = "#79B4A9","Swamp" = "#354C82"))+
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "grey", fill = NA),
        axis.line = element_line(color = "grey"),
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        strip.text = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))
Fluxalternate

## Now calculating and visualizing total richness

emergencedata <- read.csv ("EmergeBiomassTotaled.csv")

# mapping to make sure everything is classified as the correct taxa (getting rid of 
# m vs f distinctions, typos missed during data cleaning, etc)
taxa_mapping <- c(
  "mTanypodinae" = "Tanypodinae",
  "fTanypodinae" = "Tanypodinae",
  "mOrthocladiinae" = "Orthocladiinae",
  "fOrthocladiinae" = "Orthocladiinae",
  "mChironominae" = "Chironominae",
  "fChironominae" = "Chironominae",
  "mChironomidae" = "Chironomidae",
  "fChironomidae" = "Chironomidae",
  "mAnopheles" = "Anopheles",
  "fAnopheles" = "Anopheles",
  "mCulicine" = "Culicine",
  "fCulicine" = "Culicine",
  "mCulicidae" = "Culicidae",
  "fCulicidae" = "Culicidae",
  "mChironomini" = "Chironomini",
  "mTanytarsini" = "Tanytarsini",
  "mTarytarsini" = "Tanytarsini",
  "Ceratopogonoidae" = "Ceratopogonidae"
)

# setting a taxonomic hierarchy to prevent richness from becoming inflated due to 
# differences in taxonomic resolution in chironomidae
taxa_hierarchy <- c(
  "Chironominae" = 2,
  "Chironomidae" = 3
)


# data for all richness calculations besides per sample richness
resolved_data <- emergencedata %>%
  # combine taxa using mapping
  mutate(TaxaCode = ifelse(Taxa %in% names(taxa_mapping), 
                           taxa_mapping[Taxa], 
                           Taxa)) %>%
  # assign taxonomic levels with default set to 1
  mutate(TaxonomicLevel = ifelse(Taxa %in% names(taxa_hierarchy), 
                                 taxa_hierarchy[Taxa], 1)) %>%
  # resolve mixed taxonomic resolutions by retaining the finest level
  group_by(Wetland, DateCollected) %>%
  filter(TaxonomicLevel == min(TaxonomicLevel)) %>%
  ungroup()

# calculating total richness per wetland
totalrichness <- resolved_data %>%
  group_by(Wetland) %>%
  summarise(richness = n_distinct(Taxa)) %>%
  ungroup()

# total richness across all wetlands, just to have the number
totaltotalrichness <- resolved_data %>%
  summarise(richness = n_distinct(TaxaCode)) %>%
  ungroup()

# adding in metadata to total richness per wetland
totalrichnesspersite <- totalrichness %>%
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
    Wetland == "W15" ~ 106.0556,
    Wetland == "W21" ~ 163.9706,
    Wetland == "W53" ~ 188.7879,
    Wetland == "W46" ~ 224.4815,
    Wetland == "W37" ~ 142.8125,
    Wetland == "W42" ~ 325.7391,
    Wetland == "W11" ~ 345.5263,
    Wetland == "W68" ~ 201.2759,
    Wetland == "W58" ~ 329.1905,
    Wetland == "W32" ~ 136.3939,
    Wetland == "W52" ~ 284.5000
  )) %>%
  relocate(HLength, .before = 1)

Figure11Richness <- ggplot(totalrichnesspersite, aes(x = DaysInund, y = richness, 
                                       color = Type, shape = Type)) +
  geom_point(size = 4) +
  # adding confidence bands and regression lines for each wetland type
  geom_smooth(method = "lm", aes(fill = Type), alpha = 0.15, color = NA) +  
  geom_smooth(method = "lm", se = FALSE, size = 1) + 
  # adding regression equations for each wetland type
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    size = 5,
    label.x.npc = "center",
    label.y.npc = c(0.9, 0.8) # prevents overlap of two equations
  ) +
  scale_color_manual(values = c("Marsh" = "#79B4A9", "Swamp" = "#354C82")) +
  scale_fill_manual(values = c("Marsh" = "#79B4A9", "Swamp" = "#354C82")) +
  labs(x = "Average Inundation Length (Days)", 
       y = "Taxa Richness") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "grey", fill = NA),
        axis.line = element_line(color = "grey"),
        axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        strip.text = element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))
Figure11Richness

# adding both the plots together to make the complete Figure 11
library(patchwork)
Figure11Flux+Figure11Richness

