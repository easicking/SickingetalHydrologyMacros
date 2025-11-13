### Benthic trait-based community composition NMDS

# loading in packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)

# reading in the data
tot_taxa_biomass <- read.csv("BenthicBiomassTotaledwTraits.csv")

# averaging replicates
total_taxa <- tot_taxa_biomass %>%
  group_by(Wetland, Date_collected, Type, HLength, CoV, FFG) %>% 
  summarize(avg_biomass = mean(tot_biomass, na.rm = TRUE), .groups = "drop") 

# averaging biomass per site
total_biomass <- total_taxa %>%
  group_by(Wetland, FFG)%>%
  summarize(total_biomass = mean(avg_biomass, na.rm = TRUE), .groups = "drop")

# reshaping to wide format
total_nmds_data <- total_biomass %>%
  group_by(Wetland) %>% 
  pivot_wider(
    names_from = FFG,
    values_from = total_biomass,
    values_fill = 0
  ) %>%
  ungroup()

# adding back in metadata
total_nmds_data <- total_nmds_data %>%
  mutate(Veg = case_when(
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
  relocate(Veg, .before = 1) %>%
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

# running the nmds; putting the data in a dissimilarity matrix using bray
dist_matrix_totalbiomass <- vegdist(total_nmds_data[, -(1:3)], method = "bray")
nmds_totalbiomass <- metaMDS(dist_matrix_totalbiomass, k=2, trymax=1000, distance="bray", wascores=TRUE)
nmds_scores_totalbiomass <- as.data.frame(scores(nmds_totalbiomass))
nmds_scores_totalbiomass$Wetland <- total_nmds_data$Wetland
nmds_scores_totalbiomass$HLength <- total_nmds_data$HLength
nmds_scores_totalbiomass$Veg <- total_nmds_data$Veg

# plotting stress and checking if it's within reasonable limit
stressplot(nmds_totalbiomass)
nmds_totalbiomass$stress
# stress: 0.06816263

# pulling out significant traits that contribute to shaping the NMDS
envfit_scores_totalbiomass <- envfit(nmds_totalbiomass, total_nmds_data[, -(1:3)], permutations = 999)
taxa_scores_totalbiomass <- as.data.frame(scores(envfit_scores_totalbiomass, display = "vectors"))
taxa_scores_totalbiomass <- cbind(taxa_scores_totalbiomass, Taxon = rownames(taxa_scores_totalbiomass), pval = envfit_scores_totalbiomass$vectors$pvals, r2 = envfit_scores_totalbiomass$vectors$r)
sig_taxa_scores_totalbiomass <- subset(taxa_scores_totalbiomass, pval <= 0.05)

## Plotting

# loading in packages for plotting
library(ggplot2)
library(ggrepel)
# (ggrepel is used to keep labels from overlapping)

# reordering the hydroperiod lengths
nmds_scores_totalbiomass$HLength <- factor(nmds_scores_totalbiomass$HLength,
                                           levels = c("Short", "Intermediate", "Long"))
# function for making the hulls on the NMDS plots
hull_points <- function(df) {
  df[chull(df$NMDS1, df$NMDS2), ]
}
# Making Figure 4A

# making centroid for Figure 4A
VegCentroid <- 
  nmds_scores_totalbiomass %>% 
  group_by(Veg) %>% 
  summarise(axis1 = mean(NMDS1),
            axis2 = mean(NMDS2)) %>% 
  ungroup()

# making hulls
hull_points <- function(data) {
  data[chull(data$NMDS1, data$NMDS2), ]
}

hull_totaldataA <- nmds_scores_totalbiomass %>%
  group_by(Veg) %>%
  summarise(
    hull = list(hull_points(cur_data()))
  ) %>%
  unnest(hull)

hull_totaldataA <- nmds_scores_totalbiomass %>%
  group_by(Veg) %>%
  summarise(hull = list(hull_points(cur_data_all()))) %>%
  unnest(hull, names_sep = "_")

# rename the columns to avoid confusion after unnesting
hull_totaldataA <- hull_totaldataA %>%
  rename(NMDS1 = hull_NMDS1, NMDS2 = hull_NMDS2)

# nmds plot based on site totals, hulls based on vegetation
Figure4A<- ggplot() +
  # adding hulls based on vegetation type
  geom_polygon(data = hull_totaldataA,
               aes(x = NMDS1, y = NMDS2, fill = as.factor(Veg), group = Veg),
               alpha = 0.6) +
  # adding centroid for each hull, using fill for color matching
  geom_point(data = VegCentroid, 
             aes(x = axis1, y = axis2, fill = as.factor(Veg)), 
             size = 5, shape = 21) +
  # adding color
  scale_fill_manual(values = c("Marsh" = "#79B4A9",
                               "Swamp" = "#22577A"), name = "Vegetation")+
  # using the nmds scores for each site as points
  geom_point(data = nmds_scores_totalbiomass,
             aes(x = NMDS1, y = NMDS2, shape= HLength),
             size = 3) +
  # adding the taxa as dots
  geom_point(data = sig_taxa_scores_totalbiomass,
             aes(x = NMDS1, y = NMDS2),
             color = "black",
             shape = 1,
             size = 2) +
  # adding taxa point labels
  geom_text_repel(data = sig_taxa_scores_totalbiomass,
                  aes(x = NMDS1, y = NMDS2, label = Taxon),
                  color = "black",
                  size = 4, 
                  max.overlaps = 15)+
  # adding labels to the rest of the plot
  labs(x = "NMDS1", 
       y = "NMDS2",
       fill = "Vegetation",
       shape = "Hydroperiod Length")+
  theme_minimal() +
  theme(legend.position = "right")
print(Figure4A)

# making Figure 4B

# calculating centroids for adding to Figure 4B
HLengthCentroid <- 
  nmds_scores_totalbiomass %>% 
  group_by(HLength) %>% 
  summarise(axis1 = mean(NMDS1),
            axis2 = mean(NMDS2)) %>% 
  ungroup()

# making hulls for Figure 4B
hull_totaldataB <- nmds_scores_totalbiomass %>%
  group_by(HLength) %>%
  summarise(hull = list(hull_points(cur_data_all()))) %>%
  unnest(hull, names_sep = "_")

# renaming the columns to avoid confusion after unnesting (just in case)
hull_totaldataB <- hull_totaldataB %>%
  rename(NMDS1 = hull_NMDS1, NMDS2 = hull_NMDS2)

# Figure 4B
Figure4B <- ggplot() +
  # adding hulls based on hydroperiod length
  geom_polygon(data = hull_totaldataB,
               aes(x = NMDS1, y = NMDS2, fill = as.factor(HLength), group = HLength),
               alpha = 0.7) +
  # using the nmds scores for each site as points
  geom_point(data = nmds_scores_totalbiomass,
             aes(x = NMDS1, y = NMDS2, shape = Veg),
             size = 3) +
  # adding centroid for each hull
  geom_point(data = HLengthCentroid, 
             aes(x = axis1, y = axis2, fill = as.factor(HLength)), 
             size = 5, shape = 21) +
  # single color scale for both hulls and centroids
  scale_fill_manual(values = c("Short" = "#DBDEE6", 
                               "Intermediate" = "#849dc1",
                               "Long" = "#354C82"), 
                    name = "Hydroperiod Length") +
  # adding the taxa as dots
  geom_point(data = sig_taxa_scores_totalbiomass,
             aes(x = NMDS1, y = NMDS2),
             color = "black",
             shape = 1, 
             size = 2) +
  # adding labels to the taxa points
  geom_text_repel(data = sig_taxa_scores_totalbiomass,
                  aes(x = NMDS1, y = NMDS2, label = Taxon),
                  color = "black",
                  size = 4, 
                  max.overlaps = 15) +
  # adding labels for the rest of the plot
  labs(x = "NMDS1", 
       y = "NMDS2",
       shape = "Vegetation",
       fill = "Hydroperiod Length") +
  theme_minimal() +
  theme(legend.position = "right")
print(Figure4B)

# putting 4A and 4B together
library(patchwork)
Figure4A+ Figure4B

## Analysis for Figure 4A and 4B

# preparing matrix for analysis
library(dplyr)
biomass_matrix <- total_nmds_data %>%
  dplyr::select(-c(HLength, Wetland, Veg))

# type III permanova to account for any interactions
adonis2(biomass_matrix ~ HLength + Veg, data = total_nmds_data,
        permutations = 9999, method = "bray", by = "margin")

### Emergence trait-based community composition NMDS

# Packages
library(vegan) 
library(ggplot2) 
library(tidyverse) 
library(dplyr)
library(goeveg) 

# total biomass weighted community composition  per site
# reading in the data
tot_trait_biomass <- read.csv("EmergeBiomassTotaledwTraits.csv")

# averaging replicates
total_FFG <- tot_trait_biomass %>%
  group_by(Wetland, DateCollected, Type, HLength, FFG) %>% 
  summarize(avg_biomass = mean(tot_biomass, na.rm = TRUE),
            .groups = "drop")
# totaling biomass per site
total_biomass <- total_FFG %>%
  group_by(Wetland, FFG)%>%
  summarize(total_biomass = sum(avg_biomass, na.rm = TRUE),
            .groups = "drop")

# reshaping to wide format
total_nmds_data <- total_biomass %>%
  group_by(Wetland) %>% 
  pivot_wider(
    names_from = FFG,
    values_from = total_biomass,
    values_fill = 0
  ) %>%
  ungroup()

# adding back in metadata
total_nmds_data <- total_nmds_data %>%
  mutate(Veg = case_when(
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
  relocate(Veg, .before = 1) %>%
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
  relocate(HLength, .before = 1)%>%
  mutate(Days_inund = case_when(
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
  relocate(Days_inund, .before = 1)


# running the nmds; putting things in a dissimilarity matrix using bray
dist_matrix_totalbiomass <- vegdist(total_nmds_data[, -(1:4)], method = "bray")
nmds_totalbiomass <- metaMDS(dist_matrix_totalbiomass, k=2, trymax=1000, distance="bray", wascores=TRUE)
nmds_scores_totalbiomass <- as.data.frame(scores(nmds_totalbiomass))
nmds_scores_totalbiomass$Wetland <- total_nmds_data$Wetland
nmds_scores_totalbiomass$HLength <- total_nmds_data$HLength
nmds_scores_totalbiomass$Veg <- total_nmds_data$Veg
nmds_scores_totalbiomass$Days_inund <- total_nmds_data$Days_inund

# plotting stress and checking if it's within reasonable limit
stressplot(nmds_totalbiomass)
nmds_totalbiomass$stress
# stress: 0.08902729

# pulling out important taxa
envfit_scores_totalbiomass <- envfit(nmds_totalbiomass, total_nmds_data[, -(1:4)], permutations = 999)
FFG_scores_totalbiomass <- as.data.frame(scores(envfit_scores_totalbiomass, display = "vectors"))
FFG_scores_totalbiomass <- cbind(FFG_scores_totalbiomass, Taxon = rownames(FFG_scores_totalbiomass), pval = envfit_scores_totalbiomass$vectors$pvals, r2 = envfit_scores_totalbiomass$vectors$r)
sig_FFG_scores_totalbiomass <- subset(FFG_scores_totalbiomass, pval <= 0.05)

## Plotting figures 4C and 4D
# loading in packages for plotting
library(ggplot2)
library(ggrepel)
# (ggrepel is used to keep labels from overlapping)

# reordering the hydroperiod lengths
nmds_scores_totalbiomass$HLength <- factor(nmds_scores_totalbiomass$HLength, 
                                           levels = c("Short", "Intermediate", "Long"))

# function for making the hulls
hull_points <- function(df) {
  df[chull(df$NMDS1, df$NMDS2), ]
}

# making hulls

# get centroid 
HLengthCentroid <- 
  nmds_scores_totalbiomass %>% 
  group_by(HLength) %>% 
  summarise(axis1 = mean(NMDS1),
            axis2 = mean(NMDS2)) %>% 
  ungroup()


hull_totaldataD <- nmds_scores_totalbiomass %>%
  group_by(HLength) %>%
  summarise(hull = list(hull_points(cur_data_all()))) %>%
  unnest(hull, names_sep = "_")

# rename the columns to avoid confusion after unnesting
hull_totaldataD <- hull_totaldataD %>%
  rename(NMDS1 = hull_NMDS1, NMDS2 = hull_NMDS2)

# nmds plot based on site totals, hulls based on hydroperiod length
Figure4D <- ggplot() +
  # adding hulls based on hydroperiod length
  geom_polygon(data = hull_totaldataD,
               aes(x = NMDS1, y = NMDS2, fill = as.factor(HLength), group = HLength),
               alpha = 0.7) +
  # using the nmds scores for each site as points
  geom_point(data = nmds_scores_totalbiomass,
             aes(x = NMDS1, y = NMDS2, shape = Veg),
             size = 3) +
  # adding centroid for each hull, using fill for color matching
  geom_point(data = HLengthCentroid, 
             aes(x = axis1, y = axis2, fill = as.factor(HLength)), 
             size = 5, shape = 21) +
  scale_fill_manual(values = c("Short" = "#DBDEE6", "Intermediate" = "#849dc1",
                               "Long" = "#354C82"), name = "Hydroperiod Length")+
  # adding the taxa as dots
  geom_point(data = sig_FFG_scores_totalbiomass,
             aes(x = NMDS1, y = NMDS2),
             color = "black",
             shape = 1, 
             size = 2) +
  # adding labels to the taxa points
  geom_text_repel(data = sig_FFG_scores_totalbiomass,
                  aes(x = NMDS1, y = NMDS2, label = Taxon),
                  color = "black",
                  size = 4, 
                  max.overlaps = 15)+
  labs(x = "NMDS1", 
       y = "NMDS2",
       shape = "Vegetation",
       fill = "Hydroperiod Length") +
  theme_minimal() +
  theme(legend.position = "right")
print(Figure4D)

# making Figure 4C

# get centroid 
VegCentroid <- 
  nmds_scores_totalbiomass %>% 
  group_by(Veg) %>% 
  summarise(axis1 = mean(NMDS1),
            axis2 = mean(NMDS2)) %>% 
  ungroup()

# making hulls
hull_totaldataC <- nmds_scores_totalbiomass %>%
  group_by(Veg) %>%
  summarise(
    hull = list(hull_points(cur_data()))
  ) %>%
  unnest(hull)


hull_totaldataC <- nmds_scores_totalbiomass %>%
  group_by(Veg) %>%
  summarise(hull = list(hull_points(cur_data_all()))) %>%
  unnest(hull, names_sep = "_")

# rename the columns to avoid confusion after unnesting
hull_totaldataC <- hull_totaldataC %>%
  rename(NMDS1 = hull_NMDS1, NMDS2 = hull_NMDS2)

# nmds plot based on site totals, hulls based on vegetation
Figure4C <- ggplot() +
  # adding hulls based on vegetation type
  geom_polygon(data = hull_totaldataC,
               aes(x = NMDS1, y = NMDS2, fill = as.factor(Veg), group = Veg),
               alpha = 0.6) +
  # using the nmds scores for each site as points
  geom_point(data = nmds_scores_totalbiomass,
             aes(x = NMDS1, y = NMDS2, shape= HLength),
             size = 3) +
  # adding centroid for each hull, using fill for color matching
  geom_point(data = VegCentroid, 
             aes(x = axis1, y = axis2, fill = as.factor(Veg)), 
             size = 5, shape = 21) +
  scale_fill_manual(values = c("Marsh" = "#79B4A9",
                               "Swamp" = "#22577A"))+
  # adding the taxa as dots
  geom_point(data = sig_FFG_scores_totalbiomass,
             aes(x = NMDS1, y = NMDS2),
             color = "black",
             shape = 1,
             size = 2) +
  # adding taxa point labels
  geom_text_repel(data = sig_FFG_scores_totalbiomass,
                  aes(x = NMDS1, y = NMDS2, label = Taxon),
                  color = "black",
                  size = 4, 
                  max.overlaps = 15)+
  labs(x = "NMDS1", 
       y = "NMDS2",
       fill = "Vegetation",
       shape = "Hydroperiod Length")+
  theme_minimal() +
  theme(legend.position = "right")
print(Figure4C)

# plotting the plots side by side
library(patchwork)
Figure4C + Figure4D

# creating matrix for permanova analysis
library(dplyr)
biomass_matrix <- total_nmds_data %>%
  dplyr::select(-c(HLength, Wetland, Veg, Days_inund))

# permanova for differences based on vegetation and hydrology
# type III permanova 
adonis2(biomass_matrix ~ HLength + Veg, data = total_nmds_data,
        permutations = 9999, method = "bray", by = "margin")

### Emergence community across time NMDS
# Packages
library(vegan) 
library(ggplot2) 
library(tidyverse) 
library(dplyr)
library(goeveg) 

# first step is averaging the replicates in each site for each sampling date,
# to end with average taxa biomass per sampling date that can then be totaled 
# to get total taxa biomass in each site. This can be used for total biomass 
# and for how the community changes across the hydroperiod 

# reading in the data for figure 5B
taxa_biomass <- read.csv("EmergeBiomassTotaledwTraits.csv")
taxa_biomass <- taxa_biomass %>%
  rename(Biomass= tot_biomass)

# fixing dates
taxa_biomass <- taxa_biomass %>%
  mutate(Date_collected = as.Date(DateCollected, format = "%m/%d/%Y")) %>%
  mutate(Date_collected = case_when(
    Date_collected == as.Date("2023-02-21") ~ "Feb",
    Date_collected == as.Date("2024-02-21") ~ "Feb",
    Date_collected == as.Date("2023-03-06") ~ "Mar",
    Date_collected == as.Date("2023-03-05") ~ "Mar",
    Date_collected == as.Date("2023-03-21") ~ "Mar",
    Date_collected == as.Date("2023-04-05") ~ "Apr",
    Date_collected == as.Date("2023-04-20") ~ "Apr",
    Date_collected == as.Date("2023-05-03") ~ "May",
    Date_collected == as.Date("2024-05-03") ~ "May",
    Date_collected == as.Date("2023-05-12") ~ "May",
    Date_collected == as.Date("2023-05-17") ~ "May",
    Date_collected == as.Date("2023-05-31") ~ "May",
    Date_collected == as.Date("2023-06-14") ~ "Jun",
    Date_collected == as.Date("2023-06-15") ~ "Jun",
    Date_collected == as.Date("2023-06-23") ~ "Jun",
    Date_collected == as.Date("2023-06-22") ~ "Jun",
    Date_collected == as.Date("2023-07-05") ~ "Jul",
    Date_collected == as.Date("2023-07-06") ~ "Jul",
    Date_collected == as.Date("2023-07-19") ~ "Jul",
    Date_collected == as.Date("2024-07-19") ~ "Jul",
    Date_collected == as.Date("2023-08-01") ~ "Aug",
    TRUE ~ as.character(Date_collected)
  )) %>%
  relocate(Date_collected, .before = 1)

# averaging replicates, keeping in other variables
taxa_biomass <- taxa_biomass %>%
  group_by(Date_collected, HLength, Type, FFG) %>% 
  summarize(avg_biomass = mean(Biomass, na.rm = TRUE),
            .groups = "drop") 

# averaging replicates, keeping in other variables
swamp_data <- taxa_biomass %>%
  group_by(Date_collected, HLength, Type, FFG) %>%
  filter(Type == "Swamp") %>%
  ungroup()

# reshaping the data to wide format
nmds_data_s <- swamp_data%>%
  group_by(Date_collected, HLength) %>%  
  pivot_wider(
    names_from = FFG,
    values_from = avg_biomass,
    values_fill = 0
  ) %>%
  ungroup()


# running the actual nmds; putting things in a dissimilarity matrix using bray first
dist_matrix_biomass_s <- vegdist(nmds_data_s[, -(1:3)], method = "bray")
nmds_biomass_s <- metaMDS(dist_matrix_biomass_s, k=2, trymax=1000, distance="bray", wascores=TRUE)
nmds_scores_biomass_s <- as.data.frame(scores(nmds_biomass_s))
nmds_scores_biomass_s$Date_collected <- nmds_data_s$Date_collected
nmds_scores_biomass_s$HLength <- nmds_data_s$HLength
nmds_scores_biomass_s$Type <- nmds_data_s$Type

# assessing stress
stressplot(nmds_biomass_s)
nmds_biomass_s$stress
# 0.1389044

# pulling out significant trait drivers
envfit_scores_biomass_s <- envfit(nmds_biomass_s, nmds_data_s[, -(1:3)], permutations = 999)
taxa_scores_biomass_s <- as.data.frame(scores(envfit_scores_biomass_s, display = "vectors"))
taxa_scores_biomass_s <- cbind(taxa_scores_biomass_s, Taxon = rownames(taxa_scores_biomass_s), pval = envfit_scores_biomass_s$vectors$pvals, r2 = envfit_scores_biomass_s$vectors$r)
sig_taxa_scores_biomass_s <- subset(taxa_scores_biomass_s, pval <= 0.05)

# trying to plot everything!
nmds_scores_biomass_s$HLength <- factor(nmds_scores_biomass_s$HLength, 
                                        levels = c("Short", "Intermediate", "Long"))

# function for making the hulls
hull_points <- function(df) {
  df[chull(df$NMDS1, df$NMDS2), ]
}

# making hulls and centroids
DateSCentroid <- 
  nmds_scores_biomass_s %>% 
  group_by(Date_collected) %>% 
  summarise(axis1 = mean(NMDS1),
            axis2 = mean(NMDS2)) %>% 
  ungroup()

hull_data_s <- nmds_scores_biomass_s %>%
  group_by(Date_collected) %>%
  summarise(hull = list(hull_points(cur_data_all()))) %>%
  unnest(hull, names_sep = "_")

# renaming the columns to avoid confusion after unnesting
hull_data_s <- hull_data_s %>%
  rename(NMDS1 = hull_NMDS1, NMDS2 = hull_NMDS2)

# reordering months
month_order <- c("Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug")
hull_data_s$Date_collected <- factor(hull_data_s$Date_collected, levels = month_order)


# nmds with collection dates, and hulls as hydroperiod lengths
Figure5B <- ggplot() +
  # hulls based on date 
  geom_polygon(data = hull_data_s,
               aes(x = NMDS1, y = NMDS2, fill = as.factor(Date_collected), group = Date_collected),
               alpha = 0.45) +
  # adding centroid for each hull, using fill for color matching
  geom_point(data = DateSCentroid, 
             aes(x = axis1, y = axis2, fill = as.factor(Date_collected)), 
             size = 3, shape = 21) +
  scale_fill_manual(
    values = c(
      "Feb" = "#4575b4",
      "Mar" = "#74add9",
      "Apr" = "#3BA6A1",
      "May" = "#fee090",
      "Jun" = "#fdae61",
      "Jul" = "#f46d43",
      "Aug" = "#d73027"),
    name = "Month Collected"
  ) +
  # using the nmds scores for each site as points
  geom_point(data = nmds_scores_biomass_s,
             aes(x = NMDS1, y = NMDS2, shape = HLength),
             size = 2.5) +
  # adding the taxa as dots
  geom_point(data = sig_taxa_scores_biomass_s,
             aes(x = NMDS1, y = NMDS2),
             color = "black",
             shape = 1, 
             size = 2) +
  # adding labels to the taxa points
  geom_text_repel(data = sig_taxa_scores_biomass_s,
                  aes(x = NMDS1, y = NMDS2, label = Taxon),
                  color = "black",
                  size = 4, 
                  max.overlaps = 15)+
  # adding in plot labels
  labs(x = "NMDS1", 
       y = "NMDS2",
       title = "Swamp",
       shape = "Hydroperiod Length",
       fill = "Month Collected") +
  theme_minimal() +
  theme(legend.position = "right")
print(Figure5B)

# PERMANOVA
permanova_s <- adonis2(dist_matrix_biomass_s ~ HLength + Date_collected, data = nmds_data_s,
                       permutations = 9999, method = "bray", by = "margin")
permanova_s

## Figure 5A
# reading in the data for figure 5A
taxa_biomass <- read.csv("EmergeBiomassTotaledwTraits.csv")
taxa_biomass <- taxa_biomass %>%
  rename(Biomass= tot_biomass)

# simplifying dates 
taxa_biomass <- taxa_biomass %>%
  mutate(Date_collected = as.Date(DateCollected, format = "%m/%d/%Y")) %>%
  mutate(Date_collected = case_when(
    Date_collected == as.Date("2023-02-21") ~ "Feb",
    Date_collected == as.Date("2024-02-21") ~ "Feb",
    Date_collected == as.Date("2023-03-06") ~ "Mar",
    Date_collected == as.Date("2023-03-05") ~ "Mar",
    Date_collected == as.Date("2023-03-21") ~ "Mar",
    Date_collected == as.Date("2023-04-05") ~ "Apr",
    Date_collected == as.Date("2023-04-20") ~ "Apr",
    Date_collected == as.Date("2023-05-03") ~ "May",
    Date_collected == as.Date("2024-05-03") ~ "May",
    Date_collected == as.Date("2023-05-12") ~ "May",
    Date_collected == as.Date("2023-05-17") ~ "May",
    Date_collected == as.Date("2023-05-31") ~ "May",
    Date_collected == as.Date("2023-06-14") ~ "Jun",
    Date_collected == as.Date("2023-06-15") ~ "Jun",
    Date_collected == as.Date("2023-06-23") ~ "Jun",
    Date_collected == as.Date("2023-06-22") ~ "Jun",
    Date_collected == as.Date("2023-07-05") ~ "Jul",
    Date_collected == as.Date("2023-07-06") ~ "Jul",
    Date_collected == as.Date("2023-07-19") ~ "Jul",
    Date_collected == as.Date("2024-07-19") ~ "Jul",
    Date_collected == as.Date("2023-08-01") ~ "Aug",
    TRUE ~ as.character(Date_collected)
  )) %>%
  relocate(Date_collected, .before = 1)

# averaging replicates, keeping in other variables
taxa_biomass <- taxa_biomass %>%
  group_by(Date_collected, HLength, Type, FFG) %>% 
  summarize(avg_biomass = mean(Biomass, na.rm = TRUE),
            .groups = "drop") 

# averaging replicates, keeping in other variables
marsh_data <- taxa_biomass %>%
  group_by(Date_collected, HLength, Type, FFG) %>%
  filter(Type == "Marsh") %>%
  ungroup()


# reshaping the data to wide format
nmds_data_m <- marsh_data%>%
  group_by(Date_collected, HLength) %>%  
  pivot_wider(
    names_from = FFG,
    values_from = avg_biomass,
    values_fill = 0
  ) %>%
  ungroup()


# running the actual nmds; putting things in a dissimilarity matrix using bray first
dist_matrix_biomass_m <- vegdist(nmds_data_m[, -(1:3)], method = "bray")
nmds_biomass_m <- metaMDS(dist_matrix_biomass_m, k=2, trymax=1000, distance="bray", wascores=TRUE)
nmds_scores_biomass_m <- as.data.frame(scores(nmds_biomass_m))
nmds_scores_biomass_m$Date_collected <- nmds_data_m$Date_collected
nmds_scores_biomass_m$HLength <- nmds_data_m$HLength
nmds_scores_biomass_m$Type <- nmds_data_m$Type

# assessing stress
stressplot(nmds_biomass_m)
nmds_biomass_m$stress
# 0.1222988

# pulling out significant trait drivers
envfit_scores_biomass_m <- envfit(nmds_biomass_m, nmds_data_m[, -(1:3)], permutations = 999)
taxa_scores_biomass_m <- as.data.frame(scores(envfit_scores_biomass_m, display = "vectors"))
taxa_scores_biomass_m <- cbind(taxa_scores_biomass_m, Taxon = rownames(taxa_scores_biomass_m), pval = envfit_scores_biomass_m$vectors$pvals, r2 = envfit_scores_biomass_m$vectors$r)
sig_taxa_scores_biomass_m <- subset(taxa_scores_biomass_m, pval <= 0.05)


# plotting figure 5A
nmds_scores_biomass_m$HLength <- factor(nmds_scores_biomass_m$HLength, 
                                        levels = c("Short", "Intermediate", "Long"))
# function for making the hulls
hull_points <- function(df) {
  df[chull(df$NMDS1, df$NMDS2), ]
}

# making hulls and centroid
DateMCentroid <- 
  nmds_scores_biomass_m %>% 
  group_by(Date_collected) %>% 
  summarise(axis1 = mean(NMDS1),
            axis2 = mean(NMDS2)) %>% 
  ungroup()

hull_data_m <- nmds_scores_biomass_m %>%
  group_by(Date_collected) %>%
  summarise(hull = list(hull_points(cur_data_all()))) %>%
  unnest(hull, names_sep = "_")

# renaming the columns to avoid confusion after unnesting
hull_data_m <- hull_data_m %>%
  rename(NMDS1 = hull_NMDS1, NMDS2 = hull_NMDS2)

# reordering months
month_order <- c("Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug")

hull_data_m$Date_collected <- factor(hull_data_m$Date_collected, levels = month_order)


# nmds with collection dates, and hulls as hydroperiod lengths
Figure5A <- ggplot() +
  # hulls based on hydroperiod length
  geom_polygon(data = hull_data_m,
               aes(x = NMDS1, y = NMDS2, fill = as.factor(Date_collected), group = Date_collected),
               alpha = 0.45) +
  # adding centroid for each hull, using fill for color matching
  geom_point(data = DateMCentroid, 
             aes(x = axis1, y = axis2, fill = as.factor(Date_collected)), 
             size = 3, shape = 21) +
  scale_fill_manual(
    values = c(
      "Feb" = "#4575b4",
      "Mar" = "#74add9",
      "Apr" = "#3BA6A1",
      "May" = "#fee090",
      "Jun" = "#fdae61",
      "Jul" = "#f46d43",
      "Aug" = "#d73027"),
    name = "Month Collected"
  ) +
  # using the nmds scores for each site as points
  geom_point(data = nmds_scores_biomass_m,
             aes(x = NMDS1, y = NMDS2, shape = HLength),
             size = 2.5) +
  # adding the taxa as dots
  geom_point(data = sig_taxa_scores_biomass_m,
             aes(x = NMDS1, y = NMDS2),
             color = "black",
             shape = 1, 
             size = 2) +
  # adding labels to the taxa points
  geom_text_repel(data = sig_taxa_scores_biomass_m,
                  aes(x = NMDS1, y = NMDS2, label = Taxon),
                  color = "black",
                  size = 4, 
                  max.overlaps = 15)+
  
  # adding in plot labels
  labs(x = "NMDS1", 
       y = "NMDS2",
       title = "Marsh",
       shape = "Hydroperiod Length",
       fill = "Month Collected") +
  theme_minimal() +
  theme(legend.position = "right")
print(Figure5A)

# combining the two plots
library(patchwork)
Figure5A + Figure5B

# running the permanova:
permanova_m <-  adonis2(dist_matrix_biomass_m ~ HLength + Date_collected, data = nmds_data_m,
                        permutations = 9999, method = "bray", by = "margin")
permanova_m

