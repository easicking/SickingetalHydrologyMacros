### Figure 3
# loading in packages
library(dplyr) 
library(lubridate) 
library(ggplot2)
library(patchwork)

# historical hydrology

hydro_data <- read.csv("HistoricalHydrology.csv")

hydro_data <- hydro_data[-c(5:11)]

# changing column names
colnames(hydro_data) <- c("wetland","veg","date", "stage (m)") 

# changing the date column from characters to dates
hydro_data$date <- mdy(hydro_data$date)

# separating the dates into new columns (month, day, year) and add new column with
# month names
hydro_data2 <- hydro_data
hydro_data2 <- hydro_data2 %>% mutate(month = month(hydro_data$date)) 
hydro_data2 <- hydro_data2%>% mutate(month2 = month.abb[month])
hydro_data2 <- hydro_data2 %>% mutate(day = day(hydro_data$date)) 
hydro_data2 <- hydro_data2 %>% mutate(year = year(hydro_data$date))

# need to remove observations where the year is 2024
hydro_data2 <- hydro_data %>%
  mutate(month = month(date),
         month2 = month.abb[month],
         day = day(date),
         year = year(date)) %>%
  filter(year != 2024) 

# summarizing the data (mean wetland depth each month)
hydro_data_final <- hydro_data2 %>% 
  group_by(wetland = wetland,month = month2, veg = veg) %>% 
  summarize(stage=mean(`stage (m)`))


# contemporary hydrology
contemp_data <- read.csv("ContempHydrology.csv")

contemp_data <- contemp_data[-c(5:11)]

# changing column names
colnames(contemp_data) <- c("wetland","veg","date","stage (m)")

contemp_data$date <- mdy(contemp_data$date)

# separating the dates into new columns (month, day, year) and add new column with
# month names
contemp_data2 <- contemp_data
contemp_data2 <- contemp_data2 %>% mutate(month = month(contemp_data$date)) 
contemp_data2 <- contemp_data2%>% mutate(month2 = month.abb[month])
contemp_data2 <- contemp_data2 %>% mutate(day = day(contemp_data$date)) 
contemp_data2 <- contemp_data2 %>% mutate(year = year(contemp_data$date))

# summarizing the data (mean wetland depth each month)
contemp_data_final <- contemp_data2 %>% 
  group_by(wetland = wetland,month = month2, veg = veg) %>% 
  summarize(stage=mean(`stage (m)`))

# combining the plots
hydro_data_final$contempstage <-contemp_data_final$stage
combinedhydro<- hydro_data_final
combinedhydroswamp <- hydro_data_final %>%
  filter(veg == "Swamp")

combinedhydromarsh <- hydro_data_final %>%
  filter(veg == "Marsh")

combinedhydro <- combinedhydro%>%
  mutate(HLength = case_when(
    wetland == "W15" ~ "Short",
    wetland == "W21" ~ "Intermediate",
    wetland == "W53" ~ "Intermediate",
    wetland == "W46" ~ "Long",
    wetland == "W37" ~ "Short",
    wetland == "W42" ~ "Long",
    wetland == "W11" ~ "Long",
    wetland == "W68" ~ "Intermediate",
    wetland == "W58" ~ "Long",
    wetland == "W32" ~ "Short",
    wetland == "W52" ~ "Intermediate",
  )) %>%
  relocate(HLength, .before = 1)

combinedhydroswamp$wetland = factor(combinedhydroswamp$wetland, levels = c("W15", "W37", "W21", 
                                                                           "W53","W46","W42",
                                                                           "W32", "W68","W52", 
                                                                           "W11", "W58"))
combinedhydromarsh$wetland = factor(combinedhydromarsh$wetland, levels = c("W15", "W37", "W21", 
                                                                           "W53","W46","W42",
                                                                           "W32", "W68","W52", 
                                                                           "W11", "W58"))

combinedhydroswamp$month <- factor(combinedhydroswamp$month, levels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", 
                                                                      "Apr", "May", "Jun", "Jul", "Aug", "Sep"))
combinedhydromarsh$month <- factor(combinedhydromarsh$month, levels=c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", 
                                                                      "Apr", "May", "Jun", "Jul", "Aug", "Sep"))

swamp <- ggplot(combinedhydroswamp, aes(x = month, group = veg)) +
  geom_line(aes(y = stage, color = "Average Hydroperiod (1998-2023)", linetype = "Average Hydroperiod (1998-2023)"), linewidth = 1) +
  geom_line(aes(y = contempstage, color = "2023 Hydroperiod", linetype = "2023 Hydroperiod"), linewidth = 1) +
  facet_wrap(~ wetland, scales = "free_x", nrow = 1) +
  xlab('Month') +
  ylab('Stage (m)') +
  labs(title = "Swamps")+
  scale_color_manual(name = "Hydrology",
                     values = c("Average Hydroperiod (1998-2023)" = "#22577A", "2023 Hydroperiod" = "#22577A")) +
  scale_linetype_manual(name = "Hydrology",
                        values = c("Average Hydroperiod (1998-2023)" = "solid", "2023 Hydroperiod" = "dotted")) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "white", fill = NA),
        axis.line = element_line(color = "grey"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 12, margin = margin(2, 0, 2, 0))) +
  theme(text = element_text(size = 14)) +
  scale_x_discrete(labels = c("Oct" = "O", "Nov" = "N", "Dec" = "D", "Jan" = "J", 
                              "Feb" = "F", "Mar" = "M", "Apr" = "A", "May" = "M",
                              "Jun" = "J", "Jul" = "J", "Aug" = "A", "Sep" = "S"))
swamp

marsh <- ggplot(combinedhydromarsh, aes(x = month, group = veg)) +
  geom_line(aes(y = stage, color = "Average Hydroperiod (1998-2023)", linetype = "Average Hydroperiod (1998-2023)"), linewidth = 1) +
  geom_line(aes(y = contempstage, color = "2023 Hydroperiod", linetype = "2023 Hydroperiod"), linewidth = 1) +
  facet_wrap(~ wetland, scales = "free_x", nrow = 1) +
  xlab('Month') +
  ylab('Stage (m)') +
  labs(title = "Marshes")+
  scale_color_manual(name = "Hydrology",
                     values = c("Average Hydroperiod (1998-2023)" = "#5aa093", "2023 Hydroperiod" = "#5aa093")) +
  scale_linetype_manual(name = "Hydrology",
                        values = c("Average Hydroperiod (1998-2023)" = "solid", "2023 Hydroperiod" = "dotted")) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "white", fill = NA),
        axis.line = element_line(color = "grey"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 12, margin = margin(2, 0, 2, 0))) +
  theme(text = element_text(size = 14)) +
  scale_x_discrete(labels = c("Oct" = "O", "Nov" = "N", "Dec" = "D", "Jan" = "J", 
                              "Feb" = "F", "Mar" = "M", "Apr" = "A", "May" = "M",
                              "Jun" = "J", "Jul" = "J", "Aug" = "A", "Sep" = "S"))
marsh

# putting the marsh and swamp plots together
Figure3<- marsh/swamp
ggsave("Figure3Sickingetal.png", plot = Figure3, dpi = 700, width = 15, height = 6, units = "in")



hydrology<- wilcox.test(combinedhydro$stage, combinedhydro$contempstage)
hydrology

