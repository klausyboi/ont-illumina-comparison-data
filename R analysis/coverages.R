rm(list=ls())
setwd("C:/Users/thorp/Downloads")

library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)

df <- read.csv("coveragesoutput.csv")
df <- distinct(df)
ontdf <- read.csv("combined_output.csv")


average_df_ONT <- ontdf %>%
  group_by(Position, Gene) %>%
  summarise(avg_coverage = mean(Coverages, na.rm = TRUE))


df_ILL <- df %>%
  filter(str_detect(SampleName, "ILL"))

average_df_ILL <- df_ILL %>%
  group_by(Position,Gene) %>%
  summarise(avg_coverage = mean(Coverages, na.rm = TRUE))

average_df_ILL$avg_coverage <- average_df_ILL$avg_coverage*100
average_df_ILL_Mbp <- average_df_ILL %>% mutate(Position = Position / 1000000)
average_df_ONT_Mbp <- average_df_ONT %>% mutate(Position = Position / 1000000)

library(ggrepel)

gg1 <- ggplot(average_df_ONT_Mbp, aes(x = Position, y = avg_coverage, color = "ONT")) +
  geom_line() +
  geom_text_repel(data = subset(average_df_ONT_Mbp, avg_coverage < 60), 
                  aes(label = Gene),color="#5A5A5A" ,hjust = -0.1, vjust=0, size = 3, 
                  force = 10, 
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"), max.overlaps = Inf) +
  labs(y = "Coverages %", x = "Position (Mbp)") + 
  scale_color_manual(values = c("ONT" = "#FFA500")) +
  scale_x_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 1,vjust=-30),
        plot.title.position = "plot",
        legend.position = "none")



gg2 <- ggplot(average_df_ILL_Mbp, aes(x = Position, y = avg_coverage, color = "Illumina")) +
  geom_line() +
  geom_text_repel(data = subset(average_df_ILL_Mbp, avg_coverage < 60), 
                  aes(label = Gene),color="#5A5A5A" ,hjust = -0.1, vjust=0, size = 3, 
                  force = 10, 
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"), max.overlaps = Inf) +
  labs(y = "Coverages %", x = "Position (Mbp)") + 
  scale_color_manual(values = c("Illumina" = "#6AAEFF")) +
  scale_x_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 1,vjust=-30),
        plot.title.position = "plot",
        legend.position = "none")




plot_grid(gg2, gg1, align = 'v', nrow = 2, axis = 'l')


ppe_df <- df %>%
  filter(PPE != ".")


ppe_ONT_df <- ontdf %>%
  filter(PPE != ".")

ppe_ILL_df <- ppe_df %>%
  filter(str_detect(SampleName, "ILL"))

ppe_ILL_df$SampleName <- as.factor(ppe_ILL_df$SampleName)
ppe_ONT_df$SampleName <- as.factor(ppe_ONT_df$SampleName)



average_ppe_ONT <- ppe_ONT_df %>%
  group_by(Position, Gene,PPE) %>%
  summarise(avg_coverage = mean(Coverages, na.rm = TRUE))


average_ppe_ILL <- ppe_ILL_df %>%
  group_by(Position, Gene,PPE) %>%
  summarise(avg_coverage = mean(Coverages, na.rm = TRUE))
average_ppe_ILL$avg_coverage = average_ppe_ILL$avg_coverage*100
average_ppe_ILL$source <- "Illumina"
average_ppe_ONT$source <- "ONT"

average_ppe_combined <- rbind(average_ppe_ILL, average_ppe_ONT)

findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

average_ppe_combined <- average_ppe_combined %>%
  group_by(source) %>%
  mutate(outlier = ifelse(findoutlier(avg_coverage), PPE, NA))

library(ggrepel)

plot1 <-ggplot(average_ppe_combined, aes(x = "", y = avg_coverage, fill = source)) +
  geom_boxplot(position = position_dodge(width = 1), width = 1) +
  ggrepel::geom_text_repel(
    aes(label = outlier),
    size = 2.3,
    nudge_x = 0.3,
    direction = "y",
    force = 0.5
  ) +
  facet_wrap(~ source, ncol = 2, scales = "free_x") +
  labs(y = "Coverage (%)") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank()
  ) +
  scale_fill_manual(values = c("Illumina" = "lightblue", "ONT" = "orange"))


drgenes <- read.csv("coveragesoutput.csv")

stuff <- read.table("new2.txt")
rv_numbers <- stuff$V4
rv_unique <- unique(rv_numbers[grep("^Rv", rv_numbers)])

drgenes_ILL <- df_ILL[df_ILL$Gene %in% rv_unique, ]
drgenes_ONT <- ontdf[ontdf$Gene %in% rv_unique, ]


drgenes_ILL$SampleName <- as.factor(drgenes_ILL$SampleName)
drgenes_ONT$SampleName <- as.factor(drgenes_ONT$SampleName)

ggplot(drgenes_ONT, aes(x=SampleName, y=Coverages)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank())
ggplot(drgenes_ILL, aes(x=SampleName, y=Coverages)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank())


average_dr_ONT <- drgenes_ONT %>%
  group_by(Position, Gene) %>%
  summarise(avg_coverage = mean(Coverages, na.rm = TRUE))


average_dr_ILL <- drgenes_ILL %>%
  group_by(Position,Gene) %>%
  summarise(avg_coverage = mean(Coverages, na.rm = TRUE))

average_dr_ILL$source <- "Illumina"
average_dr_ONT$source <- "ONT"

average_dr_ILL$avg_coverage <- average_dr_ILL$avg_coverage*100
average_dr_combined <- rbind(average_dr_ILL, average_dr_ONT)



average_dr_combined <- average_dr_combined %>%
  group_by(source) %>%
  mutate(outlier = ifelse(findoutlier(avg_coverage), Gene, NA))


plot2<- ggplot(average_dr_combined, aes(x = "", y = avg_coverage, fill = source)) +
  geom_boxplot(position = position_dodge(width = 1), width = 1) +
  ggrepel::geom_text_repel(
    aes(label = outlier),
    size = 2.3,
    nudge_x = 0.3,
    direction = "y",
    force = 0.5
  ) +
  facet_wrap(~ source, ncol = 2, scales = "free_x") +
  labs(y = "Coverage (%)") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank()
  ) +
  scale_fill_manual(values = c("Illumina" = "lightblue", "ONT" = "orange"))

combined_data <- rbind(average_dr_combined %>% mutate(data_set = "dr"), 
                       average_ppe_combined %>% mutate(data_set = "ppe"))

find_outliers <- function(data){
  Q1 <- quantile(data, 0.25)
  Q3 <- quantile(data, 0.75)
  IQR <- Q3 - Q1
  upper_limit <- Q3 + 1.5*IQR
  lower_limit <- Q1 - 1.5*IQR
  
  return(data > upper_limit | data < lower_limit)
}

combined_data <- combined_data %>%
  group_by(data_set, source) %>%
  mutate(outlier_flag = find_outliers(avg_coverage))

outliers_data <- combined_data %>%
  filter(outlier_flag)

ggplot(combined_data, aes(x = data_set, y = avg_coverage, fill = source)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  ggrepel::geom_text_repel(
    data = outliers_data,
    aes(label = outlier),
    size = 2.3,
    force = 10,
    max.iter = 10000,  
    position = position_dodge(width = 0.8)
  ) +
  labs(y = "Coverage (%)") +
  scale_fill_manual(values = c("Illumina" = "lightblue", "ONT" = "orange")) +
  theme(panel.background = element_blank())

gene_list <- unique(average_dr_combined$Gene)

allothers <- df %>% 
  filter(PPE == "." & !(Gene %in% gene_list))


allothers_ont <- ontdf %>% 
  filter(PPE == "." & !(Gene %in% gene_list))


allOther_df_ONT <- allothers_ont %>%
  filter(str_detect(SampleName, "ONT"))
allOther_df_ILL <- allothers %>%
  filter(str_detect(SampleName, "ILL"))

allOther_average_df_ONT <- allOther_df_ONT %>%
  group_by(Position, Gene) %>%
  summarise(avg_coverage = mean(Coverages, na.rm = TRUE))

allOther_average_df_ILL <- allOther_df_ILL %>%
  group_by(Position,Gene) %>%
  summarise(avg_coverage = mean(Coverages, na.rm = TRUE))

allOther_average_df_ILL$avg_coverage <- allOther_average_df_ILL$avg_coverage*100
allOther_average_df_ONT$avg_coverage <- allOther_average_df_ONT$avg_coverage/100
allOther_average_df_ILL_Mbp <- allOther_average_df_ILL %>% mutate(Position = Position / 1000000)
allOther_average_df_ONT_Mbp <- allOther_average_df_ONT %>% mutate(Position = Position / 1000000)

allOther_average_df_ILL_Mbp$source <- "Illumina"
allOther_average_df_ONT_Mbp$source <- "ONT"
allOther_combined <- rbind(allOther_average_df_ILL_Mbp, allOther_average_df_ONT_Mbp)
allOther_combined <- allOther_combined %>%
  group_by(source) %>%
  mutate(outlier = ifelse(findoutlier(avg_coverage), Gene, NA))

total_combined <- rbind(combined_data, allOther_combined %>% mutate(data_set = "All Other"))


total_combined <- total_combined %>%
  group_by(data_set, source) %>%
  mutate(outlier_flag = find_outliers(avg_coverage))

outliers_data <- total_combined %>%
  filter(outlier_flag)

ggplot(total_combined, aes(x = data_set, y = avg_coverage, fill = source, group = interaction(data_set, source))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  ggrepel::geom_text_repel(
    data = subset(total_combined, !is.na(outlier)),
    aes(label = outlier, group = interaction(data_set, source)),
    size = 2.3,
    force = 1,
    max.iter = 10000,
    max.overlaps = 30,
    position = position_dodge(width = 0.8)
  ) +
  labs(y = "Coverage (%)") +
  scale_fill_manual(values = c("Illumina" = "lightblue", "ONT" = "orange")) +
  theme(panel.background = element_blank())




