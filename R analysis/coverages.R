rm(list=ls())
setwd("C:/Users/thorp/Downloads")

library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggrepel)

df <- read.csv("coveragesoutput.csv")
df <- distinct(df)
ontdf <- read.csv("combined_output.csv")

drgenes <- read.csv("coveragesoutput.csv")

stuff <- read.table("new2.txt")
rv_numbers <- stuff$V4
rv_unique <- unique(rv_numbers[grep("^Rv", rv_numbers)])


average_df_ONT <- ontdf %>%
  group_by(Position, Gene,PPE) %>%
  summarise(avg_coverage = mean(Coverages, na.rm = TRUE))
df_ILL <- df %>%
  filter(str_detect(SampleName, "ILL"))

average_df_ILL <- df_ILL %>%
  group_by(Position,Gene,PPE) %>%
  summarise(avg_coverage = mean(Coverages, na.rm = TRUE))

drgenes_ILL <- df_ILL[df_ILL$Gene %in% rv_unique, ]
drgenes_ONT <- ontdf[ontdf$Gene %in% rv_unique, ]

average_df_ONT <- average_df_ONT %>%
  left_join(drgenes_ONT %>%
              select(Position, Gene) %>%
              distinct() %>%
              mutate(is_drgene = TRUE),
            by = c("Position", "Gene")) %>%
  replace_na(list(is_drgene = FALSE))

average_df_ILL <- average_df_ILL %>%
  left_join(drgenes_ILL %>%
              select(Position, Gene) %>%
              distinct() %>%
              mutate(is_drgene = TRUE),
            by = c("Position", "Gene")) %>%
  replace_na(list(is_drgene = FALSE))

average_df_ILL$avg_coverage <- average_df_ILL$avg_coverage*100
average_df_ILL_Mbp <- average_df_ILL %>% mutate(Position = Position / 1000000)
average_df_ONT_Mbp <- average_df_ONT %>% mutate(Position = Position / 1000000)

average_df_ONT_Mbp <- average_df_ONT_Mbp %>%
  filter(!str_starts(Gene,"EBG"))

gg1 <- ggplot(average_df_ONT_Mbp, aes(x = Position, y = avg_coverage, color = "ONT")) +
  geom_line() +
  
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
  
  labs(y = "Coverages %", x = "Position (Mbp)") + 
  scale_color_manual(values = c("Illumina" = "#6AAEFF")) +
  scale_x_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 1,vjust=-30),
        plot.title.position = "plot",
        legend.position = "none")




plot_grid(gg2, gg1, align = 'v', nrow = 2, axis = 'l')

merged_df <- merge(average_df_ILL_Mbp, average_df_ONT_Mbp, by = c("Position", "Gene","PPE","is_drgene") )
colnames(merged_df)<- c("Position","Gene","PPE","is_drgene","avg_coverage_ILL","avg_coverage_ONT")

merged_df$dominant_tech <- ifelse(merged_df$avg_coverage_ILL > merged_df$avg_coverage_ONT, "Illumina", "ONT")
label_df <- merged_df[(merged_df$avg_coverage_ILL < 60) | (merged_df$avg_coverage_ONT < 60), ]


merged_df <- merged_df %>%
  mutate(diff = abs(avg_coverage_ILL - avg_coverage_ONT))

gg3 <- ggplot(merged_df, aes(x = avg_coverage_ILL, y = avg_coverage_ONT, color = dominant_tech)) +
  geom_point(alpha = 0.6) + 
  geom_text_repel(data = label_df, aes(label = Gene), 
                  size = 3, 
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),max.overlaps = 100) +
  scale_color_manual(values = c("Illumina" = "orange", "ONT" = "blue")) +
  labs(y = "ONT Coverages %", x = "Illumina Coverages %") + 
  theme(panel.background = element_blank(),
        plot.title = element_text(hjust = 1, vjust = -30),
        plot.title.position = "plot",
        legend.position = "top")




print(gg3)
merged_df <- merged_df %>%
  mutate(color_criteria = case_when(
    is_drgene ~ "is_drgene",
    PPE != "." ~ "PPE",
    TRUE ~ "Other"
  ))
colnames(merged_df) <- c("Position","Gene","PPE","is_drgene","avg_coverage_ILL","avg_coverage_ONT","Lower_platform","Difference","Dominant_Tech")

gg4 <- ggplot(merged_df, aes(x=avg_coverage_ILL, y=avg_coverage_ONT)) +
  geom_point(aes(color = ifelse(is_drgene, "is_drgene", ifelse(PPE != ".", "PPE", "Other"))), alpha=0.6, size=2.5) +
  geom_smooth(method="lm", se=FALSE, color="grey50") +  
  scale_color_manual(
    values = c("is_drgene" = "red", "PPE" = "blue", "Other" = "black"),
    name = "Gene Category",
    breaks = c("is_drgene", "PPE", "Other"),
    labels = c("DR Gene", "PPE Gene", "Other")
  ) +
  labs(y="ONT Coverages %", 
       x="Illumina Coverages %", 
       title="Average coverage across both platforms for all genes") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 1),
    plot.title.position = "plot",
    legend.position = "top",
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

print(gg4)






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
allOther_average_df_ONT$avg_coverage <- allOther_average_df_ONT$avg_coverage*100
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
  labs(y = "Coverage (%)") +
  scale_fill_manual(values = c("Illumina" = "lightblue", "ONT" = "orange")) +
  theme(panel.background = element_blank())

str(average_df_ILL_Mbp)
str(average_df_ONT_Mbp)

spearmansCombined <- merge(average_df_ILL_Mbp, average_df_ONT_Mbp,by = c("Position","Gene","PPE","is_drgene"))
colnames(spearmansCombined)<- c("Position","Gene","PPE","is_drgene","Average_coverage_Illumina","Average_coverage_ONT")
rho <- cor(spearmansCombined$Average_coverage_Illumina,spearmansCombined$Average_coverage_ONT,method="spearman")
print(rho)

spearmansPPE <- merge(average_ppe_ILL,average_ppe_ONT,by=c("Position","Gene","PPE"))

colnames(spearmansPPE)[colnames(spearmansPPE) == "avg_coverage.x"] <- "Average_coverage_ILL"
colnames(spearmansPPE)[colnames(spearmansPPE) == "avg_coverage.y"] <- "Average_coverage_ONT"
spearmansPPE2 <- subset(spearmansPPE,select=-c(source.x,source.y))
rho <- cor(spearmansPPE2$Average_coverage_ILL,spearmansPPE2$Average_coverage_ONT,method="spearman")
print(rho)


spearmansDR <- merge(average_dr_ILL,average_dr_ONT,by=c("Position","Gene"))
colnames(spearmansDR)[colnames(spearmansDR) == "avg_coverage.x"] <- "Average_coverage_ILL"
colnames(spearmansDR)[colnames(spearmansDR) == "avg_coverage.y"] <- "Average_coverage_ONT"
spearmansDR2 <- subset(spearmansDR,select=-c(source.x,source.y))

rho <- cor(spearmansDR2$Average_coverage_ILL,spearmansDR2$Average_coverage_ONT,method="spearman")
print(rho)

spearmansOther <- merge(allOther_average_df_ILL,allOther_average_df_ONT,by=c("Position","Gene"))

colnames(spearmansOther)[colnames(spearmansOther) == "avg_coverage.x"] <- "Average_coverage_ILL"
colnames(spearmansOther)[colnames(spearmansOther) == "avg_coverage.y"] <- "Average_coverage_ONT"
sortedSpearmansOther <- spearmansOther[order(spearmansOther$Position),]
rownames(sortedSpearmansOther) <- NULL
sortedSpearmansOther$Average_coverage_ONT <- sortedSpearmansOther$Average_coverage_ONT* 100

rho <- cor(sortedSpearmansOther$Average_coverage_ILL,sortedSpearmansOther$Average_coverage_ONT,method="spearman")
print(rho)


ratios_Other <- sortedSpearmansOther$Average_coverage_ONT / sortedSpearmansOther$Average_coverage_ILL * 100
ratios_DR <- spearmansDR2$Average_coverage_ONT / spearmansDR2$Average_coverage_ILL * 100
ratios_PPE <- spearmansPPE2$Average_coverage_ONT / spearmansPPE2$Average_coverage_ILL * 100
ratios_All <- spearmansCombined$Average_coverage_ONT / spearmansCombined$Average_coverage_Illumina * 100

median(ratios_Other)
median(ratios_DR)
median(ratios_PPE)
median(ratios_All)
min(ratios_Other)
max(ratios_Other)
min(ratios_DR)
max(ratios_DR)
min(ratios_PPE)
max(ratios_PPE)
min(ratios_All)
max(ratios_All)



snpsDF <- read.csv("Book1.csv")
snpsDF$Sample_ID <- ifelse(snpsDF$Platform == "ONT", snpsDF$Sample, lag(snpsDF$Sample))
ont_df <- snpsDF[snpsDF$Platform == "ONT", ]
illumina_df <- snpsDF[snpsDF$Platform == "Illumina", ]

colnames(ont_df)[3:6] <- paste0(colnames(ont_df)[3:6], "_ONT")
colnames(illumina_df)[3:6] <- paste0(colnames(illumina_df)[3:6], "_Illumina")

snpsDF_reformatted <- merge(ont_df, illumina_df, by = "Sample_ID")

snpsDF_reformatted <- snpsDF_reformatted[, !(colnames(snpsDF_reformatted) %in% c("Sample.x", "Platform.x", "Sample.y", "Platform.y"))]

print(snpsDF_reformatted)


snpsDF_reformatted$SNP_Ratio <- snpsDF_reformatted$no..SNPs._ONT / snpsDF_reformatted$no..SNPs._Illumina

snpsDF_reformatted$difference_SNPs <- abs(snpsDF_reformatted$no..SNPs._ONT - snpsDF_reformatted$no..SNPs._Illumina)
median_difference <- median(snpsDF_reformatted$difference_SNPs)
IQR_difference <- IQR(snpsDF_reformatted$difference_SNPs)
lower_quartile <- quantile(snpsDF_reformatted$difference_SNPs, 0.25)
upper_quartile <- quantile(snpsDF_reformatted$difference_SNPs, 0.75)

print(median_difference)
print(upper_quartile)
print(lower_quartile)

snpsDF_reformatted$ratio_SNPs <- snpsDF_reformatted$no..SNPs._ONT / snpsDF_reformatted$no..SNPs._Illumina
median_ratio <- median(snpsDF_reformatted$ratio_SNPs)


IQR_ratio <- IQR(snpsDF_reformatted$ratio_SNPs)

lower_quartile_ratio <- quantile(snpsDF_reformatted$ratio_SNPs, 0.25)
upper_quartile_ratio <- quantile(snpsDF_reformatted$ratio_SNPs, 0.75)

print(median_ratio)
print(IQR_ratio)