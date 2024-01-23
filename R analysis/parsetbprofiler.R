rm(list=ls() )
setwd("C:/Users/thorp/downloads")
library(dplyr)
##load data
tbprofilerdata <- readLines("alltogetherlol.txt")

#initialise vectors and create an empty df with specified tbprofiler columns
IDs <- c()
lineage <- c()
dr <- c()
mutations <- c()
df <- data.frame(IDs = character(), Lineage = character(), `Drug Resistance` = character(),  `Genotypic Resistance` = character(), stringsAsFactors = FALSE)

current_id <- ""
current_lin <- ""
current_dr <- ""
current_genotypic <- ""


##terate over each line in the tbprofilerdata
for (i in seq_along(tbprofilerdata)) {
  ##extract and process specific information based on the line's content
  if (grepl("ID:", tbprofilerdata[i])) {
    current_id <- gsub("ID: ","",tbprofilerdata[i])
  } else if (grepl("Strain:", tbprofilerdata[i])) {
    current_lin <- gsub("Strain: ","",tbprofilerdata[i])
  } else if (grepl("Drug-resistance:", tbprofilerdata[i])) {
    current_dr <- gsub("Drug-resistance: ","",tbprofilerdata[i])
  } else if (grepl("Resistance report", tbprofilerdata[i])) {
    i <- i + 2
    while (!grepl("Resistance variants report", tbprofilerdata[i])) {
      if (!grepl("^\\s*$", tbprofilerdata[i])) { 
        genotypic_line <- strsplit(tbprofilerdata[i], "\t")[[1]]
        drug <- genotypic_line[1]
        genotypic_resistance <- genotypic_line[2]
        mutations <- genotypic_line[3]
        ##create a new row with the extracted information
        new_row <- data.frame(IDs = current_id, Lineage = current_lin, `Drug-Resistance` = current_dr, 
                              `Genotypic Resistance` = genotypic_resistance, `Mutations` = mutations,
                              stringsAsFactors = FALSE)
        ##add the new row to the dataframe
        df <- rbind(df, new_row)
      }
      i <- i + 1
    }
  }
}
#clean the ids column from specific patterns
df$IDs <- sub("\\.\\./", "", df$IDs)
##further clean and filter the dataframe based on various criteria
df <- subset(df, !grepl("ONT.sorted.bam", IDs))
df <- na.omit(df)
df <- filter(df, Drug.Resistance != "Sensitive")
df <- filter(df, Genotypic.Resistance != "Genotypic Resistance")
df <- filter(df, IDs != "S10_ONT.sorted.bam")
df <- distinct(df)
library(dplyr)
library(stringr)

##get only pncA genes
df_pncA <- filter(df, str_detect(Mutations, fixed("pncA")))


##aggregate data to create a summary dataframe
new_df <- aggregate(df[, c("Lineage", "Drug.Resistance", "Genotypic.Resistance")], 
                    by = list(IDs = df$IDs), 
                    FUN = function(x) {
                      if(is.factor(x)) {
                        return(as.character(unique(x)))
                      } else if (is.character(x)) {
                        return(paste(unique(x), collapse = " + "))
                      } else {
                        return(toString(x))
                      }
                    })

##merge the summary dataframe with external datasets (expStats and new_df)
merged <- merge(new_df,expStats,by.x = "IDs",by.y = "filename")
mergedAll <- merge(merged,new_df, by.x="sample.code",by.y = "IDs")
##rename
colnames(mergedAll) <- c("Illumina IDs","ONT IDs","Lineage ONT","Drug Resistance ONT","Genotypic Resistance ONT","Lineage Illumina","Drug Resistance Illumina","Genotypic Resistance Illumina")
##order by ONT Ids
mergedAll <-mergedAll[order(mergedAll$`ONT IDs`),]
row.names(mergedAll) <- c()
##export data
write.csv(mergedAll,"tb_profiler_results.csv",quote=FALSE,row.names=FALSE)

##get list of resistances and filter string obtained
list_of_resistances <- df$Genotypic.Resistance
list_of_resistances <- list_of_resistances[nzchar(list_of_resistances)]
list_of_resistances <- gsub("\\t"," ",list_of_resistances)
list_of_resistances <- gsub(" R "," ",list_of_resistances)

##calculate the freq and percentage of specific resistance types
library(stringr)
total <- length(list_of_resistances)
sum(str_count(list_of_resistances,"Pyrazinamide R pncA p.Thr142Lys"))
list_of_resistances <-sort(list_of_resistances)
list_of_resistances

total_resistances <- data.frame(value = unique(list_of_resistances),
                 count = as.numeric(table(list_of_resistances)))
total_resistances
library(dplyr)
##create a dataframe to store the counts and percentages of resistances
total_resistances <- data.frame(list = list_of_resistances)
total_resistances_plus_percentage <- total_resistances %>%
  group_by(list) %>%
  summarize(count = n()) %>%
  mutate(percentage = count/139*100)

total_resistances_plus_percentage


library(tidyr)

##separate the concatenated data into individual columns
complete <-total_resistances_plus_percentage %>% separate(list, c("Drug", "Gene", "Mutation"), sep = " ")
complete$percentage <- round(complete$percentage,2)
names(complete)[5] <- "Percentage (%)"


##export the data
write.csv(complete,"total_resistances.csv",row.names = FALSE,quote = FALSE)
