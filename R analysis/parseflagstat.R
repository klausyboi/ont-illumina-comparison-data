setwd("C:/Users/thorp/downloads");
library(stringr)
library(tidyr)
library(dplyr)

flagstatIlll <- readLines("forjody/illuminaFlagstat.txt")
flagstatNano <- readLines("forjody/flagstat.txt")



  
  counter = 0
  samples <- c()
  mapped <- c()
  totalReads <- c()
  
  for (i in seq_along(flagstatIlll)) {
    
    if (grepl("==>", flagstatIlll[i])) {
      samples <- c(samples, flagstatIlll[i])
      totalReads <- c(totalReads, flagstatIlll[i+1])

      
    }else if (grepl("mapped \\(", flagstatIlll[i])) {
      mapped <- c(mapped, flagstatIlll[i])    
    }
    counter+1
  }
  for (i in seq_along(flagstatNano)) {
    
    if (grepl("==>", flagstatNano[i])) {
      samples <- c(samples, flagstatNano[i])
      totalReads <- c(totalReads, flagstatNano[i+1])
      
      
    }else if (grepl("mapped \\(", flagstatNano[i])) {
      mapped <- c(mapped, flagstatNano[i])    
    }
    counter+1
  }
  samples  
  samples <- gsub("==>","",samples)
  samples <- gsub("<==","", samples)
  samplesArray <- str_split(samples,"_",n=Inf,simplify = FALSE)
  newSamples <- c()
  
  for (i in samplesArray){
     newSamples <- c(newSamples, paste0(gsub(" ","",i[1])))
  
  #   ,"_", gsub(" ","",i[2]),"_", gsub(" ","",i[3])))}
  }
  newSamples
  newSamplesOverrule <- c()
  for (j in newSamples) {
    newString <- gsub("_NA", "", j)
    newSamplesOverrule <- c(newSamplesOverrule, newString)
  }
  newSamplesOverruled <- c()
  for(z in newSamplesOverrule){
    newString <- gsub(".fastq.gz.sam.bam.log", "", z)
    newSamplesOverruled <- c(newSamplesOverruled, newString)
  }
  newSamplesOverruled
  
  
  totalReads <- word(totalReads,1)
  mapped <- word(mapped,5)
  mapped <- gsub('\\(',"", mapped)
  
  
  nanoFrame <- data.frame(newSamplesOverruled,mapped,totalReads)
  illuFrame <- data.frame(newSamples,mapped,totalReads)

  colnames(illuFrame)[1] <- "Samples"
  colnames(illuFrame)[2] <- "% Mapped Illumina"
  colnames(illuFrame)[3] <- "Total Reads Illumina"
  colnames(nanoFrame)[1] <- "filename"
  colnames(nanoFrame)[2] <- "% Mapped Nanopore"
  colnames(nanoFrame)[3] <- "Total Reads Nanopore"
nanoFrame
illuFrame

dataTogether <- read.csv("forjody/Depth.csv",sep=",",header = TRUE)

dataCombined <- merge(dataTogether,illuFrame,by="Samples")
dataCombined <- merge(dataCombined,nanoFrame, by="filename")
colnames(dataCombined) <- c(gsub("\\."," ",colnames(dataCombined)))

samplesString <- "Sample 1, Sample 2, Sample 3, Sample 4, Sample 5, Sample 8, Sample 11, Sample 12, Sample 13, Sample 14, Sample 15, Sample 16, Sample 17, Sample 19, Sample 20, Sample 21, Sample 22, Sample 23, Sample 24, Sample 25, Sample 26, Sample 27, Sample 28, Sample 29, Sample 30, Sample 31, Sample 32, Sample 33, Sample 34, Sample 35, Sample 36, Sample 37, Sample 38, Sample 39, Sample 40, Sample 41, Sample 42, Sample 43, Sample 44, Sample 45, Sample 46, Sample 47, Sample 48, Sample 49, Sample 50, Sample 51, Sample 52, Sample 53, Sample 54, Sample 55, Sample 56, Sample 57, Sample 58, Sample 59, Sample 60, Sample 61, Sample 62, Sample 63, Sample 64, Sample 65, Sample 66, Sample 67, Sample 68, Sample 69, Sample 70, Sample 71, Sample 72"

# split the string into separate elements
samplesVector <- strsplit(samplesString, ", ")[[1]]

# create a new column with samples
dataTogether$Samples <- samplesVector
dataTogether

library(ggplot2)
library(reshape2)
df_long <- dataTogether %>% 
  pivot_longer(cols = c("Average.Depth.Nanopore", "Average.Depth.Illumina"), 
               names_to = "Sequencing Platform", 
               values_to = "Average Depth")

# Create the boxplot
ggplot(df_long, aes(x = `Sequencing Platform`, y = `Average Depth`, fill = `Sequencing Platform`)) +
  geom_boxplot() +
  labs(x = "Sequencing Platform", y = "Average Depth") +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E")) +
  theme(axis.text.x = element_text(angle=90,vjust = 0.5, hjust=1)) +
  ggtitle("Comparison of Average Depth between Nanopore and Illumina")


data <- read.csv("datacombined.csv",check.names = FALSE,header=TRUE,stringsAsFactors = FALSE)

split_dataframe <- function(df) {
  # Split the data frame by Nanopore and Illumina columns
  df_nano <- df[, c(1, 3, 5, 7, 9)]
  df_illumina <- df[, c(2, 4, 6, 8, 10)]
  
  # Rename columns
  names(df_nano) <- c("Sample", "Average Depth", "No. SNPs", "% Mapped", "Total Reads")
  names(df_illumina) <- c("Sample", "Average Depth", "No. SNPs", "% Mapped", "Total Reads")
  
  # Return a list of the split data frames
 
  return(list(df_nano, df_illumina))
}

# Split the data frame
dfs <- split_dataframe(data)

# Bind the split data frames together
result <- do.call(rbind, dfs)

# Print the result
print(result)

write.csv(dataCombined,"datacombined.csv",quote=FALSE,row.names = FALSE)



