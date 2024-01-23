setwd("C:/Users/thorp/Downloads")

library(stringr)
library(tidyr)
library(dplyr)
rm(list=ls())
##read in data
samediff <- readLines("samediff.txt")
illDiff <- readLines("illdiff.txt")
nanoDiff <- readLines("nanodiff.txt")

samples <- c()
totalReadsSame <- c()
totalReadsNano <- c()
totalReadsIll <- c()

##parse data for joint variants
for (i in seq_along(samediff)) {
  
  if (grepl("sample", samediff[i])) {
    samples <- c(samples, samediff[i])
    number <- gsub("length=","",samediff[i+1])
    totalReadsSame <- c(totalReadsSame, number)
    
    
  }
  

}
##nanopore specific
for (i in seq_along(nanoDiff)) {
  
  if (grepl("sample", nanoDiff[i])) {

    number <- gsub("length=","",nanoDiff[i+1])
    totalReadsNano <- c(totalReadsNano, number)
    
    
  }
  
  
}
##illumina specific
for (i in seq_along(illDiff)) {
  
  if (grepl("sample", illDiff[i])) {
number <- gsub("length=","",illDiff[i+1])
    totalReadsIll <- c(totalReadsIll, number)
    
    
  }
  
  
}


##process and clean
totalSamples <- c()
for (j in samples){
  newString <- gsub("_sv_calls.vcf" ,"",j)
  newString <- gsub("filename: sample","",newString)
  totalSamples <- c(totalSamples, as.integer(newString))
}
samplesSplit <- c()
samplesSplit <- str_split(totalSamples,"_",n=Inf,simplify = FALSE)
newSamples <- c()
for (z in samplesSplit){
  
  newSamples <- c(newSamples, paste0(gsub(" ","",z[2]),"_", gsub(" ","",z[3]),"_", gsub(" ","",z[4])))
}

##remove unwanted string patterns
newSamplesOverrule <- c()
for (j in newSamples) {
  newString <- gsub("_NA", "", j)
  newSamplesOverrule <- c(newSamplesOverrule, newString)
}


#samplesString <- " 1,  2,  3,  4,  5,  8,  11,  12,  13,  14,  15,  16,  17,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72"
samples <- paste( 1:72, sep = " ")

samples2 <- paste("Samples", 1:72, sep = " ")
samples2 <- samples2[-c(6,7,9,10,18)]

##create dataframes for each dataset and merge them
same_frame <- data.frame(totalSamples,totalReadsSame)
same_frame <- same_frame %>% arrange(totalSamples)
colnames(same_frame)[1]<- "Samples"

colnames(same_frame)[2] <- "Joint Variants"



nano_frame <- data.frame(totalSamples,totalReadsNano)
nano_frame <- nano_frame %>% arrange(totalSamples)

colnames(nano_frame)[1]<- "Samples"
colnames(nano_frame)[2]<- "Unique ONT Variants"

ill_frame <- data.frame(totalSamples,totalReadsIll)
ill_frame <- ill_frame %>% arrange(totalSamples)
colnames(ill_frame)[1]<- "Samples"
colnames(ill_frame)[2]<- "Unique Illumina Variants"

dataCombined <- merge(same_frame,ill_frame,by="Samples")
dataCombined <- merge(dataCombined,nano_frame, by="Samples")
dataCombined <- dataCombined %>% arrange(Samples)
dataCombined <- dataCombined[-c(6,8,9,17),]

dataCombined$`Samples` <- as.numeric(as.character(dataCombined$`Samples`))
dataCombined$`Joint Variants` <- as.numeric(as.character(dataCombined$`Joint Variants`))
dataCombined$`Unique Illumina Variants` <- as.numeric(as.character(dataCombined$`Unique Illumina Variants`))
dataCombined$`Unique ONT Variants` <- as.numeric(as.character(dataCombined$`Unique ONT Variants`))
dataCombined <- arrange(dataCombined, Samples)


dataCombined['Total'] = dataCombined['Joint Variants'] + dataCombined['Unique Illumina Variants'] + dataCombined['Unique ONT Variants']
dataCombined['Total ONT'] = dataCombined['Joint Variants']+ dataCombined['Unique ONT Variants']
dataCombined['Total Illumina'] = dataCombined['Joint Variants']+ dataCombined['Unique Illumina Variants']

##read in large dataset
sameTotal <- readLines("same_large_svs.txt")
illTotal <- readLines("illu_large_svs.txt")
nanoTotal <- readLines("nano_large_svs.txt")


##get large SV counts from text files
sameLarge = c()

for (i in seq_along(sameTotal)) {
  if (grepl("sample", sameTotal[i])) {
    number1 <- gsub(".*_sv_calls: ([0-9]+)$", "\\1", sameTotal[i])
    sameLarge <- c(sameLarge, number1)
  }
}

nanoLarge = c()
for (i in seq_along(nanoTotal)) {
  if (grepl("sample", nanoTotal[i])) {
    number1 <- gsub(".*_sv_calls: ([0-9]+)$", "\\1", nanoTotal[i])
    nanoLarge <- c(nanoLarge, number1)
  }
}
illLarge = c()
for (i in seq_along(illTotal)) {
  if (grepl("sample", illTotal[i])) {
    number1 <- gsub(".*_sv_calls: ([0-9]+)$", "\\1", illTotal[i])
    illLarge <- c(illLarge, number1)
  }
}

##create dfs for each set of large sv counts
nano_large_df <- data.frame(totalSamples,nanoLarge)
nano_large_df <- nano_large_df %>% arrange(totalSamples)
illu_large_df <- data.frame(totalSamples,illLarge)
illu_large_df <- illu_large_df %>% arrange(totalSamples)
same_large_df <- data.frame(totalSamples,sameLarge)
same_large_df <- same_large_df %>% arrange(totalSamples)

colnames(same_large_df)[1]<- "Samples"
colnames(same_large_df)[2] <- "Joint Large SVs"
colnames(nano_large_df)[1]<- "Samples"
colnames(nano_large_df)[2] <- "ONT Large SVs"
colnames(illu_large_df)[1]<- "Samples"
colnames(illu_large_df)[2] <- "Illumina Large SVs"


##merge large datasets
dataCombinedBefore <- merge(same_large_df,illu_large_df,by="Samples")
dataCombinedBefore <- merge(dataCombinedBefore,nano_large_df, by="Samples")
dataCombinedBefore <- dataCombinedBefore[-c(6,8,9,17),]
##combine all
totalCombined <- merge(dataCombinedBefore,dataCombined, by="Samples")


##export data
write.csv(totalCombined,"comparisonDelly.csv",quote = FALSE,row.names = FALSE)
