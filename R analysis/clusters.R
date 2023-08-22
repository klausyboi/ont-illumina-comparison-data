setwd("C:/Users/thorp/Downloads")
rm(list=ls())
library(cluster)


cluster_analysis <- function(threshold, distance) {
 
  clusters <- agnes(distance, diss = TRUE, method = "average")
  clusters <- as.data.frame(cutree(as.hclust(clusters), h = threshold))
  colnames(clusters) <- "cluster_id"
  
  
  clusterIDs <- unique(subset(clusters, duplicated(clusters$cluster_id))$cluster_id)
  clusters <- subset(clusters, cluster_id %in% clusterIDs)
  
  
  cluster_sizes <- table(clusters$cluster_id)
  
  total_clusters <- length(unique(clusters$cluster_id))
  median_size <- median(cluster_sizes)
  range_size <- range(cluster_sizes)
  
  return(list(total_clusters=total_clusters, median_size=median_size, range_size=range_size))
}


tbprofiler <- read.delim("tbprofiler.txt")
names(tbprofiler)[names(tbprofiler) =="sample"] <- "Sample"
matrix_data <- read.delim("matrix.tsv", row.names=1)
distance <- as.dist(matrix_data)
clusters <- agnes(distance, diss = TRUE, method = "average")
clusters <- as.data.frame(cutree(as.hclust(clusters), h = 20))
colnames(clusters) <- "cluster_id"


clusterIDs <- unique(subset(clusters, duplicated(clusters$cluster_id))$cluster_id)
clusters <- subset(clusters, cluster_id %in% clusterIDs)

clusters <- cbind(Sample = rownames(clusters), clusters)
cluster_sizes <- table(clusters$cluster_id)
hist(cluster_sizes, breaks=500, main="Distribution of Cluster Sizes",xlim=c(0,25),
     xlab="Number of Samples in Cluster", ylab="Number of Clusters", col="lightblue", border="black")
total_clusters <- length(unique(clusters$cluster_id))
print(paste("Total number of clusters:", total_clusters))
cluster_df <- as.data.frame(cluster_sizes)
names(cluster_df) <- c("Cluster_ID", "Frequency")
print(cluster_df)

median_size <- median(cluster_df$Frequency)
range_size <- range(cluster_df$Frequency)
range_diff = range_size[2] - range_size[1]
print(paste("Range of cluster sizes:", range_size[1], "to", range_size[2], "with a difference of", range_diff))

merged <- merge(clusters,tbprofiler)

merged <- merged[order(merged$cluster_id),]

write.table(merged, file = "Transmission_clusters.tsv", sep="\t", quote = F, row.names = F)




thresholds <- c(1, 5, 20)
results <- lapply(thresholds, function(t) cluster_analysis(t, distance))


for (i in 1:length(thresholds)) {
  print(paste("For threshold", thresholds[i], "SNPs:"))
  print(paste("  Total number of clusters:", results[[i]]$total_clusters))
  print(paste("  Median cluster size:", results[[i]]$median_size))
  print(paste("  Range of cluster sizes:", results[[i]]$range_size[1], "to", results[[i]]$range_size[2]))
}


#other attempt


dist_matrix <- as.dist(matrix_data)
mds_result <- cmdscale(dist_matrix)

custom_breaks <- c(2.5, 3.5, max(cluster_counts)+0.5)


hc <- hclust(dist_matrix, method="average")
clusters <- cutree(hc, h=10)

cluster_counts <- table(clusters)
hist(cluster_counts, breaks=custom_breaks, main="Distribution of Cluster Sizes", xlab="Number of Samples in Cluster", ylab="Number of Clusters", col="lightblue", border="black")

median_size <- median(cluster_sizes)
min_size <- min(cluster_sizes)
max_size <- max(cluster_sizes)

library(gplots)
heatmap.2(matrix_data, hclustfun = function(x) hclust(x, method="average"), trace="none", scale="none", key=FALSE, density.info="none", margins=c(5,5), cexRow=0.5, cexCol=0.5)
snp_diffs <- matrix_data[upper.tri(matrix_data)]
hist(snp_diffs, main="Histogram of SNP Differences", xlab="SNP Differences", breaks=50)
color_palette <- colorRampPalette(c("red", "yellow"))(256)


heatmap.2(as.matrix(matrix_data),
          main="SNP Distance Heatmap",
          trace="none",
          margins=c(8,8),
          dendrogram="both", 
          col=color_palette,
          scale="none")
