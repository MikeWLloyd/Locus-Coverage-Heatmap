source('~/GitHub/Locus-Coverage-Heatmap/missing.pattern.plot.R', chdir = TRUE)
#I forgot where I got this function from...

temp_dataset <- as.data.frame(read.table("~/GitHub/Locus-Coverage-Heatmap/example_data/matches_uce.csv", header=T, sep=","))

temp_dataset[temp_dataset==0] <- NA
holder <- temp_dataset[-c(1)]
loci <- temp_dataset$uce
named <- names(holder)

#three different plots, ordered in different ways. Use as needed. 
missing.pattern.plot(holder, x.order = T, y.order = T, clustered=F, gray.scale = F, main="UCE Loci Matches")

missing.pattern.plot(holder, clustered=T, gray.scale = F, main="UCE Loci Matches")

missing.pattern.plot(holder, x.order = T, clustered=T, gray.scale = F, main="UCE Loci Matches")


