#Code cobbled together from: 
#https://stackoverflow.com/questions/29005904/heatmap-with-categorical-variables-and-with-phylogenetic-tree-in-r
#http://www.polarmicrobes.org/merging-a-phylogenetic-tree-with-a-heatmap-in-r/

#load packages
library("ape")
library(gplots)

#----
#retrieve tree in newick format, and do a number of things to get it ready for plotting:
mytree <- read.tree("~/GitHub/Locus-Coverage-Heatmap/example_data/RAxML_bipartitions.FUN_Reduced_Set_75perc.tre")
#plot(mytree)
#identify(mytree)
#uncomment the above lines to plot tree and ID the node number to root with. 

mytree <- root(mytree, node = 79, resolve.root = FALSE)
#change the node number to the proper root value, determined above.   

mytree_brlen <- compute.brlen(mytree, method="Grafen") 
#make branches uniform length
mytree_brlen <- multi2di(mytree_brlen) 
#make tree binary

is.ultrametric(mytree_brlen)
is.rooted(mytree_brlen)
is.binary.tree(mytree_brlen)
#Is the tree all these things? If not, something is not right and the plot wont work

hc <- as.hclust(mytree_brlen) 
dend <- as.dendrogram(hc)
#turn the phylo tree to a dendrogram object, compulsory step as as.dendrogram doesn't have a method for phylo objects.

#par(mar=c(1,1,1,13))
#plot(dend, horiz=TRUE, cex.lab=0.5)
#labels(dend)
#Sanity checks for the dendrogram order. 

#----
#Get the presence/absence matrix. 
temp_dataset <- as.data.frame(read.table("~/GitHub/Locus-Coverage-Heatmap/example_data/matches_uce.csv", header=T, sep=","))
temp_dataset[temp_dataset == 0] <- NA
holder <- as.matrix(t(temp_dataset[-c(1)]))
holder <- holder[ ,order( colSums( is.na( holder ) ), decreasing = F ) ]

#force row order so that it matches the order of leaves in dendrogram
clade_order <- order.dendrogram(dend)
clade_name <- labels(dend)
clade_position <- data.frame(clade_name, clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]

new_order <- match(clade_position$clade_name, row.names(holder))
combined_ordered_matrix <- holder[new_order,]

#plot the heatmap, many options here to play with to adjust plot to your liking. 
heatmap.2(combined_ordered_matrix, Rowv=dend, Colv=F, dendrogram='row', reorderfun=reorderfun, col = 
	colorRampPalette(c("blue","green","red"))(3), na.color = 'gray90',
	key=FALSE,trace="none",tracecol=NA,density.info="none",
	cexRow=.7,cexCol=.25,srtCol=45, labCol = "",
	margins=c(1,5), lwid=c(1.5,4), lhei=c(0.3,3.5))
mtext("Presence Absence UCE", 3, line=3)