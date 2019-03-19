# installation
## Installation of the latest released version
install.packages('GOplot')
library(GOplot)
#Load the dataset
data(EC)
# Get a glimpse of the data format of the results of the functional analysis
head(EC)
head(EC$david)
# Generate the plotting object
circ <- circle_dat(EC$david, EC$genelist)
head(circ)
# Generate a simple barplot
GOBar(subset(circ, category == "BP"))
# Facet the barplot according to the categories of the terms
GOBar(circ, display = 'multiple')
# Get the Y color and facet the barplot, add a title and change the colour scale for the z-scorr
require(rPlotter)
pal = extract_colours("https://b-ssl.duitang.com/uploads/item/201605/19/20160519133736_tFaEH.jpeg",num_col = 3)
barplot(1:20, col=pal)
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', 
      zsc.col = pal)
# Generate the bubble plot with a lable threshold of 3
GOBubble(circ, labels = 3)
# Add a title, change the colour of the circles, facet the plot according to the categories and 
#change the label threshold
GOBubble(circ, title = 'Bubble plot', colour = pal, display = 'multiple', labels = 3)
# Reduce redundant terms with a gene overlap >= 0.75
reduce_circ <- reduce_overlap(circ, overlap = 0.75)
# plot it 
GOBubble(reduce_circ, labels = 2.8)
# Generate a circular visualization of the results of gene-annotation enrichment analysis
GOCircle(circ)
# Generate a circular visualization of selected terms
IDs <- c('GO:0007507', 'GO:0001568', 'GO:0001944', 'GO:0048729', 'GO:0048514', 'GO:0005886',
         'GO:0008092', 'GO:0008047')
GOCircle(circ, nsub = IDs)
# Generate a circular visualization for 10 terms
GOCircle(circ, nsub = 10)
# Define a list of genes which you think are interesting to look at. The item EC$genes of the toy
# sample contains the data frame of selected genes and their logFC. Have a look...
head(EC$genes)
# Since wo have a lot of significantly enriched processes wo selected some specific ones
EC$process
# Now it is the time to generate the binary matrix
chord <- chord_dat(circ, EC$genes, EC$process)
head(chord)
# Generate the matrix with a list of selected genes
chord <- chord_dat(data = circ, genes = EC$genes)
head(chord)
# Generat the matrix with a list of selected processes
chord <- chord_dat(data = circ, process = EC$process)
head(chord)
# Creat the plot
chord <- chord_dat(data = circ, genes = EC$genes, process = EC$process)
head(chord)
chord
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
# Display only genes which are assigned to at least three processes
GOChord(chord, limit = c(2,0), gene.order = 'logFC')
# First, we use the chord object without logFC column to creat the heatmap
chord <- chord_dat(circ, EC$genes, EC$process)
head(chord)
GOHeat(chord[,-8], nlfc = 0, fill.col = pal)
# Now we generate the heatmap with logFC valuse and user-defined colour scale
GOHeat(chord, nlfc = 1, fill.col = pal)
# Golden eye
GOCluster(circ, EC$process, clust.by = 'logFC', term.width = 2)
GOCluster(circ, EC$process, clust.by = 'term', lfc.col = c('darkgoldenrod1', 'black', 'cyan1'))
