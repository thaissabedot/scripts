#Gviz package from Bioconductor
#source: https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.pdf

#load the library on R
library(Gviz) 
library(GenomicRanges) 

#example data: A sample set of CpG island coordinates saved in a GRanges object.
data(cpgIslands) 

#Retrieve the chromosome names associated to the "cpgIslands" object
chr <- as.character(unique(seqnames(cpgIslands)))

#Retrieve the genome version related to the "cpgIslands" object
gen <- genome(cpgIslands)

#AnnotationTrack constructor: it is fairly flexible and can accomodate many different types of inputs. For instance, the start and end coordinates of the annotation features could be passed in as individual arguments start and end, as a data.frame or even as an IRanges or GRangesList object.
atrack <- AnnotationTrack(cpgIslands, name = "CpG")

#Single function that can plot all Gviz objects.
plotTracks(atrack)

#Indicate the genomic coordinates to provide some reference according to the CpG islands distance.
gtrack <- GenomeAxisTrack()

plotTracks(list(gtrack, atrack))

#Visual representation of a chromosome, with the different chromosomal staining bands indicated by color, and the centromer (if present) indicated by the shape. The necessary information to produce this visualization is stored in online data repositories, for instance at UCSC.
#If this function takes to long to load, please, load the file "itrack.rda" by doing: load("itrack.rda) and procede.
itrack <- IdeogramTrack(genome = gen, chromosome = chr)

plotTracks(list(itrack, gtrack, atrack))


#So far we have only looked at very basic annotation features and how to give a point of reference to our plots. Naturally, we also want to be able to handle more complex genomic features, such as gene models. One potential use case would be to utilize gene model information from an existing local source. Alternatively, we could dowload such data from one of the available online resources like UCSC or ENSEBML, and there are constructor functions to handle these tasks. For this example we are going to load gene model data from a stored data.frame
data(geneModels)
grtrack <- GeneRegionTrack(geneModels, genome = gen, chromosome = chr, name = "Gene Model")

plotTracks(list(itrack, gtrack, atrack, grtrack))

#plotTracks supports the "from" and "to" arguments that let us choose an arbitrary genomic range to plot, so we can zoom in or out on a particular plotting region.
plotTracks(list(itrack, gtrack, atrack, grtrack), from = 26700000, to = 26750000)

#We can also extend the view on one or both end of the plot by setting the arguments "extend.left" and "extend.right". In addition to positive or negative absolute integer values one can also provide a float value between -1 and 1 which will be interpreted as a zoom factor, i.e., a value of 0.5 will cause zooming in to half the currently displayed range.
plotTracks(list(itrack, gtrack, atrack, grtrack), extend.left = 0.5, extend.right = 1000000)


#When zooming further in it may become interesting to take a look at the actual genomic sequence at a given position, and the Gviz package provides the track class SequenceTrack that let’s you do just that. Among several other options it can draw the necessary sequence information from one of the BSgenome packages.
library(BSgenome.Hsapiens.UCSC.hg19)
strack <- SequenceTrack(Hsapiens, chromosome = chr)
plotTracks(list(itrack, gtrack, atrack, grtrack, strack), from = 26591822, to = 26591852, cex = 0.8)


#We can use the function DataTrack to add all sorts of numeric data to our genomic coordinate plots. For demonstration purposes, we can create a simple DataTrack object from randomly sampled data.
set.seed(255)
lim <- c(26700000, 26750000)
coords <- sort(c(lim[1], sample(seq(from = lim[1], to = lim[2]), 99), lim[2]))
dat <- runif(100, min = -10, max = 10)
dtrack <- DataTrack(data = dat, start = coords[-length(coords)], end = coords[-1], chromosome = chr, genome = gen, name = "Uniform")
plotTracks(list(itrack, gtrack, atrack, grtrack, dtrack), from = lim[1], to = lim[2])


#But we can also overlay multiple types of plot and compare the values across samples, for example
data(twoGroups)
dTrack <- DataTrack(twoGroups, name = "uniform")
plotTracks(dTrack, type = c("boxplot", "a", "g"))


#We can change the properties of individual track objects, such as background color, gene name, etc.
plotTracks(list(itrack, gtrack, atrack, grtrack),background.panel = "#FFFEDB", background.title = "darkblue")

#Highlighting a region
ht <- HighlightTrack(trackList = list(atrack, grtrack,dtrack), start = c(26705000, 26720000), width = 7000,chromosome = 7)
