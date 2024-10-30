# Credits -----------------------------------------------------------------

# Filter homopolymers according to a BED file
# @author: Vanessa Klauenberg
# @created: 2024-07-16
# @lastchange: 2024-08-07

# Libraries ---------------------------------------------------------------

library(GenomicRanges)

# Data --------------------------------------------------------------------

# loading the bed file
path <- "../data/xGenExomeResearchPanelV2_Exons.bed"
x <- read.delim(path, header = FALSE, 
                col.names = c("chrom", "start", "end", "gene"))
                
# load the result data from the previous execution of 01_homopolymers.r
load(file= "res.Rdata")

# Filter homopolymers -----------------------------------------------------

# convert the res file and bed file to GRanges object
res_gr <- GRanges(res)
bed_gr <- GRanges(seqnames = x$chrom,
                  ranges = IRanges(start = x$start, end = x$end))

# calculate the length of the bed file in MB for later use
range_lengths <- width(bed_gr)
lengths_mb <- range_lengths/1e6
total_length <- sum(lengths_mb)

# find the homopolymers in the bed file
overlaps <- findOverlaps(bed_gr, res_gr)
# extract the starting positions
matching_homopolymers <- res_gr[subjectHits(overlaps)]

# Save results ------------------------------------------------------------

# if there are matching positions, save the infos of the associated homopolymer
if(length(matching_homopolymers) > 0) {
  save(matching_homopolymers, file= "matching_homopolymers.Rdata")
# else just print that there are no homopolymers in the bed file
}else{
  print("No homopolymers found")
}

# save the total length for later use
save(total_length, file = "total_length.Rdata")

