# Credits -----------------------------------------------------------------

# Detect homopolymers of given length in genome assembly
# @author: Vanessa Klauenberg
# @created: 2024-07-01
# @lastchange: 2024-10-07

# Libraries ---------------------------------------------------------------

library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

# Data --------------------------------------------------------------------

# getting the length of the homopolymers (default: 5)
N = readline(prompt = "Enter the wanted homopolymere length: ");

# setting the desired version of the reference genome
v = readline(prompt = "Enter '0' for hg19(default) or '1' for hg38: ");

# Extract homopolymers ----------------------------------------------------

# convert input parameters
N = as.integer(N);
v = as.integer(v);

# only taking the hg38 if '1' is entered, otherwise it's the default hg19
if(v == 1) {
  x <- BSgenome.Hsapiens.UCSC.hg38
}else{
  x <- BSgenome.Hsapiens.UCSC.hg19
}

# only take the main chromosomes, not the alternative
main_chr <- seqlevels(x)[0:25]
# prepare an empty GRanges object for the results
res <- GRanges()
 
# loop through every chromosome				  
for (i in main_chr) {
  print(paste0("Working on chromosome: ", i))
  chromosome <- x[[i]]
  # compute the lengths of runs of equal values
  REV <- Rle(as.vector(chromosome))
  # compute the starting positions of each run
  pos <- cumsum(c(1, runLength(REV)))
  # select the runs that are >=N and consists of one of the four bases
  ix <- which(runLength(REV) >= N & runValue(REV) %in% c("A", "T", "C", "G"))
  # save the results as a GRanges object
  current_res <- GRanges(
    seqnames = Rle(i, length(ix)),
    ranges = IRanges(
      start = pos[ix],
	  end = pos[ix] + runLength(REV)[ix] - 1
    ),
    base = runValue(REV)[ix]
  )
  res <- c(res, current_res)
}

# Save results ------------------------------------------------------------

# save the results as Rdata for later use
save(res, file = "res.Rdata")
