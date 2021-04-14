## Q1
## What is the GC content of “chr22” in the “hg19” build of the human genome?
## The reference genome includes “N” bases; you will need to exclude those.

##Install packages and genome BioString

library(BiocManager)
library(BSgenome)
available.genomes() #Find appropriate genome
install("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")

## Load Chr22 into memory
hg19.chr22 = Hsapiens$chr22

## GC content
GC = letterFrequency(hg19.chr22, "GC") # number of G|C's in chr22
Total_bases = letterFrequency(hg19.chr22, "ATGC") # Total number of bases, not including 'N' or any other placeholders
GC/Total_bases # GC content as fraction


## Q2
## What is mean GC content of H3K27me3 “narrowPeak” regions from Epigenomics Roadmap from the H1 stem cell line on chr 22.

## Load data
library(AnnotationHub)
ah = AnnotationHub()
ah = subset(ah, species == "Homo sapiens")
ah1 = query(ah,c("E003", "H3K27me3", "narrowPeak")) # create new query
ah1 # check metadata to ensure correct dataset has been chosen
H3K27 = ah1[[1]] # load GRange into memory
H3K27.chr22 = keepSeqlevels(H3K27, "chr22", pruning.mode = "coarse")

## GC content of narrowPeak ranges
freq = alphabetFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg19, H3K27.chr22)) # Gets the sequences and stores base frequency
gc1 = freq[, 'C'] + freq[, 'G'] # Number of Gs and Cs in each peak
Total_bases1 = rowSums(freq) # Number of bases in each peak
gcFreq =  gc1/Total_bases1
mean(gcFreq) # The mean of the fraction of GC content in each peak


##Q3
## What is the correlation between GC content and “signalValue” of these regions (on chr22)?

H3K27.chr22.sig = as.numeric(DataFrame(H3K27.chr22)$signalValue) # Create DataFrame of signal values
cor(H3K27.chr22.sig, gcFreq) # Calculate correlation between the two vectors


## Q4
## What is the correlation between the “signalValue” of the “narrowPeak” regions and the average “fc.signal” across the same regions?

##Load data
ah2 = query(ah,c("E003", "H3K27me3", "fc.signal")) # New AH query for fc.signal
H3K27.bw = ah2[[1]]

## Import chr22 data
fc.signal22.rle = import(H3K27.bw, which = H3K27.chr22, as = "RleList")
fc.signal22.rle = fc.signal22.rle$chr22

## Calculate mean fc.signal for all narrowPeak ranges and correlation with signalValue from Q3
fc.signal22.means = aggregate(fc.signal22.rle, H3K27.chr22, FUN = mean)
cor(fc.signal22.means, H3K27.chr22.sig)


## Q5
## How many bases on chr22 have an fc.signal greater than or equal to 1?

## Get fc.signal Rle of whole of chr22
E003.fc.signal.rle = import(H3K27.bw, which = GRanges(seqnames = "chr22", ranges = IRanges(1, 10^8)), as = "RleList")
# NOTE: 10^8 was chosen as chr22 is slightly shorter in bases
E003.fc.signal.rle = E003.fc.signal.rle$chr22 # Converts from RleList to Rle containing only chr22

## Finally, calculate sum of the width of each region with fc.signal >= 1
sum(width(slice(fc.signal22a.rle, 1)))


## Q6
## Identify the regions of the genome where the signal in E003 is 0.5 or lower and the signal in E055 is 2 or higher.

## First, load both fc.signal datasets as GRanges
ah3 = query(ah,c("E055", "H3K27me3", "fc.signal")) # create new query for E005 dataset
H3K27a.bw = ah3[[1]] # Load E005 fc.signal dataset

## Import data from chr22 of E055 dataset as RleList
E055.fc.signal.rle = import(H3K27a.bw, which = GRanges(seqnames = "chr22", ranges = IRanges(1, 10^8)), as = "RleList")
E055.fc.signal.rle = E055.fc.signal.rle$chr22

## Apply filter to fc.signal values using slice function to obtain Views objects
E003.fc.signal.slice.vi = slice(E003.fc.signal.rle, upper = 0.5, includeUpper = TRUE)
E055.fc.signal.slice.vi = slice(E055.fc.signal.rle, lower = 2, includeLower = TRUE)

## Convert both to IRanges
E003.fc.signal.slice.ir = as(E003.fc.signal.slice.vi, "IRanges")
E055.fc.signal.slice.ir = as(E055.fc.signal.slice.vi, "IRanges")

## Find number of intersecting bases of the 2 IRanges
sum(width(intersect(E003.fc.signal.slice.ir, E055.fc.signal.slice.ir)))


## Q7
## What is the average observed-to-expected ratio of CpG dinucleotides for CpG Islands on chromosome 22?

## Download UCSC hg19 CpG island data from AnnotationHub
ah5 = display(ah)
#search "hg19" in genome seach bar and "CpG" in description then click "send" to store AH object in ah5
hg19.CpI.gr = ah5[[1]] # Store GRanges object
hg19.CpI.22.gr = keepSeqlevels(hg19.CpI.gr, "chr22", pruning.mode = "coarse") # Keep only chr22 ranges

## Download hg19 genom data using BSgenome
available.genomes() # View available genome to find "BSgenome.Hsapiens.UCSC.hg19"
library("BSgenome.Hsapiens.UCSC.hg19")

## Subset hg19 genome by CpI Granges to get DNAStringSet
hg19.CpI.22.DNASS = getSeq(Hsapiens, hg19.CpI.22.gr)

## average observed-to-expected ratio of CpG dinucleotides
CGFreq.CpI = dinucleotideFrequency(hg19.CpI.22.DNASS)[,7] # CG dinucleotide frequencies
CFreq.CpI = alphabetFrequency(hg19.CpI.22.DNASS)[,2] # C nucleotide frequencies
GFreq.CpI = alphabetFrequency(hg19.CpI.22.DNASS)[,3] # G nucleotide frequencies
seqLength.CpI = width(hg19.CpI.22.DNASS) # Sequence lengths
ExpFreq.CpI = (CFreq.CpI*GFreq.CpI)/seqLength.CpI # Expected frequency
mean(CGFreq.CpI/ExpFreq.CpI) # The mean observed-to-expected ratio of CpG dinucleotides for CpG Islands on chromosome 22


## Q8
## How many TATA boxes are there on chr 22 of build hg19 of the human genome (both strands)?

# Create TATA box DNAstring and reverse complement DNAstring
TATAbox = DNAString("TATAAA")
revTATAbox = reverseComplement(TATAbox)
countPattern(TATAbox, Hsapiens$chr22) + countPattern(revTATAbox, Hsapiens$chr22)


## Q9
## How many promoters of transcripts on chromosome 22 containing a coding sequence,
## contains a TATA box on the same strand as the transcript?
##
## Clarification: Use the TxDb.Hsapiens.UCSC.hg19.knownGene package to define transcripts and coding sequence.
## Here, we defined a promoter to be 900bp upstream and 100bp downstream of the transcription start site.



## Q10
## It is possible for two promoters from different transcripts to overlap, in which case the regulatory features inside the overlap
## might affect both transcripts. This happens frequently in bacteria.
##
## How many bases on chr22 are part of more than one promoter of a coding sequence?
##
## Use the TxDb.Hsapiens.UCSC.hg19.knownGene package to define transcripts and coding sequence. Here,
## we define a promoter to be 900bp upstream and 100bp downstream of the transcription start site. In this case, ignore strand in the analysis.
