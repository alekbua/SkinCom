library(plyr)
library(dplyr)
library(tidyr)
library(phyloseq)

setwd("/home/SkinCom_expt2/bowtie_counts/")

options(repr.matrix.max.rows=600, repr.matrix.max.cols=200)

# Taxa table
# Import taxonomy data with all classification fields filled in.
tax.table.all <- read.csv("../phyloseq_taxa_pilot2.csv", sep = ",")
tax.matrix.all <- as.matrix(tax.table.all)
tax.phylo.all <- tax_table(tax.matrix.all)
rownames(tax.phylo.all) <- paste0(tax.table.all$X) # Paste feature IDs from CSV file as matrix row names
tax.phylo.all <- subset(tax.phylo.all, select = -X ) # Remove extra column of row numbers
class(tax.phylo.all)
tax.phylo.all

# Read in metadata CSV.
metadata <- read.csv("../Pilot2_metadata_phyloseq.tsv", sep = "\t")
rownames(metadata) <- paste0(metadata[,1]) # Paste column "X" as row names
metadata <- subset(metadata, select = -1 ) # Remove column "X"
metadata$Concentration_f <- as.factor(metadata$Concentration)
metadata$Time_added_f <- as.factor(metadata$Time_added)

metadata.phylo <- sample_data(metadata, errorIfNULL = T)
metadata.phylo

# Loop to make counts table
file.list <- list.files(pattern = "*_output.txt") # Make a vector of relevant files names
taxa <- tax.table.all$X # Make vector "taxa" with list of taxa names
counts.table <- data.frame(taxa) # Make data frame from vector "taxa"
    
for(i in 1:length(file.list)){
    sample <- read.delim(file.list[i], h=F)
    sample.df <- data.frame(table(sample)) # Convert vector to data frame with counts per taxa
    sample.df <- sample.df[-which(sample.df$sample == "*"), ] # Remove first row, which is just "*"
    names(sample.df)[1] <- "taxa" # Change name of first column to "taxa"
    
    counts.table <- full_join(x = counts.table, y = sample.df, by = "taxa") # Transfer Freq column to counts.table
    names(counts.table)[i+1] <- substr(file.list[i], 1, 12) # Name column with the first 10 characters of file name
    }

return(counts.table)

write.csv(counts.table, "../bowtie_counts.csv", row.names=F)

# Import OTU data file.
OTU.table.abs <- read.csv("../bowtie_counts.csv", sep = ",")
OTU.matrix.abs <- data.matrix(OTU.table.abs[,-1]) # Convert to matrix. Leave first column (OTUs) out. 
rownames(OTU.matrix.abs) <- paste0(OTU.table.abs$taxa) # Paste feature IDs from CSV file as matrix row names
OTU.matrix.abs[is.na(OTU.matrix.abs)] <- 0 # Convert all NA values to 0
#head(OTU.matrix.abs)

# Importing OTU data as PhyloSeq object. 
OTU.phylo.abs = otu_table(OTU.matrix.abs, taxa_are_rows = T) # Convert from matrix to PhyloSeq OTU table

OTU.phylo.abs

physeq.abs <- phyloseq(OTU.phylo.abs, tax.phylo.all, metadata.phylo)
physeq.abs

### Noramlize for genome length

# Used following command to count # characters for each contig within .fna files used to make Bowtie2 index (Got this command off Stack overflow)
#awk '/^>/{if(NR>1){print l;} l=0; tmp=$0; next}{l+=length($0)}END{print l;}' filename.fna

# List of genome/contig sizes. Must be in order of OTU table.
genomes <- c(316264, 130867, 74540, 60174, 40920, 30891, 30824, 28427, 25311, 21563, 4941, 303791, 2649, 2450, 1386, 1324, 1160, 238118, 207960, 187988, 175212, 160455, 157637, 140763, 2560265, 2501097, 2821361, 2443604, 59661, 2499279, 4439, 4679, 8007, 17261, 24365, 6585, 2220494, 32498, 4439, 2427576, 1868883)

# Divide reads by length of genome/contig in kilobases -> RPK. Then divide each RPK by sum()/1,000,000 ('per kilobase scaling factor')
physeq_object <- physeq.abs

for (i in 1:nrow(otu_table(physeq_object))){
    otu_table(physeq_object)[i,] <- otu_table(physeq_object)[i,]/(genomes[i]/1000)
}

physeq_object <- transform_sample_counts(physeq_object, function(OTU) OTU/(sum(OTU)/1000000))

physeq.abs.norm <- physeq_object
                                                  
otu_table(physeq.abs.norm)

# Check columns sum to 1 million for RPKM values
sum(otu_table(physeq.abs.norm)[,1])
sum(otu_table(physeq.abs.norm)[,2])
sum(otu_table(physeq.abs.norm)[,3])

# Export RPKM-normalized table
write.csv(otu_table(physeq.abs.noneg), "../RPKM_counts.csv")

### "RPKM_counts.csv" was opened in Excel. RPKM average and standard deviation for each condition were calculated in Excel.
### This data was used to generate RPKM plots (code for plots available in skincom_scripts_alekbua_github.Rmd)
