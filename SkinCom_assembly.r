library(plyr)
library(dplyr)
library(phyloseq)

setwd("/home/SkinCom_expt1/bowtie_counts/")

options(repr.matrix.max.rows=600, repr.matrix.max.cols=200)

# Taxa table
# Import taxonomy data with all classification fields filled in.
tax.table.all <- read.csv("../phyloseq_taxa.csv", sep = ",")
tax.matrix.all <- as.matrix(tax.table.all)
tax.phylo.all <- tax_table(tax.matrix.all)
rownames(tax.phylo.all) <- paste0(tax.table.all$X) # Paste feature IDs from CSV file as matrix row names
tax.phylo.all <- subset(tax.phylo.all, select = -X ) # Remove extra column of row numbers
class(tax.phylo.all)
tax.phylo.all

# Read in metadata CSV.
metadata <- read.csv("../Pilot1_metadata_phyloseq.tsv", sep = "\t")
rownames(metadata) <- paste0(metadata[,1]) # Paste column "X" as row names
metadata <- subset(metadata, select = -1 ) # Remove column "X"

metadata.phylo <- sample_data(metadata, errorIfNULL = T)
head(metadata.phylo)

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
    names(counts.table)[i+1] <- substr(file.list[i], 1, 10) # Name column with the first 10 characters of file name
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

# Make PhyloSeq object.
physeq.abs = phyloseq(OTU.phylo.abs, tax.phylo.all, metadata.phylo)
physeq.abs

### Noramlize for genome length

# Used following command to count number of characters for each contig within .fna files used to make Bowtie2 index
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

# Check all columns sum to 1 million for RPKM values
sum(otu_table(physeq.abs.norm)[,1])
sum(otu_table(physeq.abs.norm)[,2])
sum(otu_table(physeq.abs.norm)[,3])

# Collapse to organism level
physeq.abs.norm.org <- tax_glom(physeq.abs.norm, taxrank=rank_names(physeq.abs.norm)[3])

# Pull out OTU table as matrix, convert to df, and save as CSV
OTU1 = as(otu_table(physeq.abs.norm), "matrix")
OTUdf = as.data.frame(OTU1)
head(OTUdf)

write.csv(OTUdf, "../bowtie2_reclass_counts_norm.csv")

library(ggplot2)
library(vegan)

setwd("/home/SkinCom_expt1/")

# Format absolute phyloseq object for diversity analysis. Remove negative control from object trimmed to species level
physeq.abs.noneg = subset_samples(physeq.abs, Condition != "Neg")
physeq.abs.norm.noneg = subset_samples(physeq.abs.norm, Condition != "Neg")

# Shannon diversity index
plot_richness(physeq.abs.noneg, x = "Condition", color = "Media", measures = "Shannon") +
    geom_point(size = 5) +
    scale_color_manual(values = c("lightgreen", "darkgreen")) +
    theme_bw() + ylab("Shannon diversity index") +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16),
         legend.text = element_text(size = 14), legend.title = element_text(size = 16),
         plot.title = element_blank())

# Data table of Shannon diversity index
alpha_diversity <- estimate_richness(physeq.abs.noneg, split = T, measures = "Shannon")

write.csv(alpha_diversity, "Shannon_diversity.csv")

# Bray-Curtis distance, all samples
bray_dist <- distance(physeq.abs.norm.noneg, method = "bray")
bray_ord <- ordinate(physeq.abs.norm.noneg, "MDS", distance = bray_dist)

shape_names <- c("square filled", "triangle filled", "circle filled", "triangle down filled", "diamond filled")

bray_plot <- plot_ordination(physeq.abs.norm.noneg, bray_ord, shape = "Condition", color = "Media") +
    scale_shape_manual("Condition", values = shape_names) +
    scale_color_manual("Media", values = c("darkgray", "black")) +
    scale_fill_manual("Media", values = c("darkgray", "black")) +
    geom_point(aes(fill = Media), size = 5) +
    theme_bw() + labs(x = "Axis 1 (61.1%)", y = "Axis 2 (31.7%)") +
    theme(axis.text = element_text(size = 14, color = "black", face = "bold"), 
          axis.title = element_text(size = 14, face = "bold"),
         legend.text = element_text(size = 14),
         legend.title = element_text(size = 14, face = "bold")) +
    theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1, "cm"),
         panel.border = element_rect(colour = "black", fill = NA, size = 1))

bray_plot

# PERMANOVA analysis to test for statistically significant differences in Bray-Curtis distance
all_permanova_media <- adonis(bray_dist ~ sample_data(physeq.abs.norm.noneg)$Media)
all_permanova_media

all_permanova_condition <- adonis(bray_dist ~ sample_data(physeq.abs.norm.noneg)$Condition + sample_data(physeq.abs.norm.noneg)$Media)
all_permanova_condition

# Beta diversity, Bray-Curtis, EM samples
# Filter for just EM data
beta_EM <- subset_samples(physeq.abs.norm.noneg, Condition == "EM")

# Beta diversity plot
bray_dist <- distance(beta_EM, method = "bray")
bray_ord <- ordinate(beta_EM, "MDS", distance = bray_dist)

bray_plot <- plot_ordination(beta_EM, bray_ord, shape = "Media") +
    geom_point(size = 5, color = "black") +
    scale_shape_manual("Media", values = c(15, 16)) +
    theme_bw() + labs(x = "Axis 1 (92.1%)", y = "Axis 2 (7.8%)") +
    theme(axis.text = element_text(size = 12, color = "black"), 
          axis.title = element_text(size = 14, face = "bold"),
         legend.text = element_text(size = 12),
         legend.title = element_text(size = 14, face = "bold")) +
    theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1, "cm"))

bray_plot

# PERMANOVA analysis for EM samples
EM_permanova <- adonis(bray_dist ~ sample_data(beta_EM)$Media)
EM_permanova
