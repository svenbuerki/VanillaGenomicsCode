setwd("~/Vanilla_Research/Genomic_heterozygosity_project/SNP_analysis")

#install and load packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install(c("gdsfmt", "SNPRelate"))

library(gdsfmt)
library(SNPRelate)


library(R.utils)


#unzip vcf file
gunzip("calls.vcf.gz", remove = FALSE)

#make an object out of vcf file
calls <- "calls.vcf"

# turn the VCF file into a less data intensive form (GDS) for easier computing
snpgdsVCF2GDS(calls,"calls.gds", method ="biallelic.only")

#look at summary of gds file
snpgdsSummary("calls.gds")

#open file
genofile <- snpgdsOpen("calls.gds")

#prune snp set to include 
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only = FALSE)

#look at snpset attributes
names(snpset)
head(snpset$chrCM028150.1) 

#get snp ids
snpset.id <- unlist(snpset)

pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = FALSE)
# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

pdf("snp_pca.pdf")
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1", pch = 19)
dev.off()




snps<-read.table("snp132_ucsc_hg19.bed",sep="\t",header=F)
colnames(snps)<-c("chr","start","end","id","score","strand")

# put the chromosomes in the good order: chr1, chr2, chr22, chrX
goodChrOrder <- paste("chr",c(1:22,"X","Y"),sep="")
snps$chr <- factor(snps$chr,levels=goodChrOrder)

# Plot the densities of snps in the bed file for each chr seperately
library(ggplot2)
snpDensity<-ggplot(snps) + 
  geom_histogram(aes(x=start),binwidth=1e6) + # pick a binwidth that is not too small 
  facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr
  ggtitle("Density of SNP-132 across hg19") +
  xlab("Position in the genome") + 
  ylab("SNP density") + 
  theme_bw() # I prefer the black and white theme

# save the plot to .png file
png("snp132_density.png",500,750)
print(snpDensity)
dev.off()
