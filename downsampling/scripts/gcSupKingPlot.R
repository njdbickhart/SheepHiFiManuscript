#!/usr/bin/env Rscript

library(ggridges)
library(ggplot2)
library(gridExtra)
library(dplyr)

outfile <- "figures/summary_taxonomic_plot.pdf"
# Should be alternating "names" and files
# ie. "ASM1" "asm1.blobtools.table.tab"
args = commandArgs(trailingOnly=TRUE)


# Input file tabs
# LEN GC KING
#alltabs <- list()
fulldata <- data.frame(LEN=NA, GC=NA, KING=NA, ASM=NA)[-1,]
fulldata$ASM <- as.character(fulldata$ASM)
fulldata$KING <- as.character(fulldata$KING)

for(l in seq(1, length(args))){
  if(l %% 2 == 0){
    next
  }
  print(args[l])
  df <- read.delim(args[l + 1], header = TRUE)
  df <- df %>% mutate(ASM = c(args[l]))
  df[df$KING == "undef","KING"] <- "no-hit"
  print(head(df))
  fulldata <- bind_rows(fulldata, df)
}

print(head(fulldata))

fulldata$ASM <- as.factor(fulldata$ASM)
fulldata$KING <- factor(fulldata$KING, levels = c("no-hit", "Viruses", "Eukaryota", "Archaea", "Bacteria"))

sumdata <- fulldata %>% group_by(ASM, KING) %>% summarize(SUM=sum(as.numeric(LEN)))
print(sumdata)

pdf(file=outfile, useDingbats=FALSE)
p3 <- ggplot(sumdata, aes(x=ASM, y=SUM, fill=ASM)) + geom_bar(stat="identity") + scale_color_brewer(palette = "Dark2") + theme_bw() + scale_fill_brewer(palette="Dark2") + theme(axis.title=element_blank(), axis.text.x = element_text(angle=45, hjust = 1)) + theme(legend.position = "none") + labs(title="Nucleotides per Kingdom") + facet_wrap(~KING, scales = "free")
p1 <- ggplot(fulldata, aes(y=KING, x=LEN, fill=ASM)) + geom_density_ridges(aes(alpha = 0.4)) + theme_bw() + theme(axis.title=element_text(size=11, face="bold"), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_brewer(palette="Dark2") + scale_x_log10(breaks=c(1000, 10000, 100000, 1000000, 6000000), limits=c(100, 6000000), labels=c("1000", "10,000", "100,000", "1,000,000", "6,000,000")) + xlab(label = "Log10 Contig Lengths (bp)") + ylab(label= "Superkingdom Taxonomic Assignment")
p2 <- ggplot(fulldata, aes(x=ASM, y=GC, fill=ASM)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, outlier.shape = NA, fill = "white") + theme_bw() + theme(axis.title=element_blank(), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + theme(legend.position = "none") + scale_fill_brewer(palette="Dark2") + labs(title="GC Content")
grid.arrange(grobs=list(p3, p2, p1), layout_matrix=rbind(c(1,1,2), c(1,1,2), c(3,3,3)))
dev.off()
