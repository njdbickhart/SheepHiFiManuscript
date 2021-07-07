#!/usr/bin/env Rscript

library(ggridges)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggridges)

# Should be alternating "names" and files
# ie. "ASM1" "asm1.blobtools.table.tab"
args = commandArgs(trailingOnly=TRUE)

first <- paste0("figures/", args[1], "_facet_cov_by_comp_scatter.pdf")
second <- paste0("figures/", args[1], "_top80_cov_by_comp.pdf")
third <- paste0("figures/", args[1], "_ridge_top95_coverage.pdf")

print(c(first, second, third))

# Input file tabs
# Contig, Coverage, Completeness, Assembly
#alltabs <- list()
fulldata <- data.frame(Contig=NA, Coverage=NA, Completeness=NA, Contamination=NA, Assembly=NA)[-1,]
fulldata$Contig <- as.character(fulldata$Contig)
fulldata$Assembly <- as.character(fulldata$Assembly)

for(l in seq(2, length(args))){
  if((l -1) %% 2 == 0){
    next
  }
  print(args[l])
  df <- read.delim(args[l + 1], header = TRUE)
  df <- df %>% mutate(Assembly = c(args[l]))
  print(head(df))
  fulldata <- bind_rows(fulldata, df)
}

fulldata$Assembly <- as.factor(fulldata$Assembly)

pdf(file=first, useDingbats = FALSE)
ggplot(data = fulldata, aes(x=Completeness, y=Coverage, colour=Assembly)) + geom_point() + scale_colour_brewer(palette='Set2') + theme_bw() + scale_y_log10() + facet_wrap(~ Assembly)
dev.off()
print("Finished first")

pdf(file=second, useDingbats = FALSE)
top80 <- fulldata[fulldata$Completeness >= 80,]
ggplot(data = top80, aes(x=Completeness, y=Coverage, colour=Assembly)) + geom_point() + scale_colour_brewer(palette='Set2') + theme_bw() + scale_y_log10() + facet_wrap(~ Assembly)
dev.off()
print("Finished second")

pdf(file=third, useDingbats = FALSE)
ggplot(data=top80[top80$Completeness >=95,], aes(y=Assembly, x=Coverage, fill=Assembly)) + geom_density_ridges() + scale_fill_brewer(palette="Dark2") + theme_bw() + scale_x_log10()
dev.off()
print("Done!")
