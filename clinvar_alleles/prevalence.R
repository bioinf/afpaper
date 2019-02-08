setwd("/mnt/2EA01BDBA01BA7FB/Working issues/WES/AF_Paper/gnomad_cleaned_corr")

library(ggplot2)

data = read.table('known_pathogenic_GNOMAD.tsv', sep='\t', header=T, quote="", 
                  stringsAsFactors = F)
head(data)
data_clean = data[sapply(data$AF, function(x) !(grepl(',', x))), ]

data_new = data_clean[data_clean$HOM == 0 & as.numeric(data_clean$AF) < 7, ]
counts = by(as.numeric(data_new$AF), data_new$GENE_NAME, sum)
ans = by(as.numeric(data_new$AN), data_new$GENE_NAME, mean)
str(counts)

totals = as.data.frame(array(counts, dim(counts), dimnames(counts)))
totals$AN = ans
colnames(totals) = c('count', 'AN')
totals$prevalence = (totals$count / totals$AN)^2

head(totals)

frequent = totals[totals$count > 3, ]
frequent

# On a larger sample

data = read.table('../gnomad_all/known_pathogenic_GNOMAD.tsv', sep='\t', header=T, quote="", 
                  stringsAsFactors = F)
head(data)
data_clean = data[sapply(data$AF, function(x) !(grepl(',', x))), ]

data_new = data_clean[data_clean$HOM == 0 & as.numeric(data_clean$AF) < 13, ]
counts = by(as.numeric(data_new$AF), data_new$GENE_NAME, sum)
ans = by(as.numeric(data_new$AN), data_new$GENE_NAME, mean)
str(counts)

totals = as.data.frame(array(counts, dim(counts), dimnames(counts)))
totals$AN = ans
colnames(totals) = c('count', 'AN')
totals$prevalence = (totals$count / totals$AN)^2

head(totals)

frequent = totals[totals$count > 3, ]
frequent
