setwd("/mnt/2EA01BDBA01BA7FB/Working issues/WES/AF_Paper/QC3/vcfResult")
library(ggplot2)

data = read.csv('selected_samples.A4.new.vcf.Method1.txt', sep='\t', header=T)
ggplot(data, aes(x=Overall.Consistency)) + geom_histogram(fill='red', col='black') +
  xlab('Heterozygous consistency (total)') + theme_bw()

ggplot(data, aes(x=Heterozygous.Consistency..CountB2A.CountA.)) + geom_histogram(fill='red', col='black') +
  xlab('Heterozygous consistency (one-way)') + ylab('log10 (pair count)') +
  theme_bw() + scale_y_log10()

data$Honest = (data$Heterozygous.Consistency..CountA2B.CountB. + 
                 data$Heterozygous.Consistency..CountB2A.CountA.) / 2

ggplot(data, aes(x=Honest)) + geom_histogram(fill='red', col='black') +
  xlab('Heterozygous consistency (mean)') + ylab('log10 (pair count)') +
  theme_bw() + scale_y_log10()

data_strange = data[data$Heterozygous.Consistency..CountB2A.CountA. >= 0.8, ]
data_good = data[!(data$SampleA %in% data_strange$SampleA), ]

write.table(unique(data_good[, 2]), file = 'clear.txt', col.names = F, 
            row.names = F, quote=F)

good_samps = read.table('final_ds.txt', header=F)
head(good_samps)

data_clean = data[data$SampleA %in% good_samps$V1 & data$SampleB %in% good_samps$V1, ]

ggplot(data_clean, aes(x=Honest)) + geom_histogram(fill='red', col='black') +
  xlab('Heterozygous consistency (mean)') + ylab('log10 (pair count)') +
  theme_bw()

mean(data_clean$Honest)

exit = FALSE
while (TRUE) {
  bad_samps = unique(data_clean[data_clean$Honest > 0.5, 'SampleA'])
  if (length(bad_samps) == 0) {
    break
  }
  tgt = as.character(bad_samps[1])
  data_clean = data_clean[data_clean$SampleA != tgt & data_clean$SampleB != tgt, ]
}

ggplot(data_clean, aes(x=Honest)) + geom_histogram(fill='red', col='black') +
  xlab('Heterozygous consistency (mean)') + ylab('log10 (pair count)') +
  theme_bw()

length(unique(data_clean$SampleA))

data_all = rbind(data, data_clean)
data_all$type = c(rep('all', nrow(data)), rep('clean', nrow(data_clean)))
ggplot(data_all, aes(x=Honest, fill=type)) + geom_histogram(col='black') +
  xlab('Heterozygous consistency (mean)') + ylab('log10 (pair count)') +
  theme_bw() + facet_wrap(~type, nrow=2) + scale_y_log10()

write.table(unique(data_clean[, 2]), file = 'cleaned_final_05.txt', col.names = F, 
            row.names = F, quote=F)
