setwd("/mnt/2EA01BDBA01BA7FB/Working issues/WES/AF_Paper/Revision")

library(ggplot2)
library(reshape2)

af_data = read.table('AF.tsv', header=T, sep='\t', stringsAsFactors = F)
head(af_data)
tt = (af_data$AF - af_data$AF_afr)^2
mean(tt)

pops = colnames(af_data)[5:10]
for (pop in pops) {
  af_data[, paste0('diff_', pop)] = sqrt((af_data$AF - af_data[, pop])^2)
}

tpl = melt(af_data[, 12:17])
ggplot(tpl, aes(x=value, fill=variable)) + geom_histogram(col='black') + 
  facet_wrap(~variable, ncol = 1) + 
  theme_bw()

m_data = data.frame(pops = colnames(af_data)[5:10],
                    F2 = colMeans(af_data[, 12:17]),
                    error = apply(af_data[, 12:17], 2, function(x) sd(x)/sqrt(nrow(af_data))))

ggplot(m_data, aes(x=pops, y=F2, fill=pops, ymin=F2-2*error, ymax=F2+2*error)) + 
  geom_bar(stat='identity', col='black') +
  geom_errorbar(width=0.6) +
  scale_y_continuous(limits=c(0, 0.15)) + 
  theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1))

pcaout = prcomp(as.matrix(t(af_data[, c(4:7, 9, 10)])))
imps = summary(pcaout)$importance

pc12 = as.data.frame(pcaout$x[, 1:2])
ggplot(pc12, aes(x=-PC1, y=PC2, col=rownames(pc12))) + geom_point(size=4) +
  theme_bw() + xlab('PC1 - 42% of variance') + ylab('PC2 - 36% of variance')

pc23 = as.data.frame(pcaout$x[, 2:3])
ggplot(pc23, aes(x=-PC2, y=PC3, col=rownames(pc23))) + geom_point(size=4) +
  theme_bw() + xlab('PC2 - 36% of variance') + ylab('PC3 - 10% of variance')

# Same on variants observed in all populations

af_data_common = af_data[rowSums(af_data[, c(4:7, 9, 10)] > 0) == 6, ]
pcaout = prcomp(as.matrix(t(af_data_common[, c(4:7, 9, 10)])))
imps = summary(pcaout)$importance

pc12 = as.data.frame(pcaout$x[, 1:2])
ggplot(pc12, aes(x=-PC1, y=PC2, col=rownames(pc12))) + geom_point(size=4) +
  theme_bw() + xlab('PC1 - 43% of variance') + ylab('PC2 - 37% of variance')

pops = colnames(af_data_common)[5:10]
for (pop in pops) {
  af_data_common[, paste0('diff_', pop)] = sqrt((af_data_common$AF - af_data_common[, pop])^2)
}

m_data = data.frame(pops = colnames(af_data_common)[5:10],
                    F2 = colMeans(af_data_common[, 12:17]),
                    error = apply(af_data_common[, 12:17], 2, function(x) sd(x)/sqrt(nrow(af_data_common))))

ggplot(m_data, aes(x=pops, y=F2, fill=pops, ymin=F2-2*error, ymax=F2+2*error)) + 
  geom_bar(stat='identity', col='black') +
  geom_errorbar(width=0.6) +
  scale_y_continuous(limits=c(0, 0.1)) + 
  theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1))
