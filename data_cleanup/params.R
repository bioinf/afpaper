setwd("/mnt/2EA01BDBA01BA7FB/Working issues/WES/AF_Paper/Revision/quality")

library(ggplot2)
library(reshape2)
library(cowplot)

dps = read.table('DP.tsv')
dps$type = c(rep('known', 311101), rep('novel', 21709))

p1 <- ggplot(dps, aes(x=type, y=V1, fill=type)) + geom_violin() +
  scale_y_continuous(limits=c(0, 50000)) + ylab('DP')

aggregate(V1~type, dps, median)

mqs = read.table('MQ.tsv')
mqs$type = c(rep('known', 311101), rep('novel', 21709))

p2 <- ggplot(mqs, aes(x=type, y=V1, fill=type)) + geom_violin(scale='width') +
  scale_y_continuous(limits=c(0, 100)) + ylab('MQ')

aggregate(V1~type, mqs, median)


vafs = read.table('VAF.tsv')
vafs$type = c(rep('known', 300000), rep('novel', 23500))

p3 <- ggplot(vafs, aes(x=type, y=V1, fill=type)) + geom_violin(scale='width') +
  scale_y_continuous(limits=c(0, 1)) + ylab('VAF')

aggregate(V1~type, vafs, median)

plot_grid(p1, p2, p3, nrow=1)

ggplot(vafs, aes(x=V1, fill=type)) + geom_density(alpha=0.4) +
  scale_x_continuous(limits=c(0, 1)) + xlab('VAF')

