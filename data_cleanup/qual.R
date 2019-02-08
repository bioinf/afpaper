setwd("/mnt/2EA01BDBA01BA7FB/Working issues/WES/AF_Paper/qual_compare")

library(ggplot2)
library(hexbin)
library(colorRamps)
mypal = matlab.like2(10000)

known = read.table('known.lst', header=F)
novel = read.table('novel.lst', header=F)

all = data.frame(qual = as.numeric(rbind(known, novel)[,1]), 
                 type = c(rep('known', nrow(known)), rep('novel', nrow(novel))))
all$logqual = log10(all$qual)

ggplot(all, aes(type, logqual)) + geom_violin()

# After filtering

known = read.table('known_cleaned.lst', header=F)
novel = read.table('novel_cleaned.lst', header=F)
known$AF = as.numeric(read.table('known_af.lst', header=F)[,1])
novel$AF = as.numeric(read.table('novel_af.lst', header=F)[,1])

all = data.frame(qual = as.numeric(rbind(known, novel)[,1]), 
                 af = as.numeric(rbind(known, novel)[,2]), 
                 type = c(rep('known', nrow(known)), rep('novel', nrow(novel))))
all$logqual = log10(all$qual)

ggplot(all, aes(type, af)) + geom_violin(scale = 'width') + scale_y_log10()

ggplot(all, aes(type, logqual)) + geom_violin()

pdf('test_free.pdf', width = 7, height=3.5)
ggplot(all, aes(x=af, y=logqual)) + geom_hex(aes(fill=log(..count..))) + 
  scale_x_log10() + facet_wrap(~type, scales = 'free') +
  scale_fill_gradientn(colours = mypal)
dev.off()

known = all[1:420187, ]
novel = all[420188:nrow(all), ]


pdf('known.pdf', width = 4, height=3.5)
ggplot(known, aes(x=af, y=logqual)) + geom_hex(aes(fill=log(..count..))) + 
  scale_x_log10() +
  scale_fill_gradientn(colours = mypal)
dev.off()

pdf('novel.pdf', width = 4, height=3.5)
ggplot(novel, aes(x=af, y=logqual)) + geom_hex(aes(fill=log(..count..))) + 
  scale_x_log10() +
  scale_fill_gradientn(colours = mypal)
dev.off()

# QD score

known = read.table('known_qd.lst', header=F)
novel = read.table('novel_qd.lst', header=F)

all = data.frame(qd = as.numeric(rbind(known, novel)[,1]), 
                 type = c(rep('known', nrow(known)), rep('novel', nrow(novel))))

ggplot(all, aes(type, qd)) + geom_violin()
