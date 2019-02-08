setwd("/mnt/2EA01BDBA01BA7FB/Working issues/WES/AF_Paper/var_counts")

library(ggplot2)
mycol1 = rgb(85, 138, 221, maxColorValue = 255)
mycol2 = rgb(255, 192, 78, maxColorValue = 255)
mycol3 = rgb(255, 121, 177, maxColorValue = 255)
mycol4 = rgb(221, 221, 221, maxColorValue = 255)
#mycol5 = rgb(100, 89, 89, maxColorValue = 255)
mycol6 = rgb(0, 189, 189, maxColorValue = 255)

counts = read.table('eff_final.tsv', sep='\t', header=T, stringsAsFactors = F)
ggplot(counts, aes(x=TYPE, y=COUNT, fill=TYPE, col=TYPE)) + geom_boxplot(lwd = 0.3) + 
  facet_wrap(~EXPTYPE) + theme_bw() +
  scale_fill_manual(values = c(mycol1, mycol2, mycol3, mycol6)) +
  scale_colour_manual(values = c(mycol1, mycol2, mycol3, mycol6)) +
  scale_y_continuous(limits = c(0, 16000)) + 
  theme(axis.text.x=element_blank())
        

aggregate(COUNT ~ TYPE + EXPTYPE, counts, median)
