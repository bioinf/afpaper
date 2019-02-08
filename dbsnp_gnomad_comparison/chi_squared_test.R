





df <- as_data_frame(fread('SnpEffects_correct.tsv', header = F))
summ <- df %>% group_by(V2, V3) %>% summarise(observed=n())
summ['expected'] = ifelse(summ['V3']=='known', 459495/508381, (508381-459495)/508381)
df_list <- split(summ, summ$V2)
stat_test <- lapply(df_list, function(x) chisq.test(x = x$observed, p = x$expected))
stat_test <- do.call('rbind', stat_test)
stat_data <- data.frame(stat_test) %>% select(statistic, parameter, p.value)
fwrite(stat_data, 'Chi_squared_test.tsv', row.names = T, sep='\t')
