library(sleuth)
data_tab <- read.table("~/miniProject_elias_crum/sleuth_table.txt",header=TRUE,stringsAsFactors=FALSE)
data_tab
so <- sleuth_prep(data_tab)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
library(dplyr)
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval) 
sleuth_significant
write.table(dplyr::select(sleuth_significant, target_id, test_stat, pval, qval), 
            file="EF999921_sleuth_results.txt",quote = FALSE,row.names = FALSE)
