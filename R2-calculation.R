library(data.table)
library(writexl)
direction = 'CUD_to_2016EA'

directions_list = list.files(getwd())
tmp <- strsplit(directions_list, split=".", fixed=TRUE)
directions_list <- unlist(lapply(tmp, head, 1)) 

for (direction in directions_list){
	print(direction)
	dat = read.csv(paste0('./', direction, '.csv'))
	my_dat = subset(dat, select=c('SNP','effect_allele.exposure','other_allele.exposure','eaf.exposure','beta.exposure','se.exposure','pval.exposure', 'beta.outcome','se.outcome', 'pval.outcome', 'samplesize.exposure'))
	names(my_dat) = c('SNP','A1','A2','EAF','beta.exposure','se.exposure','pval.exposure', 'beta.outcome','se.outcome', 'pval.outcome', 'N')
	y = my_dat
	write_xlsx(y, paste0('./', direction, '.xlsx'))
}

y$MAF = ifelse(y$EAF>=0.5, 1-y$EAF, y$EAF)
y$r2_2 = 2 * (abs(y$beta.exposure)) ** 2 * y$MAF * (1 - y$MAF) 
y$Fstatistic = (y$beta.exposure/y$se.exposure)^2
n = y$N
k = nrow(y)
y$sum_r2 = sum(y$r2)
y$Fstatistic_overall = ((n - k - 1)/k) * (y$sum_r2/(1-y$sum_r2))
write_xlsx(y, paste0('./', direction, '.xlsx'))