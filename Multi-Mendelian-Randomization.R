##### Multivariable MR #####
##### Multivariable MR #####
library(TwoSampleMR)
library(MendelianRandomization) 
library(MVMR)
library(data.table)
library(writexl)
library(readxl)
library(stringr)

pre_path =  "/mnt/data1/user/chendongze/Project/EA_CUD_Project/DataFormat_For_LDSC/disease_"

exposure_trait_1 = 'CUD'
exposure_trait_2= 'Intelligence'  
outcome_trait = '2016EA'  


exp_1 = fread(paste0(pre_path, exposure_trait_1, '.txt'))
exp_1$BETA = log(exp_1$OR)
exp_1_raw = exp_1
exp_2 = fread(paste0(pre_path, exposure_trait_2, '.txt'))
exp_2$BETA = log(exp_2$OR)
exp_2_raw = exp_2
outcome = fread(paste0(pre_path, outcome_trait, '.txt'))
outcome$BETA = log(outcome$OR)
outcome_raw = outcome

exp_1_ivs = fread(paste0('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/MDD_GI_Project/', exposure_trait_1, '.csv'))
exp_2_ivs = fread(paste0('/home/dongzechen/GWAS/04-twosampleMR/Manual_IVs/MDD_GI_Project/', exposure_trait_2, '.csv'))

ivs = Reduce(union,  list(exp_1_ivs$SNP, exp_2_ivs$SNP))  
ivs = data.frame(SNP=ivs); nrow(ivs)
ivs_independent = clump_data(ivs, clump_kb = 10000, clump_r2 = 0.001, pop = "EUR")
nrow(ivs_independent)


# check for proxy
common_snp = intersect(as.character(ivs_independent$SNP), exp_1_raw$SNP); length(common_snp)
need_proxy_snp = setdiff(as.character(ivs_independent$SNP), common_snp)
write_xlsx(data.frame(need_proxy_snp),'snps_need_proxy.xlsx')

temp_snp = fread('./snps.txt')
write.table(data.frame(intersect(temp_snp$SNP, exp_1_raw$SNP)),'snps_have_proxy.txt',col.names=TRUE,row.names = FALSE, sep=",",quote=FALSE)

proxy_snp_1 = c('rs4356032'='rs4261101', 'rs3852788'='rs11643192','rs7231748'='rs12958048','rs3742786'='rs200343480',
				'rs10809510'='rs10959913','rs1549212'='rs4869056','rs3821709'='rs16851483','rs9910745'='rs12940622',
				'rs10483389'='rs11847697','rs13180906'='rs116755193','rs2163544'='rs62099069','rs4494619'='rs11663393',
				'rs12967135'='rs6567160','rs6752378'='rs10182181','rs10423928'='rs2287019','rs7248181'='rs3810291',
				'rs7859780'='rs4740619','rs9957264'='rs7243357','rs1412239'='rs10968576','rs2181375'='rs11165643',
				'rs8010984'='rs10132280','rs12548773'='rs17405819','rs6700902'='rs11583200','rs6714241'='rs11688816',
				'rs29942'='rs29941','rs2430307'='rs2245368','rs2192497'='rs1016287','rs1949772'='rs2365389',
				'rs1307813'='rs12885454','rs9579083'='rs12016871','rs12269171'='rs7899106') 
	
SNPS = c(common_snp, names(proxy_snp_1)); length(SNPS)
exp_1 = exp_1_raw[exp_1_raw$SNP %in% SNPS, ]
exp_1 = subset(exp_1, select=c('SNP', 'BETA', 'SE', 'P'))


if (length(proxy_snp) != 0){
	outcome_dat$SNP = str_replace_all(outcome_dat$SNP, proxy_snp)
	}

exp_1 = exp_1_raw[exp_1_raw$SNP %in% as.character(ivs_independent$SNP), ]
exp_1 = subset(exp_1, select=c('SNP', 'BETA', 'SE', 'P'))
exp_2 = exp_2_raw[exp_2_raw$SNP %in% as.character(ivs_independent$SNP), ]
exp_2 = subset(exp_2, select=c('SNP', 'BETA', 'SE', 'P'))

# exp_2 = extract_outcome_data(snps = as.character(ivs_independent$SNP), outcomes = 'ieu-b-35') 
# exp_2 = subset(exp_2,select=c('SNP','beta.outcome','se.outcome', 'pval.outcome'))
# names(exp_2) = c('SNP', 'BETA', 'SE', 'P')

outcome = outcome_raw[outcome_raw$SNP %in% as.character(ivs_independent$SNP), ]
outcome = subset(outcome, select=c('SNP', 'BETA', 'SE', 'P'))

dat = merge(exp_1, exp_2, by='SNP')
dat = merge(dat, outcome, by='SNP')
nrow(dat)  # common identify in exposure and outcome SNP
dat = dat[dat$P > 5e-8, ]   
nrow(dat)
Num_ivs = nrow(dat)  

trait_1_specific_num = length(intersect(dat$SNP, exp_1_ivs$SNP))
trait_2_specific_num = length(intersect(dat$SNP, exp_2_ivs$SNP))


rawdat_mvmr1 = subset(dat, select=c('BETA.x','BETA.y', 'SE.x', 'SE.y', 'BETA', 'SE', 'SNP'))
names(rawdat_mvmr1) = c('exp1_beta','exp2_beta', 'exp1_se', 'exp2_se', 'outcome_beta', 'outcome_se', 'SNP')
rawdat_mvmr1 = na.omit(rawdat_mvmr1)

write_xlsx(rawdat_mvmr1, paste0('./Project/MDD_GI_Project/',exposure_trait_1,'-',exposure_trait_2, '_to_', outcome_trait,'_IVs.xlsx'), col_names = TRUE)  # 保存工具变量

F.data <- format_mvmr(BXGs = rawdat_mvmr1[,c(1,2)],
					BYG = rawdat_mvmr1[,5],
					seBXGs = rawdat_mvmr1[,c(3,4)],
					seBYG = rawdat_mvmr1[,6],
					RSID = rawdat_mvmr1[,7])

# head(F.data)
# mvmrcovmatrix<-matrix(c(1, 0.45, 0.45,1), nrow = 2, ncol = 2)  
# Xcovmat<-phenocov_mvmr(mvmrcovmatrix,F.data[,3:4])

sres <- strength_mvmr(r_input = F.data, gencov = 0)
# sres <- strength_mvmr(r_input = F.data, gencov = Xcovmat)

pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
# pres <- pleiotropy_mvmr(r_input = F.data, gencov = Xcovmat)

res <- ivw_mvmr(r_input = F.data)
# res <- qhet_mvmr(F.data, mvmrcovmatrix, CI = T, iterations = 1000)  

out_df = cbind(data.frame(confounders_mediators=c(exposure_trait_1,exposure_trait_2), outcome = outcome_trait,  Raw_Trait_NumIVs=c(length(exp_1_ivs$SNP), length(exp_2_ivs$SNP)), 
				Final_Trait_NumIVs=c(trait_1_specific_num, trait_2_specific_num), NumIVs =Num_ivs), as.data.frame(res)[,c(1,2,4)],
				conditional_F_statistic=c(sres$exposure1, sres$exposure2), Qstat=c(pres$Qstat,pres$Qstat), Qpval=c(pres$Qpval, pres$Qpval))


if (file.exists('./MDD_GI_MVMR_Result.xlsx')){
	old_dat = read_excel('./MDD_GI_MVMR_Result.xlsx')
	new_dat = rbind(old_dat, "", out_df)
	write_xlsx(new_dat,'./MDD_GI_MVMR_Result.xlsx') 
} else {
	write_xlsx(out_df,'./MDD_GI_MVMR_Result.xlsx') 
} 
