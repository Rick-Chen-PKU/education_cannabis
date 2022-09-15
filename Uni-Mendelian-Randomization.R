######  Univariate Mendelian randomization analysis  ######
library(TwoSampleMR)
library(data.table)
library(MendelianRandomization) 
library(writexl)
library(readxl)
library(simex)
library(MRPRESSO)
library(mr.raps)
library(stringr)


path = "/mnt/data1/user/chendongze/Project/EA_CUD_Project/DataFormat_For_LDSC/disease_"

Two_sample_MR = function(exposure_list, outcome_trait){ 
	# Read outcome data
	disease_path = paste0( path, outcome_trait,c('.txt'))
	print("Reading outcome data")
	outcome = fread(disease_path)   
	outcome$BETA=log(as.numeric(outcome$OR))
	outcome$Trait = outcome_trait
	outcome1 = outcome
	
	for (exposure_trait in exposure_list){
		# Read exposure data
		exposure = fread(paste0("./Manual_IVs/", exposure_trait, ".csv")) 
		if ('EAF' %in% colnames(exposure)){
			exp_dat = read_exposure_data(
				filename = paste0("./Manual_IVs/", exposure_trait, ".csv"),
				sep = ",",
				phenotype_col = "Trait",
				snp_col = "SNP",
				eaf_col = "EAF",
				beta_col = "BETA",
				se_col = "SE",
				effect_allele_col = "A1",
				other_allele_col = "A2",
				pval_col = "P",
				samplesize_col = "N")    
		} else {
			exp_dat = read_exposure_data(
				filename = paste0("./Manual_IVs/", exposure_trait, ".csv"),
				sep = ",",
				phenotype_col = "Trait",
				snp_col = "SNP",
				beta_col = "BETA",
				se_col = "SE",
				effect_allele_col = "A1",
				other_allele_col = "A2",
				pval_col = "P",
				samplesize_col = "N")
		}

		# Check if need proxy SNP
		outcome = outcome1
		common_snp = intersect(exp_dat$SNP, outcome$SNP); length(common_snp)
		need_proxy_snp = setdiff(exp_dat$SNP, common_snp)
		print(paste0("Needed SNPs: ", paste(need_proxy_snp,collapse=", ")))
		
		if (length(need_proxy_snp) != 0){
			# proxy_snp = c("rs12517460"='rs12187824' , "rs4502776"='rs34426618') # for intelligence → EA  
			# proxy_snp = c("rs2014265"='rs5886709')  # for ADHD → EA
			# proxy_snp = c("rs55964653"='rs28813180')  # for CigDay → EA
			# proxy_snp = c("rs17564493"='rs192818565', "rs7140993" = "rs8005528")  # for EA → CUD
			# proxy_snp = c("rs62252499"='rs113642272', "rs144447022" = "rs3130834")  # for EverTakenCannabis → EA
			# proxy_snp = c("rs62056778"='rs192818565')  # for EA → 2019CUD
			proxy_snp = c("rs62056778"='rs192818565')  # for EA → EverTakenCannabi
		} else {
			proxy_snp = c()
		}
		
		SNPS = c(exp_dat$SNP, names(proxy_snp))
		if ('EAF' %in% colnames(outcome)){
			outcome = subset(outcome, select=c('SNP','A1', 'A2', 'BETA', 'SE', 'P', 'N', 'EAF', 'Trait'))
			print("Writing the outcome to csv")
			outcome = outcome[outcome$SNP %in% SNPS, ]
			write.csv(outcome,"./outcome.csv",row.names = FALSE,quote = FALSE)   
			print("Clumping the outcome, take some time")
			outcome_dat <- read_outcome_data(
			snps = SNPS,
			filename = "./outcome.csv",
			sep = ",",
			phenotype_col= "Trait",
			snp_col = "SNP",
			eaf_col = "EAF",
			beta_col = "BETA",
			se_col = "SE",
			effect_allele_col = "A1",
			other_allele_col = "A2",
			pval_col = "P",
			samplesize_col = "N")
		} else {
			outcome = subset(outcome, select=c('SNP','A1', 'A2', 'BETA', 'SE', 'P', 'N', 'Trait'))
			print("Writing the outcome to csv")
			outcome = outcome[outcome$SNP %in% SNPS, ]
			write.csv(outcome,"./outcome.csv",row.names = FALSE,quote = FALSE)   
			print("Clumping the outcome, take some time")
			outcome_dat <- read_outcome_data(
			snps = SNPS,
			filename = "./outcome.csv",
			sep = ",",
			phenotype_col= "Trait",
			snp_col = "SNP",
			beta_col = "BETA",
			se_col = "SE",
			effect_allele_col = "A1",
			other_allele_col = "A2",
			pval_col = "P",
			samplesize_col = "N")
		}
		
		if (length(proxy_snp) != 0){
			outcome_dat$SNP = str_replace_all(outcome_dat$SNP, proxy_snp)
		}
		
		dat = harmonise_data(exp_dat, outcome_dat, action=2); nrow(dat) 
		dat = dat[dat$ambiguous == 'FALSE', ]; nrow(dat)  # remove ambiguous SNPs
		dat = steiger_filtering(dat) # steiger_filtering 
		dat = dat[dat$steiger_dir == 'TRUE', ]; nrow(dat)
		write.table(dat, paste0('./Project/EA_CUD_Project/', exposure_trait, '_to_', outcome_trait, '_instruments.csv'), col.names=TRUE,row.names = FALSE, sep=",",quote=FALSE)
		
		# Read IVs
		dat = data.frame(fread(paste0('./Project/EA_CUD_Project/', exposure_trait, '_to_', outcome_trait, '_instruments.csv'))) 
		# potential_snp = c('rs11783093', 'rs1392816', 'rs7783012')  # for CUD exposure
		# potential_snp = c('rs2875907', 'rs9919557', 'rs10499')  # for LCU exposure
		# potential_snp = c('rs10061788','rs1008078','rs11191193','rs112634398','rs11690172','rs11712056','rs12646808','rs12671937',
						# 'rs12682297','rs12772375','rs12969294','rs12987662','rs13294439','rs13402908','rs1402025','rs148734725','rs17167170',
						# 'rs192818565','rs2431108','rs2456973','rs2837992','rs324886','rs34305371','rs35761247','rs4500960','rs572016',
						# 'rs62263923','rs7306755','rs7767938','rs895606','rs9320913','rs9537821')  # for 2016EA exposure
		# potential_snp = c('rs10061788','rs1008078','rs11191193','rs112634398','rs12646808','rs12671937',
						# 'rs12682297','rs12772375','rs12969294','rs13402908','rs17167170',
						# 'rs192818565','rs2456973','rs324886','rs34305371','rs4500960',
						# 'rs62263923','rs7306755','rs895606','rs9320913')  # for 2016EA exposure				
		# potential_snp = c('rs1004787','rs10750025','rs11030084','rs11692435','rs11940694','rs1229984','rs1260326','rs13024996','rs13107325',
							# 'rs1713676','rs17177078','rs2472297','rs281379','rs28680958','rs28929474','rs3748034','rs3803800','rs4092465',
							# 'rs4916723','rs6460047','rs6787172','rs6951574','rs77165542','rs823114','rs828867')
		# dat =  dat[dat$SNP %in% setdiff(dat$SNP,potential_snp), ]; nrow(dat)
		MRInputObject <- mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse =  dat$se.outcome, 
		exposure = exposure_trait, outcome = outcome_trait, snps = dat$SNP, effect_allele = dat$effect_allele.exposure, other_allele = dat$other_allele.exposure, eaf = dat$eaf.exposure)
		# IVW
		IVWObject1 <- MendelianRandomization::mr_ivw(MRInputObject,model= "default",robust = FALSE, penalized = FALSE, correl = FALSE, weights ="simple", psi = 0,distribution = "normal",alpha = 0.05)
		# IVWObject1 
		vec1 = c(IVWObject1@Exposure, IVWObject1@Outcome, IVWObject1@SNPs, 'IVW raw', round(IVWObject1@Estimate,3), round(exp(IVWObject1@Estimate),3), 
				paste0("(", as.character(round(IVWObject1@CILower,3)),",", as.character(round(IVWObject1@CIUpper,3)), ")"),
				paste0("(", as.character(round(exp(IVWObject1@CILower),3)),",", as.character(round(exp(IVWObject1@CIUpper),3)), ")"), round(IVWObject1@StdError,3), IVWObject1@Pvalue,  round(IVWObject1@Heter.Stat[1],2), IVWObject1@Heter.Stat[2],
				round((IVWObject1@Heter.Stat[1] - IVWObject1@SNPs + 1 )/IVWObject1@Heter.Stat[1], 3) )
		IVWObject2 <- MendelianRandomization::mr_ivw(MRInputObject,model= "default",robust = TRUE, penalized = TRUE, correl = FALSE, weights ="simple", psi = 0,distribution = "normal",alpha = 0.05)
		# IVWObject2 
		vec2 = 	c(IVWObject2@Exposure, IVWObject2@Outcome, IVWObject2@SNPs, 'IVW robust', round(IVWObject2@Estimate,3), round(exp(IVWObject2@Estimate),3),
				paste0("(", as.character(round(IVWObject2@CILower,3)),",", as.character(round(IVWObject2@CIUpper,3)), ")"),
				paste0("(", as.character(round(exp(IVWObject2@CILower),3)),",", as.character(round(exp(IVWObject2@CIUpper),3)), ")"), round(IVWObject2@StdError,3), IVWObject2@Pvalue, NA, NA, NA )
		
		# Median-based
		MedianObject1 <-MendelianRandomization::mr_median(MRInputObject,weighting = "simple",distribution ="normal",alpha = 0.05,iterations = 10000,seed = 314159265)  
		# MedianObject1 
		vec3 = 	c(MedianObject1@Exposure, MedianObject1@Outcome, MedianObject1@SNPs, 'Simple median', round(MedianObject1@Estimate,3), round(exp(MedianObject1@Estimate),3),
				paste0("(", as.character(round(MedianObject1@CILower,3)),",", as.character(round(MedianObject1@CIUpper,3)), ")"),
				paste0("(", as.character(round(exp(MedianObject1@CILower),3)),",", as.character(round(exp(MedianObject1@CIUpper),3)), ")"), round(MedianObject1@StdError,3), MedianObject1@Pvalue, NA, NA, NA )
		MedianObject2 <-MendelianRandomization::mr_median(MRInputObject,weighting = "weighted",distribution ="normal",alpha = 0.05,iterations = 10000,seed = 314159265) 
		# MedianObject2
		vec4 = 	c(MedianObject2@Exposure, MedianObject2@Outcome, MedianObject2@SNPs, 'Weighted median', round(MedianObject2@Estimate,3), round(exp(MedianObject2@Estimate),3),
				paste0("(", as.character(round(MedianObject2@CILower,3)),",", as.character(round(MedianObject2@CIUpper,3)), ")"),
				paste0("(", as.character(round(exp(MedianObject2@CILower),3)),",", as.character(round(exp(MedianObject2@CIUpper),3)), ")"), round(MedianObject2@StdError,3), MedianObject2@Pvalue, NA, NA, NA )
		
		# Mode-based
		MBEObject <- mr_mbe(MRInputObject, weighting = "weighted", stderror = "delta", phi = 1, seed = 314159265, iterations = 10000, distribution = "normal", alpha = 0.05)
		# MBEObject
		vec5 = 	c(MBEObject@Exposure, MBEObject@Outcome, MBEObject@SNPs, 'Weighted mode', round(MBEObject@Estimate,3), round(exp(MBEObject@Estimate),3),
				paste0("(", as.character(round(MBEObject@CILower,3)),",", as.character(round(MBEObject@CIUpper,3)), ")"), 
				paste0("(", as.character(round(exp(MBEObject@CILower),3)),",", as.character(round(exp(MBEObject@CIUpper),3)), ")"), round(MBEObject@StdError,3), MBEObject@Pvalue, NA, NA, NA )
		
		# Debiased inverse-variance weighted method
		DIVWObject = mr_divw(MRInputObject, over.dispersion = TRUE, alpha = 0.05, diagnostics = FALSE)
		# DIVWObject
		vec6 = 	c(DIVWObject@Exposure, DIVWObject@Outcome, DIVWObject@SNPs, 'DIVW', round(DIVWObject@Estimate,3), round(exp(DIVWObject@Estimate),3),
				paste0("(", as.character(round(DIVWObject@CILower,3)),",", as.character(round(DIVWObject@CIUpper,3)), ")"),
				paste0("(", as.character(round(exp(DIVWObject@CILower),3)),",", as.character(round(exp(DIVWObject@CIUpper),3)), ")"), round(DIVWObject@StdError,3), DIVWObject@Pvalue, NA, NA, NA )
		# MR-RAPS & MR-PRESSO
		res2 = mr.raps.mle(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, over.dispersion=TRUE, diagnostics=FALSE)
		res3 = mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", 
						OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 2000,  SignifThreshold = 0.05)
		vec7 = c(exposure_trait, outcome_trait, IVWObject1@SNPs, 'MR-RAPS', round(res2$beta.hat,3), round(exp(res2$beta.hat),3), 
					paste0('(', as.character(round(res2$beta.hat-1.96*res2$beta.se, 2)) , ',', as.character(round(res2$beta.hat+1.96*res2$beta.se, 2)), ')'), 
					paste0('(', as.character(round(exp(res2$beta.hat-1.96*res2$beta.se), 2)) , ',', as.character(round(exp(res2$beta.hat+1.96*res2$beta.se), 2)), ')'), 
					round(res2$beta.se,3), as.numeric(res2$beta.p.value), NA, NA, NA)
	
		vec81 = c(exposure_trait, outcome_trait, IVWObject1@SNPs,'MR-PRESSO:raw', round(res3$`Main MR results`$`Causal Estimate`[1],3), round(exp(res3$`Main MR results`$`Causal Estimate`[1]),2), 
					paste0('(', as.character(round(res3$`Main MR results`$`Causal Estimate`[1]-1.96*res3$`Main MR results`$`Sd`[1], 2)) , ',', as.character(round(res3$`Main MR results`$`Causal Estimate`[1]+1.96*res3$`Main MR results`$`Sd`[1], 2)), ')'), 
					paste0('(', as.character(round(exp(res3$`Main MR results`$`Causal Estimate`[1]-1.96*res3$`Main MR results`$`Sd`[1]), 2)) , ',', as.character(round(exp(res3$`Main MR results`$`Causal Estimate`[1]+1.96*res3$`Main MR results`$`Sd`[1]), 2)), ')'), 
					round(res3$`Main MR results`$`Sd`[1], 3), res3$`Main MR results`$`P-value`[1], NA, NA, NA)
					
		vec82 = c(exposure_trait, outcome_trait, IVWObject1@SNPs - length(res3$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),'MR-PRESSO:Outlier-corrected', 
					round(res3$`Main MR results`$`Causal Estimate`[2],3), round(exp(res3$`Main MR results`$`Causal Estimate`[2]),2), 
					paste0('(', as.character(round(res3$`Main MR results`$`Causal Estimate`[2]-1.96*res3$`Main MR results`$`Sd`[2], 2)) , ',', as.character(round(res3$`Main MR results`$`Causal Estimate`[2]+1.96*res3$`Main MR results`$`Sd`[2], 2)), ')'), 
					paste0('(', as.character(round(exp(res3$`Main MR results`$`Causal Estimate`[2]-1.96*res3$`Main MR results`$`Sd`[2]), 2)) , ',', as.character(round(exp(res3$`Main MR results`$`Causal Estimate`[2]+1.96*res3$`Main MR results`$`Sd`[2]), 2)), ')'), 
					round(res3$`Main MR results`$`Sd`[2], 3), res3$`Main MR results`$`P-value`[2], NA, NA, NA)
					
		
		out_df = rbind(vec1, vec2, vec3, vec4, vec5, vec6, vec7, vec81, vec82)
		out_df = data.frame(out_df)
		names(out_df) = c('Exposure', 'Outcome', 'N_snp', 'Method', 'beta', 'OR', 'beta_CI', 'OR_CI', 'SE', 'p_value', 'Q_Stat', 'Q_pvalue', 'I_square')
		
		# Sensitivity analyses
		sen1 = mr_heterogeneity(dat, method_list=c("mr_egger_regression", "mr_ivw_fe"))   # heterogeneity test
		sen2 = mr_pleiotropy_test(dat)  # Horizontal pleiotropy test, egger_intercept test
		out_df$egger_intercept_p_value = sen2$pval
	
		
		if (file.exists('./RemoveSNP_MR_Result.xlsx')){
            old_dat = read_excel('./RemoveSNP_MR_Result.xlsx')
            new_dat = rbind(old_dat,"",out_df)
            write_xlsx(new_dat,'./RemoveSNP_MR_Result.xlsx') 
        } else {
            write_xlsx(out_df,'./RemoveSNP_MR_Result.xlsx') 
        } 

		}
}
	
exposure_list  = c('CUD', 'LCU')
outcome_list = c('2016EA')

for (outcome_trait in outcome_list){
	Two_sample_MR(exposure_list, outcome_trait)
}
