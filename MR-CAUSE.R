##### MR-CAUSE #####
##### MR-CAUSE #####

## Step 1: Format Data for CAUSE
library(data.table)
library(readr)
library(dplyr)
library(cause)

X1 <- fread('/mnt/data1/user/chendongze/Project/EA_CUD_Project/DataFormat_For_LDSC/disease_CUD.txt')  # exposure
X2 <- fread('/mnt/data1/user/chendongze/Project/EA_CUD_Project/DataFormat_For_LDSC/disease_2016EA.txt')   # outcome

X1 = data.frame(X1)
X2 = data.frame(X2)
X1$BETA = log(X1$OR)
X2$BETA = log(X2$OR)

X <- gwas_merge(X1, X2, snp_name_cols = c("SNP", "SNP"), 
                       beta_hat_cols = c("BETA", "BETA"), 
                       se_cols = c("SE", "SE"), 
                       A1_cols = c("A1", "A1"), 
                       A2_cols = c("A2", "A2"),
					   pval_cols = c("P", "P"))
head(X)

## Step 2: Calculate nuisance parameters
set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)   

## Step 3: LD Pruning
r2_thresh = 0.01
pval_thresh = 1e-5
X_clump <- X %>% rename(rsid = snp, pval = p1) %>% ieugwasr::ld_clump(dat = ., clump_r2 = r2_thresh, clump_p = pval_thresh, plink_bin = genetics.binaRies::get_plink_binary(), pop = "EUR")
top_vars <- X_clump$rsid
print(length(top_vars))

## Step 4: Fit CAUSE
res <- cause(X=X, variants = top_vars, param_ests = params)  
res$loos[[2]]
loo::pareto_k_table(res$loos[[2]])
res$loos[[3]]
loo::pareto_k_table(res$loos[[3]])

## Step 5: Look at Results
class(res)
names(res)
res$elpd
summary(res, ci_size=0.95)















