if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!"TCGAbiolinks" %in% installed.packages())
  BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(survival)
library(survminer)


cancer_name <- "kirc"
setwd("D:/C_backup/Data/")
query <- GDCquery(project = paste0("TCGA-",toupper(cancer_name)),
                  data.category = "Transcriptome Profiling",
                  workflow.type = "HTSeq - FPKM-UQ",
                  sample.type = "Solid Tissue Normal")

GDCdownload(query)
data <- GDCprepare(query)
patients <- data$patient

expr <- data@assays@data$`HTSeq - FPKM-UQ`
expr_gene_list <- data@rowRanges$ensembl_gene_id

setwd("D:/Age_Project/MSB/code")

inc_gene <- read.table("./requirements/KIRC_inc.txt", header = FALSE)$V1
dec_gene <- read.table("./requirements/KIRC_dec.txt", header = FALSE)$V1

clinical <- read.table("./requirements/kirc_tcga_clinical_data.tsv", sep = "\t", header = TRUE)
clinical <- clinical[sapply(patients, function(arg) which(clinical$Patient.ID == arg)),]
age <- data$days_to_birth/-365.25
status <- factor(ifelse(is.na(data$days_to_death), 0, 1))
surv_time <- clinical$Overall.Survival..Months.




#### Cox regression ####
inc_gene_idx <- sapply(inc_gene, function(arg) which(expr_gene_list == arg))
inc_gexp <- expr[inc_gene_idx,]
inc_gexp <- apply(inc_gexp,1,function(e) (e-min(e))/(max(e)-min(e)))
inc_gexp <- apply(t(inc_gexp),2,mean)
inc_gexp <- scale(inc_gexp, center = T)


dec_gene_idx <- sapply(dec_gene, function(arg) which(expr_gene_list == arg))
dec_gexp <- expr[dec_gene_idx,]
dec_gexp <- apply(dec_gexp,1,function(e) (e-min(e))/(max(e)-min(e)))
dec_gexp <- apply(t(dec_gexp),2,mean)
dec_gexp <- scale(dec_gexp, center = T)


surv_df <- cbind(surv_time,status, age, inc_gexp, dec_gexp)
colnames(surv_df) <- c("time", "status", "age", "inc_gexp", "dec_gexp")
surv_df <- as.data.frame(surv_df)

# with age
cox_result_age <- coxph(Surv(time,status) ~ age, data = surv_df)
summary(cox_result_age)

# with increased genes
cox_result_inc <- coxph(Surv(time,status) ~ inc_gexp , data = surv_df)
summary(cox_result_inc)

# with decreased genes
cox_result_dec <- coxph(Surv(time,status) ~ dec_gexp , data = surv_df)
summary(cox_result_dec)


#### Kaplan-Meier estimator ####
inc_group <- ifelse(surv_df$inc_gexp >= median(surv_df$inc_gexp), "High","Low")
dec_group <- ifelse(surv_df$dec_gexp >= median(surv_df$dec_gexp), "High","Low")

surv_df <- data.frame(surv_df, inc_group = inc_group, dec_group = dec_group)

# with increased genes
survfit_inc <- survfit(Surv(time,status)~inc_group, data = surv_df)
ggsurvplot(survfit_inc, data = surv_df)

logrank_diff_inc <- survdiff(Surv(time,status)~inc_group, rho = 0, data = surv_df)
logrank_p_inc <- pchisq(logrank_diff_inc$chisq, length(logrank_diff_inc$n)-1, lower.tail = F)
print(logrank_p_inc)


# with decreased genes
survfit_dec <- survfit(Surv(time,status)~dec_group, data = surv_df)
ggsurvplot(survfit_dec, data = surv_df)

logrank_diff_dec <- survdiff(Surv(time,status)~dec_group, rho = 0, data = surv_df)
logrank_p_dec <- pchisq(logrank_diff_dec$chisq, length(logrank_diff_dec$n)-1, lower.tail = F)
print(logrank_p_dec)


