if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!"TCGAbiolinks" %in% installed.packages())
  BiocManager::install("TCGAbiolinks")

protein_coding_ensembl <- read.table("./requirements/protein_coding_ensembl.txt", sep="\t", header=FALSE)

cancer_name <- "kirc"
query <- GDCquery(project = paste0("TCGA-",toupper(cancer_name)),
                  data.category = "Transcriptome Profiling",
                  workflow.type = "HTSeq - FPKM-UQ",
                  sample.type = "Solid Tissue Normal")

GDCdownload(query)
data <- GDCprepare(query)

patients <- data$patient
age <- data$days_to_birth/-365.25


expr <- data@assays@data$`HTSeq - FPKM-UQ`
expr_gene_list <- data@rowRanges$ensembl_gene_id
protein_coding_gene_idx <- which(expr_gene_list %in% protein_coding_ensembl)
expr_protein_coding <- expr[protein_coding_gene_idx,]
zero_expr_samples <- apply(expr_protein_coding, 1, function(expr) length(which(expr == 0)))

ifelse(length(patients) > 30/0.8,
       valid_idx <- which(zero_expr_samples < 30),
       valid_idx <- which(zero_expr_samples < 0.8*length(patients)))

expr_protein_coding <- expr_protein_coding[valid_idx,]



### Identifying aging-related gene ####
# Results may vary due to bootstrapping
aging_gene_idx <- NULL
for(i in 1:nrow(expr_protein_coding)){
  test <- cbind(expr_protein_coding[i,], age)
  test <- as.data.frame(test)
  colnames(test)[1:2] <- c("exp","age")
  model1 <- lm(formula = exp ~ age, data = test)
  lm_summary <- summary(model1)
  coeff_mat <- lm_summary$coefficients
  if((coeff_mat[2,4] != "NaN") & (coeff_mat[2,4] < 0.05)){
    aging_gene_idx <- c(aging_gene_idx, i)
  }
}

pca_dt <- prcomp(t(expr_protein_coding),
                 center = T,
                 retx =T)
pca_dt$x <- pca_dt$x[,1:10]
pca_dt$x <- pca_dt$x[,apply(pca_dt$x,2,function(e) cor.test(e,age)$p.value)>0.05]
ncol(pca_dt$x)

increase_age_gene_idx <- NULL
decrease_age_gene_idx <- NULL
pb <- progress_bar$new(total = length(aging_gene_idx))
for(i in aging_gene_idx){
  pb$tick()
  test <- cbind(expr_protein_coding[i,], age, pca_dt$x)
  test <- as.data.frame(test)
  colnames(test)[1:2] <- c("exp","age")
  model2 <- lm(formula = exp ~ ., data = test)
  lm_summary <- summary(model2)
  coeff_mat <- lm_summary$coefficients
  if((coeff_mat[2,4] != "NaN") & (coeff_mat[2,4] < 0.05)){
    if(coeff_mat[2,1] > 0){
      increase_age_gene_idx <- c(increase_age_gene_idx, i)
    }
    else{
      decrease_age_gene_idx <- c(decrease_age_gene_idx, i)
    }
  }
}

candidate_gene_idx <- c(increase_age_gene_idx, decrease_age_gene_idx)

n <- 100
gene_bootstrap_count <- c(rep(0,nrow(expr_protein_coding)))
for(iter in 1:n){
set.seed(iter + 985791)
bootstrap_idx <- sample(1:ncol(expr_protein_coding),
                   round(0.7*ncol(expr_protein_coding)),
                   replace = T)
bootstrap_expr <- expr_protein_coding[,bootstrap_idx]
bootstrap_pca_dt <- prcomp(t(bootstrap_expr),
                      center = T,
                      retx =T)
bootstrap_pca_dt$x <- bootstrap_pca_dt$x[,1:10]
bootstrap_pca_dt$x <- bootstrap_pca_dt$x[,apply(bootstrap_pca_dt$x,2,function(e) cor.test(e,age[bootstrap_idx])$p.value)>0.05]
for(i in candidate_gene_idx){
  test <- cbind(bootstrap_expr[i,], age[bootstrap_idx], bootstrap_pca_dt$x)
  test <- as.data.frame(test)
  colnames(test)[1:2] <- c("exp","age")
  model1 <- lm(formula = exp ~ age, data = test)
  lm_summary <- summary(model1)
  coeff_mat <- lm_summary$coefficients
  if((coeff_mat[2,4] != "NaN") & (coeff_mat[2,4] < 0.05)){
    model2 <- lm(formula = exp ~ ., data = test)
    lm_summary <- summary(model2)
    coeff_mat <- lm_summary$coefficients
    if((coeff_mat[2,4] != "NaN") & (coeff_mat[2,4] < 0.05)){
      gene_bootstrap_count[i] <- gene_bootstrap_count[i] + 1
    }
  }
}
}
increase_age_gene_idx_valid <- intersect(increase_age_gene_idx, which(gene_bootstrap_count >= 50))
decrease_age_gene_idx_valid <- intersect(decrease_age_gene_idx, which(gene_bootstrap_count >= 50))


increase_age_gene <- expr_gene_list[protein_coding_gene_idx][valid_idx][increase_age_gene_idx_valid]
decrease_age_gene <- expr_gene_list[protein_coding_gene_idx][valid_idx][decrease_age_gene_idx_valid]


