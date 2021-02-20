if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")

library(DESeq2)
library(TCGAbiolinks)

protein_coding_ensembl <- read.table("D:/Age_Project/MSB/code/requirements/protein_coding_ensembl.txt", sep="\t", header=FALSE)$V1


cancer <- "kirc"
query <- GDCquery(project = paste0("TCGA-",toupper(cancer)),
                  data.category = "Transcriptome Profiling",
                  workflow.type = "HTSeq - Counts",
                  sample.type = "Solid Tissue Normal")
GDCdownload(query)
data_normal <- GDCprepare(query)
expr_gene_list <- data_normal@rowRanges$ensembl_gene_id
protein_coding_gene_idx <- which(expr_gene_list %in% protein_coding_ensembl)

patients <- data_normal$patient
COUNTS_normal <- data_normal@assays@data$`HTSeq - Counts`

tumor_barcode = query$results[[1]]$cases[sapply(patients, function(arg) which(query$results[[1]]$cases.submitter_id == arg))]


query <- GDCquery(project = paste0("TCGA-",toupper(cancer)),
                  data.category = "Transcriptome Profiling",
                  workflow.type = "HTSeq - Counts",
                  sample.type = "Primary Tumor",
                  barcode = tumor_barcode)
GDCdownload(query)
data_tumor <- GDCprepare(query)
COUNTS_tumor <- data_tumor@assays@data$`HTSeq - Counts`



#### Finding cancer DEG ####

cts <- cbind(COUNTS_normal[protein_coding_gene_idx,], COUNTS_tumor[protein_coding_gene_idx,])
colnames(cts) <- c(data_normal$sample, data_tumor$sample)
coldata <- data.frame(condition = c(rep("Normal", ncol(COUNTS_normal)), rep("Tumor", ncol(COUNTS_tumor))),
                      patient = c(patients,patients))
rownames(coldata) <- colnames(cts)
dds <- DESeqDataSetFromMatrix(countData=cts, colData = coldata, design = ~ condition + patient)
dds <- DESeq(dds)
res <- results(dds,alpha = 0.01,lfcThreshold = 2,contrast = c("condition","Tumor","Normal"))



##### 
query <- GDCquery(project = paste0("TCGA-",toupper(cancer)),
                  data.category = "Transcriptome Profiling",
                  workflow.type = "HTSeq - FPKM-UQ",
                  sample.type = "Solid Tissue Normal")
GDCdownload(query)
data_normal <- GDCprepare(query)


expr_normal <- data_normal@assays@data$`HTSeq - FPKM-UQ`


query <- GDCquery(project = paste0("TCGA-",toupper(cancer)),
                  data.category = "Transcriptome Profiling",
                  workflow.type = "HTSeq - FPKM-UQ",
                  sample.type = "Primary Tumor",
                  barcode = tumor_barcode)
GDCdownload(query)
data_tumor <- GDCprepare(query)
expr_tumor <- data_tumor@assays@data$`HTSeq - FPKM-UQ`

expr_tumor <- expr_tumor[,sapply(patients, function(arg) which(data_tumor$patient == arg))]



inc_gene <- read.table("D:/Age_Project/MSB/code/requirements/KIRC_inc.txt", header = FALSE)$V1
dec_gene <- read.table("D:/Age_Project/MSB/code/requirements/KIRC_dec.txt", header = FALSE)$V1

inc_gene_idx <- sapply(inc_gene, function(arg) which(expr_gene_list == arg))
inc_gexp <- expr_normal[inc_gene_idx,]
inc_gexp <- apply(inc_gexp,1,function(e) (e-min(e))/(max(e)-min(e)))
inc_gexp <- apply(t(inc_gexp),2,mean)
inc_gexp <- scale(inc_gexp, center = T)


dec_gene_idx <- sapply(dec_gene, function(arg) which(expr_gene_list == arg))
dec_gexp <- expr_normal[dec_gene_idx,]
dec_gexp <- apply(dec_gexp,1,function(e) (e-min(e))/(max(e)-min(e)))
dec_gexp <- apply(t(dec_gexp),2,mean)
dec_gexp <- scale(dec_gexp, center = T)

up_deg <- expr_gene_list[protein_coding_gene_idx][which(res$padj < 0.01 & res$log2FoldChange > 2)]
down_deg <- expr_gene_list[protein_coding_gene_idx][which(res$padj < 0.01 & res$log2FoldChange < -2)]


up_gene_idx <- sapply(up_deg, function(arg) which(expr_gene_list == arg))
up_gexp <- expr_tumor[up_gene_idx,]
up_gexp <- apply(up_gexp,1,function(e) (e-min(e))/(max(e)-min(e)))
up_gexp <- apply(t(up_gexp),2,mean)
up_gexp <- scale(up_gexp, center = T)


down_gene_idx <- sapply(down_deg, function(arg) which(expr_gene_list == arg))
down_gexp <- expr_tumor[down_gene_idx,]
down_gexp <- apply(down_gexp,1,function(e) (e-min(e))/(max(e)-min(e)))
down_gexp <- apply(t(down_gexp),2,mean)
down_gexp <- scale(down_gexp, center = T)



clinical <- read.table("D:/Age_Project/MSB/code/requirements/kirc_tcga_clinical_data.tsv", sep = "\t", header = TRUE)
clinical <- clinical[sapply(patients, function(arg) which(clinical$Patient.ID == arg)),]
age <- data_normal$days_to_birth/-365.25
status <- factor(ifelse(is.na(data_normal$days_to_death), 0, 1))
surv_time <- clinical$Overall.Survival..Months.



my_data <- data.frame(surv_time = surv_time,
                      age = age,
                      status = status,
                      inc_gexp = inc_gexp,
                      dec_gexp = dec_gexp,
                      up_gexp = up_gexp,
                      down_gexp = down_gexp)

library(caret)
library(ROCR)
library(reshape2)

plan <- c("status ~ down_gexp",
          "status ~ up_gexp",
          "status ~ dec_gexp",
          "status ~ inc_gexp",
          "status ~ up_gexp + inc_gexp")

k <- 5
n <- 1000
id <- rep_len(1:5,length(pat))
auc_df <- NULL
data_list <- list()
pb <- progress_bar$new(total = n)
for(j in 1:n){
  pb$tick()
  temp_data <- rbind(my_data[which(my_data$status==0),][sample(1:length(which(my_data$status==0))),],
                     my_data[which(my_data$status==1),][sample(1:length(which(my_data$status==1))),])
  sample_data <- data.frame(temp_data,id=id)
  data_list[[j]] <- sample_data
  
  auc_vec <-NULL
  for(form in plan){
    temp <- 0
    for(i in 1:k){
      Train <- sample_data[which(sample_data$id != i),]
      Test <- sample_data[which(sample_data$id == i),]
      model <- glm(form, family = "binomial",data = Train)
      prob <- predict(model, newdata = Test, type = "response")
      pred <- prediction(prob, Test$status)
      auc <- performance(pred, measure = "auc")
      #print(auc@y.values[1])
      temp <- temp + as.numeric(auc@y.values[1])
    }
    temp <- temp/k
    auc_vec <- c(auc_vec,temp)
  }
  auc_df <- rbind(auc_df,auc_vec)
}

rownames(auc_df) <- 1:n
colnames(auc_df) <- c("Downregulated in tumor tissue",
                      "Upregulated in tumor tissue",
                      "Decreasing with age in normal tissue",
                      "Increasing with age in normal tissue",
                      "Upregulated in tumor tissue + Increasing with age in normal tissue")
apply(auc_df,2,mean)

ggplot(data = melt(auc_df), aes(x=Var2, y= value)) + theme_bw() + 
  geom_boxplot(aes(fill = Var2), width = 0.35) +
  ylab("AUC") + xlab(element_blank()) + 
  guides(fill=guide_legend(title="Input Variables")) +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 15))
