#This code allows calculating abundance and thermal stability scores
#for the TPP analysis of SARS-CoV-2 infection dynamics 


#Install required packages
library(tidyverse)
library(vsn)
library(limma)
library(fdrtool)

#Load sample mapping file (contains information of what is in each TMT channel in each sample)
sample_mapping <- read_delim('Sample_mapping.csv',delim = ',')

#Make sure that protein files are in MS_data folder (data can be downloaded from PRIDE - files ending with '_proteins.txt')
results_files <- dir('MS_data')

#Import protein files, perform vsn normalization per temperature and calculate log2 fold-changes to uninfected sample
prot_melt <- bind_rows(lapply(results_files, function(MS_file){
  
  #Import protein file, remove contaminants and reverse database hits and filter for proteins with at least 2 unique peptides
  prot_tab <- read_delim(file.path('MS_data',MS_file), delim = "\t") %>%
    filter(!grepl("##", protein_id),
           qupm > 1) %>%
    group_by(gene_name) %>%
    filter(qupm == max(qupm), top3 == max(top3)) %>%
    ungroup
  
  #Perform vsn normalization per time point (first 8 labels correspond to first temperature of the plex, and last 8 the second one)
  signal_sum_mat <- as.matrix(prot_tab %>%
                                dplyr::select(matches("signal_sum_")))

  signal_sum_mat[is.infinite(log2(signal_sum_mat)) | is.nan(log2(signal_sum_mat))] <- NA

  vsn_fit_1 <- vsn2(signal_sum_mat[,1:8])
  vsn_norm_mat_1 <- as.data.frame(predict(vsn_fit_1, signal_sum_mat[,1:8]))
  vsn_norm_mat_1 <- vsn_norm_mat_1 - vsn_norm_mat_1[,1]
  vsn_fit_2 <- vsn2(signal_sum_mat[,9:16])
  vsn_norm_mat_2 <- as.data.frame(predict(vsn_fit_2, signal_sum_mat[,9:16]))
  vsn_norm_mat_2 <- vsn_norm_mat_2 - vsn_norm_mat_2[,1]
  vsn_norm_mat <- cbind(vsn_norm_mat_1,vsn_norm_mat_2)
  vsn_norm_mat$gene_name <- prot_tab$gene_name

  vsn_norm_mat %>%
    tbl_df %>%
    gather(key,value,-gene_name) %>%
    mutate(MS_sample_number = sub('_merged_results.*','',sub('.S[0-9]{4}','',sub(sub('S[0-9]{4}.*','',MS_file),'',MS_file)))) %>%
    dplyr::rename(log2fc = value)

}))

#Merging data from all experiments and annotating each data point
prot_merged <- prot_melt %>%
  left_join(sample_mapping %>%
              mutate(key = paste0('signal_sum_',TMT_label)) %>%
              dplyr::select(Time_point,Replicate,MS_sample_number,key,Temperature,Sample),
            by = c('key' = 'key','MS_sample_number' = 'MS_sample_number')) %>%
  dplyr::select(-key) %>%
  group_by(gene_name, Time_point, Sample) %>%
  mutate(rep_abun = sum(Temperature == 40),
         rep_stab = n()) %>%
  ungroup() %>%
  filter(Time_point != 0)

#Calculating abundance scores
abun_scores_raw <- prot_merged %>% 
  filter(rep_abun >= 2,
         Temperature %in% sort(unique(sample_mapping$Temperature))[1:2]) %>% 
  group_by(gene_name,Time_point,Replicate,Sample) %>% 
  mutate(abun_score = mean(log2fc, na.rm = TRUE)) %>% 
  ungroup() %>% 
  dplyr::select(gene_name,Time_point,Replicate,Sample,abun_score) %>%
  distinct()

abun_scores <- abun_scores_raw %>% 
  mutate(condition = paste0('t',Time_point,'_sample',Sample,'_rep',Replicate)) %>% 
  dplyr::select(gene_name,condition,abun_score) %>%
  spread(condition,abun_score)

abun_scores <- as.data.frame(abun_scores)
rownames(abun_scores) <- abun_scores$gene_name
abun_scores$gene_name <- NULL

#Calculating thermal stability scores and weights for limma analysis
stab_scores_raw <- prot_merged %>% 
  filter(rep_abun >= 2, 
         rep_stab >= 10) %>%
  left_join(abun_scores_raw) %>% 
  group_by(gene_name,Time_point,Replicate,Sample) %>% 
  mutate(stab_score = sum(log2fc-abun_score, na.rm = TRUE),
         stab_score_weight = n()) %>% 
  ungroup() %>% 
  dplyr::select(gene_name,Time_point,Replicate,Sample,stab_score,stab_score_weight) %>% 
  distinct() %>% 
  filter(stab_score_weight > 4)

stab_scores <- stab_scores_raw %>% 
  mutate(condition = paste0('t',Time_point,'_sample',Sample,'_rep',Replicate)) %>% 
  dplyr::select(gene_name,condition,stab_score) %>%
  spread(condition,stab_score)

stab_scores <- as.data.frame(stab_scores)
rownames(stab_scores) <- stab_scores$gene_name
stab_scores$gene_name <- NULL

stab_weights <- stab_scores_raw %>% 
  mutate(condition = paste0('t',Time_point,'_sample',Sample,'_rep',Replicate)) %>% 
  dplyr::select(gene_name,condition,stab_score_weight) %>%
  spread(condition,stab_score_weight)

stab_weights <- as.data.frame(stab_weights)
rownames(stab_weights) <- stab_weights$gene_name
stab_weights$gene_name <- NULL

#limma analysis and FDR tool correction
limma_results <- bind_rows(lapply(paste0(rep(paste0('t',unique(sample_mapping$Time_point)[-1]),each = length(unique(sample_mapping$Sample))),'_sample',unique(sample_mapping$Sample)), function(time_point) {
  abun_scores_temp <- abun_scores[,grep(time_point,colnames(abun_scores)),drop = FALSE]
  abun_scores_temp <- abun_scores_temp[rowSums(!is.na(abun_scores_temp)) > 1,,drop = FALSE]
  stab_scores_temp <- stab_scores[,grep(time_point,colnames(stab_scores)),drop = FALSE]
  stab_weights_temp <- stab_weights[,grep(time_point,colnames(stab_weights)),drop = FALSE]
  stab_scores_temp <- stab_scores_temp[rowSums(!is.na(stab_scores_temp)) > 1,,drop = FALSE]
  stab_weights_temp <- stab_weights_temp[rownames(stab_scores_temp),]
  
  metadata <- data.frame(ID = colnames(abun_scores_temp),
                         time_point = gsub("_.+","",colnames(abun_scores_temp)),
                         sample = gsub("_.+","",gsub(".+_sample","",colnames(abun_scores_temp))),
                         rep = gsub(".+_","",colnames(abun_scores_temp)))
  rownames(metadata) <- metadata$ID
  
  abundance.dataE <- ExpressionSet(assayData = as.matrix(abun_scores_temp),
                                   phenoData = AnnotatedDataFrame(metadata))
  stability.dataE <- ExpressionSet(assayData = as.matrix(stab_scores_temp),
                                   phenoData = AnnotatedDataFrame(metadata))
  
  abun_comparison <- eBayes(lmFit(abundance.dataE))
  stab_comparision <- eBayes(lmFit(stability.dataE,weights = stab_weights_temp^2))
  
  abun_res <- limma::topTable(abun_comparison, sort.by = "t",  
                              coef = 1, number = Inf)
  stab_res <- limma::topTable(stab_comparision, sort.by = "t",  
                              coef = 1, number = Inf)
  
  abun_res$gene_name <- rownames(abun_res)
  abun_res$comparison <- 'abundance'
  abun_res$time_point <- unique(metadata$time_point)
  abun_res$sample <- unique(metadata$sample)
  stab_res$gene_name <- rownames(stab_res)
  stab_res$comparison <- 'stability'
  stab_res$time_point <- unique(metadata$time_point)
  stab_res$sample <- unique(metadata$sample)
  
  abun_res$FDR_tool <- fdrtool(abun_res$t,verbose = FALSE, plot = FALSE)$qval
  stab_res$FDR_tool <- fdrtool(stab_res$t,verbose = FALSE, plot = FALSE)$qval
  
  rbind(abun_res,stab_res) %>% 
    tbl_df() %>% 
    dplyr::select(gene_name,time_point,sample,comparison,logFC,P.Value,adj.P.Val,FDR_tool) %>% 
    dplyr::rename(p_value = P.Value,
                  FDR_limma = adj.P.Val)
}))

#Final table with all abundance/thermal stability scores and annotation if significant hit
scaled_limma_results <- limma_results %>% 
  mutate(time_point = as.numeric(sub('t','',time_point))) %>% 
  group_by(comparison) %>% 
  mutate(z_score = scale(logFC)) %>% 
  ungroup() %>% 
  mutate(hit = abs(z_score) > qnorm(0.975) & FDR_tool < 0.05)


