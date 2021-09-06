library(tidyverse)
library(data.table)
library(rlist)
library(parallel)

nunique <- function(x) length(unique(x))

### функция для расчета 
get.TF.from.ghdb <- function(genehancer_id){
  x_TF <- filter(GHTF600, GHid %in% genehancer_id)
  x_GHTF <- bind_rows(lapply(genehancer_id, function(x) x_TF[x_TF$GHid %in% x,]))
  table(x_GHTF$TF)
} # подаем сразу вектор значений



                                        ## upload genehancer melted DB and GHTF
genehancerDB <- as_data_frame(fread('/hdd20tb/marat/shizo/genehancerDB.melted.csv', drop = 1))
GHTF <- fread('/hdd20tb/marat/shizo/GeneHancer_TFBSs.txt')

genehancerDB600 <- genehancerDB[genehancerDB$score >= 600,]
genehancerDB <- genehancerDB[genehancerDB$genehancer_id %in% GHTF$GHid,]
nunique(genehancerDB$genehancer_id) #275062
quantile(genehancerDB$score, probs = 0.9)

GHTF600 <- GHTF %>% filter(GHid %in% genehancerDB600$genehancer_id)



                                              ## upload control genes
#final_control <- fread('/hdd20tb/marat/shizo/final_df_bgn_symbol.csv')
#final_control600 <- final_control[final_control$SYMBOL %in% genehancerDB600$connected_genes]
#nrow(final_control600) # 9880 genes

#cntrl <- fread('/hdd20tb/marat/shizo/df_cntrl_symbol.csv')
#cntrl600 <- cntrl %>% filter(SYMBOL %in% genehancerDB600$connected_genes)

shapiro_cntrl <- fread('shapired_df_symbol.csv') #10755
shapiro_cntrl600 <- shapiro_cntrl %>% filter(SYMBOL %in% genehancerDB600$connected_genes) #10548




# upload Brain Shizo Specific Genes 
bssg <- as_data_frame(fread("/hdd20tb/marat/shizo/gene_symbol_wilcox_rin_regression_bn.csv", drop = 1))
nunique(bssg$gene) # 5620

bssg600 <- bssg[bssg$gene %in% genehancerDB600$connected_genes,]
nunique(bssg600$gene) # 5439


pv_total = list()
cell_types = unique(bssg600$Acronym)

for(cell in cell_types){
  #getting genehancer_id list from gene list
  cell_ghdb600 <-genehancerDB600 %>% filter(connected_genes %in% 
                                           bssg600[bssg600$Acronym == cell,]$gene)
  cell_GHTF600 <- get.TF.from.ghdb(cell_ghdb600$genehancer_id)
  #random sampling
  random_cell <- replicate(1000, as_tibble(final_control600$SYMBOL, name = NULL) %>% 
                                sample_n(length(bssg600[bssg600$Acronym == cell,]$gene)), simplify=F)
  random_cell_ghdb600 <- mclapply(random_cell, function(x)
    filter(genehancerDB600, connected_genes %in% x$value), mc.cores=25)
  random_cell_GHTF600 <- unlist(mclapply(random_cell_ghdb600, function(x)
    get.TF.from.ghdb(x$genehancer_id), mc.cores=25))
  random_cell_GHTF600 <- split(unname(random_cell_GHTF600), names(random_cell_GHTF600))
  random_cell_GHTF600 <- random_cell_GHTF600[names(random_cell_GHTF600) %in% names(cell_GHTF600)]
  #permutitaion test
  p_val = c()
  for(tf in 1:length(cell_GHTF600)){
    p_val[tf] <- sum(random_cell_GHTF600[[tf]] >= cell_GHTF600[tf]) / 1001
  }
  names(p_val) <- names(cell_GHTF600)
  pv_total3 = c(pv_total3, list(p_val))
}
names(pv_total) = cell_types



                                # p-adjust and getiing significant TFs
padj_total = lapply(pv_total, function(x) p.adjust(x, method = "BH"))


padj_sign = lapply(padj_total, function(x) x[x <= 0.05])

save(padj_sign, file="padj_sign_df_final_cntrl.RData")





rm(pv_total)
#



####make df from padj list
sapply(padj_sign, length)

padj_sign = padj_sign[sapply(padj_sign, length) != 0]


padj_sign_df = data.frame(cell_type = as.character(), TF= as.character(), p_adj=as.numeric())
rm(i)
for (i in 1:length(padj_sign)){
  df = data.frame(cell_type = names(padj_sign)[i], 
                              TF = names(padj_sign[[i]]), 
                              p_adj = padj_sign[[i]]) 

  padj_sign_df = bind_rows(padj_sign_df, df)
}

rownames(padj_sign_df) = 1:nrow(padj_sign_df)

write.table(padj_sign_df, 'result_tables/TFbs by final_control.csv', sep=';')



