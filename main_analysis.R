library(tidyverse)
library(data.table)
library(rlist)
library(parallel)

nunique <- function(x) length(unique(x))

### Function for countiong of TF's bibding sites 
get.TF.from.ghdb <- function(genehancer_id){
  x_TF <- filter(GHTF600, GHid %in% genehancer_id)
  x_GHTF <- bind_rows(lapply(genehancer_id, function(x) x_TF[x_TF$GHid %in% x,]))
  table(x_GHTF$TF)
} 



                                        ## upload genehancer melted DB and GHTF
# enhancer-gene association dataset (output of 'gff.attributes_tidy.R' script)
genehancerDB <- as_data_frame(fread('genehancerDB.melted.csv', drop = 1))
# enhancer TF's binding sites dataset
GHTF <- fread('GeneHancer_TFBSs.txt')

# Select enhancers with known TF's bs only
genehancerDB <- genehancerDB[genehancerDB$genehancer_id %in% GHTF$GHid,]

                                              ## upload control genes
# Upload control genes for your task and select genes with known enhancers only
#control_genes <- fread('your_control_genes.csv')
#final_control <- control_genes[control_genes$SYMBOL %in% genehancerDB$connected_genes]

                                              ## upload of testing DEGs
hsbg <- as_data_frame(fread("DEGs.csv", drop = 1)) # hsbg = human-specific brain genes
# select DEGs with kniwn enhancers only
hsbg <- hsbg[hsbg$gene %in% genehancerDB$connected_genes,]

                                              ## main analysis
pv_total = list()
cell_types = unique(hsbg$Acronym)

for(cell in cell_types){
  #getting genehancer_id list from gene list
  cell_hsbg <-genehancerDB %>% filter(connected_genes %in% 
                                           hsbg[hsbg$Acronym == cell,]$gene)
  cell_GHTF <- get.TF.from.ghdb(cell_hsbg$genehancer_id)
  #random sampling
  random_cell <- replicate(1000, as_tibble(final_control$SYMBOL, name = NULL) %>% 
                                sample_n(length(hsbg[hsbg$Acronym == cell,]$gene)), simplify=F)
  random_cell_ghdb <- mclapply(random_cell, function(x)
    filter(genehancerDB, connected_genes %in% x$value), mc.cores=25)
  random_cell_GHTF <- unlist(mclapply(random_cell_ghdb, function(x)
    get.TF.from.ghdb(x$genehancer_id), mc.cores=25))
  random_cell_GHTF <- split(unname(random_cell_GHTF), names(random_cell_GHTF))
  random_cell_GHTF <- random_cell_GHTF[names(random_cell_GHTF) %in% names(cell_GHTF)]
  #permutitaion test
  p_val = c()
  for(tf in 1:length(cell_GHTF)){
    p_val[tf] <- sum(random_cell_GHTF[[tf]] >= cell_GHTF[tf]) / 1000
  }
  names(p_val) <- names(cell_GHTF)
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

for (i in 1:length(padj_sign)){
  df = data.frame(cell_type = names(padj_sign)[i], 
                              TF = names(padj_sign[[i]]), 
                              p_adj = padj_sign[[i]]) 

  padj_sign_df = bind_rows(padj_sign_df, df)
}

rownames(padj_sign_df) = 1:nrow(padj_sign_df)

write.table(padj_sign_df, 'result_tables/TFbs by final_control.csv', sep=';')



