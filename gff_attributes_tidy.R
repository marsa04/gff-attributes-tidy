library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(parallel)


#upload .gff file
genehancerDB <- fread('/Users/maratsabirov/Desktop/HumanBrain/genehancerV5by1stJune2021/genehancer.gff')

# function for parsing 
getAttributeField1 <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}

#get_genehancer_id
db.test.genehancer_id <- getAttributeField1(genehancerDB$attributes, field = "genehancer_id")
db.test.genehancer_id <- unlist(lapply(db.test.genehancer_id, function(x) setNames(x,"genehancer_id")))
db.test.genehancer_id[1:10]

#get_score
db.test1 <- gsub(genehancerDB$attributes, pattern = "genehancer_id=", replacement = "",fixed = T)
db.test2 <- strsplit(db.test1, split = ";connected_gene=", fixed = TRUE)
db.test.score <- sapply(db.test2, function(x) getAttributeField1(x, field = "score"))
db.test.score.1 <- sapply(db.test.score, function(x) x[!is.na(x)])

#get_connected_gene
db.test.connected_gene <- strsplit(db.test1, split = ";score=", fixed = TRUE)
db.test.connected_gene.1 <- sapply(db.test.connected_gene, function(x) 
  getAttributeField1(x, field = "connected_gene"))
db.test.connected_gene.2 <- sapply(db.test.connected_gene.1, function(x) x[!is.na(x)])

# get_DB
i <- 0
for (i in 1:length(db.test.score.1)) {
  db.test.score.1[[i]] <- setNames(db.test.score.1[[i]], db.test.connected_gene.2[[i]])}
j <- 0
for (j in 1:length(db.test.score.1)){
  db.test.score.1[[j]] <- c(db.test.genehancer_id[j], db.test.score.1[[j]])}

#combine all generates attributes in one dataset
genehancerDB.melted <- bind_rows(mclapply(db.test.score.1, 
                                        function(x) gather(as.data.frame(t(x)), 
                                                           key = "connected_genes", 
                                                           value = "score", -V1), mc.cores=detectCores()))

colnames(genehancerDB.melted)[1] <- "genehancer_id"
genehancerDB.melted$score <- as.numeric(genehancerDB.melted$score)
str(genehancerDB.melted)

# save generated dataset 
write.table(genehancerDB.melted, file = "/Users/maratsabirov/Desktop/HumanBrain/genehancerV5by1stJune2021/genehancerDB.melted.csv", sep = ",")













