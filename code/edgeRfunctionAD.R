#edgeR function only for AD

load("/Users/stavroula/edgeR/All_Counts/HIP_Counts.rdata")
load("/Users/stavroula/edgeR/All_Counts/FWM_Counts.rdata")
load("/Users/stavroula/edgeR/All_Counts/PCx_Counts.rdata")
load("/Users/stavroula/edgeR/All_Counts/TCx_Counts.rdata")

HIP_counts <- subset(HIP_counts, select = -c(GeneID, Length))
FWM_counts <- subset(FWM_counts, select = -c(GeneID, Length))

#treatments by treatment.ipynb AD only
HIP_treatment = c(3,3,2,2,2,3,3,3,2,3,2,3,2,2,2,3,2,3,3,2,3,2,3,2,3,3,2,2,2,3,3,2,2,2,2,2,2,3,2,3,2,2,3,2,3,2,3,2,2)
FWM_treatment = c(3,2,2,2,3,2,3,3,2,3,2,2,2,2,3,2,3,3,2,3,2,3,2,3,3,2,2,2,3,2,2,2,2,3,3,2,3,2,2,3,2,3,3,2,3,2)
PCx_treatment = c(3,2,2,2,3,2,3,3,2,3,2,2,2,2,3,2,3,3,2,2,3,2,3,3,2,2,2,3,2,2,2,2,3,3,2,3,2,2,3,2,3,3,2,3,2,2)
TCx_treatment = c(3,3,2,2,3,2,3,2,3,3,2,3,2,2,2,2,3,2,3,3,2,3,2,3,2,3,3,2,2,2,3,2,2,2,2,2,3,2,3,3,2,3,2,3,3,2,3,2,2,2)

HIP_list = read.csv("/Users/stavroula/edgeR/24-7-2022/OnlyADcsv/HIP_AD_Only.csv")
FWM_list = read.csv("/Users/stavroula/edgeR/24-7-2022/OnlyADcsv/FWM_AD.csv")
PCx_list = read.csv("/Users/stavroula/edgeR/24-7-2022/OnlyADcsv/PCx_AD.csv")
TCx_list = read.csv("/Users/stavroula/edgeR/24-7-2022/OnlyADcsv/TCx_AD.csv")

edgeRfunctionADonly <- function(p,q,r) {
  a = r$donor_name
  b = lapply(a, function(x) strsplit(x, ".", fixed=TRUE)[[1]][[3]])
  c = r$structure_acronym
  d = as.character(r$rna_well)
  e = do.call(rbind.data.frame, b)
  e$structure = c
  e$rna = d
  f = apply(e,1,paste,collapse="_")
  f = paste0(f, ".sorted.bam")
  #f = paste0("X", f)
  p = p[, f]
  countMatrix <- p

  #Data Cleaning
  cpm_log <- cpm(countMatrix, log = TRUE)
  median_log2_cpm <- apply(cpm_log, 1, median)
  hist(median_log2_cpm)
  expr_cutoff <- -1
  abline(v = expr_cutoff, col = "red", lwd = 3)
  sum(median_log2_cpm > expr_cutoff)
  data_clean <- countMatrix[median_log2_cpm > expr_cutoff, ]
  treatment <- q
  
  #DGEList
  d <- DGEList(counts = data_clean, group = treatment)
  
  #CalcNormFactors
  cnf <- calcNormFactors(d)
  
  #design matrix definition
  designMatrix <- model.matrix(~treatment)
  
  #EstimateDisp
  estDisp <- estimateDisp(cnf, designMatrix)
  
  #glmQLFit
  fit <- glmQLFit(estDisp, designMatrix)
  
  #glmQLFTest
  results <- glmQLFTest(fit)
  
  #topTags
  tags <- topTags(results, 'all')
  
  #save topTags full list
  write.csv(tags, "/Users/stavroula/edgeR/24-7-2022/8-9-2022/tagsADonly.csv", row.names = TRUE)

  print(summary(decideTestsDGE(results, adjust.method="none", p.value=0.001)))
  
  #Plots edgeR, p=0.001
  a <- tags
  a <- a$table
  b <- data.frame(a$logFC, a$PValue)
  rownames(b) <- rownames(a)
  EnhancedVolcano(b,
                  lab = rownames(b),
                  x = 'a.logFC',
                  y = 'a.PValue',
                  title = 'Some title',
                  pCutoff = 0.001,
                  FCcutoff = 0.1,
                  pointSize = 3.0,
                  labSize = 6.0)
}

#example calling function:
#edgeRfunctionADonly(HIP_counts, HIP_treatment, HIP_list)