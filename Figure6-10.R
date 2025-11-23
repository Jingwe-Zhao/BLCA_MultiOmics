######bulk RNA-Seq
library(ggrisk)
library(survival)
library(rms)
library(patchwork)
library(dplyr)
library(GEOquery)

#GSE31684----
gse <- getGEO("GSE31684", destdir = "D:/DATA/BLCA/GSE31684/")
expr_matrix <- exprs(gse[[1]])
pheno_data <- pData(gse[[1]])

colnames(pheno_data)
table(pheno_data$`last known status:ch1`)
meta <- pheno_data[, c("last known status:ch1","survival.months:ch1"), drop = FALSE]


ribofl<-c('ACP1','FLAD1','BLVRB','ACP2','ACP5','RFK','ENPP1','ENPP3')

expr_final<-as.data.frame(expr_final)
expr3 <- subset(expr_final, rownames(expr_final) %in% ribofl)
expr3<-as.data.frame(t(expr3))


expr3$time=meta$`survival.months:ch1`
expr3$status=meta$`last known status:ch1`

table(expr3$status)

expr3$status <- ifelse(expr3$status %in% c('DOC','DOD'), 1, 0)
table(expr3$status)

expr3$time<-as.numeric(expr3$time)

sum(is.na(expr3))  
sapply(expr3, function(x) sum(is.na(x))) 

expr3_clean <- na.omit(expr3)

sum(is.na(expr3_clean)) 


fit2 <- coxph(Surv(time, status)~., expr3_clean)
fit2

ggrisk(fit2,
       color.A = c(low = '#6395C7',high = '#E06EAD'),
       color.B = c(code.0 ='#6395C7',code.1 = '#E06EAD'), 
       color.C = c(low = "#6395C7",median = 'white',high = '#E06EAD'))


#GSE48256----
gse <- getGEO("GSE48256", destdir = "D:/DATA/BLCA/GSE48256/")
expr_matrix <- exprs(gse[[1]])
pheno_data <- pData(gse[[1]])

colnames(pheno_data)
table(pheno_data$`last known status:ch1`)
meta <- pheno_data[, c("last known status:ch1","survival.months:ch1"), drop = FALSE]


ribofl<-c('ACP1','FLAD1','BLVRB','ACP2','ACP5','RFK','ENPP1','ENPP3')

expr_final<-as.data.frame(expr_final)
expr3 <- subset(expr_final, rownames(expr_final) %in% ribofl)
expr3<-as.data.frame(t(expr3))


expr3$time=meta$`survival.months:ch1`
expr3$status=meta$`last known status:ch1`

table(expr3$status)

expr3$status <- ifelse(expr3$status %in% c('DOC','DOD'), 1, 0)
table(expr3$status)

expr3$time<-as.numeric(expr3$time)

sum(is.na(expr3))  
sapply(expr3, function(x) sum(is.na(x))) 

expr3_clean <- na.omit(expr3)

sum(is.na(expr3_clean)) 


fit2 <- coxph(Surv(time, status)~., expr3_clean)
fit2

ggrisk(fit2,
       color.A = c(low = '#6395C7',high = '#E06EAD'),
       color.B = c(code.0 ='#6395C7',code.1 = '#E06EAD'), 
       color.C = c(low = "#6395C7",median = 'white',high = '#E06EAD'))


#GSE32894----
gse <- getGEO("GSE32894", destdir = "D:/DATA/BLCA/GSE32894/")
expr_matrix <- exprs(gse[[1]])
pheno_data <- pData(gse[[1]])

colnames(pheno_data)
table(pheno_data$`last known status:ch1`)
meta <- pheno_data[, c("last known status:ch1","survival.months:ch1"), drop = FALSE]


ribofl<-c('ACP1','FLAD1','BLVRB','ACP2','ACP5','RFK','ENPP1','ENPP3')

expr_final<-as.data.frame(expr_final)
expr3 <- subset(expr_final, rownames(expr_final) %in% ribofl)
expr3<-as.data.frame(t(expr3))


expr3$time=meta$`survival.months:ch1`
expr3$status=meta$`last known status:ch1`

table(expr3$status)

expr3$status <- ifelse(expr3$status %in% c('DOC','DOD'), 1, 0)
table(expr3$status)

expr3$time<-as.numeric(expr3$time)

sum(is.na(expr3))  
sapply(expr3, function(x) sum(is.na(x))) 

expr3_clean <- na.omit(expr3)

sum(is.na(expr3_clean)) 


fit2 <- coxph(Surv(time, status)~., expr3_clean)
fit2

ggrisk(fit2,
       color.A = c(low = '#6395C7',high = '#E06EAD'),
       color.B = c(code.0 ='#6395C7',code.1 = '#E06EAD'), 
       color.C = c(low = "#6395C7",median = 'white',high = '#E06EAD'))



##GSVA-----
library(GSVA)

expr_final2 <- as.data.frame(t(expr_final))
expr_final2$SampleID <- rownames(expr_final2)

pheno_data2<-pheno_data[,1, drop = FALSE]
pheno_data2$SampleID <- rownames(pheno_data2)


merged_df <- merge(expr_final2, pheno_data2, by = "SampleID", all = TRUE)
rownames(merged_df)<-merged_df$SampleID
merged_df<-merged_df[,-1]

dim(merged_df)
colnames(merged_df)[1:20] 


merged_df$stage <- substr(merged_df$title, nchar(merged_df$title)-1, nchar(merged_df$title))
merged_df$title<-NULL
table(merged_df$stage)

save(merged_df,file = 'D:/DATA/WUST/BLCA/视黄素返修/data/GSE31684/merged_df.rdata')

meta<-merged_df[,'stage', drop = FALSE]
meta$CB<-row.names(meta)

dim(merged_df)
colnames(merged_df)[1:20] 

merged_df$stage<-NULL

merged_df2<-(t(merged_df))
class(merged_df2)
merged_df2<-as.matrix(merged_df2)

merged_df2[1:5, 1:5]

ribofl<-c('ACP1','FLAD1','BLVRB','ACP2','ACP5','RFK','ENPP1','ENPP3')
ribofl_gset <- list(Riboflavin_Related = ribofl)

gsvaPar <- ssgseaParam(exprData = merged_df2, 
                       geneSets = ribofl_gset,
                       normalize = TRUE)
gsva_data <- gsva(gsvaPar, verbose = FALSE)


#merge GSVA and stage
ES<-as.data.frame(t(gsva_data))
ES$CB=rownames(ES)

score<-merge(ES,meta,by='CB')

COL2<-c("Ta" = "#7BAFDE", "T1" = "#FDB462", "T2" = "#BEAED4", "T3" = "#7570B3", "T4" = "#1965B0")

score$stage <- factor(score$stage, levels = c("Ta", "T1", "T2", "T3", "T4"))

ggplot(score, aes(x = stage, y = Riboflavin_Related)) +     
  geom_boxplot(aes(fill = stage), color = "black", width = 0.7, outlier.shape = NA) +  
  geom_jitter(width = 0.2, size = 1, alpha = 0.6, color = "grey20") +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.8, color = "black") +  
  stat_compare_means(method = "anova", label.y =1.96, size = 6) +  
  scale_fill_manual(values = COL2) +    
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        plot.title = element_text(size = 16, color = "black"),
        text = element_text(size = 14, color = "black")) +
  ylab('Riboflavin metabolism Score') + 
  xlab('Stage')+ggtitle('GSE31684')


do_mr_online_one_gene <- function(gene,token,outcome,pval){
  
  #'@note This function only support the one by one gene input
  
  library(clusterProfiler)
  library(forestploter)
  library(TwoSampleMR)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  
  print("进行代理设定")
  Sys.setenv(OPENGWAS_JWT = token)
  ieugwasr::get_opengwas_jwt()
  ieugwasr::user()
  
  print("获取并筛选暴露的工具变量")
  gene <- gene[!duplicated(gene)]
  df <- clusterProfiler::bitr(geneID = gene,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
  id <- df$ENSEMBL;id <- id[!duplicated(id)]
  exposure_id <- paste0('eqtl-a-',id)
  exposure_dat <- extract_instruments(outcomes = exposure_id,
                                      p1 = 5e-08, 
                                      clump = TRUE,
                                      r2 = 0.001,
                                      p2 = 5e-08,
                                      kb = 10000)
  R2a = 2*exposure_dat$beta.exposure*exposure_dat$beta.exposure*exposure_dat$eaf.exposure*(1 - exposure_dat$eaf.exposure)
  R2b = 2*exposure_dat$se.exposure*exposure_dat$se.exposure*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1 - exposure_dat$eaf.exposure)
  R2 = R2a/(R2a+R2b)
  exposure_dat$F_statistics <- R2*(exposure_dat$samplesize.exposure-2)/(1-R2)
  exposure_dat <- exposure_dat[exposure_dat$F_statistics > 10,]
  exposure_dat$exposure <- substr(exposure_dat$exposure,15,41)
  exposure_dat$ENSEMBL <- exposure_dat$exposure
  exposure_dat <- left_join(exposure_dat,df, by = "ENSEMBL")
  exposure_dat$exposure <- id
  exposure_dat$SYMBOL <- gene
  if (nrow(exposure_dat) < 1) {
    stop("基因有效eQTL少于1个 无法进行MR分析")
  }
  
  print("获取结局的工具变量并进行harmonise")
  outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP,outcomes = outcome)
  exposure_dat <- exposure_dat[exposure_dat$SNP %in% outcome_dat$SNP,]
  
  print("执行MR分析并保存数据")
  harmonised_dat <- harmonise_data(exposure_dat, outcome_dat)
  mr_modified <- function(dat = harmonised_dat, prop_var_explained = T){
    mr_res <- mr(dat)
    pve <- dat %>% 
      dplyr::select(id.exposure, beta.exposure, se.exposure, samplesize.exposure) %>% 
      dplyr::group_by(id.exposure) %>% 
      dplyr::summarise(pve = sum((beta.exposure^2)/(beta.exposure^2 + samplesize.exposure*se.exposure^2)))
    if(prop_var_explained)
    {
      mr_res <- mr_res %>% 
        dplyr::left_join(pve, by = "id.exposure")
    }
    return(mr_res)
  }
  mr_res <- mr_modified(harmonised_dat, prop_var_explained = TRUE)
  save(exposure_dat,outcome_dat,mr_res,harmonised_dat,file ='res_MR.rdata')
  
  print("所有SNP-Gene火山图绘制")
  volcano_plot <- function(.data,
                           number_comparasion = 1,
                           title = "eQTL",
                           legend.position = "none") {
    p_thershold <- 0.05/number_comparasion
    p <- .data %>%
      mutate(y = -log10(pval),
             label = ifelse(pval < p_thershold, exposure, NA)) %>%
      ggplot(aes(x = b, y = y)) +
      geom_point(aes(size = pve), alpha = 0.5, color = "#0072b5") +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_hline(yintercept = -log10(p_thershold), linetype = 2) +
      theme_classic() +
      theme(panel.grid = element_blank(),
            legend.title = element_text(size = 6.5),
            legend.text = element_text(size = 6.5),
            legend.position = legend.position) +
      labs(x = "ln(OR)",
           y = parse(text = "-log[10]*(italic(P)-value)"),
           title = title) + 
      scale_size(name = "PVE",
                 breaks = c(0.2*1:3)) +
      ggrepel::geom_label_repel(aes(label = label), size = 3)
    plot(p)
  }
  mr_res %>%
    dplyr::filter(method %in% c("Wald ratio","Inverse variance weighted",'MR Egger','Weighted median','Simple mode','Weighted mode')) %>% 
    volcano_plot(number_comparasion = 1)
  mr_res <- mr_res %>% 
    dplyr::mutate(SYMBOL = target)
  if (!exists('result')) {
    dir.create('result/')
  } else {print('结果文件夹已经存在')}
  ggsave(filename = 'result/1.volcanal.pdf',width = 6.5,height = 6.5)
  
  print("靶基因敏感性分析")
  harmonised_dat1 <- harmonised_dat[harmonised_dat$exposure == id,]
  harmonised_dat1$SNP <- paste0(harmonised_dat1$exposure,'_',harmonised_dat1$SNP)
  mr_res1 <- mr_res
  mr_scatter_plot(mr_results = mr_res1, dat = harmonised_dat1)
  ggsave(filename = 'result/2.scatter.pdf',width = 6.5,height = 6.5)
  res_single <- mr_singlesnp(harmonised_dat1)
  mr_forest_plot(res_single)
  ggsave(filename = 'result/3.forest.pdf',width = 6.5,height = 6.5)
  mr_funnel_plot(singlesnp_results = res_single)
  ggsave(filename = 'result/4.funnel.pdf',width = 6.5,height = 6.5)
  mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(harmonised_dat1))
  ggsave(filename = 'result/5.leaveout.pdf',width = 7.5,height = 7.5)
  
  # data <- table1
  # data$HR <- round(data$or,4)
  # data$HR.95L <- round(data$or_lci95,4)
  # data$HR.95H <- round(data$or_uci95,4)
  # data$pvalue <- round(as.numeric(data$`P value`),3)
  # data$pvalue = ifelse(data$pvalue<0.001,"<0.001",data$pvalue)
  # data <- data %>% 
  #   select(Gene,SNP,HR,HR.95L,HR.95H,pvalue)
  
  print("森林图绘制")
  table1 <- mr_res %>% 
    dplyr::filter(pval < threshold & method %in% c("Wald ratio","Inverse variance weighted",'MR Egger','Weighted median','Simple mode','Weighted mode')) %>% 
    dplyr::left_join(exposure_dat, by = "SYMBOL")
  table1 <- table1 %>% 
    generate_odds_ratios() %>% 
    mutate(`OR (95% CI)` = sprintf("%.3f (%.3f, %.3f)",or,or_lci95, or_uci95),
           `P value` = scales::scientific(pval),
           `PVE` = paste0(sprintf("%.2f",100 * pve),"%"),
           `F statistics` = sprintf("%.2f",F_statistics)) %>% 
    dplyr::select(Gene = SYMBOL, `ENSEMBL ID` = ENSEMBL,
                  SNP, `Effect allele` = effect_allele.exposure, 
                  `OR (95% CI)`, `P value`, 
                  PVE, `F statistics`,SNP,or,or_lci95,or_uci95)
  mr_res2 <- mr_res[mr_res$exposure %in% table1$Gene,]
  result_or <- generate_odds_ratios(mr_res2)
  mydata <- result_or
  mydata$` ` <- paste(rep(" ", 20), collapse = " ")
  mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.3f (%.3f - %.3f)",mydata$or, mydata$or_lci95, mydata$or_uci95))
  mydata$or <- round(x = mydata$or,3)
  mydata$or_uci95 <- round(x = mydata$or_uci95,3)
  mydata$or_lci95 <- round(x = mydata$or_lci95,3)
  mydata$pval <- round(x = mydata$pval,3)
  pdf(file = 'result/6.forest_all.pdf',width = 8.5,height = 7.5)
  forest(mydata[,c(4,2,3,6,14,15,9,13)],
         est = mydata$or,
         lower = round(mydata$or_lci95,3), 
         upper = round(mydata$or_uci95,3),
         sizes =0.3,
         ci_column =8 ,
         ref_line = 1,
         xlim = c(0, 4),
  )
  dev.off()
}



library(coloc)
library(tidyverse)

pqtl_locus_data <- vroom("./data/GCST90274765.tsv")
disease_locus_data <- dat

# 定义所使用的关键参数
GENE_CHR <- 9
GENE_START_POS <- 34689570  
GENE_END_POS <- 34691276    
CIS_WINDOW_BP <- 1000000 
COLOC_WINDOW_KB <- 500   
COLOC_WINDOW_BP <- COLOC_WINDOW_KB * 1000
N_CASES <- 145 + 2234;N_CONTROLS <- 351740 + 169588  
s_value <- N_CASES / (N_CASES + N_CONTROLS)

# 定义cis区域边界位置
cis_start <- max(1, GENE_START_POS - CIS_WINDOW_BP)
cis_end <- GENE_END_POS + CIS_WINDOW_BP
print(paste("已定义 'cis' 区域: chr", GENE_CHR, ":", cis_start, "-", cis_end))

# 从pQTL数据中筛选出所有cis-SNP
cis_snps_data <- pqtl_locus_data[
  pqtl_locus_data$chromosome == GENE_CHR &
    pqtl_locus_data$base_pair_location >= cis_start &
    pqtl_locus_data$base_pair_location <= cis_end &
    !is.na(pqtl_locus_data$p_value), 
]
if(nrow(cis_snps_data) == 0) {
  stop(paste("错误：在您定义的 cis 区域 (chr 9) 内没有找到任何 SNP。请检查基因坐标。"))
}

# 找到p值最小的 SNP
lead_snp_row <- cis_snps_data[which.min(cis_snps_data$p_value), ]
lead_snp_rsid <- lead_snp_row$rsid
lead_chr <- lead_snp_row$chromosome
lead_pos <- lead_snp_row$base_pair_location
print(paste("已找到 *cis-Lead SNP*: ", lead_snp_rsid, " (P:", lead_snp_row$p_value, ") at chr", lead_chr, ":", lead_pos))

# 以该SNP为界定义上下游区域
start_pos <- lead_pos - window_bp
end_pos <- lead_pos + window_bp
if(start_pos < 0) { start_pos <- 1 }
print(paste("已定义分析窗口: chr", lead_chr, ":", start_pos, "-", end_pos))

# 寻找pQTL和disease中处于该区域的SNP
coord_map <- pqtl_locus_data %>%
  dplyr::select(rsid, chromosome, base_pair_location) %>%
  dplyr::filter(!is.na(rsid)) %>%
  dplyr::distinct(rsid, .keep_all = TRUE)
disease_data_fixed <- inner_join(disease_locus_data, coord_map, by = c("rsids" = "rsid"))
print("正在筛选 pQTL 窗口...")
pqtl_locus_filtered <- pqtl_locus_data %>%
  dplyr::filter(chromosome == lead_chr & base_pair_location >= start_pos & base_pair_location <= end_pos)
print(paste("pQTL 窗口中找到", nrow(pqtl_locus_filtered), "个 SNPs"))
print("正在筛选 Disease 窗口...")
disease_locus_filtered <- disease_data_fixed %>%
  dplyr::filter(chromosome == lead_chr & base_pair_location >= start_pos & base_pair_location <= end_pos)
print(paste("Disease 窗口中找到", nrow(disease_locus_filtered), "个 SNPs"))
print("正在准备和合并数据...")
data1_pre <- data.frame(
  snp_id = pqtl_locus_filtered$rsid,
  beta = pqtl_locus_filtered$beta,
  varbeta = (pqtl_locus_filtered$standard_error)^2,
  N = pqtl_locus_filtered$n,
  MAF = pqtl_locus_filtered$effect_allele_frequency
)
data1_pre <- subset(data1_pre, !is.na(snp_id) & !is.na(varbeta) & !is.na(MAF) & !is.na(N))
data2_pre <- data.frame(
  snp_id = disease_locus_filtered$rsids, 
  beta = disease_locus_filtered$beta,
  varbeta = (disease_locus_filtered$sebeta)^2, 
  MAF = disease_locus_filtered$af_alt        
)
data2_pre <- subset(data2_pre, !is.na(snp_id) & !is.na(varbeta) & !is.na(MAF))
merged_data <- merge(data1_pre, data2_pre, by = "snp_id")

# 预处理及过滤SNP
if (any(duplicated(merged_data$snp_id))) {
  print("警告：存在重复的SNP ID，正在移除...")
  merged_data <- merged_data[!duplicated(merged_data$snp_id), ]
}
if (nrow(merged_data) < 50) { 
  stop(paste("错误：共有SNP数量太少 (", nrow(merged_data), ")。无法进行共定位。"))
}
D1 <- list(
  beta = merged_data$beta.x,
  varbeta = merged_data$varbeta.x,
  N = merged_data$N,
  MAF = merged_data$MAF.x,
  type = "quant",
  snp = merged_data$snp_id
)
D2 <- list(
  beta = merged_data$beta.y,
  varbeta = merged_data$varbeta.y,
  N = N_CASES + N_CONTROLS, 
  s = N_CASES / (N_CASES + N_CONTROLS), 
  MAF = merged_data$MAF.y,
  type = "cc",
  snp = merged_data$snp_id
)

# 执行coloc分析
res_coloc_locus <- coloc.abf(dataset1 = D1, dataset2 = D2)
print(res_coloc_locus)

# 保存关键数据
dat_list <- list(
  results_list = results_list,
  disease_locus_data = disease_locus_data,
  pqtl_locus_data = pqtl_locus_data
)
qs::qsave(dat_list,file = "data/all_data.qs")