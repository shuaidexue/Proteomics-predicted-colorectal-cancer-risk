
###identify colorectal cancer-related differential expression protein in discovery stage
library(ggplot2)
library(reshape2)
library(dplyr)
library(readr)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(impute)
library(reshape2)
library(preprocessCore)
library(data.table)
library(ggsci)
library(openxlsx)

data.quant=read.csv("CRC_protein_quantification_matrix.csv",header=T) 
sample=read_tsv(file = "sample.txt")
sample$ID=paste(sample$Batch,sample$`TMT label`-125,sep="_")

library(plyr)
data.quant <- rename(data.quant,c("Protein.accession"="Protein accession"))
data.quant <- rename(data.quant,c("Gene.name"="Gene name"))
data.quant <- rename(data.quant,c("Protein.description"="Protein description"))

limma.filter=function(data.quant,s,fdr=0.01,fc=1.5){
  data.diff=subset(data.quant,select = c("Protein accession","Gene name","Protein description",s$ID))
  MA=data.diff[rowSums(!is.na(data.diff[s$ID]))>nrow(s)/3,s$ID]
  fit <- lmFit(MA)
  fit <- eBayes(fit)
  out=topTable(fit,n=Inf)
  data.diff=merge(data.diff,subset(out,select=c("logFC","P.Value","adj.P.Val")),by="row.names")
  out.data=volcano.plot(data.diff,s,fdr,fc)
  return(list(data=out.data$data,plot=out.data$plot))
}

volcano.plot=function(data.diff,s,fdr=0.01,fc=1.5,scale.v="none"){
  data.diff$Regulation="unchanged"
  data.diff$Regulation[data.diff$adj.P.Val<fdr&data.diff$logFC>log2(fc)]='up-regulated'
  data.diff$Regulation[data.diff$adj.P.Val<fdr&data.diff$logFC< -log2(fc)]='down-regulated'
  data.diff.filter=subset(data.diff,Regulation!="unchanged")
  #--------- heatmap plot --------------
  p1=pheatmap(data.diff.filter[rowSums(!is.na(data.diff.filter[s$ID]))>nrow(s)/2,s$ID],border_color = 'white',show_colnames = F,silent = T,
              breaks = seq(-2,2,length.out=100),show_rownames = F,cluster_cols = F,scale = scale.v,
              annotation_col = annot_sample,cluster_rows = T,width = 8,height = 10)
  data.diff.s=as.data.frame(table(data.diff$Regulation))
  colnames(data.diff.s)=c("Regulation","Number of proteins")
  data.diff$`Gene name`[data.diff$`Gene name`=="---"]=data.diff$`Protein accession`[data.diff$`Gene name`=="---"]
  data.diff$label=NA
  data.diff$label[data.diff$Regulation=="up-regulated"]=data.diff$`Gene name`[data.diff$Regulation=="up-regulated"]
  maxy=max(-log10(data.diff$adj.P.Val))*1.1
  p2=ggplot(data.diff,aes(x=logFC,y=-log10(adj.P.Val),color=Regulation))+
    geom_point(size=1)+theme_classic2()+xlim(-3,3)+
    scale_color_manual(values = c(unchanged="grey",`up-regulated`="red",`down-regulated`="steelblue"))+
    annotation_custom(grob = tableGrob(data.diff.s,rows = NULL,theme = ttheme_default(base_size = 10)),
                      xmin=3,ymin=maxy*0.9)+ylim(0,maxy)+
    geom_hline(yintercept = -log10(fdr),linetype=2)+
    geom_vline(xintercept = c(-log2(fc),log2(fc)),linetype=2)+
    geom_text_repel(aes(label=label))
  p=ggarrange(p1$gtable,p2,nrow = 1,labels = "AUTO")
  return(list(plot=p,data=data.diff))
}

TL=brewer.pal(3,"Dark2")
names(TL)=c("Left colon","Right colon","Rectum")
annot_col=list(`Tumor location`=TL)

#------- before ----------
sample=subset(sample,sample$`Tumor location`=='Rectum'|sample$`Tumor location`=='Right colon'|sample$`Tumor location`=='Left colon')
cor.data=cor(data.quant[sample$ID[sample$`TMT label`!=129]],use="pairwise.complete.obs")
annot_sample=as.data.frame(subset(sample,select=c("Batch","TMT label","Sex","Age","Tumor location")))
annot_sample$Batch=as.numeric(sub(annot_sample$Batch,pattern = "S0",replacement = ""))
rownames(annot_sample)=sample$ID
cor.tree=pheatmap(cor.data,breaks = seq(-1,1,length.out=100),annotation_col = annot_sample,show_rownames = F,show_colnames = F,border_color = 'white',treeheight_col = 0,cutree_rows = 2,cutree_cols = 2,silent = T,annotation_colors = annot_col)
data.pca=impute.knn(data=as.matrix(data.quant[rowSums(!is.na(data.quant[sample$ID[sample$`TMT label`!=129]]))>75,sample$ID[sample$`TMT label`!=129]]))
pca = prcomp(data.pca$data,scale = TRUE,center = T)
pca.d=as.data.frame(pca$rotation)
pca.d=merge(annot_sample,pca.d[1:2],by="row.names")
fig5.b=ggplot(pca.d,aes(PC1,PC2,color=`Tumor location`))+geom_point(size=3)+theme_classic2()+geom_text_repel(aes(label=Row.names))+scale_color_brewer(type="qual",palette = 2)
fig5=ggarrange(cor.tree$gtable,fig5.b,nrow = 1,labels = 'AUTO')
ggsave(plot = fig5,file="Figure. Pearson correlation and PCA of all sample.pdf",width = 15,height = 6)

#------- after ----------
sample.class=cutree(cor.tree$tree_row,2)
sample.filter=names(sample.class[sample.class==1])
sample.class.d=as.data.frame(sample.class)
sample.class.d$id=row.names(sample.class.d)
write.csv(sample.class.d,"sample.class.d.csv",row.names=F)

cor.data.new=cor(data.quant[sample.filter],use="pairwise.complete.obs")
fig6.a=pheatmap(cor.data.new,breaks = seq(-1,1,length.out=100),annotation_col = annot_sample,show_rownames = F,show_colnames = F,border_color = 'white',treeheight_col = 0,silent = T,annotation_colors = annot_col)
data.pca=impute.knn(data=as.matrix(data.quant[rowSums(!is.na(data.quant[sample.filter]))>75,sample.filter]))
pca = prcomp(data.pca$data,scale = TRUE,center = T)
pca.d=as.data.frame(pca$rotation)
pca.d=merge(annot_sample,pca.d[1:2],by="row.names")
fig6.b=ggplot(pca.d,aes(PC1,PC2,color=`Tumor location`))+geom_point(size=3)+theme_classic2()+geom_text_repel(aes(label=Row.names))+scale_color_brewer(type="qual",palette = 2)
fig6=ggarrange(fig6.a$gtable,fig6.b,nrow = 1,labels = 'AUTO')
ggsave(plot = fig6,file="Figure. Pearson correlation and PCA of QC passed sample.pdf",width = 15,height = 6)

sample.del=names(sample.class[sample.class==2])
fwrite(subset(sample,ID %in% sample.del),row.names = F,sep="\t",file="Delete.sample.xls")
ggplot(subset(sample,ID %in% sample.del),aes(Batch,as.factor(`TMT label`)))+geom_tile(fill='steelblue')+
  theme_bw()+theme(axis.text.x=element_text(angle = 60,hjust=1))


###abundance difference protein
diff.list=list()
s=subset(sample,ID %in% sample.filter)
s$Type="Cancer"

for(i in unique(s$Type)){
  temp.sample=subset(s,Type==i)
  out.data=limma.filter(data.quant,temp.sample)
  write.table(subset(out.data$data,Regulation!="unchanged",
                     select=c("Protein accession","Gene name",
                              "Protein description","logFC",
                              "adj.P.Val","Regulation")),
              sep="\t",quote = F,row.names =F)
  diff.list[[i]]=out.data$data
}

write.table(out.data$data,file = "Result.protein.difference.deprsssion.csv",sep=",",row.names = F)


###validate differential expression protein in UKBB

id=read.csv("protein-id-coding.csv",head=T)

library(openxlsx)
pro <- subset(out.data,adj.P.Val<0.05)
proid=merge(id,pro,by.x='protein',by.y='Protein')
aa=proid[,c("coding")]

library(dplyr)
library(survival)
library(survminer)

load("E:\\UKB-protein-data\\update-ukb-3000pro\\52231\\crc_Qcancer_52231.Rdata")
load("E:\\UKB-protein-data\\update-ukb-3000pro\\olink_data_66354_3000.Rdata")
library(plyr)
head(data)
mydata<-data[,c("f.eid.66354","protein_id","result")]
mydata$id=mydata$f.eid.66354

individual_data<-merg
individual_data$f.eid.66354=individual_data$f.eid
individual_data=individual_data[,c("f.eid.66354","CRC_OS_timefinal","CRC_OS","sex","age.rec")]
individual_data=individual_data[individual_data$f.eid.66354 %in% data$f.eid.66354, ]

result_data <- data.frame(
  protein_id = numeric(0),
  events_number = numeric(0),
  total_number = numeric(0),
  HR = numeric(0),
  LCI = numeric(0),
  HCI = numeric(0),
  P = numeric(0),
  stringsAsFactors = FALSE
)
protein_id <- aa

for (i in 1:253) {
  current_protein_id <- protein_id[i]
  subset_data <- subset(mydata, mydata$protein_id == current_protein_id)
  merge_data <- merge(subset_data, individual_data, by.x = "f.eid.66354", by.y = "f.eid.66354",all.y=T)
  merge_data=as.data.frame(apply(merge_data,2,function(x) as.numeric(as.character(x))))
  df1 <- sapply(merge_data, function(x){
    x[is.na(x)] <- mean(x, na.rm=T)
    x
  })
  merge_data=as.data.frame(df1)
  merge_data=merge_data[!duplicated(merge_data$f.eid.66354), ]
  res.cox <- coxph(Surv(CRC_OS_timefinal, CRC_OS == 1) ~ result, data = merge_data)
  result <- data.frame(
    protein_id = current_protein_id,
    HR = exp(coef(res.cox)["result"]),
    LCI = exp(confint(res.cox)["result", 1]),
    HCI = exp(confint(res.cox)["result", 2]),
    P = summary(res.cox)$coefficients["result", "Pr(>|z|)"],
    events_number = summary(res.cox)$nevent,
    total_number = summary(res.cox)$n,
    stringsAsFactors = FALSE
  )
  result_data <- rbind(result_data, result)
}
result_data$beta=log(result_data$HR)
result_data$se = abs(log(result_data$HR)/qnorm(result_data$P/2)
write.table(result_data,file = "Res.protein.crc.ukb.52231.csv",sep=",",row.names = F)


###identify protein passed validation
data=read.csv("Result.protein.difference.deprsssion.csv",sep=",",row.names = F)

merg=merge(result_data,data,by.x='Gene name',by.y='protein')
dat=subset(merg,(log2FC>0&beta>0&P<0.05)|(log2FC<0&beta<0&P<0.05))
write.csv(dat,"validated.88pro_ukpro_crc_cox.52231.csv",row.names=F)


##LASSO-Cox proportional hazards model further prioritized protein biomarkers as predictors
library(caret)
library(pROC)
library(survival)

load("UKB.crc.88pro.Rdata")
data3=data2
library(glmnet)
library(survival)

y <- Surv(data3$CRC_OS_timefinal, data3$CRC_OS == 1)
x <- data.matrix(data3[, 4:91])

#LASSO
set.seed(12)
lasso <- glmnet(x, y, family = "cox", alpha = 1)
print(lasso)

tiff(file="fig1.lasso-52231.tiff",units='cm',width=18,height=18,res=300,
     compression = "lzw")
plot(lasso, xvar = "lambda", label = TRUE) 
dev.off()

fitCV<- cv.glmnet(x, y, family = "cox", type.measure = "deviance", nfolds = 5)

tiff(file="fig2.lasso-52231.tiff",units='cm',width=18,height=18,res=300,
     compression = "lzw")
plot(fitCV)
dev.off()

fitCV$lambda.1se
coef(fitCV, s = "lambda.1se")

fitCV$lambda.min
coef(fitCV, s = "lambda.min")


##construct protein score
library(openxlsx)
id <- read.xlsx("data.xlsx",sheet='keep2-lasso-1se')
aa=id[,c("coding")]
load("t_olink_data_66354_2923pro_52682.Rdata")

all_data=all_protein[,c("id.66354",aa[1])]
for (i in 2:15) {
  da=all_protein[,c("id.66354",aa[i])]
  all_data=merge(all_data,da,by.x="id.66354",by.y="id.66354")
}

df1 <- sapply(all_data, function(x){
  x[is.na(x)] <- mean(x, na.rm=T)
  x
})
df1=as.data.frame(df1)
all_data=df1

merg=merg33
datt=merg[,c("f.eid","CRC_OS_timefinal","CRC_OS")]

data3=merge(datt,all_data,by.x="f.eid",by.y="id.66354")

model <- coxph(Surv(CRC_OS_timefinal, CRC_OS) ~.,data=data3[,-1]);
summary(model)
print(summary(model)$coef)

predict <- predict(model,type='risk',newdata=data3)
data3$Proteinscore_15pro=predict


##construct QCancerscore
library(caret)
library(pROC)
library(survival)
load("crc_Qcancer_52231_proS_prs.Rdata")
merg=merg33
fold_pre <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~sex+age.rec+Townse_index_imp_median+bmi_imp_median+ethnic+smoke+alcohol+familyhistory_crc+diabetes+uc+
                    bowel_polyps+breast_cancer+uterine_cancer+ovarian_cancer+cervical_cancer+lung_cancer+blood_cancers+oral_cancers+region,data=merg)
fold_pre
fold_predict <- predict(fold_pre,type='risk',newdata=merg)
merg$QCancerscore=fold_predict

fold_pre <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~sex+age.rec+Townse_index_imp_median+bmi_imp_median+ethnic+smoke+alcohol+familyhistory_crc+diabetes+uc+
                    bowel_polyps+breast_cancer+uterine_cancer+ovarian_cancer+cervical_cancer+lung_cancer+blood_cancers+oral_cancers,data=merg)
fold_pre
fold_predict <- predict(fold_pre,type='risk',newdata=merg)
merg$QCancerscore_noregion=fold_predict

merg33=merg
save(merg33,file="crc_Qcancer_52231_proS_prs.Rdata")


#########construct PRS
library(bigsnpr)
library(bigstatsr)
info <- readRDS("map_hm3_plus.rds")		

xx <- snp_readBGI("ukb_c1_22_map_hm3pus52231.bgen.bgi")
snp_id <- with(xx, paste(chromosome, position, allele1, allele2, sep = "_"))
snp_readBGEN("ukb_c1_22_map_hm3pus52231.bgen","ukb_c1_22_map_hm3pus52231",list(snp_id))

obj.bigSNP <- snp_attach("ukb_c1_22_map_hm3pus52231.rds")

G   <- obj.bigSNP$genotypes
(NCORES <- 10)
big_counts(G, ind.col = 1:8)

CHR <- obj.bigSNP$map$chromosome
CHR=as.numeric(CHR)

POS <- obj.bigSNP$map$physical.pos
library(bigreadr)

fam <- bigreadr::fread2("ukb_map_hm3pus52231.sample")[-1,]
fam$CRC_OS=as.numeric(fam$CRC_OS)
y   <- fam$CRC_OS


load("eacrc_rs.Rdata")
library(plyr)
data <- rename(data,c("standard_error"="beta_se"))
data <- rename(data,c("other_allele"="a0"))
data <- rename(data,c("effect_allele"="a1"))
data <- rename(data,c("chromosome"="chr"))
data <- rename(data,c("base_pair_location"="pos"))
sumstats=data

sumstats$n_eff <- 4 / (1 / 100204 + 1 / 154587) 
sumstats$chr=as.numeric(sumstats$chr)

set.seed(10)
ind.val <- sample(nrow(G), 36561)
ind.test <- setdiff(rows_along(G), ind.val)
ind.all=1:nrow(G)

##Matching variants between genotype data and summary statistics
map <- setNames(obj.bigSNP$map[-2], c("chr", "rsid", "pos", "a0", "a1"))
map=map[,-c(6:7)]

library(dplyr)
map=map[,-c(6)]
df_beta <- snp_match(sumstats, map, join_by_pos = F)

###Computing LDpred2 scores genome-wide
library(RhpcBLASctl)
POS2 <- snp_asGeneticPos(CHR, POS, dir = "genetic-maps", ncores = NCORES,type = c("hapmap"))

ind.row <- rows_along(G)
maf <- snp_MAF(G, ind.row = ind.row, ind.col = df_beta$`_NUM_ID_`, ncores = NCORES)
maf_thr <- 1 / sqrt(length(ind.row))  # threshold I like to use
df_beta2 <- df_beta[maf > maf_thr, ]
df_beta=df_beta2

##QC
af_ref <- big_colstats(G, ind.col = df_beta$`_NUM_ID_`, ncores = NCORES)$sum / (2 * nrow(G))
sd_ref <- sqrt(2 * af_ref * (1 - af_ref))
sd_ss <- with(df_beta, 2 / sqrt(n_eff * beta_se^2 + beta^2))
is_bad <-
  sd_ss < (0.5 * sd_ref) | sd_ss > (sd_ref + 0.1) | 
  sd_ss < 0.05 | sd_ref < 0.05

df_beta2 <- df_beta[!is_bad, ]
df_beta=df_beta2

tmp <- tempfile(tmpdir = "./tmp-data-eacrc-52231")

for (chr in 1:22) {
  
  print(chr)
  ind.chr <- which(df_beta$chr == chr)
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  
  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                   infos.pos = POS2[ind.chr2], ncores = NCORES)
  
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}
file.size(corr$sbk) / 1024^3  # file size in GB

(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))
h2_est <- ldsc[["h2"]]

##LDpred2(-grid): grid of models
(h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4))								
(p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))							
(params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))

set.seed(10)
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)

pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])
params$score <- apply(pred_grid[ind.val, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  summary(glm(y[ind.val] ~ x, family = "binomial"))$coef["x", 3]

library(dplyr)
params %>%
  mutate(sparsity = colMeans(beta_grid == 0), id = row_number()) %>%
  arrange(desc(score)) %>%
  mutate_at(c("score", "sparsity"), round, digits = 3) %>%
  slice(1:10)
 
best_beta_grid <- params %>%
  mutate(id = row_number()) %>%
  arrange(desc(score)) %>%
  slice(1) %>%
  print() %>% 
  pull(id) %>% 
  beta_grid[, .]
 
pred_grid_best <- big_prodVec(G, best_beta_grid, ind.row = ind.test,
                    ind.col = df_beta[["_NUM_ID_"]])
pcor(pred_grid_best, y[ind.test], NULL)

##LDpred2-auto: automatic model
coef_shrink <- 0.5 

set.seed(10) 
multi_auto <- snp_ldpred2_auto(
  corr, df_beta, h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.2, length.out = 30), ncores = NCORES,
  use_MLE = FALSE, 
  allow_jump_sign = FALSE, shrink_corr = coef_shrink)
str(multi_auto, max.level = 1)

str(multi_auto[[1]], max.level = 1)

library(ggplot2)
auto <- multi_auto[[1]]  # first chain

(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
pred_auto <- big_prodVec(G, beta_auto, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
pcor(pred_auto, y[ind.test], NULL)

all_pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.all, ind.col = df_beta[["_NUM_ID_"]])
all_pred_grid_best <- big_prodVec(G, best_beta_grid, ind.row = ind.all,
                              ind.col = df_beta[["_NUM_ID_"]])
all_pred_auto <- big_prodVec(G, beta_auto, ind.row = ind.all, ind.col = df_beta[["_NUM_ID_"]])

dd=fam
dd$all_pred_grid_best=all_pred_grid_best
dd$all_pred_auto=all_pred_auto
prs=dd
save(prs,file="genome-wide_crcprs_ldpred2.Rdata")

##construct combinedscore
load("crc_Qcancer_52231_proS_prs.Rdata")
data=merg33
data$QCancerscore_noregion_scale=scale(data$QCancerscore_noregion)
data$all_pred_grid_best_scale=scale(data$all_pred_grid_best)
data$Proteinscore_15pro_scale=scale(data$Proteinscore_15pro)

#combined
library(caret)
library(pROC)
library(survival)

fold_pre <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~QCancerscore_noregion_scale+all_pred_grid_best_scale+Proteinscore_15pro_scale,data=data)
fold_pre
fold_predict <- predict(fold_pre,type='risk',newdata=data)
data$combinedscore=fold_predict
merg33=data
save(merg33,file="crc_Qcancer_52231_proS_prs.Rdata")


##estimate HR
load("crc_Qcancer_52231_proS_prs.Rdata")
library(survival)
tt=merg33
cox2 <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~scale(tt$QCancerscore),data=tt)
cc1=broom::tidy(cox2, exponentiate = T, conf.int = T)

cox2 <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~scale(tt$QCancerscore_noregion),data=tt)
cc2=broom::tidy(cox2, exponentiate = T, conf.int = T)

cox2 <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~scale(tt$Proteinscore_15pro),data=tt)
cc3=broom::tidy(cox2, exponentiate = T, conf.int = T)

cox2 <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~scale(tt$all_pred_grid_best),data=tt)
cc4=broom::tidy(cox2, exponentiate = T, conf.int = T)

cox2 <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~scale(tt$all_pred_auto),data=tt)
cc5=broom::tidy(cox2, exponentiate = T, conf.int = T)

all_resu=rbind(cc1,cc2,cc3,cc4,cc5)
all_resu$case=table(tt$CRC_OS)['1']
all_resu$control=table(tt$CRC_OS)['0']
write.table(all_resu, file = "score-perSD-crc-cox-result.csv",sep=",",row.names = F)


##randomly split data into training and validation cohorts with a 7:3 ratio
load("crc_Qcancer_52231.Rdata")
library(caret)
data=merg
set.seed(2)
dl = createDataPartition(data$CRC_OS,p = 0.7, list = F)
train = data[dl,]
test = data[-dl,]
head(test)
table(test$CRC_OS)
save(train,file="train.crc_Qcancer_52231.Rdata")
save(test,file="test.crc_Qcancer_52231.Rdata")


###estimate AUC in train dataset with 5-fold cross-validation
load("train.crc_Qcancer_52231.Rdata")
load("test.crc_Qcancer_52231.Rdata")
load("crc_Qcancer_52231_proS_prs.Rdata")

library(pROC)
library(caret)
library(survival)
data3=merg33[merg33$f.eid %in% train$f.eid, ]
data4=merg33[merg33$f.eid %in% test$f.eid, ]

#
max=0  
num=0 
set.seed(222)
folds <- createFolds(y=data3$CRC_OS,k=5)
auc_value<-as.numeric()
for(i in 1:5){
  test<- data3[ folds[[i]],] 
  train <- data3[-folds[[i]],] 
  model<- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~QCancerscore_noregion,data=train,x=T)
  fold_predict <- predict(model,type='lp',newdata=test)
  auc_value<- append(auc_value,as.numeric(auc(as.numeric(unlist(test[,c("CRC_OS")])),fold_predict)))
}
num<-which.max(auc_value)

fold_test1 <- data3[folds[[num]],]   
fold_train1 <- data3[-folds[[num]],]
fold_pre1 <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~QCancerscore_noregion,data=fold_train1)##Qcancer
fold_predict1 <- predict(fold_pre1,type='lp',newdata=fold_test1)
roc_curve1 <- roc(as.numeric(fold_test1[,c("CRC_OS")]),fold_predict1)##Qcancer
round(ci(roc_curve1),2)##95%CI


##
max=0  
num=0 
set.seed(222)
folds <- createFolds(y=data3$CRC_OS,k=5)

auc_value<-as.numeric()
for(i in 1:5){
  test<- data3[ folds[[i]],] 
  train <- data3[-folds[[i]],] 
  model<- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~all_pred_grid_best,data=train,x=T)
  fold_predict <- predict(model,type='lp',newdata=test)
  auc_value<- append(auc_value,as.numeric(auc(as.numeric(unlist(test[,c("CRC_OS")])),fold_predict)))
}
num<-which.max(auc_value)

fold_test2 <- data3[folds[[num]],]   
fold_train2 <- data3[-folds[[num]],]
fold_pre2 <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~all_pred_grid_best,data=fold_train2)##PRS
fold_predict2 <- predict(fold_pre2,type='lp',newdata=fold_test2)
roc_curve2 <- roc(as.numeric(fold_test2[,c("CRC_OS")]),fold_predict2)##PRS
round(ci(roc_curve2),2)##95%CI

##
max=0  
num=0 
auc_value<-as.numeric()
set.seed(222)
folds <- createFolds(y=data3$CRC_OS,k=5)
for(i in 1:5){
  test<- data3[ folds[[i]],] 
  train <- data3[-folds[[i]],] 
  model<- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~Proteinscore_15pro,data=train,x=T)
  fold_predict <- predict(model,type='lp',newdata=test)
  auc_value<- append(auc_value,as.numeric(auc(as.numeric(unlist(test[,c("CRC_OS")])),fold_predict)))
}
num<-which.max(auc_value)

fold_test3 <- data3[folds[[num]],]   
fold_train3 <- data3[-folds[[num]],]
fold_pre3 <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~Proteinscore_15pro,data=fold_train3)##protein
fold_predict3 <- predict(fold_pre3,type='lp',newdata=fold_test3)
roc_curve3 <- roc(as.numeric(fold_test3[,c("CRC_OS")]),fold_predict3)##protein
round(ci(roc_curve3),2)##95%CI

##
max=0  
num=0 
auc_value<-as.numeric()
set.seed(222)
folds <- createFolds(y=data3$CRC_OS,k=5)

for(i in 1:5){
  test<- data3[ folds[[i]],]
  train <- data3[-folds[[i]],]
  model<- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~QCancerscore_noregion+all_pred_grid_best+Proteinscore_15pro,data=train,x=T)
  fold_predict <- predict(model,type='lp',newdata=test)
  auc_value<- append(auc_value,as.numeric(auc(as.numeric(unlist(test[,c("CRC_OS")])),fold_predict)))
}
num<-which.max(auc_value)
fold_test4 <- data3[folds[[num]],]   
fold_train4 <- data3[-folds[[num]],]
fold_pre4 <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~QCancerscore_noregion+all_pred_grid_best+Proteinscore_15pro,data=fold_train4)##combined
fold_predict4 <- predict(fold_pre4,type='lp',newdata=fold_test4)
roc_curve4 <- roc(as.numeric(fold_test4[,c("CRC_OS")]),fold_predict4)##combined
round(ci(roc_curve4),2)##95%CI

#
set.seed(1)
roc.test(roc_curve1,roc_curve2,method='bootstrap',boot.n=500)#PRS vs Qcancer
set.seed(1)
roc.test(roc_curve1,roc_curve3,method='bootstrap',boot.n=500)#proS vs Qcancer
set.seed(1)
roc.test(roc_curve1,roc_curve4,method='bootstrap',boot.n=500)#combined vs Qcancer


pdf(file="AUC.crc-QCancerscore_noregion.all_pred_grid_best.Proteinscore_15pro-52231.pdf",width=6,height=6)
plot(roc_curve1,col="#0072B5FF",title=FALSE,lwd=1.5,
     legacy.axes=TRUE 
     ,xlim=c(1,-0.1),ylim=c(0,1)
)
plot(roc_curve2,col="#458B00",title=FALSE,lwd=1.5,add=T,
     legacy.axes=TRUE 
     ,xlim=c(1,-0.1),ylim=c(0,1)
)     
plot(roc_curve3,col="#E18727FF",title=FALSE,lwd=1.5,add=T,
     legacy.axes=TRUE 
     ,xlim=c(1,-0.1),ylim=c(0,1)
)
plot(roc_curve4,col="#990000",title=FALSE,lwd=1.5,add=T,
     legacy.axes=TRUE 
     ,xlim=c(1,-0.1),ylim=c(0,1)
)

legend("bottomright",
       c(paste("QCancer-S: 0.71 (0.66-0.76)"),
         paste("PRS: 0.74 (0.69-0.80)"),
         paste("Pro-S: 0.66 (0.61-0.71)"),
         paste("Combined: 0.79 (0.75-0.84)")
       ),
       col=c("#0072B5FF","#458B00","#E18727FF","#990000"),lwd=1.5,cex=1)
dev.off()


###estimate AUC in test dataset
load("train.crc_Qcancer_52231.Rdata")
load("test.crc_Qcancer_52231.Rdata")
load("crc_Qcancer_52231_proS_prs.Rdata")

library(pROC)
library(caret)
library(survival)
data3=merg33[merg33$f.eid %in% train$f.eid, ]
data4=merg33[merg33$f.eid %in% test$f.eid, ]

max=0  
num=0 
auc_value<-as.numeric()
model<- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~QCancerscore_noregion
  ,data=data3,x=T)
fold_predict <- predict(model,type='lp',newdata=data4)
auc_value<- append(auc_value,as.numeric(auc(as.numeric(unlist(data4[,c("CRC_OS")])),fold_predict)))
print(auc_value)
roc_curve1 <- roc(as.numeric(data4[,c("CRC_OS")]),fold_predict)##Qcancer
round(ci(roc_curve1),2)##95%CI


auc_value<-as.numeric()
model<- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~all_pred_grid_best
  ,data=data3,x=T)
fold_predict <- predict(model,type='lp',newdata=data4)
auc_value<- append(auc_value,as.numeric(auc(as.numeric(unlist(data4[,c("CRC_OS")])),fold_predict)))
print(auc_value)
roc_curve2 <- roc(as.numeric(data4[,c("CRC_OS")]),fold_predict)##PRS
round(ci(roc_curve2),2)##95%CI

auc_value<-as.numeric()
model<- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~Proteinscore_15pro
              ,data=data3,x=T)
fold_predict <- predict(model,type='lp',newdata=data4)
auc_value<- append(auc_value,as.numeric(auc(as.numeric(unlist(data4[,c("CRC_OS")])),fold_predict)))
print(auc_value)
roc_curve3 <- roc(as.numeric(data4[,c("CRC_OS")]),fold_predict)##ProS
round(ci(roc_curve3),2)##95%CI

auc_value<-as.numeric()
model<- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~QCancerscore_noregion+all_pred_grid_best+Proteinscore_15pro
  ,data=data3,x=T)
fold_predict <- predict(model,type='lp',newdata=data4)
auc_value<- append(auc_value,as.numeric(auc(as.numeric(unlist(data4[,c("CRC_OS")])),fold_predict)))
print(auc_value)
roc_curve4 <- roc(as.numeric(data4[,c("CRC_OS")]),fold_predict)##combined
round(ci(roc_curve4),2)##95%CI

##
set.seed(1)
roc.test(roc_curve1,roc_curve2,method='bootstrap',boot.n=500)#PRS vs Qcancer
set.seed(1)
roc.test(roc_curve1,roc_curve3,method='bootstrap',boot.n=500)#ProS vs Qcancer
set.seed(1)
roc.test(roc_curve1,roc_curve4,method='bootstrap',boot.n=500)#Com vs Qcancer


###calibration
library(pec)
library(rms)
library(survival)
library(survminer)

load("train.crc_Qcancer_52231.Rdata")
load("test.crc_Qcancer_52231.Rdata")
load("crc_Qcancer_52231_proS_prs.Rdata")

library(pROC)
library(caret)
library(survival)
data1=merg33[merg33$f.eid %in% train$f.eid, ]
data2=merg33[merg33$f.eid %in% test$f.eid, ]


cox1 <- coxph(Surv(CRC_OS_timefinal,CRC_OS==1)~QCancerscore_noregion,
              data=data1,x =TRUE,y =TRUE)
cox2 <- coxph(Surv(CRC_OS_timefinal,CRC_OS==1)~all_pred_grid_best,
              data=data1,x =TRUE,y =TRUE)
cox3 <- coxph(Surv(CRC_OS_timefinal,CRC_OS==1)~Proteinscore_15pro,
              data=data1,x =TRUE,y =TRUE)
cox4 <- coxph(Surv(CRC_OS_timefinal,CRC_OS==1)~QCancerscore_noregion+all_pred_grid_best+Proteinscore_15pro,
              data=data1,x =TRUE,y =TRUE)

library(riskRegression)
cox_score <- Score(list("QCancer-S"=cox1,
                        "PRS"=cox2,
                        "ProS"=cox3,
                        "Combined"=cox4),
                   formula=Surv(CRC_OS_timefinal,CRC_OS)~1,
                   data=data1,
                   plots="cal"
)

tiff(file="train.calibration-52231.tiff",units='cm',width=14,height=13,res=300,
     compression = "lzw")
plotCalibration(cox_score,method="quantile",cens.method="local",xlim = c(0,0.06), ylim = c(0, 0.06),
                auc.in.legend=F,
                brier.in.legend=F,
                legend=T,
                lty=c(1,2,4,3),
                type="l",
                xlab = "Predicted Risk",
                ylab = "Observerd Risk",
                col=c("#0072B5FF","#458B00","#E18727FF","#990000"))
dev.off()

cox_score <- Score(list("QCancer-S"=cox1,
                        "PRS"=cox2,
                        "ProS"=cox3,
                        "Combined"=cox4),
                   formula=Surv(CRC_OS_timefinal,CRC_OS)~1,
                   data=data2,
                   plots="cal"
)

tiff(file="test.calibration-52231.tiff",units='cm',width=14,height=13,res=300,
     compression = "lzw")
plotCalibration(cox_score,method="quantile",cens.method="local",xlim = c(0,0.06), ylim = c(0, 0.06),
                auc.in.legend=F,
                brier.in.legend=F,
                legend=T,
                lty=c(1,2,4,3),
                type="l",
                xlab = "Predicted Risk",
                ylab = "Observerd Risk",
                col=c("#0072B5FF","#458B00","#E18727FF","#990000"))
dev.off()

###nomogram
library(regplot)
library(rms)
load("crc_Qcancer_52231_proS_prs.Rdata")
data1=merg33
coxfit <- cph(Surv(CRC_OS_timefinal,CRC_OS) ~ QCancerscore_noregion+all_pred_grid_best+Proteinscore_15pro,
              data = data1, x=T,y=T,surv = T
)

dd <- datadist(data1)
options(datadist = "dd")

f2 <- psm(Surv(CRC_OS_timefinal,CRC_OS==1) ~QCancerscore_noregion+all_pred_grid_best+Proteinscore_15pro
          , data =  data1, dist='lognormal') 

med <- Quantile(f2)
surv <- Survival(f2)

nom <- nomogram(f2, fun=list(function(x) surv(365.25*5, x),
                             function(x) surv(365.25*10, x),
                             function(x) surv(365.25*15, x)),
                lp=F,
                funlabel=c("5 years CRC-free",
                           "10 years CRC-free",
                           "15 years CRC-free" ))

pdf(file="nomogram-52231.pdf",width=13,height=6)

plot(nom, xfrac=.2,col.grid=c("#BC3C29FF","#0072B5FF"))
dev.off()

##Kaplan-Meier curves
library(ggplot2)
library("survival")
library(survminer)
load("train.crc_Qcancer_52231.Rdata")
load("test.crc_Qcancer_52231.Rdata")
load("crc_Qcancer_52231_proS_prs.Rdata")
data1=merg33[merg33$f.eid %in% train$f.eid, ]

fit=survfit(Surv(CRC_OS_timefinal/365.25,CRC_OS)~combinedscore_4.3, data=data1)
fit

tiff(file="train.combinedscore_4.3_KM-52231.tiff",units='cm',width=15,height=15,res=300,
     compression = "lzw")
ggsurvplot(fit,
           fun = "cumhaz", 
           conf.int = TRUE,
           pval = TRUE,
           pval.method = TRUE,
           pval.size=4,
           pval.method.size=4,
           pval.coord=c(2.5,0.02),
           pval.method.coord=c(0.1,0.02),
           legend.labs = c("Low","Medium","High"),
           xlim=c(0,15.5),
           break.x.by=3, 
           risk.table = T,
           legend.title = "",
           xlab="Follow-up time (years)",
           ylab="Cumulative hazard of CRC",
           palette = c("green4","#2f5688","#CC0000"),
           ggtheme = theme_bw()
)
dev.off()

data1=merg33[merg33$f.eid %in% test$f.eid, ]
fit=survfit(Surv(CRC_OS_timefinal/365.25,CRC_OS)~combinedscore_4.3, data=data1)
fit

tiff(file="test.combinedscore_4.3_KM-52231.tiff",units='cm',width=15,height=15,res=300,
     compression = "lzw")
ggsurvplot(fit,
           fun = "cumhaz", 
           conf.int = TRUE,
           pval = TRUE,
           pval.method = TRUE,
           pval.size=4,
           pval.method.size=4,
           pval.coord=c(2.5,0.02),
           pval.method.coord=c(0.1,0.02),
           legend.labs = c("Low","Medium","High"),
           xlim=c(0,15.5),
           break.x.by=3, 
           risk.table = T,
           legend.title = "",
           xlab="Follow-up time (years)",
           ylab="Cumulative hazard of CRC",
           palette = c("green4","#2f5688","#CC0000"),
           ggtheme = theme_bw()
)
dev.off()

###DCA analysis
library(rmda)
library(ggDCA)
library(ggplot2)
library(ggprism)
library(rms)
library(caret)
library(survival)

load("crc_Qcancer_52231_proS_prs.Rdata")
data4=merg33

model1<- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~QCancerscore_noregion
              ,data=data4,x=T)
model2<- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~QCancerscore_noregion+all_pred_grid_best+Proteinscore_15pro
              ,data=data4,x=T)
dca_cph <- dca(model1, model2, model.names = c("QCancer-S", "Combined"))
save(dca_cph,file="dca_cph_result_52231.Rdata")

###plot
data=dca_cph
pdf(file="crc-DCA-52231.pdf",width=6,height=5.5)
ggplot(data,
       linetype =F,
       lwd = 0.4)+theme_classic()+  
  theme_prism(base_size =10)+
  theme(legend.position=c(0.8,0.6))+
  scale_x_continuous(
    limits = c(0, 0.1),
    guide = "prism_minor"
  ) +
  scale_y_continuous(
    limits = c(-0.001, 0.015),
    guide = "prism_minor"
  )+
  scale_colour_prism(
    palette = "candy_bright",
    name = "Cylinders",
    label = c("QCancerscore", "Combined","ALL","None")
  )+labs(x = "Risk Threshold", y = "Net Benefit",title = "")+scale_colour_manual(values = c("#0072B5FF","#990000","black","gray67"))

dev.off()

###RAP analysis
library(boot)
library(survival)
library(foreign)
library(rms)
load("crc_Qcancer_52231_proS_prs.Rdata")
data=merg33

RAP <- function(data,indices){
  dat <- data[indices,]
  cox <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~relevel(all_pred_grid_best_3, ref="Q2")+relevel(Proteinscore_15pro_3, ref="Q2")+age.rec+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=dat)
  cox[["coefficients"]][[1]]/cox[["coefficients"]][["age.rec"]]
  
}
set.seed(1)
results <- boot(data=data, statistic=RAP, R=500)
print(results)

###
RAP <- function(data,indices){
  dat <- data[indices,]
  cox <- coxph(Surv(CRC_OS_timefinal,CRC_OS) ~relevel(combinedscore_3, ref="Q2")+age.rec+sex+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data=dat)
  
  cox[["coefficients"]][[1]]/cox[["coefficients"]][["age.rec"]]
  
}
set.seed(1)
results <- boot(data=data, statistic=RAP, R=500)
print(results)


###10-year cumulative risk
##10-year cumulative risk across combinedscore (three groups)
##combinedscore_Q1
load("crc_Qcancer_52231_proS_prs.Rdata")
data3=merg33
data=subset(data3,combinedscore_3=='Q1')

result=data.frame(matrix(nrow = 0, ncol = 1))
aa=unique(data$age.rec)
for (i in aa) {
  dd=subset(data,age.rec==i)
  as.data.frame(table(dd$CRC_OS)["1"]/(sum(dd$CRC_OS_timefinal)/365.25))
  result=rbind(result, as.data.frame(table(dd$CRC_OS)["1"]/(sum(dd$CRC_OS_timefinal)/365.25)))
}
result$age=aa

library(plyr)
result<- rename(result,c('table(dd$CRC_OS)["1"]/(sum(dd$CRC_OS_timefinal)/365.25)'="age_specific_annual_incidence_rate"))
result$age_specific_annual_incidence_rate=ifelse(is.na(result$age_specific_annual_incidence_rate),0,result$age_specific_annual_incidence_rate)
result = result[order(result[,"age"]),]
write.table(result,file = "combinedscore_3_Q1-age-specific annual incidence rate_52231.csv",sep=",",row.names = F)

dat=read.csv("combinedscore_3_Q1-age-specific annual incidence rate_52231.csv",header=T)
dat$age=as.numeric(dat$age)
head(dat)

aa= c(39:61)
result=data.frame(matrix(nrow = 0, ncol = 1))

for (i in aa) {
  dd=subset(dat,age>=i&age<10+i)
  result=rbind(result, as.data.frame(sum(dd$age_specific_annual_incidence_rate)))
} 
result$age=aa
library(plyr)
result<- rename(result,c("sum(dd$age_specific_annual_incidence_rate)"="incidence_rate_10year"))

result$risk_10year=1-exp(-result$incidence_rate_10year)
result$risk_10year_100=100*result$risk_10year
result$group='combinedscore_3_Q1'
write.table(result,file = "combinedscore_3_Q1_age-specific 10-year cumulative incidence rate and risk_52231.csv",sep=",",row.names = F)

##
##combinedscore_Q2
load("crc_Qcancer_52231_proS_prs.Rdata")
data3=merg33
data=subset(data3,combinedscore_3=='Q2')

result=data.frame(matrix(nrow = 0, ncol = 1))
aa=unique(data$age.rec)
for (i in aa) {
  dd=subset(data,age.rec==i)
  as.data.frame(table(dd$CRC_OS)["1"]/(sum(dd$CRC_OS_timefinal)/365.25))
  result=rbind(result, as.data.frame(table(dd$CRC_OS)["1"]/(sum(dd$CRC_OS_timefinal)/365.25)))
}
result$age=aa

library(plyr)
result<- rename(result,c('table(dd$CRC_OS)["1"]/(sum(dd$CRC_OS_timefinal)/365.25)'="age_specific_annual_incidence_rate"))
result$age_specific_annual_incidence_rate=ifelse(is.na(result$age_specific_annual_incidence_rate),0,result$age_specific_annual_incidence_rate)
result = result[order(result[,"age"]),]
write.table(result,file = "combinedscore_3_Q2-age-specific annual incidence rate_52231.csv",sep=",",row.names = F)

dat=read.csv("combinedscore_3_Q2-age-specific annual incidence rate_52231.csv",header=T)
dat$age=as.numeric(dat$age)
head(dat)

aa= c(40:61)
result=data.frame(matrix(nrow = 0, ncol = 1))

for (i in aa) {
  dd=subset(dat,age>=i&age<10+i)
  result=rbind(result, as.data.frame(sum(dd$age_specific_annual_incidence_rate)))
} 
result$age=aa
library(plyr)
result<- rename(result,c("sum(dd$age_specific_annual_incidence_rate)"="incidence_rate_10year"))

result$risk_10year=1-exp(-result$incidence_rate_10year)
result$risk_10year_100=100*result$risk_10year
result$group='combinedscore_3_Q2'
write.table(result,file = "combinedscore_3_Q2_age-specific 10-year cumulative incidence rate and risk_52231.csv",sep=",",row.names = F)

##
##combinedscore_Q3
load("crc_Qcancer_52231_proS_prs.Rdata")
data3=merg33
data=subset(data3,combinedscore_3=='Q3')

result=data.frame(matrix(nrow = 0, ncol = 1))
aa=unique(data$age.rec)
for (i in aa) {
  dd=subset(data,age.rec==i)
  as.data.frame(table(dd$CRC_OS)["1"]/(sum(dd$CRC_OS_timefinal)/365.25))
  result=rbind(result, as.data.frame(table(dd$CRC_OS)["1"]/(sum(dd$CRC_OS_timefinal)/365.25)))
}
result$age=aa


library(plyr)
result<- rename(result,c('table(dd$CRC_OS)["1"]/(sum(dd$CRC_OS_timefinal)/365.25)'="age_specific_annual_incidence_rate"))
result$age_specific_annual_incidence_rate=ifelse(is.na(result$age_specific_annual_incidence_rate),0,result$age_specific_annual_incidence_rate)
result = result[order(result[,"age"]),]

write.table(result,file = "combinedscore_3_Q3-age-specific annual incidence rate_52231.csv",sep=",",row.names = F)


dat=read.csv("combinedscore_3_Q3-age-specific annual incidence rate_52231.csv",header=T)
dat$age=as.numeric(dat$age)
head(dat)

aa= c(40:61)
result=data.frame(matrix(nrow = 0, ncol = 1))

for (i in aa) {
  dd=subset(dat,age>=i&age<10+i)
  result=rbind(result, as.data.frame(sum(dd$age_specific_annual_incidence_rate)))
} 
result$age=aa
library(plyr)
result<- rename(result,c("sum(dd$age_specific_annual_incidence_rate)"="incidence_rate_10year"))

result$risk_10year=1-exp(-result$incidence_rate_10year)
result$risk_10year_100=100*result$risk_10year
result$group='combinedscore_3_Q3'
write.table(result,file = "combinedscore_3_Q3_age-specific 10-year cumulative incidence rate and risk_52231.csv",sep=",",row.names = F)


####
dat1=read.csv("age-specific 10-year cumulative incidence rate and risk_52231.csv",header=T)
dat2=read.csv("combinedscore_3_Q1_age-specific 10-year cumulative incidence rate and risk_52231.csv",header=T)
dat3=read.csv("combinedscore_3_Q2_age-specific 10-year cumulative incidence rate and risk_52231.csv",header=T)
dat4=read.csv("combinedscore_3_Q3_age-specific 10-year cumulative incidence rate and risk_52231.csv",header=T)

dd=rbind(dat1,dat2,dat3,dat4)
table(dd$group)
write.table(dd,file = "summary.combinedscore_3.age-specific 10-year cumulative incidence rate and risk_52231.csv",sep=",",row.names = F)

##plot
library(ggplot2)
dat=read.csv("summary.combinedscore_3.age-specific 10-year cumulative incidence rate and risk_52231.csv",header=T)

tiff(file="combinedscore_3-10-year-risk-52231.tiff",units='cm',width=15,height=11,res=300,
     compression = "lzw")
ggplot(dat, aes(age, risk_10year_100, linetype = Group, color = Group))+
  geom_point(size =0)+
  geom_line(size =0.65)+ 
  labs(
    x = "Age (years)",
    y = "10-year cumulative risk (%)")+

  scale_x_continuous(
    limits = c(40, 60)
  ) +
  scale_linetype_manual(values=c("solid", "solid", "solid","solid"))+
  theme_bw() +scale_color_manual(values = c('#fdbd10',"#ec1c24","#458B00",'#0066b2'))+geom_hline(aes(yintercept=0.867478772135688), size = 0.4,                                                                                             
                                                                                                 color = "black", linetype = "dashed")
dev.off()

##the calculation of 10-year cumulative risk for PRS (three groups) and ProS(three groups) was consistent with combinedscore 

