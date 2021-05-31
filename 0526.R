.libPaths( c( .libPaths(), "D:/R_BOOKS/R/library") )
setwd('D:/science/course/Bioinformatics')

####DATA READ & PREPROCESSING####
library(readr)
Metabolome <- read_delim("D:/science/course/Bioinformatics/Metabolome.txt", 
                         "\t", escape_double = FALSE, locale = locale(encoding = "GB2312"), 
                         trim_ws = TRUE)
Meta <- Metabolome[,6:101]
Meta <- as.data.frame(lapply(Meta, as.numeric))
#1 Normalization
NormMeta <- scale(Meta,center = T,scale = T)
NormMeta_rownames <-Metabolome$Compounds
NormMeta_colnames <-colnames(Metabolome)[6:101]
NormMeta <- as.data.frame(NormMeta,row.names = NormMeta_rownames)
colnames(NormMeta) <- NormMeta_colnames
####Analysis#####
#2.1 calculate LOG2(Fold change) 
NormMeta$H_mean <- rowMeans(NormMeta_m[,1:10])
NormMeta$M_mean <- rowMeans(NormMeta_m[,69:96])
NormMeta$S_mean <- rowMeans(NormMeta_m[,47:68])
NormMeta$F_mean <- rowMeans(NormMeta_m[,11:46])
NormMeta$M_VS_H <- log2(abs(NormMeta$M_mean/NormMeta$H_mean))
NormMeta$S_VS_H <- log2(abs(NormMeta$S_mean/NormMeta$H_mean))
NormMeta$F_VS_H <- log2(abs(NormMeta$F_mean/NormMeta$H_mean))
NormMeta <- as.data.frame(NormMeta)
#2.2 calculate P-value
library(car)
#注：采用shapiro test发现不满足正态性假设且Levene test发现不满足方差齐性假设
j=1
  for (j in 1:431) {
    NormMeta$M_VS_H_Pvalue[j] <- wilcox.test(as.numeric(NormMeta[j,1:10]),
                                              as.numeric(NormMeta[j,69:96]))$p.value  
    
  }

j=1
for (j in 1:431) {
  NormMeta$S_VS_H_Pvalue[j] <- wilcox.test(as.numeric(NormMeta[j,1:10]),
                                           as.numeric(NormMeta[j,47:68]))$p.value  
  
}

j=1
for (j in 1:431) {
  NormMeta$F_VS_H_Pvalue[j] <- wilcox.test(as.numeric(NormMeta[j,1:10]),
                                           as.numeric(NormMeta[j,11:46]))$p.value  
  
}

#2.3 Heat map
library(pheatmap)
library(ggplot2)
NormMeta <- NormMeta[order(-abs(NormMeta$F_VS_H)),]
NormMeta_m <- as.matrix(NormMeta[,101:103])
NormMeta_m <- NormMeta_m[1:50,]
heatmap1 <- pheatmap(NormMeta_m,border='black',
                     scale = "none",
                     labels_col = c('Mild',
                                    'Severe',
                                    'Fatal'),
                     angle_col  = 0,
                     cluster_rows = F,
                     cluster_cols = F,
                     legend_breaks = c(-3,0,3),
                     cellwidth = 45,
                     cellheight = 12,
                     display_numbers = T,
                     legend = F,
                     fontsize_row = 12
                     )

NormMeta_m2 <- as.matrix(NormMeta[,104:106])
NormMeta_m2 <- NormMeta_m2[1:50,]
heatmap2 <- pheatmap(NormMeta_m2,border='black',
                     scale = "none",
                     labels_col = c('Mild',
                                    'Severe',
                                    'Fatal'),
                     angle_col  = 0,
                     cluster_rows = F,
                     cluster_cols = F,
                     cellwidth = 45,
                     cellheight = 12,
                     display_numbers = matrix(
                       ifelse(NormMeta_m2<0.001,'***','*'),nrow(NormMeta_m2)
                     ),
                     number_color = 'black',
                     fontsize_number = 8,
                    legend = F,
                     fontsize_row = 12,
                    labels_row = ''
)


p1 <- ggplotify::as.ggplot(heatmap1)
p2 <- ggplotify::as.ggplot(heatmap2)
cowplot::plot_grid(p2,p1,hjust = 1,vjust = 0.5,
                   labels = c('Fig 1.1','Fig 1.2'),
                   align = 'h',rel_widths = c(1.3,1.3))->p3
p3 <- p3+theme(text = element_text(family = 'Times New Roman'),
               plot.margin = unit(c(1,1,1,1),'cm')
  
)
tiff(file = 'Heatmap.tiff', width = 15, height = 45, units = "in", 
     compression = "lzw",res = 200)
p3
dev.off()
library(dplyr)
library(tidyverse)
stat_dat <- Metabolome[,6:101]
stat_dat <- as.data.frame(t(stat_dat))
colnames(stat_dat) <- Metabolome$Compounds
stat_dat$Samples <-rownames(stat_dat) 
stat_dat$Label <- case_when(str_detect(stat_dat$Samples,'H')==T~0,
                            str_detect(stat_dat$Samples,'M')==T~1,
                            str_detect(stat_dat$Samples,'S')==T~2,
                            str_detect(stat_dat$Samples,'F')==T~3)
stat_dat1 <- data.frame(stat_dat1$Samples,stat_dat1$Label,stat_dat1[,1:431])
colnames(stat_dat1)[1:2] <- c('Samples','Label')
write.csv(stat_dat1,'Statistical analysis.csv')
#enrichment analysis
mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-c("L-Aspartic Acid","O-Acetyl-L-serine","Ginkgoic acid","2-Nonanone","L-Tyrosine","L-Lysine","L-Phenylalanine","L-Proline","3-Hydroxy-3-Methylpentane-1,5-Dioic Acid","Serotonin","Urea","Phthalic Acid","4-Nitrophenol","Inositol","D-Fructose","D-Mannose","L-Ascorbate","6-Aminocaproic Acid","Homovanillic Acid","Pantetheine","6-Phosphogluconic Acid Trisodium Salt","Phosphoric Acid","6-Methylmercaptopurine","2-Hydroxycinnamic acid","2-Pyrrolidinone","2-(3,4-dimethoxyphenyl)ethanamine","Isonicotinic acid","1,2-Dichloroethane","2-Picoline","2,2,2-Trichloroethanol","Methylcysteine","2-Acetylfuran","p-Cymene","DL-Leucine","Dibutyl phthalate","FURFURYL ALCOHOL","1-Phenylethanol","Choline chloride","Phosphocholine")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "blood", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "net", "png", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "png", 72, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_1_", "net", "png", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_1_", "png", 72, width=NA)
#pathway analysis
library('MetaboAnalystR')
mSet<-InitDataObjects("conc", "pathora", FALSE)
cmpd.vec<-c("L-Cystine","L-Threonine","L-Arginine","L-Aspartic Acid","L-Citrulline","L-Glutamic Acid","L-Isoleucine","L-Pyroglutamic Acid","L-Serine","Allantoin","Dopamine","Hexanoyl Glycine","L-Asparagine Anhydrous","L-Carnosine","L-Cysteine","L-Glutamine","L-Homocitrulline","N-Acetyl-L-Glutamic Acid","N-Acetyl-L-Leucine","N-Acetylneuraminic Acid","N-Isovaleroylglycine","N-Propionylglycine","O-Phospho-L-Serine","S-Sulfo-L-Cysteine","Trans-4-Hydroxy-L-Proline","Benzoic Acid","Terephthalic Acid","(3,4-Dimethoxyphenyl) Acetic Acid","2-Methoxybenzoic Acid","3,4-Dihydroxybenzeneacetic Acid","3-Hydroxyanthranilic Acid","P–Hydroxyphenyl Acetic Acid","Phenethylamine","Dulcitol","1-Methylxanthine","2,6-Dihydroxypurine","5-Methylcytosine","Adenosine 5'-Monophosphate","Hypoxanthine","Thymidine","Β-Nicotinamide Mononucleotide","Epinephrine","L-Thyroxine","Pyrroloquinoline Quinone","Citric Acid","L-Malic Acid","Succinic Acid","Cis-Aconitic Acid","Pyruvic Acid","Citramalic Acid","Tryptamine","1,5-Anhydro-D-Glucitol","D-Trehalose","D-Erythronolactone","D-Arabinose","L-Fucose","L-Rhamnose","D-Glucoronic Acid","D-Glyceric Acid","L-Gulonic-Γ-Lactone","Orotic Acid","Pantothenate","3-Indolebutyric Acid","Indole-3-Carboxaldehyde","2,6-Diaminooimelic Acid","2-Aminoethanesulfinic Acid","2-Aminoethanesulfonic Acid","2-Hydroxybutanoic Acid","2-Hydroxyisocaproic Acid","3-Hydroxy-3-Methyl Butyric Acid","3-Hydroxybutyrate","3-Methylcrotonyl Glycine","4-Hydroxy-2-Oxoglutaric Acid","Adipic Acid","Creatine","Dodecanedioic Aicd","Kinic Acid","L-3-Phenyllactic Acid","L-Lactic Acid","L-Dihydroorotic Acid","Malonicacid","Mandelic Acid","Methylmalonic Acid","N-Γ-Acetyl-N-2-Formyl-5-Methoxykynurenamine","Phenyllactate (Pla)","Phenylpyruvic Acid","Pyrrole-2-Carboxylic Acid","Shikimic Acid","Subericacid","Neopterin","4-(Aminomethyl)-5-(Hydroxymethyl)-2-Methylpyridin-3-Ol","Ureidoisobutyric Acid","N-Acetylglycine","Thiamine Triphosphate","N-Acetylthreonine","Hypoxanthine-9-β-D-Arabinofuranoside","D-Sedoheptuiose 7-Phosphate","Oxoadipic Acid","3-Methoxy-4-Hydroxyphenylethyleneglycol Sulfate","1-Methyluric Acid","Aminomalonic Acid","3-Methyl-2-Oxobutanoic Acid","Allysine","N-Acetylmethionine","Ribulose-5-Phosphate","Argininosuccinic acid","2-Deoxyribose 1-Phosphate","N-Acetylglucosamine 1-Phosphate","Estrone","(3-Methoxy-4-hydroxyphenyl)ethylene glycol sulfate","Xanthosine","estrone 3-sulfate","UDP-glucose","3'-Sialyllactose","1-Methylguanine","DL-3,4-Dihydroxyphenyl glycol","dihydrotachysterol","Indoleacrylic acid","2-Methylbenzoic acid","Hydroxyphenyllactic acid","INDOLE-3-CARBINOL","DL-o-Tyrosine","L-Erythrulose","N-lactoyl-phenylalanine","Cyclamic acid","D-(+)-Malic acid","deoxyguanosine 5'-monophosphate (dGMP)","5-Methoxytryptophol","Uridine triphosphate (UTP)","6β-hydroxytestosterone","Carbamoyl phosphate","O-Acetyl-L-serine","Indoxylsulfuric acid","Salicyluric acid","Formononetin","Ellagic acid","N-Oleoyl Glycine","Ginkgoic acid","Hydroquinone","Anthranilic acid","Succinic anhydride","3-Methylsalicylic acid","Indoleacetaldehyde","D-Xylulose 5-phosphate","Vanillylmandelic acid","Gamma-Glu-Leu","O-Desmethylnaproxen","Ethylsalicylate","Maltol","2-nonanol","1-Naphthol","2-n-Pentylfuran","2-Methyl-5-nitroimidazole-1-ethanol","Aspirin","Sorbic acid","Phenoxyacetic acid","4-Ethylphenol","4-tert-Octylphenol","PHENYL-BETA-D-GLUCOPYRANOSIDE","2-(Methylthio)ethanol","4-Hydroxybenzyl alcohol","2-Nonanone","Propylparaben","2,4-Di-tert-butylphenol","Methylparaben","Taurodeoxycholic acid sodium salt","Xylose","4-Hydroxyhippurate","octadecanedioate (C18)","N-acetyl-beta-alanine","Scyllo inositol","TRP-GLU","N-acetylornithine","N-Amidino-L-Aspartate","Dethiobiotin","Mono-Methyl Glutarate","Phenylpyruvate","D-Homocysteine","Glycine","L-Tyrosine","L-Lysine","L-Ornithine","L-Alanine","L-Histidine","L-Leucine","L-Methionine","L-Phenylalanine","L-Proline","L-Tryptophan","L-Valine","(5-L-Glutamyl)-L-Amino Acid","3-Chloro-L-Tyrosine","3-Hydroxy-3-Methylpentane-1,5-Dioic Acid","5-Hydroxy-L-Tryptophan","5-Oxoproline","Asp-Phe","Betaine","Glutathione Reducedform","L-Cystathionine","L-Dopa","L-Saccharopine","L-Theanine","Methionine Sulfoxide","N,N-Dimethylglycine","N6-Acetyl-L-Lysine","N-Acetyl-L-Tyrosine","N-Acetylmannosamine","N-Glycyl-L-Leucine","Nα-Acetyl-L-Arginine","Phenylacetyl-L-Glutamine","Phe-Phe","S-(5-Adenosy)-L-Homocysteine","Serotonin","Urea","L-Alanyl-L-Lysine","N-Acetylhistamine","Methyleugenol","P-Coumaric Acid","Syringic Acid","4-Ethylbenzoic Acid","4-Hydroxybenzoic Acid","Phthalic Acid","1,4-Dihydro-1-Methyl-4-Oxo-3-Pyridinecarboxamide","4-Pyridoxic Acid","Theobromine","Acetylcholine Chloride","Choline","Putrescine","Diethanolamine","4-Nitrophenol","Inositol","Myoinositol","Vanillin","1,7-Dimethylxanthine","1-Methyladenine","1-Methylhistidine","3'-Aenylic Acid","5,6-Dihydro-5-Methyluracil","5-Methyluridine","7-Methylxanthine","Adenine","Cytidine-5-Monophosphate","Cytosine","Guanosine Monophosphate","Inosine","Purine","Thymine","Uracil","Uridine","3,3',5-Triiodo-L-Thyronine","Norepinephrine","D-Arabitol","D-Mannitol","D-Fructose","D-Mannose","D-Melezitose","Lactose","Maltose","N-Acetyl-D-Glucosamine","Raffinose","Gluconic Acid","L-Ascorbate","Nicotinamide","4-Oxoretinol","Nicotinic Acid","Nicotinuric Acid","Pantothenol","Riboflavin","Trigonelline","3-Indolepropionic Acid","Indole-3-Acetic Acid","Methyl Indole-3-Acetate","3,4,5-Trimethoxycinnamic Acid","4-Acetamidobutyric Acid","4-Guanidinobutyric Acid","5-Aminovaleric Acid","6-Aminocaproic Acid","7-Methyluric Acid","Azelaic Acid","Chlorogenic Acid","Creatinine","Dl-2-Aminooctanoic Acid","D-Pipecolinic Acid","Guanidineacetic Acid","Hippuric Acid","Homovanillic Acid","Hydrocinnamic Acid","Kinurenine","Kynurenic Acid","L-Homoserine","Maleic Acid","N'-Formylkynurenine","P-Aminohippuric Acid","Cinnamic Acid","Uric Acid","4-Hydroxybenzaldehyde","Taurocholic Acid Sodium Salt Hydrate","2-(Formylamino)Benzoic Acid","Acetaminophen Glucuronide","Cyclic Amp","N6-Succinyl Adenosine","3-Hydroxyhippuric Acid","2-(Dimethylamino)Guanosine","7-Methylguanine","Β-Pseudouridine","Isoxanthopterin","H-Homoarg-Oh","Pantetheine","4-Hydroxyretinoic Acid","L-Pipecolic Acid","5'-Deoxy-5'-(Methylthio) Adenosine","Sarcosine","3-Aminoisobutanoic Acid","Cys-Gly","Trans-Citridic Acid","6-Phosphogluconic Acid Trisodium Salt","Phosphoric Acid","2-Aminoadipic Acid","Hordenine","6-Methylmercaptopurine","Ergothioneine","N-Formylmethionine","Methoxyindoleacetic Acid","Sn-Glycero-3-Phosphocholine","3-Hydroxykynurenine","P-Aminobenzoate","Indole","Isoquinoline","Piperidine","5-Aminolevulinate","2-Methylguanosine","Glucosamine","2-Hydroxycinnamic acid","trans-3-Indoleacrylic acid","DL-Stachydrine","L-Norleucine","N-Acetyl-L-alanine","Indole-3-acetamide","Triethyl phosphate","δ-Valerolactam","2-Pyrrolidinone","Androsterone","Dihydrouracil","5-amino-1-[3,4-dihydroxy-5-(hydroxymethyl)oxolan-2-yl]imidazole-4-carboxamide","N1-Acetylspermine","Asp-Phe methyl ester","Pterine","N-Acetylphenylalanine","N-Alpha-acetyllysine","Ritalinic acid","Dihydroactinidiolide","2'-Hydroxy-5'-methylacetophenone","Spermidine","2-(3,4-dimethoxyphenyl)ethanamine","N-Acetyl-L-Histidine","Isonicotinic acid","Furfural","L-phenylalanyl-L-proline","3-Hydroxy-DL-kynurenine","Methyl cinnamate","Benzylcinnamate","1,2-Dichloroethane","m-Cresol","Cyclohexylamine","2-Picoline","Butylamine","6-Methyl-5-hepten-2-one","2,2,2-Trichloroethanol","Hydroxyacetone","Methylcysteine","2-Acetylfuran","Benzophenone","TranexamicAcid","Triethylamine","2,5-Dimethyl pyrazine","2,6-Di-tert-butyl-4-methylphenol","(E)-2-Octen-1-ol","Naproxen","2-Pentadecanone","Methyl dihydrojasmonate","Pyrazine","Methylisobutyrate","N-Nitrosodiethylamine","3,5-Dinitrosalicylic acid","Barbituric acid","2-Aminophenol","1,3-diisopropylbenzene","p-Cymene","Terpinolene","gamma-Terpinene","alpha-Terpinene","DL-Leucine","Glycyrrhetinic acid","(+)-borneol","p-Menth-1-en-4-ol","Norambreinolide","2-Undecanone","Undecanal","2-Pentyl-3-phenyl-2-propenal","PYRROLIDINE","Octanal","Pyrene","2-(Methylthio)benzothiazole","5-METHYLFURFURAL","N-ETHYLACETAMIDE","1-PENTADECANOL","2-Methyl-1-propanol","1-Aminopropan-2-ol","Dibutyl phthalate","Pulegone","FURFURYL ALCOHOL","1-Phenylethanol","Nα-Acetyl-L-glutamine","Choline chloride","glycylphenylalanine","L-Glutamic acid","proline betaine","Bisphenol propane","N-methylalanine","tryptophan betaine","O-Succinyl-L-Homoserine","Cytidine 5'-Diphosphocholine","N-Methyl-L-Glutamate","N-Methyl-D-Aspartic Acid","N-Alpha-Acetyl-L-Asparagine","Phosphocholine","2,4-Dihydroxypteridine","L-Methionine Sulfoximine","2'-Deoxycytidine 5'-Diphosphate","Lumichrome","L-Tryptophanamide")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
mSet<-PlotPathSummary(mSet, F, "path_view_0_", "png", 72, width=NA)
#network analysis
mSet<-InitDataObjects("conc", "network", FALSE)
mSet<-SetOrganism(mSet, "hsa")
cmpdListFile<-"tmp.csv"
cmpdList<-readChar(cmpdListFile, file.info(cmpdListFile)$size)
mSet<-PerformCmpdMapping(mSet, cmpdList, "hsa", "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-GetNetworkGeneMappingResultTable(mSet)
mSet<-GetNetworkGeneMappingResultTable(mSet)
mSet<-GetNetworkGeneMappingResultTable(mSet)
mSet<-PrepareNetworkData(mSet);
mSet<-PrepareKeggQueryJson(mSet)
mSet<-PerformKOEnrichAnalysis_KO01100(mSet, "pathway", "network_enrichment_pathway_0");
mSet<-SearchNetDB(mSet, "pheno", "metabo_phenotypes", FALSE, 0.5)
mSet<-CreateGraph(mSet)
####Biomarker analysis####
#Decision tree
tree_data <- t(NormMeta[,1:96])
stat_dat1[,3:433] <- tree_data[,1:431]
data1 <- stat_dat1
data1 <- data1[,-1]
data1$Label <- as.factor(data1$Label)
data11 <- data1[]
library(rpart)
set.seed(20210517) #随机抽样设置种子
train<-sample(nrow(data1),0.7*nrow(data1)) #抽样函数，第一个参数为向量，nrow()返回行数 后面的是抽样参数前
tdata<-data1[train,] #根据抽样参数列选择样本，都好逗号是选择行
vdata<-data1[-train,] #删除抽样行
colnames(tdata) <- c('Label',Metabolome$Compounds)
colnames(vdata) <- c('Label',Metabolome$Compounds)
tdata$Label <- case_when(str_detect(tdata$Label,'0')==T~'Healthy adults',
                         str_detect(tdata$Label,'1')==T~'Mild cases',
                         str_detect(tdata$Label,'2')==T~'Severe cases',
                         str_detect(tdata$Label,'3')==T~'Fatal cases')
vdata$Label <- case_when(str_detect(vdata$Label,'0')==T~'Healthy adults',
                         str_detect(vdata$Label,'1')==T~'Mild cases',
                         str_detect(vdata$Label,'2')==T~'Severe cases',
                         str_detect(vdata$Label,'3')==T~'Fatal cases')
library("pROC")
dt <- rpart(formula = Label~.,data = tdata)
predtree<-predict(dt,newdata=vdata,type="class")
predtree <- as.numeric(predtree)
unique(vdata1$Label)
printcp(dt1)
#
tdata1 <- tdata[,c('Label',row.names(tmp))]
vdata1 <- vdata[,c('Label',row.names(tmp))]
dt1 <- rpart(formula = Label~.,data = tdata1)
predtree<-predict(dt1,newdata=vdata1,type="class")
#predtree1 <- as.numeric(predtree1)
roc(vdata1$Label,predtree1, plot=TRUE, levels=c(2,3),
    print.thres=TRUE, print.auc=TRUE)
roc(vdata1$Label,predtree1, plot=TRUE, levels=c(1,2),
    print.thres=TRUE, print.auc=TRUE)
roc(vdata1$Label,predtree1, plot=TRUE, levels=c(1,2),
    print.thres=TRUE, print.auc=TRUE)
roc(vdata1$Label,predtree1, plot=TRUE, levels=c(0,2),
    print.thres=TRUE, print.auc=TRUE)
printcp(dt1)
#decision tree
library(rpart.plot)
rpart.plot(dt,branch=1,type=4, fallen.leaves=T,
           cex=0.6, sub="The decision tree")

tmp <- NormMeta[abs(NormMeta$M_VS_H)>1&
                  abs(NormMeta$F_VS_H)>1&
                  abs(NormMeta$S_VS_H)>1&
                  NormMeta$M_VS_H_Pvalue<0.05&
                  NormMeta$S_VS_H_Pvalue<0.05&
                  NormMeta$F_VS_H_Pvalue<0.05,]
#SVM
mSet<-InitDataObjects("pktable", "roc", FALSE)
mSet<-Read.TextData(mSet, "Replacing_with_your_file_path", "rowu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ContainMissing(mSet)
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckData(mSet)
mSet<-ContainMissing(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "MeanCenter", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
mSet<-SetAnalysisMode(mSet, "univ")
mSet<-PrepareROCData(mSet)
mSet<-CalculateFeatureRanking(mSet)
mSet<-Perform.UnivROC(mSet, "Putrescine", 0, "png", 72, F, T, "closest.topleft", F, "sp", 0.2)
mSet<-PlotRocUnivBoxPlot(mSet, "Putrescine", 0, "png", 72, T, FALSE)
mSet<-Perform.UnivROC(mSet, "Phosphoric acid", 1, "png", 72, F, T, "closest.topleft", F, "sp", 0.2)
mSet<-PlotRocUnivBoxPlot(mSet, "Phosphoric acid", 1, "png", 72, T, FALSE)
mSet<-Perform.UnivROC(mSet, "4-Ethylphenol", 2, "png", 72, F, T, "closest.topleft", F, "sp", 0.2)
mSet<-PlotRocUnivBoxPlot(mSet, "4-Ethylphenol", 2, "png", 72, T, FALSE)
mSet<-Perform.UnivROC(mSet, "Pyroglutamic acid", 3, "png", 72, F, T, "closest.topleft", F, "sp", 0.2)
mSet<-PlotRocUnivBoxPlot(mSet, "Pyroglutamic acid", 3, "png", 72, T, FALSE)
mSet<-Perform.UnivROC(mSet, "Pyroglutamic acid", 4, "png", 72, T, T, "closest.topleft", F, "sp", 0.2)
mSet<-PlotRocUnivBoxPlot(mSet, "Pyroglutamic acid", 4, "png", 72, T, FALSE)
mSet<-SetAnalysisMode(mSet, "explore")
mSet<-PrepareROCData(mSet)
mSet<-PerformCV.explore(mSet, "svm", "svm", 2)
mSet<-PlotProbView(mSet, "cls_prob_0_", "png", 72, -1, 0, 0)
mSet<-PlotImpVars(mSet, "cls_imp_0_", "png", 72, -1, "freq", 15);
mSet<-PlotAccuracy(mSet, "cls_accu_0_", "png", 72)
mSet<-PlotROC(mSet, "cls_roc_0_", "png", 72, 0, "threshold", 0, 0, "fpr", 0.5)
mSet<-PlotROC(mSet, "cls_roc_1_", "png", 72, 0, "threshold", 1, 0, "fpr", 0.2)
mSet<-PlotProbView(mSet, "cls_prob_1_", "png", 72, 1, 0, 0)
mSet<-SetAnalysisMode(mSet, "test")
mSet<-PrepareROCData(mSet)
mSet<-UpdateKmeans(mSet, 5)
mSet<-SetAnalysisMode(mSet, "test")
mSet<-PrepareROCData(mSet)
selected.cmpds <- c("Putrescine", "Phosphoric acid", "4-Ethylphenol", "L-Aspartic acid", "Neopterin");
selected.smpls <- c()
mSet<-SetCustomData(mSet, selected.cmpds, selected.smpls)
mSet<-PerformCV.test(mSet, "svm", 2)
mSet<-PlotROC(mSet, "cls_test_roc_0_", "png", 72, 0, "threshold", 0, 0, "fpr", 0.5)
mSet<-PlotProbView(mSet, "cls_test_prob_0_", "png", 72, -1, 0, 0)
mSet<-PlotTestAccuracy(mSet, "cls_test_accu_0_", "png", 72)
mSet<-GetAccuracyInfo(mSet)
#neuron network
library(lattice)
library(caret)
library(e1071)
library(nnet)
tdata1$Label<-relevel(as.factor(tdata1$Label),ref = "Healthy adults")
vdata1$Label<-relevel(as.factor(vdata1$Label),ref = "Healthy adults")
mult.model<-multinom(Label~.,data=tdata1,Hess = T)
summary(mult.model)
z <- summary(mult.model)$coefficients/summary(mult.model)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1))*2
p
exp(coef(mult.model))
pre_logistic<-predict(mult.model,newdata = vdata1) 
table(vdata1$Label,pre_logistic)
conMat4<-confusionMatrix(factor(pre_logistic),factor(vdata1$Label))
#OPLS-DA
BiocManager::install('ropls')
library(ropls)
fit31 <- opls(x=tdata1[,2:40],y=tdata1$Label,orthoI=0)
adata <- rbind(tdata1,vdata1)
fit32 <- opls(x=adata[,2:40],y=adata$Label,orthoI=0)
vipvn <- getVipVn(fit32)
vipvn_select <- vipvn[vipvn>1]
head(vipvn_select)











