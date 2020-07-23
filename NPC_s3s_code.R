# 1. Head -----------------------------------------------------------------

output_directory <- "/media/data/Josh/NPC_s3s_paper"
setwd(output_directory)
date_time <- gsub("-", "", Sys.Date())


# 2. Load files -----------------------------------------------------------

# read file to match ENS IDs with HUGO gene names
HUGO_ENS_match <- read.table("HUGO_ENS_matching.txt", header = TRUE, sep = "\t", as.is = TRUE)
table(duplicated(HUGO_ENS_match$Approved.symbol)) # no duplicates here
table(duplicated(HUGO_ENS_match$Ensembl.ID.supplied.by.Ensembl.)) # 4209 duplicates
hugo_name_lookup <- HUGO_ENS_match$Approved.symbol
names(hugo_name_lookup) <- HUGO_ENS_match$Ensembl.ID.supplied.by.Ensembl.

# read study files
count_fil <- read.table("counts_ENS.txt", header = TRUE, sep = "\t", as.is = TRUE)    # read in counts matrix
rownames(count_fil) <- count_fil$X
count_fil$X <- NULL
colnames(count_fil) <- sub("X", "", colnames(count_fil))

tdata_fil <- read.table("technical_data.txt", header = TRUE, sep = "\t", as.is = TRUE) # read in technical annotation file
rownames(tdata_fil) <- tdata_fil$sample_ID
identical(rownames(tdata_fil), colnames(count_fil))

cdata <- read.table("clinical_data.txt", header = TRUE, sep = "\t", as.is = TRUE) # read in clinical annotation file
tcdata <- merge(tdata_fil,cdata,"patient_index", all=T)
tcdata <- tcdata[!is.na(tcdata$sample_ID),]
rownames(tcdata) <- tcdata$sample_ID



# 3. PCA ---------------------------------------------------------------------

# define levels, colors, and symbols
cell_type_levels <- c("PNS", "NAT", "DYS", "TUM", "TME", "CEL", "SQM", "SQL")
cell_type_levels_text <- c("normal", "normal-adjacent to tumor", "dysplastic", "tumor", "microenvironment", "C666-1 cell line", "squamous", "squamous-inflamed")
cell_type_symbols <- c(1,1,15,15,3,11,19,10)
cell_type_colors <- c("black","green3","hotpink","black","black","chocolate4","blue","blue")
batch_levels <- c("A1", "A2", "A3", "A4", "A5", "A6", "H1")
batch_colors <- c("green3", "lightblue3", "hotpink3",  "lavenderblush3", "orange2", "lightgoldenrod4", "royalblue2")


#   3.1 PCA for nasopharyngeal cells ------------------------------------------------------

ns_libraries <- rownames(tdata_fil)[!tdata_fil$cell_type %in% c("SQM", "SQL")]
count_ns <- count_fil[,ns_libraries]
tdata_ns <- tdata_fil[ns_libraries,]

# get tpm values, log2 transform and center
tpm_ns <- 1E6 * sweep(count_ns, 2, colSums(count_ns), "/")
tpm_ns_high <- tpm_ns[rowMeans(tpm_ns)>15,] # filter for high-expressing genes
dim(tpm_ns_high)
tpm_ns_high_log2 <- log2(tpm_ns_high + 1)
tpm_ns_high_log2cen <- tpm_ns_high_log2 - rowMeans(tpm_ns_high_log2)

# singular value decomposition
svd_ns <- svd(tpm_ns_high_log2cen)
# percent variance explained
plot(svd_ns$d^2/sum(svd_ns$d^2),ylab="Percent Variance Explained",xlab="Principal Component", col =" royalblue",type ="o", pch=19)
percent <- round(100*(svd_ns$d^2/sum(svd_ns$d^2)),digits=2)
# define levels for cell type
table(tdata_ns$cell_type)
tdata_ns$cell_type_factor <- factor(tdata_ns$cell_type, levels= cell_type_levels)

# PCA plot
pdf(file=paste("plots/",date_time, "_PCA_nasopharyngeal_text.pdf", sep=""),height=5,width=7)
plot(svd_ns$v[,1],svd_ns$v[,2],ylab=paste("PC2 (", percent[2], "% of variance)", sep=""), xlab=paste("PC1 (", percent[1], "% of variance)", sep=""),
     pch = cell_type_symbols[as.numeric(tdata_ns$cell_type_factor)], 
     col=cell_type_colors[as.numeric(tdata_ns$cell_type_factor)],
     xlim=c(-0.22, 0.27), ylim=c(-0.13, 0.22))
with(tdata_ns, legend("bottomright", pch = cell_type_symbols[1:6], col=cell_type_colors[1:6], legend = cell_type_levels_text[1:6], title="  Cell type  ", cex=0.8))
dev.off()


pdf(file=paste("plots/",date_time, "_PCA_nasopharyngeal.pdf", sep=""),height=6,width=7.5)
# by batch
plot(svd_ns$v[,1],svd_ns$v[,2],ylab=paste("PC2 (", percent[2], "% of variance)", sep=""), xlab=paste("PC1 (", percent[1], "% of variance)", sep=""),
     pch = cell_type_symbols[as.numeric(tdata_ns$cell_type_factor)], 
     col=batch_colors[as.numeric(tdata_ns$batch_factor)],
     xlim=c(-0.22, 0.24), ylim=c(-0.13, 0.22),
     main = "Library Batch")
with(tdata_ns, legend("bottomright", pch = 15, col=batch_colors, legend = batch_levels, title="  Batch  "))

# by age
redbluepal <- colorRampPalette(c('red','yellow','green','blue'))
plot(svd_ns$v[,1],svd_ns$v[,2],ylab=paste("PC2 (", percent[2], "% of variance)", sep=""), xlab=paste("PC1 (", percent[1], "% of variance)", sep=""),
     pch = cell_type_symbols[as.numeric(tdata_ns$cell_type_factor)], 
     col = redbluepal(25)[as.numeric(tdata_ns$block_age)+1],
     xlim=c(-0.22, 0.24), ylim=c(-0.13, 0.22),
     main = "Age of Sample")
legend("bottomright", title="Age (Years)", legend=c(0,6,12,18,24), col=redbluepal(25)[c(1,7,13,19,25)], pch=20)
dev.off()


# 3.2 PCA for panendoscopy samples ------------------------------------------------

pan_libraries <- rownames(tdata_fil)[tdata_fil$cell_type %in% c("PNS", "SQM", "SQL")]
count_pan <- count_fil[,pan_libraries]
tdata_pan <- tdata_fil[pan_libraries,]

tpm_pan <- 1E6 * sweep(count_pan, 2, colSums(count_pan), "/")
tpm_pan_high <- tpm_pan[rowMeans(tpm_pan)>15,]
dim(tpm_pan_high)
tpm_pan_high_log2 <- log2(tpm_pan_high + 1)
tpm_pan_high_log2cen <- tpm_pan_high_log2 - rowMeans(tpm_pan_high_log2)

svd_pan <- svd(tpm_pan_high_log2cen)
names(svd_pan) 
plot(svd_pan$d^2/sum(svd_pan$d^2),ylab="Percent Variance Explained",xlab="Principal Component", col =" royalblue",type ="o", pch=19)
percent <- round(100*(svd_pan$d^2/sum(svd_pan$d^2)),digits=2)

table(tdata_pan$cell_type_epi)
cell_type_epi_levels <- c("BOT", "PIR", "PAL", "TONLP", "TONLR", "PNS")
tdata_pan$cell_type_factor <- factor(tdata_pan$cell_type_epi, levels=cell_type_epi_levels)

pdf(file=paste("plots/",date_time, "_PCA_upperairway.pdf", sep=""),height=5,width=7)
plot(svd_pan$v[,1],svd_pan$v[,2],ylab=paste("PC2 (", percent[2], "% of variance)", sep=""), xlab=paste("PC1 (", percent[1], "% of variance)", sep=""),
     pch = c(11,8,2,9,3,0)[as.numeric(tdata_pan$cell_type_factor)], 
     xlim = c(-0.3, 0.6), ylim = c(-0.7, 0.4))
with(tdata_pan, legend("bottomright", pch = c(11,8,2,9,3,0), col="black", 
                     legend = c("Base of Tongue", "Piriform Sinus", "Palate", 
                                "Tonsil (Lymphoid-Poor)", "Tonsil (Lymphoid-Rich)", "Nasopharynx"), 
                     title = expression(bold("Location")), cex=1))
dev.off()


# 4. Define groups, find DE Genes along the spectrum, perform GSEA ------------------------------------------------------
de_function <- function(de_counts, de_pdata, group_header)
{
  library(DESeq2)
  library(BiocParallel)
  dds <- DESeqDataSetFromMatrix(de_counts, colData = de_pdata, design = ~ temp_group)
  dds2 <- DESeq(dds, fitType = "local", parallel = T, BPPARAM = MulticoreParam(workers = 38)) # eats a lot of memory!)
  res <- results(dds2)
  res <- as.data.frame(res)
  
  res$ENS_id <- rownames(res)
  res$ENS_id <- gsub("\\..*","",res$ENS_id)
  res$symbol <- NA
  res$symbol[res$ENS_id %in% names(hugo_name_lookup)] <- hugo_name_lookup[res$ENS_id[res$ENS_id %in% names(hugo_name_lookup)]]
  res$symbol_or_ENS_id <- res$ENS_id
  res$symbol_or_ENS_id[res$ENS_id %in% names(hugo_name_lookup)] <- hugo_name_lookup[res$ENS_id[res$ENS_id %in% names(hugo_name_lookup)]]
  
  res <- res[order(res$padj),]
  
  library(WriteXLS)
  WriteXLS(res, ExcelFileName=paste("output/", date_time, "_", group_header, "_deseq.xlsx", sep=""), row.names=TRUE)
  return(res)
}

fgsea_function <- function(de_res, de_header)
{
  library(tidyverse)
  ranks_res <- de_res %>% 
    dplyr::select(symbol, stat) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(symbol) %>% 
    summarize(stat=mean(stat))
  ranks_res <- deframe(ranks_res)
  
  library(fgsea)
  fg_res_go <- fgsea(pathways=path_go, stats=ranks_res, nperm=10000)
  fg_res_go$gene_ratio <- lengths(fg_res_go$leadingEdge) / fg_res_go$size
    fg_res_go_tidy <- fg_res_go %>%
    as_tibble() %>%
    arrange(NES)
  fg_res_go_tidy$GO_type <- NA
  fg_res_go_tidy$GO_type[fg_res_go_tidy$pathway %in% names(path_go_bioprocess)] <- "Biological Process"
  fg_res_go_tidy$GO_type[fg_res_go_tidy$pathway %in% names(path_go_cellcomp)] <- "Cellular Component"
  fg_res_go_tidy$GO_type[fg_res_go_tidy$pathway %in% names(path_go_molfx)] <- "Molecular Function"
  fg_res_go_tidy <- fg_res_go_tidy[,c(ncol(fg_res_go_tidy),1:(ncol(fg_res_go_tidy)-1))]
  
  # # plot top 20 pathways
  topPathwaysUp <- fg_res_go[head(order(-NES), n=15), pathway]
  topPathwaysDown <- fg_res_go[head(order(NES), n=15), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  pdf(file=paste("output/",date_time, "_", de_header, "_fgsea_top.pdf", sep=""),height=8,width=12)
  plotGseaTable(path_go[topPathways], ranks_res, fg_res_go, gseaParam = 0.5)
  dev.off()
  
  library(data.table)
  fwrite(fg_res_go_tidy, file=paste("output/", date_time, "_", de_header, "_fgsea_GO.txt", sep=""), 
         sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
  fg_res_go_tidy <- fg_res_go_tidy[fg_res_go_tidy$padj < 0.05,]
  fwrite(fg_res_go_tidy, file=paste("output/", date_time, "_", de_header, "_fgsea_GO_sig.txt", sep=""), 
         sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
  
    
  return(fg_res_go)
}

# define groups
samples_PNS <- rownames(tdata_ns)[tdata_ns$cell_type == "PNS"]
samples_NAT <- rownames(tdata_ns)[tdata_ns$cell_type == "NAT"]
samples_DYS <- rownames(tdata_ns)[tdata_ns$cell_type == "DYS"]
samples_TUM <- rownames(tdata_ns)[tdata_ns$cell_type == "TUM"]
samples_CEL <- rownames(tdata_ns)[tdata_ns$cell_type == "CEL"]
samples_TME <- rownames(tdata_fil)[tdata_fil$cell_type == "TME"]

# pathways
library(fgsea)
path_hallmark <- gmtPathways("/media/data/Josh/NPC_s3s_all/gmt/h.all.v6.2.symbols.gmt")
path_hallmark_descriptions <- read.table("/media/data/Josh/NPC_s3s_all/gmt/hallmark_descriptions.csv", header = TRUE, sep = "\t", as.is = TRUE)
path_hallmark_immune <- path_hallmark_descriptions$Hallmark.Name[path_hallmark_descriptions$Process.Category == "immune"]
path_hallmark_proliferation <- path_hallmark_descriptions$Hallmark.Name[path_hallmark_descriptions$Process.Category == "proliferation"]
path_hallmark_additional <- c("TNFA_SIGNALING_VIA_NFKB")

path_allmsdb <- gmtPathways("/media/data/Josh/NPC_s3s_all/gmt/msigdb.v6.2.symbols.gmt")
path_go <- gmtPathways("/media/data/Josh/NPC_s3s_all/gmt/c5.all.v6.2.symbols.gmt")
path_go_bioprocess <- gmtPathways("/media/data/Josh/NPC_s3s_all/gmt/c5.bp.v6.2.symbols.gmt")
path_go_cellcomp <- gmtPathways("/media/data/Josh/NPC_s3s_all/gmt/c5.cc.v6.2.symbols.gmt")
path_go_molfx <- gmtPathways("/media/data/Josh/NPC_s3s_all/gmt/c5.mf.v6.2.symbols.gmt")


#   4.1 Normal adjacent (NAT) vs Panendoscopy normal (PNS) ---------------------------------------------
grp_header <- "normal-adj_vs_normal"
grp_A <- samples_NAT
grp_B <- samples_PNS
count_temp <- count_ns[,c(grp_A, grp_B)]
anno_temp <- tdata_ns[c(grp_A,grp_B),]
identical(colnames(count_temp), rownames(anno_temp))
anno_temp$temp_group <- NA
anno_temp$temp_group[rownames(anno_temp) %in% grp_A] <- "A"
anno_temp$temp_group[rownames(anno_temp) %in% grp_B] <- "B"
anno_temp$temp_group <- factor(anno_temp$temp_group, levels=c("A","B"))
anno_temp$temp_group <- relevel(anno_temp$temp_group, ref="B") # B is the reference

res_NAT_vs_PNS <- de_function(count_temp, anno_temp, grp_header) 
library(tidyverse)
ranks_NAT_vs_PNS <- res_NAT_vs_PNS %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))
ranks_NAT_vs_PNS <- deframe(ranks_NAT_vs_PNS)

# FGSEA Hallmark
fg_NAT_vs_PNS_hallmark <- fgsea(pathways=path_hallmark, stats=ranks_NAT_vs_PNS, nperm=1000)
fg_NAT_vs_PNS_hallmark_tidy <- fg_NAT_vs_PNS_hallmark %>%
  as_tibble() %>%
  arrange(desc(ES))
fwrite(fg_NAT_vs_PNS_hallmark_tidy, file=paste("output/", date_time, "_", grp_header, "_fgsea_HALLMARK.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
fwrite(fg_NAT_vs_PNS_hallmark_tidy[fg_NAT_vs_PNS_hallmark_tidy$padj < 0.05,], file=paste("output/", date_time, "_", grp_header, "_fgsea_HALLMARK_sig.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)

# FGSEA GO
fg_NAT_vs_PNS_go <- fgsea_function(res_NAT_vs_PNS, grp_header)

# volcano genes
vp <- res_NAT_vs_PNS
fp <- fg_NAT_vs_PNS_go
vp_response_to_interferon_type1 <- unlist(fp[fp$pathway=="GO_RESPONSE_TO_TYPE_I_INTERFERON","leadingEdge"][[1]])
vp_virus_response <- unlist(fp[fp$pathway=="GO_DEFENSE_RESPONSE_TO_VIRUS","leadingEdge"][[1]])
vp_innate <- unlist(fp[fp$pathway=="GO_INNATE_IMMUNE_RESPONSE","leadingEdge"][[1]])
vp_spry <- c("SPRY1", "SPRY2")
vp_text <- c("SPRY1", "SPRY2", "CXCL9", "IFITM2", "ISG15", "IFI6", "IFIT3", "CCL11", "CXCL10", "LCN2", "STAT1", "APOL1", "C1QB") # "IFITM1"


pdf(file=paste("plots/",date_time, "_vol_NATvsPNS.pdf", sep=""),height=5,width=6.25)
with(vp, plot(log2FoldChange, -log10(padj), pch=20, main=" ", col=alpha("black", 0.2), xlim=c(-8,8)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_innate,], points(log2FoldChange, -log10(padj), pch=20, cex=1.2, col=alpha("yellow", 0.6)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_virus_response,], points(log2FoldChange, -log10(padj), pch=20, cex=1.2, col=alpha("orange", 0.8)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_response_to_interferon_type1,], points(log2FoldChange, -log10(padj), pch=20, cex=1.2, col=alpha("red", 0.6)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_spry,], points(log2FoldChange, -log10(padj), pch=20, cex=1.8, col=alpha("purple", 0.6)))

vp_select <- vp[vp$symbol %in% vp_text,]
library(basicPlotteR)
addTextLabels(vp_select$log2FoldChange, -log10(vp_select$padj), vp_select$symbol, 
              cex.label = 0.8, 
              col.label="black", col.background = "white")

legend("topleft", pch = c(16,16,16), col=c("yellow","orange","red"), 
         legend = c("innate immune response", "defense response to virus", "response to type I interferon"), 
         title = expression(bold("Pathway")), cex=0.8)
dev.off()

#   4.2 Dysplasia (DYS) vs Normal-adjacent (NAT) ---------------------------------------------
grp_header <- "dysplasia_vs_normal-adj"
grp_A <- samples_DYS
grp_B <- samples_NAT
count_temp <- count_ns[,c(grp_A, grp_B)]
anno_temp <- tdata_ns[c(grp_A,grp_B),]
identical(colnames(count_temp), rownames(anno_temp))
anno_temp$temp_group <- NA
anno_temp$temp_group[rownames(anno_temp) %in% grp_A] <- "A"
anno_temp$temp_group[rownames(anno_temp) %in% grp_B] <- "B"
anno_temp$temp_group <- factor(anno_temp$temp_group, levels=c("A","B"))
anno_temp$temp_group <- relevel(anno_temp$temp_group, ref="B") # B is the reference

res_DYS_vs_NAT <- de_function(count_temp, anno_temp, grp_header) 
library(tidyverse)
ranks_DYS_vs_NAT <- res_DYS_vs_NAT %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))
ranks_DYS_vs_NAT <- deframe(ranks_DYS_vs_NAT)

# FGSEA Hallmark
library(fgsea)
fg_DYS_vs_NAT_hallmark <- fgsea(pathways=path_hallmark, stats=ranks_DYS_vs_NAT, nperm=1000)
fg_DYS_vs_NAT_hallmark_tidy <- fg_DYS_vs_NAT_hallmark %>%
  as_tibble() %>%
  arrange(desc(ES))
fwrite(fg_DYS_vs_NAT_hallmark_tidy, file=paste("output/", date_time, "_", grp_header, "_fgsea_HALLMARK.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
fwrite(fg_DYS_vs_NAT_hallmark_tidy[fg_DYS_vs_NAT_hallmark_tidy$padj < 0.05,], file=paste("output/", date_time, "_", grp_header, "_fgsea_HALLMARK_sig.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)

# FGSEA GO
library(fgsea)
fg_DYS_vs_NAT_go <- fgsea(pathways=path_go, stats=ranks_DYS_vs_NAT, nperm=1000)
fg_DYS_vs_NAT_go_tidy <- fg_DYS_vs_NAT_go %>%
  as_tibble() %>%
  arrange(desc(NES))

#   4.3 Tumor (TUM) vs Dysplasia (DYS) ---------------------------------------------
grp_header <- "tumor_vs_dysplasia"
grp_A <- samples_TUM
grp_B <- samples_DYS
count_temp <- count_ns[,c(grp_A, grp_B)]
anno_temp <- tdata_ns[c(grp_A,grp_B),]
identical(colnames(count_temp), rownames(anno_temp))
anno_temp$temp_group <- NA
anno_temp$temp_group[rownames(anno_temp) %in% grp_A] <- "A"
anno_temp$temp_group[rownames(anno_temp) %in% grp_B] <- "B"
anno_temp$temp_group <- factor(anno_temp$temp_group, levels=c("A","B"))
anno_temp$temp_group <- relevel(anno_temp$temp_group, ref="B") # B is the reference

res_TUM_vs_DYS <- de_function(count_temp, anno_temp, grp_header) 
library(tidyverse)
ranks_TUM_vs_DYS <- res_TUM_vs_DYS %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))
ranks_TUM_vs_DYS <- deframe(ranks_TUM_vs_DYS)

# FGSEA Hallmark
fg_TUM_vs_DYS_hallmark <- fgsea(pathways=path_hallmark, stats=ranks_TUM_vs_DYS, nperm=1000)
fg_TUM_vs_DYS_hallmark_tidy <- fg_TUM_vs_DYS_hallmark %>%
  as_tibble() %>%
  arrange(desc(ES))
fwrite(fg_TUM_vs_DYS_hallmark_tidy, file=paste("output/", date_time, "_", grp_header, "_fgsea_HALLMARK.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
fwrite(fg_TUM_vs_DYS_hallmark_tidy[fg_TUM_vs_DYS_hallmark_tidy$padj < 0.05,], file=paste("output/", date_time, "_", grp_header, "_fgsea_HALLMARK_sig.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)


# FGSEA GO
fg_TUM_vs_DYS_go <- fgsea(pathways=path_go, stats=ranks_TUM_vs_DYS, nperm=1000)
fg_TUM_vs_DYS_go_tidy <- fg_TUM_vs_DYS_go %>%
  as_tibble() %>%
  arrange(desc(NES))


#   4.4 Tumor (TUM) vs Normal (PNS) ---------------------------------------------
grp_header <- "tumor_vs_normal"
grp_A <- samples_TUM
grp_B <- samples_PNS
count_temp <- count_ns[,c(grp_A, grp_B)]
anno_temp <- tdata_ns[c(grp_A,grp_B),]
identical(colnames(count_temp), rownames(anno_temp))
anno_temp$temp_group <- NA
anno_temp$temp_group[rownames(anno_temp) %in% grp_A] <- "A"
anno_temp$temp_group[rownames(anno_temp) %in% grp_B] <- "B"
anno_temp$temp_group <- factor(anno_temp$temp_group, levels=c("A","B"))
anno_temp$temp_group <- relevel(anno_temp$temp_group, ref="B") # B is the reference

res_TUM_vs_PNS <- de_function(count_temp, anno_temp, grp_header) 
library(tidyverse)
ranks_TUM_vs_PNS <- res_TUM_vs_PNS %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))
ranks_TUM_vs_PNS <- deframe(ranks_TUM_vs_PNS)

# FGSEA Hallmark
library(fgsea)
fg_TUM_vs_PNS_hallmark <- fgsea(pathways=path_hallmark, stats=ranks_TUM_vs_PNS, nperm=1000)
fg_TUM_vs_PNS_hallmark$gene_ratio <- lengths(fg_TUM_vs_PNS_hallmark$leadingEdge) / fg_TUM_vs_PNS_hallmark$size

fg_TUM_vs_PNS_hallmark_tidy <- fg_TUM_vs_PNS_hallmark %>%
  as_tibble() %>%
  arrange(desc(ES))
fwrite(fg_TUM_vs_PNS_hallmark_tidy, file=paste("output/", date_time, "_", grp_header, "_fgsea_HALLMARK.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
fwrite(fg_TUM_vs_PNS_hallmark_tidy[fg_TUM_vs_PNS_hallmark_tidy$padj < 0.05,], file=paste("output/", date_time, "_", grp_header, "_fgsea_HALLMARK_sig.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)


fg_TUM_vs_PNS_hallmark_select <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_G2M_CHECKPOINT",
  # "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_MYC_TARGETS_V2",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_KRAS_SIGNALING_UP",
  # "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  # "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  # "HALLMARK_COMPLEMENT",
  # "HALLMARK_APOPTOSIS",
  # "HALLMARK_UV_RESPONSE_UP",
  "HALLMARK_HYPOXIA",
  "HALLMARK_XENOBIOTIC_METABOLISM",
  "HALLMARK_ANDROGEN_RESPONSE",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_FATTY_ACID_METABOLISM",
  # "HALLMARK_HEME_METABOLISM",
  # "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_PEROXISOME")

library(stringr)
fg_TUM_vs_PNS_hallmark_dp <- fg_TUM_vs_PNS_hallmark[fg_TUM_vs_PNS_hallmark$pathway %in% fg_TUM_vs_PNS_hallmark_select,]
fg_TUM_vs_PNS_hallmark_dp$pathway <- str_sub(fg_TUM_vs_PNS_hallmark_dp$pathway, start =10)
fg_TUM_vs_PNS_hallmark_dp$pathway <- str_replace_all(fg_TUM_vs_PNS_hallmark_dp$pathway, "_", " ")

pdf(file=paste("plots/",date_time, "_hallmark_dp_TUMvsPNS.pdf", sep=""),height=5,width=6)
library(ggplot2)
ggplot(fg_TUM_vs_PNS_hallmark_dp) +
  geom_point(aes(NES, reorder(pathway, NES), size = gene_ratio, color=padj)) +
  scale_color_gradient(low="red", high = "blue") +
  ylab("Hallmark Process") + 
  xlab("Normalized Enrichment Score") +
  labs(color="p-adj", size = "gene ratio") 
dev.off()

# FGSEA GO
library(fgsea)
fg_TUM_vs_PNS_go <- fgsea_function(res_TUM_vs_PNS, grp_header)

# volcano genes
vp <- res_TUM_vs_PNS
fp <- fg_TUM_vs_PNS_go

vp_dna_replication <- unlist(fp[fp$pathway=="GO_DNA_REPLICATION","leadingEdge"][[1]])
vp_leu_migration <- unlist(fp[fp$pathway=="GO_LEUKOCYTE_MIGRATION","leadingEdge"][[1]])
vp_virus_response <- unlist(fp[fp$pathway=="GO_RESPONSE_TO_VIRUS","leadingEdge"][[1]])
vp_cellcycle <- unlist(fp[fp$pathway=="GO_CELL_CYCLE_CHECKPOINT","leadingEdge"][[1]])
vp_cilia_org <- unlist(fp[fp$pathway=="GO_CILIUM_ORGANIZATION","leadingEdge"][[1]])

vp_text <- c("TOP2A","CHEK1","MCM7", "CDK1", "MKI67", "TP53", "E2F1", 
             "VCAM1", "IFI44L", "CXCL9", "CXCL10","CXCL11", "IFIT3", "CCL11")
library (calibrate)
pdf(file=paste("plots/",date_time, "_vol_TUMvsPNS.pdf", sep=""),height=5,width=6.25)
with(vp, plot(log2FoldChange, -log10(padj), pch=20, col=alpha("grey", 0.2), xlim=c(-8,9.5), ylim=c(0,28)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_cilia_org,], points(log2FoldChange, -log10(padj), pch=20, cex=1.3, col=alpha("chocolate4", 0.5)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_leu_migration,], points(log2FoldChange, -log10(padj), pch=20, cex=1.3, col=alpha("magenta", 0.5)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_dna_replication,], points(log2FoldChange, -log10(padj), pch=20, cex=1.3, col=alpha("green", 0.5)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_virus_response,], points(log2FoldChange, -log10(padj), pch=20, cex=1.3, col=alpha("orange", 0.5)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_cellcycle,], points(log2FoldChange, -log10(padj), pch=20, cex=1.3, col=alpha("blue", 0.5)))

vp_select <- vp[vp$symbol %in% vp_text,]
library(basicPlotteR)
addTextLabels(vp_select$log2FoldChange, -log10(vp_select$padj), vp_select$symbol, 
              cex.label = 0.8, 
              col.label="black", col.background = "white")

# with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_text,], textxy(log2FoldChange, -log10(padj), labs=symbol, cex=.6))
legend("topleft", pch = c(16,16,16), col=c("green","blue","orange", "magenta", "chocolate4"), bg="white",
       legend = c("dna replication", "cell cycle checkpoint", "response to virus", "leukocyte migration", "cilium organization"), 
       title = expression(bold("Pathway")), cex=0.9)
dev.off()


vp_nfkb_noncanonical <- c("CD40", "RELB", "NFKB2", "TNFRSF1B")
vp_nfkb_canonical <- c("NFKBIA", "NFKB1", "TRAF2", "RELA", "TNFRSF1A")
vp_text <- c(vp_nfkb_canonical, vp_nfkb_noncanonical)

pdf(file=paste("plots/",date_time, "_vol_nfkb.pdf", sep=""),height=5,width=6)
with(vp, plot(log2FoldChange, -log10(padj), pch=20, main=" ", col=alpha("grey", 0.2), xlim=c(-2,3.5), ylim=c(0,6)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_nfkb_canonical,], points(log2FoldChange, -log10(padj), pch=20, cex=1.4, col=alpha("green", 1)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_nfkb_noncanonical,], points(log2FoldChange, -log10(padj), pch=20, cex=1.4, col=alpha("blue", 1)))
abline(h=1.301, col="red", lty = 2, lwd = 2)
text(-1.6, 1.6, "p adj = 0.05", cex = .9)

vp_select <- vp[vp$symbol %in% vp_text,]
library(basicPlotteR)
addTextLabels(vp_select$log2FoldChange, -log10(vp_select$padj), vp_select$symbol, 
              cex.label = 0.8, 
              col.label="black", col.background = "white")
legend("topleft", pch = c(16,16), col=c("green","blue"), bg="white",
       legend = c("Canonical", "Non-canonical"), 
       title = expression(bold("     NF-kB Pathway Genes     ")), cex=1)
dev.off()


fg_TUM_vs_PNS_select <- c(
  "GO_DNA_REPLICATION",
  "GO_LEUKOCYTE_MIGRATION",
  "GO_LYMPHOCYTE_ACTIVATION",
  "GO_RESPONSE_TO_VIRUS",
  "GO_CELL_CYCLE_G1_S_PHASE_TRANSITION",
  "GO_RESPONSE_TO_INTERFERON_GAMMA",
  "GO_CELL_CYCLE_CHECKPOINT",
  "GO_IMMUNE_RESPONSE",
  "GO_MITOTIC_NUCLEAR_DIVISION",
  "GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
  "GO_INNATE_IMMUNE_RESPONSE",
  "GO_RESPONSE_TO_IONIZING RADIATION",
  "GO_REGULATION_OF_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION",
  "GO_MITOTIC_CELL_CYCLE",
  "GO_REGULATION_OF_LEUKOCYTE_CHEMOTAXIS",
  "GO_RESPONSE_TO_TYPE_I_INTERFERON",
  "GO_NIK_NF_KAPPAB_SIGNALING",
  "GO_CILIUM_ORGANIZATION",
  "GO_AXONEMAL_DYNEIN_COMPLEX_ASSEMBLY",
  "GO_CILIUM_MOVEMENT",
  "GO_SMOOTHENED_SIGNALING_PATHWAY",
  "GO_CELLULAR_GLUCURONIDATION",
  "GO_GLUTATHIONE")


library(stringr)
fg_TUM_vs_PNS_go_dp <- fg_TUM_vs_PNS_go[fg_TUM_vs_PNS_go$pathway %in% fg_TUM_vs_PNS_select,]
fg_TUM_vs_PNS_go_dp$pathway <- str_sub(fg_TUM_vs_PNS_go_dp$pathway, start =4)
fg_TUM_vs_PNS_go_dp$pathway <- str_replace_all(fg_TUM_vs_PNS_go_dp$pathway, "_", " ")

pdf(file=paste("plots/",date_time, "_gsea_dp_TUMvsPNS.pdf", sep=""),height=5,width=8.5)
library(ggplot2)
ggplot(fg_TUM_vs_PNS_go_dp) +
  geom_point(aes(NES, reorder(pathway, NES), size = gene_ratio, color=pval)) +
  scale_color_gradient(low="red", high = "blue") +
  ylab("GO Biological Process") + 
  xlab("Normalized Enrichment Score") +
  labs(color="p-adj", size = "gene ratio")
dev.off()

#   4.5 Hallmark plots -----
plot_hallmark_fx <- function(fgp, fgp_head)
{
  fgp$pathway <- gsub("HALLMARK_", "", fgp$pathway)
  fgp$color_padj <- NA
  fgp$color_padj[fgp$padj < 0.05] <- "yes"
  fgp$color_padj[fgp$padj >= 0.05] <- "no"
  fgp$color_padj <- factor(fgp$color_padj, levels=c("yes","no"))
  fgp$color_highlight <- NA
  fgp$color_highlight[fgp$padj < 0.05] <- 1
  fgp$color_highlight[fgp$padj >= 0.05] <- 2
  fgp$color_highlight[fgp$pathway %in% path_hallmark_immune] <- 3
  fgp$color_highlight[fgp$pathway %in% path_hallmark_proliferation] <- 4
  fgp$color_highlight[fgp$pathway %in% path_hallmark_additional] <- 5
  color_bars <- c("gray80", "gray90", "orangered3", "blue")
  color_axis <- c("black", "black", "orangered3", "blue", "purple")
  
  library(ggplot2)
  p <- ggplot(fgp, aes(reorder(pathway,ES), ES)) 
  fgplot <- p + 
    geom_col(aes(fill=color_padj), 
             col=color_bars[rev(fgp$color_highlight)]
    ) +
    scale_fill_manual("padj < 0.05", labels=c("yes","no"), values=color_bars) +
    # labs(x="Pathway", y="Enrichment Score", title=paste("Hallmark ", fgp_head, sep="")) + 
    labs(x="Pathway", y="Enrichment Score") + 
        coord_flip() +
    theme_minimal() +
    theme(axis.text.y = element_text(color=color_axis[rev(fgp$color_highlight)]))
  return(fgplot)
}

pdf(file=paste("plots/",date_time, "_hallmark_NATvsPNS.pdf", sep=""),height=9,width=6)
plot_hallmark_fx(fg_NAT_vs_PNS_hallmark_tidy, "NAT vs PNS")
dev.off()

pdf(file=paste("plots/",date_time, "_hallmark_DYSvsNAT.pdf", sep=""),height=9,width=6)
plot_hallmark_fx(fg_DYS_vs_NAT_hallmark_tidy, "DYS vs NAT")
dev.off()

pdf(file=paste("plots/",date_time, "_hallmark_TUMvsDYS.pdf", sep=""),height=9,width=6)
plot_hallmark_fx(fg_TUM_vs_DYS_hallmark_tidy, "TUM vs DYS")
dev.off()

pdf(file=paste("plots/",date_time, "_hallmark_TUMvsPNS.pdf", sep=""),height=9,width=6)
plot_hallmark_fx(fg_TUM_vs_PNS_hallmark_tidy, "TUM vs PNS")
dev.off()



#   4.6 Heatmap combining Hallmark genes ------------------------------------
# get a combined df
fg_NAT_vs_PNS_hallmark_tidy$group <- "NAT_vs_PNS"
fg_DYS_vs_NAT_hallmark_tidy$group <- "DYS_vs_NAT"
fg_TUM_vs_DYS_hallmark_tidy$group <- "TUM_vs_DYS"
fg_hallmark_df <- list(fg_NAT_vs_PNS_hallmark_tidy, fg_DYS_vs_NAT_hallmark_tidy, fg_TUM_vs_DYS_hallmark_tidy)
fg_hallmark_combined <- do.call(rbind, fg_hallmark_df)

# filter
fg_hallmark_combined <- fg_hallmark_combined[fg_hallmark_combined$padj < 0.05,]
dim(fg_hallmark_combined)
names(fg_hallmark_combined$leadingEdge) <- fg_hallmark_combined$pathway
hallmark_genelist <- fg_hallmark_combined$leadingEdge
# combine based on names
hallmark_genelist <- sapply(unique(names(hallmark_genelist)), function(x) unname(unlist(hallmark_genelist[names(hallmark_genelist)==x])), simplify=FALSE)
for (i in c(1:length(hallmark_genelist)))
  {hallmark_genelist[[i]] <- unique(hallmark_genelist[[i]])}
length(hallmark_genelist)

hallmark_genes_all <- hallmark_genelist[[1]]
for (i in c(2:length(hallmark_genelist)))
  {hallmark_genes_all <- c(hallmark_genes_all, hallmark_genelist[[i]])}
hallmark_genes_all <- unique(hallmark_genes_all)
length(hallmark_genes_all)

cases_path <- c(samples_PNS, samples_NAT, samples_DYS, samples_TUM)

count_hallmark <- count_ns
ens_temp <- rownames(count_hallmark)
length(unique(ens_temp))
ens_temp[1:57820] <- gsub("\\..*","",ens_temp[1:57820])
length(unique(ens_temp))

ens_temp[ens_temp %in% names(hugo_name_lookup)] <- hugo_name_lookup[ens_temp[ens_temp %in% names(hugo_name_lookup)]]
length(ens_temp)
table(duplicated(ens_temp))
rownames(count_hallmark) <- ens_temp

tpm_hallmark <- 1E6 * sweep(count_hallmark, 2, colSums(count_hallmark), "/")
tpm_hallmark <- tpm_hallmark[hallmark_genes_all,cases_path]
tpm_hallmark_lg2cen <- log2(tpm_hallmark +1)
tpm_hallmark_lg2cen <- t(apply(tpm_hallmark_lg2cen, 1, function(x) (x - mean(x)) / sd(x)))
# tpm_hallmark_lg2cen <- tpm_hallmark_lg2cen - rowMeans(tpm_hallmark_lg2cen)
range(tpm_hallmark_lg2cen)
sd(tpm_hallmark_lg2cen)

anno_hallmark <- tdata_ns[cases_path,c("cell_type","batch")]
anno_hallmark$batch <- NULL
library(dplyr)
anno_hallmark$cell_type <- recode(anno_hallmark$cell_type, 'PNS' = "normal", 'NAT' = "normal-adjacent",
                                  'DYS' = "dysplasia", 'TUM' = "tumor")

identical(rownames(anno_hallmark), colnames(tpm_hallmark_lg2cen))

colorlist <- list(cell_type=c("normal"="green3", "normal-adjacent"="khaki2", "dysplasia"="orange", "tumor"="indianred1"),
                  batch=c("A1"="darkolivegreen1","A2"="darkolivegreen2","A3"="darkolivegreen3","A4"="darkolivegreen4",
                          "A5"="darkgreen","A6"="darkolivegreen", "B1"="lightskyblue1", "B2"="lightskyblue",
                          "B3"="lightskyblue3", "H1"="mistyrose3", "M1"="chocolate4", "P1"="orange")
                  )
breaksList = seq(-4, 4, by = 0.1)

library(pheatmap)
library(RColorBrewer)
pdf(file=paste("plots/",date_time, "_hallmark_heatmap.pdf", sep=""),height=6,width=6.5)
hm_hallmark_genes <- pheatmap(
  tpm_hallmark_lg2cen,
  annotation_col = anno_hallmark,
  annotation_colors = colorlist,
  show_colnames = F,
  show_rownames = F,
  fontsize_col = 8,
  fontsize_row = 8,
  cluster_rows=T,
  cluster_cols=F,
  col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(breaksList)),
  annotation_names_col = FALSE,
  breaks = breaksList,
  border_color="black"
)
hm_hallmark_genes
dev.off()



# 5. Panendoscopy -----------------------------------------------------

samples_PNS <- rownames(tdata_pan)[tdata_pan$cell_type == "PNS"]
samples_SQM <- rownames(tdata_pan)[tdata_pan$cell_type == "SQM"]
samples_SQMSQL <- rownames(tdata_pan)[tdata_pan$cell_type %in% c("SQM", "SQL")]

grp_header <- "nasopharynx_vs_squamous"
grp_A <- samples_PNS
grp_B <- samples_SQM
count_temp <- count_pan[,c(grp_A, grp_B)]
anno_temp <- tdata_pan[c(grp_A,grp_B),]
identical(colnames(count_temp), rownames(anno_temp))
anno_temp$temp_group <- NA
anno_temp$temp_group[rownames(anno_temp) %in% grp_A] <- "A"
anno_temp$temp_group[rownames(anno_temp) %in% grp_B] <- "B"
anno_temp$temp_group <- factor(anno_temp$temp_group, levels=c("A","B"))
anno_temp$temp_group <- relevel(anno_temp$temp_group, ref="B") # B is the reference

res_PNS_vs_SQM <- de_function(count_temp, anno_temp, grp_header) 
fg_PNS_vs_SQM_go <- fgsea_function(res_PNS_vs_SQM, grp_header)


fg_PNS_vs_SQM_select <- c(
"GO_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM",
"GO_KERATINOCYTE_DIFFERENTIATION",
"GO_TRANSLATIONAL_INITIATION",
"GO_KERATINIZATION",
"GO_CHROMOSOME_SEGREGATION",
"GO_CELLULAR_RESPIRATION",
"GO_MITOTIC_NUCLEAR_DIVISION",
"GO_CELL_CYCLE_G1_S_PHASE_TRANSITION",
"GO_DNA_REPLICATION",
"GO_CILIUM_MORPHOGENESIS",
"GO_CILIUM_ORGANIZATION",
"GO_AXONEMAL_DYNEIN_COMPLEX_ASSEMBLY",
"GO_SMOOTHENED_SIGNALING_PATHWAY",
"GO_NIK_NF_KAPPAB_SIGNALING",
"GO_LYMPHOCYTE_CHEMOTAXIS",
"GO_CALCIUM_ION_TRANSMEMBRANE_TRANSPORT")

fg_PNS_vs_SQM_go_dp <- fg_PNS_vs_SQM_go[fg_PNS_vs_SQM_go$pathway %in% fg_PNS_vs_SQM_select,]
fg_PNS_vs_SQM_go_dp$pathway <- str_sub(fg_PNS_vs_SQM_go_dp$pathway, start =4)
fg_PNS_vs_SQM_go_dp$pathway <- str_replace_all(fg_PNS_vs_SQM_go_dp$pathway, "_", " ")

pdf(file=paste("plots/",date_time, "_gsea_dp_PNSvsSQM.pdf", sep=""),height=5,width=7.8)
library(ggplot2)
ggplot(fg_PNS_vs_SQM_go_dp) +
  geom_point(aes(NES, reorder(pathway, NES), size = gene_ratio, color=padj)) +
  scale_color_gradient(low="red", high = "blue") +
  ylab("GO Biological Process") + 
  xlab("Normalized Enrichment Score") +
  labs(color="p-adj", size = "gene ratio") +
  ggtitle(" ") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
dev.off()

vp <- res_PNS_vs_SQM
fp <- fg_PNS_vs_SQM_go
vp_keratinocytediff <- unlist(fp[fp$pathway=="GO_KERATINOCYTE_DIFFERENTIATION","leadingEdge"][[1]])
vp_keratinfilament <- unlist(fp[fp$pathway=="GO_KERATIN_FILAMENT","leadingEdge"][[1]])
vp_epidermisdiff <- unlist(fp[fp$pathway=="GO_EPIDERMIS_DEVELOPMENT","leadingEdge"][[1]])
vp_ciliumorg <- unlist(fp[fp$pathway=="GO_CILIUM_ORGANIZATION","leadingEdge"][[1]])
vp_smoothened <- unlist(fp[fp$pathway=="GO_SMOOTHENED_SIGNALING_PATHWAY","leadingEdge"][[1]])
vp_text <- c("S100A8","S100A14", "S100A9", "KRT4", "KRT13", "SPRR1A", "SPRR2",
             "FOXJ1", "DNAAF1", "HYDIN", "TEKT2", "RSPH4A", "PTCH1", "STK36")

pdf(file=paste("plots/",date_time, "_vol_PNSvsSQM.pdf", sep=""),height=4.8,width=5.4)

with(vp, plot(log2FoldChange, -log10(padj), pch=20, 
              main=" ", 
              col=alpha("grey", 0.6),
              xlim = c(-8, 9), ylim=c(0,120)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_epidermisdiff,], points(log2FoldChange, -log10(padj), pch=20, cex=1.5, col=alpha("blue", 0.4)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_keratinfilament,], points(log2FoldChange, -log10(padj), pch=20, cex=1.5, col=alpha("purple", 0.6)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_ciliumorg,], points(log2FoldChange, -log10(padj), pch=20, cex=1.5, col=alpha("green3", 0.4)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_smoothened,], points(log2FoldChange, -log10(padj), pch=20, cex=1.5, col=alpha("red", 0.6)))
legend("topright", pch = c(16,16,16,16), col=c("blue","purple","green","red"), 
       legend = c("epidermis differentiation", "keratin filament", "cilia organization", "smoothened pathway"), 
       title = expression(bold("Pathway")), cex=0.8)
vp_select <- vp[vp$symbol %in% vp_text,]
library(basicPlotteR)
addTextLabels(vp_select$log2FoldChange, -log10(vp_select$padj), vp_select$symbol, 
              cex.label = 0.7, 
              col.label="black", col.background = "white")
dev.off()



# 5.1 Panendoscopy plots --------------------------------------------------

cases_panendo <- c(samples_PNS, samples_SQM)
count_panendo <- count_fil[,cases_panendo]
ens_temp <- rownames(count_panendo)
length(unique(ens_temp))
ens_temp[1:57820] <- gsub("\\..*","",ens_temp[1:57820])
length(unique(ens_temp))

ens_temp[ens_temp %in% names(hugo_name_lookup)] <- hugo_name_lookup[ens_temp[ens_temp %in% names(hugo_name_lookup)]]
length(ens_temp)
table(duplicated(ens_temp))
rownames(count_panendo) <- ens_temp

tpm_panendo <- 1E6 * sweep(count_panendo, 2, colSums(count_panendo), "/")
tpm_panendo <- tpm_panendo[,cases_panendo]
tpm_panendo_lg2cen <- log2(tpm_panendo +1)
tpm_panendo_lg2cen <- t(apply(tpm_panendo_lg2cen, 1, function(x) (x - mean(x)) / sd(x)))
anno_panendo <- tdata_fil[cases_panendo,]
tdata_ns$cell_type_factor <- factor(tdata_ns$cell_type, levels= cell_type_levels)
identical(rownames(anno_panendo),colnames(tpm_panendo))
# TRUE


g_panendo <- c("CX3CL1", "CCL20", "SAA1")
res_barchart <- g_panendo
pdf(file=paste("plots/",date_time, "_box_panendo.pdf", sep=""),height=2,width=3.2)
for (i in c(1:length(res_barchart))){
library(tidyr)
untidy <- as.data.frame(t(tpm_panendo[res_barchart[i],]))
untidy$cases <- row.names(untidy)
untidy$cell_type_factor <- anno_panendo[row.names(untidy),'cell_type']
tidied <- untidy %>%
  gather(res_barchart[i], key="gene", value="tpm")
library(ggplot2)
p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor)) +
  geom_boxplot() + 
  labs(title = "", x = "", y = "tpm", fill = "") +
  scale_fill_manual(labels = c("Nasopharyngeal", "Squamous"), 
                    values=c("darkseagreen2", "chartreuse4")) +
  theme_classic() +
  theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12))
print(p)}
dev.off()



# 5.2 Panendoscopy Deconvolution ------------------------------------------
cases_panendo <- c(samples_PNS, samples_SQM)
decon_xcell <- read.table("xCell_deconvoluted.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
decon_xcell <- t(decon_xcell)
rownames(decon_xcell) <- substring(rownames(decon_xcell), 2)

decon_xcell_panendo <- as.data.frame(t(decon_xcell[cases_panendo,]))
res_barchart <- rownames(decon_xcell_panendo)

pdf(file=paste("plots/",date_time, "_box_decon_panendo.pdf", sep=""),height=2,width=3.5)
for (i in c(1:length(res_barchart))){
  library(tidyr)
  untidy <- as.data.frame(t(decon_xcell_panendo[res_barchart[i],]))
  untidy$cases <- row.names(untidy)
  untidy$cell_type_factor <- anno_panendo[row.names(untidy),'cell_type']
  tidied <- untidy %>%
    gather(res_barchart[i], key="gene", value="tpm")
  t_temp <- t.test(tidied$tpm[tidied$cell_type_factor == "PNS"], tidied$tpm[tidied$cell_type_factor == "SQM"],
         var.equal=TRUE, conf.level=0.95)
  library(ggplot2)
  p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor)) +
    geom_boxplot() + 
    labs(title = "", x = paste("p = ", signif(t_temp$p.value, 3)), y = "score", fill = "") +
    scale_fill_manual(labels = c("Nasopharyngeal", "Squamous"), 
                      values=c("darkseagreen2", "chartreuse4")) +
    theme_classic() +
    theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12))
  print(p)}
dev.off()


# 6. Plot Progression -----------------------------------------------------

cases_prog <- c(samples_PNS, samples_NAT, samples_DYS, samples_TUM, samples_TME)
count_prog <- count_fil[,cases_prog]
ens_temp <- rownames(count_prog)
length(unique(ens_temp))
ens_temp[1:57820] <- gsub("\\..*","",ens_temp[1:57820])
length(unique(ens_temp))

ens_temp[ens_temp %in% names(hugo_name_lookup)] <- hugo_name_lookup[ens_temp[ens_temp %in% names(hugo_name_lookup)]]
length(ens_temp)
table(duplicated(ens_temp))
rownames(count_prog) <- ens_temp

tpm_prog <- 1E6 * sweep(count_prog, 2, colSums(count_prog), "/")
tpm_prog <- tpm_prog[,cases_prog]
tpm_prog_lg2cen <- log2(tpm_prog +1)
tpm_prog_lg2cen <- t(apply(tpm_prog_lg2cen, 1, function(x) (x - mean(x)) / sd(x)))
anno_prog <- tdata_fil[cases_prog,]
anno_prog$cell_type_factor <- factor(anno_prog$cell_type, levels= cell_type_levels)
identical(rownames(anno_prog),colnames(tpm_prog))
# TRUE


# Progression plots
library(cowplot)
res_barchart <- c("CCL20", "CX3CL1", "SAA1")
pdf(file=paste("plots/",date_time, "_box_immune.pdf", sep=""),height=2.5,width=8)
myplots <- list()
for (i in c(1:length(res_barchart))){
  library(tidyr)
  untidy <- as.data.frame(t(tpm_prog[res_barchart[i],]))
  untidy$cases <- row.names(untidy)
  untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
  tidied <- untidy %>%
    gather(res_barchart[i], key="gene", value="tpm")
  library(ggplot2)
  p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor)) +
    geom_boxplot() + 
    labs(title = "", x = "", y = "", fill = "") +
    scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                      values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
    theme_classic() +
    theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=14), 
          legend.position = "none")
  myplots[[i]] <- p  # add each plot into plot list
}
cowplot::plot_grid(plotlist = myplots, ncol=4)
dev.off()


res_barchart <- "LMP-1"
pdf(file=paste("plots/",date_time, "_box_LMP1.pdf", sep=""),height=3,width=4)
for (i in c(1:length(res_barchart))){
  library(tidyr)
  untidy <- as.data.frame(t(tpm_prog[res_barchart[i],]))
  untidy$cases <- row.names(untidy)
  untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
  tidied <- untidy %>%
    gather(res_barchart[i], key="gene", value="tpm")
  library(ggplot2)
  p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor)) +
    geom_boxplot() + 
    labs(title = "", x = "", y = "tpm", fill = "") +
    scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                      values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
    theme_classic() +
    theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12))
  print(p)
}
dev.off()

g_spry <- c("SPRY1", "SPRY2")
pdf(file=paste("plots/",date_time, "_box_SPRY.pdf", sep=""),height=2.5,width=5)
res_barchart <- g_spry
library(tidyr)
untidy <- as.data.frame(t(tpm_prog[res_barchart,]))
untidy$cases <- row.names(untidy)
untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
tidied <- untidy %>%
  gather(res_barchart, key="gene", value="tpm")
library(ggplot2)
p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor))
p + geom_boxplot() + 
  labs(title = "", x = "", y = "tpm", fill = "") +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                    values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
  theme_classic() +
  theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12))
dev.off()


g_fgf <- c("FGF1", "FGF2")
pdf(file=paste("plots/",date_time, "_box_FGF.pdf", sep=""),height=2.5,width=5)
res_barchart <- g_fgf
library(tidyr)
untidy <- as.data.frame(t(tpm_prog[res_barchart,]))
untidy$cases <- row.names(untidy)
untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
tidied <- untidy %>%
  gather(res_barchart, key="gene", value="tpm")

# range of tpm for FGF1 & FGF2
range(tidied$tpm[tidied$cell_type_factor == "PNS" & tidied$gene =="FGF1"])
table(tidied$tpm[tidied$cell_type_factor == "TUM" & tidied$gene =="FGF1"] > 10.5)
range(tidied$tpm[tidied$cell_type_factor == "PNS" & tidied$gene =="FGF2"])
table(tidied$tpm[tidied$cell_type_factor == "TUM" & tidied$gene =="FGF2"] > 19.5)

library(ggplot2)
p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor))
p + geom_boxplot() + 
  labs(title = "", x = "", y = "tpm", fill = "") +
  scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                    values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
  theme_classic() +
  theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12))
dev.off()



g_fgfr <- c("FGFR1","FGFR2","FGFR3")
pdf(file=paste("plots/",date_time, "_box_FGFR.pdf", sep=""),height=2.5,width=5)
res_barchart <- g_fgfr
library(tidyr)
untidy <- as.data.frame(t(tpm_prog[res_barchart,]))
untidy$cases <- row.names(untidy)
untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
tidied <- untidy %>%
  gather(res_barchart, key="gene", value="tpm")
library(ggplot2)
p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor))
p + geom_boxplot() + 
  labs(title = "", x = "", y = "tpm", fill = "") +
  scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                    values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
  theme_classic() +
  theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12), legend.position = "none")
dev.off()

g_egf <- c("EGF", "TGFA", "BTC", "AREG", "EREG", "EPGN", "IGF1")
pdf(file=paste("plots/",date_time, "_box_EGF_A.pdf", sep=""),height=2.5,width=11)
res_barchart <- g_egf
library(tidyr)
untidy <- as.data.frame(t(tpm_prog[res_barchart,]))
untidy$cases <- row.names(untidy)
untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
tidied <- untidy %>%
  gather(res_barchart, key="gene", value="tpm")
library(ggplot2)
p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor))
p + geom_boxplot() + 
  labs(title = "", x = "", y = "tpm", fill = "") +
  scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                    values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
  theme_classic() +
  theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12), legend.position = "none")
dev.off()

g_egf <- c("HBEGF", "IGF2")
pdf(file=paste("plots/",date_time, "_box_EGF_B.pdf", sep=""),height=2.5,width=3.5)
res_barchart <- g_egf
library(tidyr)
untidy <- as.data.frame(t(tpm_prog[res_barchart,]))
untidy$cases <- row.names(untidy)
untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
tidied <- untidy %>%
  gather(res_barchart, key="gene", value="tpm")
library(ggplot2)
p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor))
p + geom_boxplot() + 
  labs(title = "", x = "", y = "tpm", fill = "") +
  scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                    values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
  theme_classic() +
  theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12), legend.position = "none")
dev.off()


g_egfr <- c("EGFR", "ERBB2", "ERBB3", "IGF1R")
pdf(file=paste("plots/",date_time, "_box_EGFR.pdf", sep=""),height=2.5,width=6.5)
res_barchart <- g_egfr
library(tidyr)
untidy <- as.data.frame(t(tpm_prog[res_barchart,]))
untidy$cases <- row.names(untidy)
untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
tidied <- untidy %>%
  gather(res_barchart, key="gene", value="tpm")
library(ggplot2)
p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor))
p + geom_boxplot() + 
  labs(title = "", x = "", y = "tpm", fill = "") +
  scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                    values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
  theme_classic() +
  theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12), legend.position = "none")
dev.off()

g_vegf <- c("VEGFA", "VEGFB")
pdf(file=paste("plots/",date_time, "_box_VEGF.pdf", sep=""),height=2.5,width=3.5)
res_barchart <- g_vegf
library(tidyr)
untidy <- as.data.frame(t(tpm_prog[res_barchart,]))
untidy$cases <- row.names(untidy)
untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
tidied <- untidy %>%
  gather(res_barchart, key="gene", value="tpm")
library(ggplot2)
p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor))
p + geom_boxplot() + 
  labs(title = "", x = "", y = "tpm", fill = "") +
  scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                    values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
  theme_classic() +
  theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12), legend.position = "none")
dev.off()


g_vegfr <- c("FLT1", "KDR", "FLT4")
pdf(file=paste("plots/",date_time, "_box_VEGFR.pdf", sep=""),height=2.5,width=5)
res_barchart <- g_vegfr
library(tidyr)
untidy <- as.data.frame(t(tpm_prog[res_barchart,]))
untidy$cases <- row.names(untidy)
untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
tidied <- untidy %>%
  gather(res_barchart, key="gene", value="tpm")
library(ggplot2)
p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor))
p + geom_boxplot() + 
  labs(title = "", x = "", y = "tpm", fill = "") +
  scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                    values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
  theme_classic() +
  theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12), legend.position = "none")
dev.off()

g_pdgf <- c("PDGFA", "PDGFB")
pdf(file=paste("plots/",date_time, "_box_PDGF.pdf", sep=""),height=2.5,width=3.5)
res_barchart <- g_pdgf
library(tidyr)
untidy <- as.data.frame(t(tpm_prog[res_barchart,]))
untidy$cases <- row.names(untidy)
untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
tidied <- untidy %>%
  gather(res_barchart, key="gene", value="tpm")
library(ggplot2)
p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor))
p + geom_boxplot() + 
  labs(title = "", x = "", y = "tpm", fill = "") +
  scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                    values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
  theme_classic() +
  theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12), legend.position = "none")
dev.off()

g_pdgfr <- c("PDGFRA", "PDGFRB")
pdf(file=paste("plots/",date_time, "_box_PDGFR.pdf", sep=""),height=2.5,width=3.5)
res_barchart <- g_pdgfr
library(tidyr)
untidy <- as.data.frame(t(tpm_prog[res_barchart,]))
untidy$cases <- row.names(untidy)
untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
tidied <- untidy %>%
  gather(res_barchart, key="gene", value="tpm")
library(ggplot2)
p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor))
p + geom_boxplot() + 
  labs(title = "", x = "", y = "tpm", fill = "") +
  scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                    values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
  theme_classic() +
  theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12), legend.position = "none")
dev.off()


res_barchart <- c("LGALS1", "LGALS3")
pdf(file=paste("plots/",date_time, "_box_GAL.pdf", sep=""),height=2.5,width=5)
library(tidyr)
untidy <- as.data.frame(t(tpm_prog[res_barchart,]))
untidy$cases <- row.names(untidy)
untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
tidied <- untidy %>%
  gather(res_barchart, key="gene", value="tpm")
library(ggplot2)
p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor))
p + geom_boxplot() + 
  labs(title = "", x = "", y = "tpm", fill = "") +
  scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                    values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
  theme_classic() +
  theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12))
dev.off()


res_barchart <- c("RELA", "NFKB1", "TNFRSF1A", "TRAF2", "NFKBIA")
pdf(file=paste("plots/",date_time, "_box_nfkb_canonical.pdf", sep=""),height=2.5,width=10)
myplots <- list()
for (i in c(1:length(res_barchart))){
  library(tidyr)
  untidy <- as.data.frame(t(tpm_prog[res_barchart[i],]))
  untidy$cases <- row.names(untidy)
  untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
  tidied <- untidy %>%
    gather(res_barchart[i], key="gene", value="tpm")
  library(ggplot2)
  p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor)) +
    geom_boxplot() + 
    labs(title = "", x = "", y = "", fill = "") +
    scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                      values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
    theme_classic() +
    theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=15), 
          legend.position = "none")
  myplots[[i]] <- p  # add each plot into plot list
}
cowplot::plot_grid(plotlist = myplots, ncol=5)
dev.off()

res_barchart <- c("RELB", "NFKB2", "TNFRSF1B", "CD40")
pdf(file=paste("plots/",date_time, "_box_nfkb_noncanonical.pdf", sep=""),height=2.5,width=10)
myplots <- list()
for (i in c(1:length(res_barchart))){
  library(tidyr)
  untidy <- as.data.frame(t(tpm_prog[res_barchart[i],]))
  untidy$cases <- row.names(untidy)
  untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
  tidied <- untidy %>%
    gather(res_barchart[i], key="gene", value="tpm")
  library(ggplot2)
  p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor)) +
    geom_boxplot() + 
    labs(title = "", x = "", y = "", fill = "") +
    scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                      values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
    theme_classic() +
    theme(axis.text.y=element_text(size=12), axis.text.x=element_text(size=15), 
          legend.position = "none")
  myplots[[i]] <- p  # add each plot into plot list
}
cowplot::plot_grid(plotlist = myplots, ncol=5)
dev.off()



# 7. Correlation Plots for TUM samples -----------------------------------------------

tpm_tum <- tpm_ns[,samples_TUM]
ens_temp <- rownames(tpm_tum)
length(unique(ens_temp))
ens_temp[1:57820] <- gsub("\\..*","",ens_temp[1:57820])
length(unique(ens_temp))

ens_temp[ens_temp %in% names(hugo_name_lookup)] <- hugo_name_lookup[ens_temp[ens_temp %in% names(hugo_name_lookup)]]
length(ens_temp)
table(duplicated(ens_temp))
rownames(tpm_tum) <- ens_temp
t_tpm_tum <- as.data.frame(t(tpm_tum))

pdf(file=paste("plots/",date_time, "_LMP1_correlation_titles.pdf", sep=""),height=5,width=7.5)
par(mfrow=c(2,3))

c <- cor.test(log(t_tpm_tum$'LMP-1'+1),log(t_tpm_tum$NFKB1+1), method="pearson")
cs <- paste("r = ", signif(c$estimate, digits = 3), ",", "p = ", signif(c$p.value, digits = 3))
plot(log(t_tpm_tum$'LMP-1'+1),log(t_tpm_tum$NFKB1 +1), xlab = paste("LMP1 log2(tpm) \n", cs), 
     ylab = "NFKB1 log2(tpm)" , pch=20, main = "NFKB1")
abline(lm(log(t_tpm_tum$NFKB1 +1) ~ log(t_tpm_tum$'LMP-1'+1)), col = "red")

c <- cor.test(log(t_tpm_tum$'LMP-1'+1),log(t_tpm_tum$NFKB2+1), method="pearson")
cs <- paste("r = ", signif(c$estimate, digits = 3), ",", "p = ", signif(c$p.value, digits = 3))
plot(log(t_tpm_tum$'LMP-1'+1),log(t_tpm_tum$NFKB2 +1), xlab = paste("LMP1 log2(tpm) \n", cs), 
     ylab = "NFKB2 log2(tpm)" , pch=20, main = "NFKB2")
abline(lm(log(t_tpm_tum$NFKB2 +1) ~ log(t_tpm_tum$'LMP-1'+1)), col = "red")

c <- cor.test(log(t_tpm_tum$'LMP-1'+1),log(t_tpm_tum$RELA+1), method="pearson")
cs <- paste("r = ", signif(c$estimate, digits = 3), ",", "p = ", signif(c$p.value, digits = 3))
plot(log(t_tpm_tum$'LMP-1'+1),log(t_tpm_tum$RELA +1), xlab = paste("LMP1 log2(tpm) \n", cs), 
     ylab = "RELA log2(tpm)" , pch=20, main = "RELA")
abline(lm(log(t_tpm_tum$RELA +1) ~ log(t_tpm_tum$'LMP-1'+1)), col = "red")

c <- cor.test(log(t_tpm_tum$'LMP-1'+1),log(t_tpm_tum$RELB+1), method="pearson")
cs <- paste("r = ", signif(c$estimate, digits = 3), ",", "p = ", signif(c$p.value, digits = 3))
plot(log(t_tpm_tum$'LMP-1'+1),log(t_tpm_tum$RELB +1), xlab = paste("LMP1 log2(tpm) \n", cs), 
     ylab = "RELB log2(tpm)" , pch=20, main = "RELB")
abline(lm(log(t_tpm_tum$RELB +1) ~ log(t_tpm_tum$'LMP-1'+1)), col = "red")

c <- cor.test(log(t_tpm_tum$'LMP-1'+1),log(t_tpm_tum$NFKBIA+1), method="pearson")
cs <- paste("r = ", signif(c$estimate, digits = 3), ",", "p = ", signif(c$p.value, digits = 3))
plot(log(t_tpm_tum$'LMP-1'+1),log(t_tpm_tum$NFKBIA +1), xlab = paste("LMP1 log2(tpm) \n", cs), 
     ylab = "NFKBIA log2(tpm)" , pch=20, main = "NFKBIA" )
abline(lm(log(t_tpm_tum$NFKBIA +1) ~ log(t_tpm_tum$'LMP-1'+1)), col = "red")
dev.off()

# 8. Deconvoluted samples -----------------------------------------------

# read in deconvoluted profiles
decon_cib <- read.table("CIBERSORTx_deconvoluted.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)

decon_xcell <- read.table("xCell_deconvoluted.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
decon_xcell <- t(decon_xcell)
rownames(decon_xcell) <- substring(rownames(decon_xcell), 2)


# 8.1 Deconvoluted progression --------------------------------------------

cases_prog <- c(samples_PNS, samples_NAT, samples_DYS, samples_TUM, samples_TME)
decon_xcell_t <- t(decon_xcell)
decon_xcell_prog <- decon_xcell_t[,cases_prog]
identical(rownames(anno_prog),colnames(decon_xcell_prog))
# TRUE
decon_xcell_prog <- as.data.frame(t(decon_xcell[cases_prog,]))


#       Progression plots for patient samples only -------------------------------------------------------

pdf(file=paste("plots/",date_time, "_box_decon_progression.pdf", sep=""),height=2,width=3.5)
res_barchart <- rownames(decon_xcell_prog)
for (i in c(1:length(res_barchart))){
  library(tidyr)
  untidy <- as.data.frame(t(decon_xcell_prog[res_barchart[i],]))
  untidy$cases <- row.names(untidy)
  untidy$cell_type_factor <- anno_prog[row.names(untidy),'cell_type_factor']
    tidied <- untidy %>%
    gather(res_barchart[i], key="gene", value="tpm")
  library(ggplot2)
  p <- ggplot(tidied, aes(gene,tpm,fill = cell_type_factor)) + 
    geom_boxplot() + 
    labs(title = "", x = "", y = "score", fill = "") +
    scale_fill_manual(labels = c("Normal", "Normal-adjacent", "Dysplastic", "Tumor", "Microenvironment"), 
                      values=c("darkseagreen2", "khaki3","gold3","palevioletred2","paleturquoise3")) +
    # scale_y_continuous(limits = c(0, 100)) +
    theme_classic() +
    theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12))
  print(p)}
dev.off()



# 8.2. TME correlates -----------------------------------------------

# remove any duplicates
tdata_paired_tme <- tdata_ns[tdata_ns$cell_type == "TME",]
tdata_paired_tme <- tdata_paired_tme[duplicated(tdata_paired_tme$patient_index) == FALSE,]
tdata_paired_tme_patient_index <- tdata_paired_tme$patient_index
unique(tdata_paired_tme_patient_index)

tdata_paired_tum <- tdata_ns[tdata_ns$cell_type == "TUM",]
tdata_paired_tum <- tdata_paired_tum[duplicated(tdata_paired_tum$patient_index) == FALSE,]
tdata_paired_tum_patient_index <- tdata_paired_tum$patient_index

# get those with both TME and TUM
tdata_paired_tumtme_patient_index <- intersect(tdata_paired_tme_patient_index, tdata_paired_tum_patient_index)
# sort both correctly
tdata_paired_tme <- tdata_paired_tme[tdata_paired_tme$patient_index %in% tdata_paired_tumtme_patient_index,]
tdata_paired_tme <- tdata_paired_tme[order(tdata_paired_tme$patient_index),]

tdata_paired_tum <- tdata_paired_tum[tdata_paired_tum$patient_index %in% tdata_paired_tumtme_patient_index,]
tdata_paired_tum <- tdata_paired_tum[order(tdata_paired_tum$patient_index),]
identical(tdata_paired_tme$patient_index, tdata_paired_tum$patient_index)

paired_tme_cases <- rownames(tdata_paired_tme)
paired_tum_cases <- rownames(tdata_paired_tum)

tpm_paired_tme <- as.data.frame(t(tpm_prog[,paired_tme_cases]))
tpm_paired_tum <- as.data.frame(t(tpm_prog[,paired_tum_cases]))
decon_xcell_tme_paired <- as.data.frame(decon_xcell[paired_tme_cases,])
decon_cib_tme_paired <- as.data.frame(decon_cib[paired_tme_cases,])
dim(decon_cib_tme_paired)

# plot

pdf(file=paste("plots/",date_time, "_corr_CXCLs_CCL20.pdf", sep=""),height=3.5,width=3.5)

x_temp <- tpm_paired_tum$CXCL11
y_temp <- tpm_paired_tme$CXCR3
cor_temp <- cor.test(log2(x_temp+1),log2(y_temp+1), method="pearson")
plot(log2(x_temp+1),log2(y_temp +1), xlab = "Tumor CXCL11 log2(tpm)", ylab = "TME CXCR3 log2(tpm)", pch=20,
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2)
abline(lm(log2(y_temp +1) ~ log2(x_temp+1)), col = "red")
cor.test(log2(x_temp+1),log2(y_temp+1), method="pearson")

x_temp <- tpm_paired_tum$CCL20
y_temp <- tpm_paired_tme$CCR6
cor_temp <- cor.test(log2(x_temp+1),log2(y_temp+1), method="pearson")
plot(log2(x_temp+1),log2(y_temp +1), xlab = "Tumor CCL20 log2(tpm)", ylab = "TME CCR6 log2(tpm)", pch=20,
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2)
abline(lm(log2(y_temp +1) ~ log2(x_temp+1)), col = "red")

dev.off()


pdf(file=paste("plots/",date_time, "_corr_CCL20_cib.pdf", sep=""),height=3.5,width=3.5)
x_temp <- tpm_paired_tum$`CCL20`
z_temp <- decon_cib_tme_paired
colnames(z_temp)[1]
for (i in colnames(z_temp)) {
  y_temp <- z_temp[,i]
  cor_temp <- cor.test(log(x_temp+1),log(y_temp+1), method="pearson")
  plot(log(x_temp+1),log(y_temp +1), xlab = "Tumor CCL20 log2(tpm)", ylab = i, pch=20, font.lab = 2,
       sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)))
  abline(lm(log(y_temp +1) ~ log(x_temp+1)), col = "red")
}
dev.off()

pdf(file=paste("plots/",date_time, "_corr_CCL20_xcell.pdf", sep=""),height=3.5,width=3.5)
x_temp <- tpm_paired_tum$`CCL20`
z_temp <- decon_xcell_tme_paired
colnames(z_temp)[1]
for (i in colnames(z_temp)) {
  y_temp <- z_temp[,i]
  cor_temp <- cor.test(log(x_temp+1),log(y_temp+1), method="pearson")
  plot(log(x_temp+1),log(y_temp +1), xlab = "tumor CCL20 tpm", ylab = i, pch=20, font.lab=2,
       sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)))
  abline(lm(log(y_temp +1) ~ log(x_temp+1)), col = "red")
}
dev.off()




# 8.3 TME Heatmap ----------------------------------------------------

decon_cib_tme <- as.data.frame(decon_cib[paired_tme_cases,])
dim(decon_cib_tme)
decon_cib_tme_hm <- as.data.frame(t(decon_cib_tme[,c(1:22)]))


##  get LMP-1 status of the paired tumor
tcdata_tme <- tcdata[paired_tme_cases,]
tcdata_tme$LMP1_tumor <- NA
tcdata_tme$LMP1_tumor <- as.numeric(tpm_paired_tum$`LMP-1`)
tcdata_tme$LMP1_tumor_scale <- scale(tcdata_tme$LMP1_tumor, center = TRUE, scale = TRUE)
tcdata_tme$'LMP1 status' <- NA
tcdata_tme$'LMP1 status'[tcdata_tme$LMP1_tumor > 0 & tcdata_tme$LMP1_tumor < 5.7] <- "low"
tcdata_tme$'LMP1 status'[tcdata_tme$LMP1_tumor == 0] <- "negative"
tcdata_tme$'LMP1 status'[tcdata_tme$LMP1_tumor >= 5.7] <- "high"
tcdata_tme$'LMP1 status'

tcdata_tme$'Any recurrence' <- tcdata_tme$Recurrence_Any
tcdata_tme$'Histology pattern' <- tcdata_tme$Histo_Pattern

tcdata_tme$M <- droplevels(tcdata_tme$M) # remove X levels
tcdata_tme$Stage <- recode(tcdata_tme$Stage, '2'="II", '3'="III", '4'="IV")

anno_cib_tme <- tcdata_tme[,c("Any recurrence", "Stage", "M", "N", "T", "Histology pattern", "LMP1 status"), drop=FALSE]
colorlist <- list(type = c("tumor"="#EB6E80", "normal"="#008F95"), 
                  T=c("1"="#e6e6e6", "2"="#b3b3b3", "3"="#808080", "4"="#4d4d4d", "NA"="#ffffff"),
                  N=c("0"="#e6e6e6", "1"="#b3b3b3", "2"="#808080", "3"="#4d4d4d", "NA"="#ffffff"),
                  M=c("0"="#e6e6e6", "1"="#4d4d4d", "NA"="#ffffff"),
                  Stage=c("I"="#e6e6e6", "II"="#b3b3b3", "III"="#808080", "IV"="#4d4d4d", "NA"="#ffffff"),
                  'Any recurrence'=c("yes"="#29a329", "no"="#d6f5d6", "NA"="#ffffff"),
                  'Histology pattern' = c("Schmincke"="skyblue1", "Regaud"="skyblue3","Undetermined"="seashell3"),
                  'LMP1 status'=c("negative"="seashell3", "low"="plum1", "high"= "orchid3"))


pdf(file=paste(date_time, "tme_hm.pdf", sep=""),height=6,width=7)

library(pheatmap)
hm_cib_tme <- pheatmap(
  decon_cib_tme_hm,
  annotation_col = anno_cib_tme,
  annotation_colors = colorlist,
  show_colnames = T,
  show_rownames = T,
  fontsize_col = 8,
  fontsize_row = 8,
  cluster_rows=F,
  cluster_cols=T,
  clustering_method = "complete"
)
print(hm_cib_tme)
dev.off()


# Save --------------------------------------------------------------------
setwd(output_directory)
save.image("NPC_s3s_code_environment.Rdata")

