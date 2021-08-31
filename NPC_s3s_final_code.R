# 1. Head -----------------------------------------------------------------

output_directory <- "/Users/joshuatay/Documents/Research/NPC_s3s_paper/final_code"
setwd(output_directory)
# load("NPC_s3s_paper_final.Rdata")
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
tpm_fil <- 1E6 * sweep(count_fil, 2, colSums(count_fil), "/")

tdata_fil <- read.table("technical_data.txt", header = TRUE, sep = "\t", as.is = TRUE) # read in technical annotation file
rownames(tdata_fil) <- tdata_fil$sample_ID
identical(rownames(tdata_fil), colnames(count_fil))

stats_align <- read.table("stats_star_featurecounts.txt", header = TRUE, sep = "\t", as.is = TRUE)   
rownames(stats_align) <- stats_align$X
stats_align <- stats_align[colnames(count_fil),]
identical(rownames(stats_align), colnames(count_fil))

cdata <- read.table("clinical_data.txt", header = TRUE, sep = "\t", as.is = TRUE) # read in clinical annotation file
tcdata <- merge(tdata_fil,cdata,"patient_index", all=T)
tcdata <- tcdata[!is.na(tcdata$sample_ID),]
rownames(tcdata) <- tcdata$sample_ID

# 3. Clustering ---------------------------------------------------------------------

# define colors and symbols
cell_type_levels <- c("PNS", "NAT", "DYS", "TUM", "TME", "CEL", "SQM", "SQL")
cell_type_levels_text <- c("normal", "normal-adjacent to tumor", "dysplastic", "tumor", "microenvironment", "C666-1 cell line", "squamous", "squamous-inflamed")
cell_type_symbols <- c(1,1,15,15,3,11,19,10)
cell_type_colors <- c("black","green3","hotpink","black","black","chocolate4","blue","blue")
batch_levels <- c("A1", "A2", "A3", "A4", "A5", "A6", "H1")
batch_colors <- c("green3", "lightblue3", "hotpink3",  "lavenderblush3", "orange2", "lightgoldenrod4", "royalblue2")



#   3.1 PCA for nasopharyngeal cells only ------------------------------------------------------

ns_libraries <- rownames(tdata_fil)[!tdata_fil$cell_type %in% c("SQM", "SQL")]
count_ns <- count_fil[,ns_libraries]
tdata_ns <- tdata_fil[ns_libraries,]

tpm_ns <- 1E6 * sweep(count_ns, 2, colSums(count_ns), "/")
tpm_ns_high <- tpm_ns[rowMeans(tpm_ns)>15,]
dim(tpm_ns_high)
tpm_ns_high_log2 <- log2(tpm_ns_high + 1)
tpm_ns_high_log2cen <- tpm_ns_high_log2 - rowMeans(tpm_ns_high_log2)

svd_ns <- svd(tpm_ns_high_log2cen)
names(svd_ns) # Look at the percent variance explained
pdf(file=paste("plots/",date_time, "_PCA_nasopharyngeal_variance.pdf", sep=""),height=5,width=6.25)
plot(svd_ns$d^2/sum(svd_ns$d^2),ylab="Percent Variance Explained",xlab="Principal Component", col =" royalblue",type ="o", pch=19)
dev.off()
percent <- round(100*(svd_ns$d^2/sum(svd_ns$d^2)),digits=2)

table(tdata_ns$cell_type)
tdata_ns$cell_type_factor <- factor(tdata_ns$cell_type, levels= cell_type_levels)

pdf(file=paste(date_time, "_PCA_nasopharyngeal_text.pdf", sep=""),height=5,width=7)
plot(svd_ns$v[,1],svd_ns$v[,2],ylab=paste("PC2 (", percent[2], "% of variance)", sep=""), xlab=paste("PC1 (", percent[1], "% of variance)", sep=""),
     pch = cell_type_symbols[as.numeric(tdata_ns$cell_type_factor)], 
     col=cell_type_colors[as.numeric(tdata_ns$cell_type_factor)],
     xlim=c(-0.22, 0.27), ylim=c(-0.13, 0.22),
     main = paste("PCA by Cell Type, n = ", length(rownames(tdata_ns)), sep="")
     )
with(tdata_ns, legend("bottomright", pch = cell_type_symbols[1:6], col=cell_type_colors[1:6], legend = cell_type_levels_text[1:6], title="  Cell type  ", cex=0.8))
dev.off()


# by batch
tdata_ns$batch_factor <- factor(tdata_ns$batch, levels= batch_levels)
plot(svd_ns$v[,1],svd_ns$v[,2],ylab=paste("PC2 (", percent[2], "% of variance)", sep=""), xlab=paste("PC1 (", percent[1], "% of variance)", sep=""),
     pch = cell_type_symbols[as.numeric(tdata_ns$cell_type_factor)], 
     col=batch_colors[as.numeric(tdata_ns$batch_factor)],
     xlim=c(-0.22, 0.24), ylim=c(-0.13, 0.22),
     main = "Library Batch")
with(tdata_ns, legend("bottomright", pch = 15, col=batch_colors, legend = batch_levels, title="  Batch  "))

# plot by age
redbluepal <- colorRampPalette(c('red','yellow','green','blue'))
plot(svd_ns$v[,1],svd_ns$v[,2],ylab=paste("PC2 (", percent[2], "% of variance)", sep=""), xlab=paste("PC1 (", percent[1], "% of variance)", sep=""),
     pch = cell_type_symbols[as.numeric(tdata_ns$cell_type_factor)], 
     col = redbluepal(25)[as.numeric(tdata_ns$block_age)+1],
     xlim=c(-0.22, 0.24), ylim=c(-0.13, 0.22),
     main = "Age of Sample")
legend("bottomright", title="Age (Years)", legend=c(0,6,12,18,24), col=redbluepal(25)[c(1,7,13,19,25)], pch=20)

#   3.2 UMAP for unique primary tumors ------------------------------------------------------

# select unique primary tumor libraries based on percentage mappability
tum_libraries <- rownames(tdata_fil)[tdata_fil$cell_type %in% c("TUM") & tdata_fil$recurrence_case == "no"] # exclude recurrent cases
stats_align_tum <- stats_align[tum_libraries,]
stats_align_tum <- stats_align_tum[order(stats_align_tum$STAR_mqc.generalstats.star.uniquely_mapped_percent, decreasing =TRUE), ] # order based on mapping percentage
tdata_tum_unq <- tdata_fil[rownames(stats_align_tum),]
identical(rownames(stats_align_tum), rownames(tdata_tum_unq))
length(unique(tdata_tum_unq$patient_index))
samples_TUM_UNQ <- rownames(tdata_tum_unq)[duplicated(tdata_tum_unq$patient_index) == FALSE]

count_tum_unq <- count_fil[,samples_TUM_UNQ]
tdata_tum_unq <- tdata_fil[samples_TUM_UNQ,]

tpm_tum_unq <- 1E6 * sweep(count_tum_unq, 2, colSums(count_tum_unq), "/")
tum_hm <- tpm_tum_unq[rowMeans(tpm_tum_unq)>10,]
dim(tum_hm)
library(uwot)
set.seed(1234)
tum_um <- umap(t(tum_hm), 
               # n_neighbors = 15,
               metric = "correlation")
plot(tum_um[,1], tum_um[,2])
tum_um_df <- as.data.frame(tum_um)
rownames(tum_um_df) <- rownames(tdata_tum_unq)
samples_TUM_UNQ_cluster1 <- rownames(tum_um_df)[tum_um_df$V1 < 1.4]
samples_TUM_UNQ_cluster2 <- rownames(tum_um_df)[tum_um_df$V1 > 1.4]

tcdata_TUM_UNQ <- tcdata[rownames(tdata_tum_unq),]
dim(tcdata_TUM_UNQ)
tcdata_TUM_UNQ$cluster <- NA
tcdata_TUM_UNQ$cluster[tcdata_TUM_UNQ$sample_ID %in% samples_TUM_UNQ_cluster1] <- 1
tcdata_TUM_UNQ$cluster[tcdata_TUM_UNQ$sample_ID %in% samples_TUM_UNQ_cluster2] <- 2
table(tcdata_TUM_UNQ$cluster, tcdata_TUM_UNQ$Recurrence_Any)

identical(rownames(tcdata_TUM_UNQ), colnames(tpm_tum_unq))
identical(rownames(tcdata_TUM_UNQ), colnames(tum_hm))

png(file=paste(date_time, "_", "TUM_UNIQUE_umap.png", sep=""), res = 300,  height=5, width=6, unit = "in")
plot(tum_um[,1], tum_um[,2],
     pch = 15, 
     cex = 1.2, # xlim=c(-0.5, 0.2), ylim=c(-0.2, 0.6),
     col=c("black", "red")[as.factor(tcdata_TUM_UNQ$Recurrence_Any)],
     main = paste("Microdissected tumor regions, n =", length(rownames(tcdata_TUM_UNQ))),
     ylab = "UMAP2", xlab = "UMAP1"
)
legend("topleft", pch = 15, pt.cex = 1.2,  col=c("black", "red"), 
       legend = c("No", "Yes"), title="Recurrence", cex=0.9)
dev.off()

# 3.3 PCA for panendoscopy ------------------------------------------------

pan_libraries <- rownames(tdata_fil)[tdata_fil$cell_type %in% c("PNS", "SQM", "SQL")]
count_pan <- count_fil[,pan_libraries]
tdata_pan <- tdata_fil[pan_libraries,]

tpm_pan <- 1E6 * sweep(count_pan, 2, colSums(count_pan), "/")
tpm_pan_high <- tpm_pan[rowMeans(tpm_pan)>15,]
dim(tpm_pan_high)
tpm_pan_high_log2 <- log2(tpm_pan_high + 1)
tpm_pan_high_log2cen <- tpm_pan_high_log2 - rowMeans(tpm_pan_high_log2)

svd_pan <- svd(tpm_pan_high_log2cen)
names(svd_pan) # Look at the percent variance explained
plot(svd_pan$d^2/sum(svd_pan$d^2),ylab="Percent Variance Explained",xlab="Principal Component", col =" royalblue",type ="o", pch=19)
percent <- round(100*(svd_pan$d^2/sum(svd_pan$d^2)),digits=2)

table(tdata_pan$cell_type_epi)
cell_type_epi_levels <- c("BOT", "PIR", "PAL", "TONLP", "TONLR", "PNS")
tdata_pan$cell_type_factor <- factor(tdata_pan$cell_type_epi, levels=cell_type_epi_levels)

pdf(file=paste(date_time, "_PCA_upperairway.pdf", sep=""),height=5,width=7)
plot(svd_pan$v[,1],svd_pan$v[,2],ylab=paste("PC2 (", percent[2], "% of variance)", sep=""), xlab=paste("PC1 (", percent[1], "% of variance)", sep=""),
     pch = c(11,8,2,9,3,0)[as.numeric(tdata_pan$cell_type_factor)], 
     xlim = c(-0.3, 0.6), ylim = c(-0.7, 0.4))
     # main = "Gene Expression of Upper Aiway Epithelium (n = 23)")
with(tdata_pan, legend("bottomright", pch = c(11,8,2,9,3,0), col="black", 
                     legend = c("Base of Tongue", "Piriform Sinus", "Palate", 
                                "Tonsil (Lymphoid-Poor)", "Tonsil (Lymphoid-Rich)", "Nasopharynx"), 
                     title = expression(bold("Location")), cex=1))
dev.off()

# 4. DE-gene and FGSEA functions ------------------------------------------------------
de_function <- function(de_counts, de_pdata, group_header)
{
  library(DESeq2)
  library(BiocParallel)
  dds <- DESeqDataSetFromMatrix(de_counts, colData = de_pdata, design = ~ temp_group)
  dds2 <- DESeq(dds, fitType = "local", parallel = T, BPPARAM = MulticoreParam(workers = 15)) 
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
  WriteXLS(res, ExcelFileName=paste(date_time, "_", group_header, "_deseq.xlsx", sep=""), row.names=TRUE)
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
  pdf(file=paste(date_time, "_", de_header, "_fgsea_top.pdf", sep=""),height=8,width=12)
  plotGseaTable(path_go[topPathways], ranks_res, fg_res_go, gseaParam = 0.5)
  dev.off()
  
  library(data.table)
  fwrite(fg_res_go_tidy, file=paste(date_time, "_", de_header, "_fgsea_GO.txt", sep=""), 
         sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
  fg_res_go_tidy <- fg_res_go_tidy[fg_res_go_tidy$padj < 0.05,]
  fwrite(fg_res_go_tidy, file=paste(date_time, "_", de_header, "_fgsea_GO_sig.txt", sep=""), 
         sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)

  return(fg_res_go)
}


# 5. Define groups, load gene sets --------------------------------------------------------
# gene sets can be downloaded from http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

samples_PNS <- rownames(tdata_ns)[tdata_ns$cell_type == "PNS"]
samples_NAT <- rownames(tdata_ns)[tdata_ns$cell_type == "NAT"]
samples_NAT_UNQ <- rownames(tdata_ns)[tdata_ns$cell_type == "NAT" & tdata_ns$sample_ID != "10204_A35_E1_NAT"] # A35_E1_NAT has a duplicate, remove the sample with lower mappability.
samples_DYS <- rownames(tdata_ns)[tdata_ns$cell_type == "DYS"]
samples_TUM <- rownames(tdata_ns)[tdata_ns$cell_type == "TUM"]
samples_CEL <- rownames(tdata_ns)[tdata_ns$cell_type == "CEL"]
samples_TME <- rownames(tdata_ns)[tdata_fil$cell_type == "TME"]
samples_TME_UNQ <- rownames(tdata_ns)[tdata_ns$cell_type == "TME" & tdata_ns$sample_ID != "10258_A34_S1_TME"] # remove one duplicate


# pathways 
library(fgsea)
path_hallmark <- gmtPathways("gmt/h.all.v6.2.symbols.gmt")
path_hallmark_descriptions <- read.table("gmt/hallmark_descriptions.csv", header = TRUE, sep = "\t", as.is = TRUE)
path_hallmark_immune <- path_hallmark_descriptions$Hallmark.Name[path_hallmark_descriptions$Process.Category == "immune"]
path_hallmark_proliferation <- path_hallmark_descriptions$Hallmark.Name[path_hallmark_descriptions$Process.Category == "proliferation"]
path_hallmark_additional <- c("TNFA_SIGNALING_VIA_NFKB")

path_allmsdb <- gmtPathways("gmt/msigdb.v6.2.symbols.gmt")
path_go <- gmtPathways("gmt/c5.all.v6.2.symbols.gmt")
path_go_bioprocess <- gmtPathways("gmt/c5.bp.v6.2.symbols.gmt")
path_go_cellcomp <- gmtPathways("gmt/c5.cc.v6.2.symbols.gmt")
path_go_molfx <- gmtPathways("gmt/c5.mf.v6.2.symbols.gmt")


#   5.1 Compare Normal-adjacent VS Normal (NAT_vs_PNS) ---------------------------------------------
grp_header <- "NAT_vs_PNS"
grp_A <- samples_NAT_UNQ
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

library(fgsea)
library(data.table)
# FGSEA Hallmark
fg_NAT_vs_PNS_hallmark <- fgsea(pathways=path_hallmark, stats=ranks_NAT_vs_PNS, nperm=1000)
fg_NAT_vs_PNS_hallmark_tidy <- fg_NAT_vs_PNS_hallmark %>%
  as_tibble() %>%
  arrange(desc(ES))
fwrite(fg_NAT_vs_PNS_hallmark_tidy, file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
fwrite(fg_NAT_vs_PNS_hallmark_tidy[fg_NAT_vs_PNS_hallmark_tidy$padj < 0.05,], file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK_sig.txt", sep=""), 
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


pdf(file=paste(date_time, "_vol_NATvsPNS.pdf", sep=""),height=5,width=6.25)
with(vp, plot(log2FoldChange, -log10(padj), pch=20, main=" ", col=alpha("black", 0.2), xlim=c(-8,8)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_innate,], points(log2FoldChange, -log10(padj), pch=20, cex=1.2, col=alpha("yellow", 0.6)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_virus_response,], points(log2FoldChange, -log10(padj), pch=20, cex=1.2, col=alpha("orange", 0.8)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_response_to_interferon_type1,], points(log2FoldChange, -log10(padj), pch=20, cex=1.2, col=alpha("red", 0.6)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_spry,], points(log2FoldChange, -log10(padj), pch=20, cex=1.2, col=alpha("black", 0.9)))
vp_select <- vp[vp$symbol %in% vp_text,]
library(basicPlotteR)
addTextLabels(vp_select$log2FoldChange, -log10(vp_select$padj), vp_select$symbol, 
              cex.label = 0.8, 
              col.label="black", col.background = "white")
legend("topleft", pch = c(16,16,16), col=c("yellow","orange","red"), 
         legend = c("innate immune response", "defense response to virus", "response to type I interferon"), 
         title = expression(bold("Pathway")), cex=0.8)
dev.off()

#   5.2 Compare Dysplasia VS Normal-adjacent (DYS_vs_NAT) ---------------------------------------------
grp_header <- "DYS_vs_NAT"
grp_A <- samples_DYS
grp_B <- samples_NAT_UNQ
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
fwrite(fg_DYS_vs_NAT_hallmark_tidy, file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
fwrite(fg_DYS_vs_NAT_hallmark_tidy[fg_DYS_vs_NAT_hallmark_tidy$padj < 0.05,], file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK_sig.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)


# FGSEA GO
library(fgsea)
fg_DYS_vs_NAT_go <- fgsea_function(res_DYS_vs_NAT, grp_header)


#   5.3 Compare Tumor VS Dysplasia (TUM_vs_DYS) ---------------------------------------------
grp_header <- "TUM_vs_DYS"
grp_A <- samples_TUM_UNQ
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
fwrite(fg_TUM_vs_DYS_hallmark_tidy, file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
fwrite(fg_TUM_vs_DYS_hallmark_tidy[fg_TUM_vs_DYS_hallmark_tidy$padj < 0.05,], file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK_sig.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)


# FGSEA GO
library(fgsea)
fg_TUM_vs_DYS_go <- fgsea_function(res_TUM_vs_DYS, grp_header)


#   5.4 Compare Tumor VS Normal (TUM_vs_PNS) -------------------------------------------------------------------------

grp_header <- "TUM_vs_PNS"
grp_A <- samples_TUM_UNQ
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
fwrite(fg_TUM_vs_PNS_hallmark_tidy, file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
fwrite(fg_TUM_vs_PNS_hallmark_tidy[fg_TUM_vs_PNS_hallmark_tidy$padj < 0.05,], file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK_sig.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)

library(stringr)
fg_TUM_vs_PNS_hallmark_dp <- fg_TUM_vs_PNS_hallmark[fg_TUM_vs_PNS_hallmark$padj < 0.01,]
fg_TUM_vs_PNS_hallmark_dp$pathway <- str_sub(fg_TUM_vs_PNS_hallmark_dp$pathway, start =10)
fg_TUM_vs_PNS_hallmark_dp$pathway <- str_replace_all(fg_TUM_vs_PNS_hallmark_dp$pathway, "_", " ")

pdf(file=paste(date_time, "_hallmark_dp_TUMvsPNS.pdf", sep=""),height=5,width=6)
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

png(file=paste(date_time, "_vol_TUMvsPNS.png", sep=""),res = 300, height=5, width=6.25, unit ="in")
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
legend("topleft", pch = c(16,16,16), col=c("green","blue","orange", "magenta", "chocolate4"), bg="white",
       legend = c("dna replication", "cell cycle checkpoint", "response to virus", "leukocyte migration", "cilium organization"), 
       title = expression(bold("Pathway")), cex=0.9)
dev.off()

vp_nfkb_noncanonical <- c("CD40", "RELB", "NFKB2", "TNFRSF1B")
vp_nfkb_canonical <- c("NFKBIA", "NFKB1", "TRAF2", "RELA", "TNFRSF1A")
vp_text <- c(vp_nfkb_canonical, vp_nfkb_noncanonical)

png(file=paste(date_time, "_vol_nfkb.png", sep=""),res=300,height=5,width=6,units="in")
with(vp, plot(log2FoldChange, -log10(padj), pch=20, main=" ", col=alpha("grey", 0.2), xlim=c(-2,3.5), ylim=c(0,6)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_nfkb_canonical,], points(log2FoldChange, -log10(padj), pch=20, cex=1.4, col=alpha("green", 1)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_nfkb_noncanonical,], points(log2FoldChange, -log10(padj), pch=20, cex=1.4, col=alpha("blue", 1)))
# with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_nfkb_atypical,], points(log2FoldChange, -log10(padj), pch=20, cex=1.4, col=alpha("chocolate4", 1)))
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


vp_text <- c("FGF1", "FGF2", "EGF", "IGF1", "IGF2", "VEGFA", "VEGFB", "PDGFA", "PDGFB")
png(file=paste(date_time, "_vol_rtk_ligands.pdf", sep=""),res=300, height=5,width=6, units = "in")
with(vp, plot(log2FoldChange, -log10(padj), pch=20, main=" ", col=alpha("grey", 0.2), xlim=c(-2,3.5), ylim=c(0,6)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_text,], points(log2FoldChange, -log10(padj), pch=20, cex=1.4, col=alpha("black", 1)))
abline(h=1.301, col="red", lty = 2, lwd = 2)
text(-1.6, 1.6, "p adj = 0.05", cex = .9)

vp_select <- vp[vp$symbol %in% vp_text,]
library(basicPlotteR)
addTextLabels(vp_select$log2FoldChange, -log10(vp_select$padj), vp_select$symbol, 
              cex.label = 0.8, 
              col.label="black", col.background = "white")
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

pdf(file=paste(date_time, "_gsea_dp_TUMvsPNS.pdf", sep=""),height=5,width=8.5)
library(ggplot2)
ggplot(fg_TUM_vs_PNS_go_dp) +
  geom_point(aes(NES, reorder(pathway, NES), size = gene_ratio, color=pval)) +
  scale_color_gradient(low="red", high = "blue") +
  ylab("GO Biological Process") + 
  xlab("Normalized Enrichment Score") +
  labs(color="p-adj", size = "gene ratio")
  # ggtitle("GSEA TUM vs PNS")
dev.off()


#   5.5 Compare Tumor VS Microenvironment (TUM_vs_TME) ---------------------------------------------
grp_header <- "TUM_vs_TME"
grp_A <- samples_TUM_UNQ
grp_B <- samples_TME_UNQ
count_temp <- count_ns[,c(grp_A, grp_B)]
anno_temp <- tdata_ns[c(grp_A,grp_B),]
identical(colnames(count_temp), rownames(anno_temp))
anno_temp$temp_group <- NA
anno_temp$temp_group[rownames(anno_temp) %in% grp_A] <- "A"
anno_temp$temp_group[rownames(anno_temp) %in% grp_B] <- "B"
anno_temp$temp_group <- factor(anno_temp$temp_group, levels=c("A","B"))
anno_temp$temp_group <- relevel(anno_temp$temp_group, ref="B") # B is the reference

res_TUM_vs_TME <- de_function(count_temp, anno_temp, grp_header) 
library(tidyverse)
ranks_TUM_vs_TME <- res_TUM_vs_TME %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))
ranks_TUM_vs_TME <- deframe(ranks_TUM_vs_TME)

# FGSEA Hallmark
fg_TUM_vs_TME_hallmark <- fgsea(pathways=path_hallmark, stats=ranks_TUM_vs_TME, nperm=1000)
fg_TUM_vs_TME_hallmark_tidy <- fg_TUM_vs_TME_hallmark %>%
  as_tibble() %>%
  arrange(desc(ES))
fwrite(fg_TUM_vs_TME_hallmark_tidy, file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
fwrite(fg_TUM_vs_TME_hallmark_tidy[fg_TUM_vs_TME_hallmark_tidy$padj < 0.05,], file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK_sig.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)

# FGSEA GO
library(fgsea)
fg_TUM_vs_TME_go <- fgsea_function(res_TUM_vs_TME, grp_header)


#   5.6 Compare Normal VS Microenvironment (PNS_vs_TME) ---------------------------------------------
grp_header <- "PNS_vs_TME"
grp_A <- samples_PNS
grp_B <- samples_TME_UNQ
count_temp <- count_ns[,c(grp_A, grp_B)]
anno_temp <- tdata_ns[c(grp_A,grp_B),]
identical(colnames(count_temp), rownames(anno_temp))
anno_temp$temp_group <- NA
anno_temp$temp_group[rownames(anno_temp) %in% grp_A] <- "A"
anno_temp$temp_group[rownames(anno_temp) %in% grp_B] <- "B"
anno_temp$temp_group <- factor(anno_temp$temp_group, levels=c("A","B"))
anno_temp$temp_group <- relevel(anno_temp$temp_group, ref="B") # B is the reference

res_PNS_vs_TME <- de_function(count_temp, anno_temp, grp_header) 
library(tidyverse)
ranks_PNS_vs_TME <- res_PNS_vs_TME %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))
ranks_PNS_vs_TME <- deframe(ranks_PNS_vs_TME)

# FGSEA Hallmark
fg_PNS_vs_TME_hallmark <- fgsea(pathways=path_hallmark, stats=ranks_PNS_vs_TME, nperm=1000)
fg_PNS_vs_TME_hallmark_tidy <- fg_PNS_vs_TME_hallmark %>%
  as_tibble() %>%
  arrange(desc(ES))
fwrite(fg_PNS_vs_TME_hallmark_tidy, file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
fwrite(fg_PNS_vs_TME_hallmark_tidy[fg_PNS_vs_TME_hallmark_tidy$padj < 0.05,], file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK_sig.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)


# FGSEA GO
library(fgsea)
fg_PNS_vs_TME_go <- fgsea_function(res_PNS_vs_TME, grp_header)



#   5.7 Comparing Cluster 1 vs Cluster 2 tumors (cluster1_vs_cluster2)  ---------------------------------------------
grp_header <- "cluster1_vs_cluster2"
grp_A <- samples_TUM_UNQ_cluster1
grp_B <- samples_TUM_UNQ_cluster2
count_temp <- count_ns[,c(grp_A, grp_B)]
anno_temp <- tdata_ns[c(grp_A,grp_B),]
identical(colnames(count_temp), rownames(anno_temp))
anno_temp$temp_group <- NA
anno_temp$temp_group[rownames(anno_temp) %in% grp_A] <- "A"
anno_temp$temp_group[rownames(anno_temp) %in% grp_B] <- "B"
anno_temp$temp_group <- factor(anno_temp$temp_group, levels=c("A","B"))
anno_temp$temp_group <- relevel(anno_temp$temp_group, ref="B") # B is the reference

res_CLUS1_vs_CLUS2 <- de_function(count_temp, anno_temp, grp_header) 
library(tidyverse)
ranks_CLUS1_vs_CLUS2 <- res_CLUS1_vs_CLUS2 %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))
ranks_CLUS1_vs_CLUS2 <- deframe(ranks_CLUS1_vs_CLUS2)

library(fgsea)
# FGSEA Hallmark
fg_CLUS1_vs_CLUS2_hallmark <- fgsea(pathways=path_hallmark, stats=ranks_CLUS1_vs_CLUS2, nperm=1000)
fg_CLUS1_vs_CLUS2_hallmark_tidy <- fg_CLUS1_vs_CLUS2_hallmark %>%
  as_tibble() %>%
  arrange(desc(ES))
fwrite(fg_CLUS1_vs_CLUS2_hallmark_tidy, file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)
fwrite(fg_CLUS1_vs_CLUS2_hallmark_tidy[fg_CLUS1_vs_CLUS2_hallmark_tidy$padj < 0.05,], file=paste(date_time, "_", grp_header, "_fgsea_HALLMARK_sig.txt", sep=""), 
       sep = "\t", sep2 = c("",",",""), row.names=FALSE,col.names=TRUE,verbose = TRUE)

# FGSEA GO
fg_CLUS1_vs_CLUS2_go <- fgsea_function(res_CLUS1_vs_CLUS2, grp_header)

# volcano genes
vp <- res_CLUS1_vs_CLUS2
fp <- fg_CLUS1_vs_CLUS2_go
# CLUS2
vp_humoral_response <- unlist(fp[fp$pathway=="GO_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_CIRCULATING_IMMUNOGLOBULIN","leadingEdge"][[1]])
vp_collagen_org <- unlist(fp[fp$pathway=="GO_COLLAGEN_FIBRIL_ORGANIZATION","leadingEdge"][[1]])
# CLUS1
vp_dna_pack <- unlist(fp[fp$pathway=="GO_DNA_PACKAGING","leadingEdge"][[1]])
vp_chromatin_assem <- unlist(fp[fp$pathway=="GO_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY","leadingEdge"][[1]])

vp_text <- c("HIST1H1C", "SMARCC2", "NCAPG2",
             "MMP11", "COL3A1", "IGKC", "IGLC2", "IGLL5") # "IFITM1"
png(file=paste(date_time, "_vol_CLUS1vs2.png", sep=""), res = 300,  height=5, width=6, unit = "in")
with(vp, plot(log2FoldChange, -log10(padj), pch=20, main=" ", col=alpha("black", 0.2), xlim=c(-8,8)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_collagen_org,], points(log2FoldChange, -log10(padj), pch=20, cex=1.2, col=alpha("orange", 0.4)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_humoral_response,], points(log2FoldChange, -log10(padj), pch=20, cex=1.2, col=alpha("yellow", 0.4)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_dna_pack,], points(log2FoldChange, -log10(padj), pch=20, cex=1.2, col=alpha("blue", 0.4)))
with(subset(vp, -log10(padj)>0 & abs(log2FoldChange)>0)[vp$symbol %in% vp_chromatin_assem,], points(log2FoldChange, -log10(padj), pch=20, cex=1.2, col=alpha("green", 0.4)))
vp_select <- vp[vp$symbol %in% vp_text,]
library(basicPlotteR)
addTextLabels(vp_select$log2FoldChange, -log10(vp_select$padj), vp_select$symbol, 
              cex.label = 0.8, 
              col.label="black", col.background = "white")
legend("topright", pch = 16, col=c("yellow","orange","blue", "green"), cex = 0.7,
       legend = c("collagen fibril organization", "humoral immune response (circulating Ig)",
                  "DNA packaging", "chromatin assembly / disassembly"), 
       title = expression(bold("Pathway")))
dev.off()

#heatmap differentially expressed genes

hm_genes <- res_CLUS1_vs_CLUS2$ENS_id[res_CLUS1_vs_CLUS2$padj < 0.05]
hm_genes <- hm_genes[!is.na(hm_genes)]
clus1vs2_hm <- tpm_tum_unq[hm_genes,]
breaksList = seq(-6, 6, by = 0.1)
library(pheatmap)
library(RColorBrewer)
png(file=paste("plots/", date_time, "_hm_unq_tum.png", sep=""), res = 300,  height=4, width=5, unit = "in")
y <- pheatmap(clus1vs2_hm,
              clustering_method = "ward",
              scale = "row",
              show_colnames = T, show_rownames = F, fontsize_col = 5,
              # labels_col = tdata_tum_unq$patient_index_num,
              col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(breaksList)),
              cluster_rows = T)
dev.off()

#   6.1 LMP-1 distribution -----------------------------------

# get tpm for tumor samples only
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


# distribution of LMP1
pdf(file=paste(date_time, "_dot_LMP1.pdf", sep=""),height=3,width=1.5)
ggplot(t_tpm_tum, aes(x = "", y = log2(t_tpm_tum$`LMP-1` + 1))) +
  geom_jitter(size = 0.7) +
  labs(title = "", x = "LMP1", y = "log2(tpm)", fill = "") +
  theme_classic() +
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=11)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = "crossbar", width = 0.5, color="red", size = 0.2)
dev.off()

table(t_tpm_tum$`LMP-1` == 0)  # 70 pos, 29 neg
median((t_tpm_tum$`LMP-1`))


#   6.2 Hallmark plots -----
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
             col=color_bars[fgp$color_highlight]
    ) +
    scale_fill_manual("padj < 0.05", labels=c("yes","no"), values=color_bars) +
    # labs(x="Pathway", y="Enrichment Score", title=paste("Hallmark ", fgp_head, sep="")) + 
    labs(x="Pathway", y="Enrichment Score") + 
        coord_flip() +
    theme_minimal() +
    theme(axis.text.y = element_text(color=color_axis[rev(fgp$color_highlight)]))
  return(fgplot)
}

pdf(file=paste(date_time, "_hallmark_NATvsPNS.pdf", sep=""),height=9,width=6)
plot_hallmark_fx(fg_NAT_vs_PNS_hallmark_tidy, "NAT vs PNS")
dev.off()

pdf(file=paste(date_time, "_hallmark_DYSvsNAT.pdf", sep=""),height=9,width=6)
plot_hallmark_fx(fg_DYS_vs_NAT_hallmark_tidy, "DYS vs NAT")
dev.off()

pdf(file=paste(date_time, "_hallmark_TUMvsDYS.pdf", sep=""),height=9,width=6)
plot_hallmark_fx(fg_TUM_vs_DYS_hallmark_tidy, "TUM vs DYS")
dev.off()

pdf(file=paste(date_time, "_hallmark_TUMvsPNS.pdf", sep=""),height=9,width=6)
plot_hallmark_fx(fg_TUM_vs_PNS_hallmark_tidy, "TUM vs PNS")
dev.off()



#   6.3 Heatmap combining Hallmark genes ------------------------------------
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
pdf(file=paste(date_time, "_hallmark_heatmap.pdf", sep=""),height=6,width=6.5)
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
  # gaps_col=c(26,500),
  # scale="row",
  col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(breaksList)),
  annotation_names_col = FALSE,
  breaks = breaksList,
  border_color="black"
  # main="2559 Hallmark Leading Genes (n = 125)"
)
hm_hallmark_genes
dev.off()

# 7. Panendoscopy -----------------------------------------------------

samples_PNS <- rownames(tdata_pan)[tdata_pan$cell_type == "PNS"]
samples_SQM <- rownames(tdata_pan)[tdata_pan$cell_type == "SQM"]
samples_SQMSQL <- rownames(tdata_pan)[tdata_pan$cell_type %in% c("SQM", "SQL")]

grp_header <- "PNS_vs_SQM"
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

pdf(file=paste(date_time, "_gsea_dp_PNSvsSQM.pdf", sep=""),height=5,width=7.8)
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

png(file=paste(date_time, "_vol_PNSvsSQM.png", sep=""),res=300, height=4.8,width=5.4, units = "in")

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



# 7.1 Panendoscopy boxplots --------------------------------------------------

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
pdf(file=paste(date_time, "_box_panendo.pdf", sep=""),height=2,width=3.2)
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


# 7.3 Panendoscopy deconvolution ------------------------------------------
cases_panendo <- c(samples_PNS, samples_SQM)
decon_xcell <- read.table("xCell_20191226_decon_mix_xCell_0211010120.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
decon_xcell <- t(decon_xcell)
rownames(decon_xcell) <- substring(rownames(decon_xcell), 2)

decon_xcell_panendo <- as.data.frame(t(decon_xcell[cases_panendo,]))
res_barchart <- rownames(decon_xcell_panendo)

pdf(file=paste(date_time, "_box_decon_panendo.pdf", sep=""),height=2,width=3.5)
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


# 8. Plot progression across spectrum-----------------------------------------------------

# Define patient samples
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

g_spry <- c("SPRY1", "SPRY2")
g_fgfr <- c("FGFR1","FGFR2","FGFR3")
g_egf <- c("EGF", "TGFA", "BTC", "AREG", "EREG", "EPGN", "IGF1")
g_egf <- c("HBEGF", "IGF2")
g_egfr <- c("EGFR", "ERBB2", "ERBB3", "IGF1R")
g_vegf <- c("VEGFA", "VEGFB")
g_vegfr <- c("FLT1", "KDR", "FLT4")
g_pdgf <- c("PDGFA", "PDGFB")
g_pdgfr <- c("PDGFRA", "PDGFRB")
g_gal <- c("LGALS1", "LGALS3")
g_nfkb1 <- c("RELA", "NFKB1", "TNFRSF1A", "TRAF2", "NFKBIA")
g_nfkb2 <- c("RELB", "NFKB2", "TNFRSF1B", "CD40", "LMP-1")
g_nfkbia <- "NFKBIA"
g_immune <- c("CCL20", "CX3CL1", "SAA1")
g_lmp1 <- "LMP-1"

res_barchart <- g_lmp1  # replace res_barchart with genes as needed
pdf(file=paste("plots/",date_time, "_box_LMP1.pdf", sep=""),height=2.5,width=8)
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


# 9. Deconvoluted samples -----------------------------------------------

# export tpm+prog matrix for decon
decon_mix <- tpm_fil
ens_temp <- rownames(decon_mix)
length(unique(ens_temp))
ens_temp[1:57820] <- gsub("\\..*","",ens_temp[1:57820])
length(unique(ens_temp))

ens_temp[ens_temp %in% names(hugo_name_lookup)] <- hugo_name_lookup[ens_temp[ens_temp %in% names(hugo_name_lookup)]]
length(ens_temp)
table(duplicated(ens_temp))
rownames(decon_mix) <- ens_temp

decon_mix <- cbind(rownames(decon_mix),decon_mix)
colnames(decon_mix)[1] <- "GeneSymbol"
write.table(decon_mix, file=paste(date_time, "_decon_mix.txt", sep=""), 
            sep = "\t", quote = FALSE, row.names=FALSE, col.names=TRUE)

tdata_tme <- tdata_ns[tdata_ns$cell_type == "TME",]
tdata_tme <- tdata_tme[duplicated(tdata_tme$patient_index) == FALSE,]
samples_TME_unique <- rownames(tdata_tme)
tpm_tme <- as.data.frame(t(tpm_prog[,samples_TME_unique]))


# read in deconvoluted profiles
decon_cib <- read.table("CIBERSORTx_Job4_Adjusted.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
decon_cib$Macrophage.total <- rowSums(cbind(decon_cib$Macrophages.M0, decon_cib$Macrophages.M1, decon_cib$Macrophages.M2))

decon_xcell <- read.table("xCell_20191226_decon_mix_xCell_0211010120.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
decon_xcell <- t(decon_xcell)
rownames(decon_xcell) <- substring(rownames(decon_xcell), 2)

# 9.1 Deconvoluted TUM vs PNS --------------------------

grp_header <- "decon_TUM_VS_PNS"
grp_A <- samples_TUM_UNQ
grp_B <- samples_PNS

decon_temp <- as.data.frame(decon_xcell[c(grp_A, grp_B),])
# View(decon_temp)
decon_temp$group_status <- NA
decon_temp[row.names(decon_temp) %in% grp_A,]$group_status <- "Tumor"
decon_temp[row.names(decon_temp) %in% grp_B,]$group_status <- "Normal"

library(ggplot2)
pdf(file=paste(date_time, "_xcell_TUM_vs_PNS.pdf", sep=""),height=6,width=6)
for (j in c(1:67))
{
  cell_type <- colnames(decon_temp[,j, drop=FALSE])
  print(
    ggplot(aes(y=decon_temp[,j], x=decon_temp$group_status), data=decon_temp) +
      geom_boxplot() +
      labs(title = paste("xcell", cell_type, " "), x = "Region", y = cell_type) # +
    # theme(plot.title = element_text(hjust = 0.5))  # to center title
  )
}
dev.off()

# t-test
ttest_header1 <- t(matrix(c(grp_header, "xCell")))
ttest_header2 <- t(matrix(c("Cell_type", "Group1_mean", "Group2_mean", "p-value", "Bonferroni_corrected_p-value")))
write.table(ttest_header1, file=paste(date_time, grp_header, "_xcell_ttest.txt", sep=""), append=F, sep="\t", row.names=F, col.names=F, quote=F)
write.table(ttest_header2, file=paste(date_time, grp_header, "_xcell_ttest.txt", sep=""), append=T, sep="\t", row.names=F, col.names=F, quote=F)

test_data <- decon_temp
for (i in 1:(ncol(test_data)-1))
{
  # print(colnames(test_data)[i])
  ttest_result <- t.test(test_data[,i][test_data$group_status == "Tumor"], test_data[,i][test_data$group_status == "Normal"], var.equal = TRUE, paired = FALSE)
  # print(ttest_result)
  ttest_list <- t(matrix(c(colnames(test_data)[i], ttest_result$estimate, ttest_result$p.value, ttest_result$p.value *65)))
  write.table(ttest_list, file=paste(date_time, grp_header, "_xcell_ttest.txt", sep=""), append=T, sep="\t", row.names=F, col.names=F, quote=F)
}


# cibersort
decon_temp <- as.data.frame(decon_cib[c(grp_A, grp_B),])
# View(decon_temp)
decon_temp$group_status <- NA
decon_temp[row.names(decon_temp) %in% grp_A,]$group_status <- "Tumor"
decon_temp[row.names(decon_temp) %in% grp_B,]$group_status <- "Normal"

library(ggplot2)
pdf(file=paste(date_time, "_cib_TUM_vs_PNS.pdf", sep=""),height=6,width=6)
for (j in c(1:22))
{
  cell_type <- colnames(decon_temp[,j, drop=FALSE])
  print(
    ggplot(aes(y=decon_temp[,j], x=decon_temp$group_status), data=decon_temp) +
      geom_boxplot() +
      labs(title = paste("cib", cell_type, " "), x = "Region", y = cell_type) # +
    # theme(plot.title = element_text(hjust = 0.5))  # to center title
  )
}
dev.off()

# t-test
ttest_header1 <- t(matrix(c(grp_header, "CIBERSORT")))
ttest_header2 <- t(matrix(c("Cell_type", "Group1_mean", "Group2_mean", "p-value", "Bonferroni_corrected_p-value")))
write.table(ttest_header1, file=paste(date_time, grp_header, "cib_ttest.txt", sep="_"), append=F, sep="\t", row.names=F, col.names=F, quote=F)
write.table(ttest_header2, file=paste(date_time, grp_header, "cib_ttest.txt", sep="_"), append=T, sep="\t", row.names=F, col.names=F, quote=F)

test_data <- decon_temp
for (i in 1:22)
{
  # print(colnames(test_data)[i])
  ttest_result <- t.test(test_data[,i][test_data$group_status == "Tumor"], test_data[,i][test_data$group_status == "Normal"], var.equal = TRUE, paired = FALSE)
  # print(ttest_result)
  ttest_list <- t(matrix(c(colnames(test_data)[i], ttest_result$estimate, ttest_result$p.value, ttest_result$p.value *22)))
  write.table(ttest_list, file=paste(date_time, grp_header, "cib_ttest.txt", sep="_"), append=T, sep="\t", row.names=F, col.names=F, quote=F)
}



# 9.2 Deconvoluted TUM vs DYS --------------------------

grp_header <- "decon_TUM_VS_DYS"
grp_A <- samples_TUM_UNQ
grp_B <- samples_DYS

decon_temp <- as.data.frame(decon_xcell[c(grp_A, grp_B),])
# View(decon_temp)
decon_temp$group_status <- NA
decon_temp[row.names(decon_temp) %in% grp_A,]$group_status <- "Tumor"
decon_temp[row.names(decon_temp) %in% grp_B,]$group_status <- "Dysplasia"

library(ggplot2)
pdf(file=paste(date_time, "_xcell_TUM_vs_DYS.pdf", sep=""),height=6,width=6)
for (j in c(1:67))
{
  cell_type <- colnames(decon_temp[,j, drop=FALSE])
  print(
    ggplot(aes(y=decon_temp[,j], x=decon_temp$group_status), data=decon_temp) +
      geom_boxplot() +
      labs(title = paste("xcell", cell_type, " "), x = "Region", y = cell_type) # +
    # theme(plot.title = element_text(hjust = 0.5))  # to center title
  )
}
dev.off()

# t-test
ttest_header1 <- t(matrix(c(grp_header, "xCell")))
ttest_header2 <- t(matrix(c("Cell_type", "Group1_mean", "Group2_mean", "p-value", "Bonferroni_corrected_p-value")))
write.table(ttest_header1, file=paste(date_time, grp_header, "_xcell_ttest.txt", sep=""), append=F, sep="\t", row.names=F, col.names=F, quote=F)
write.table(ttest_header2, file=paste(date_time, grp_header, "_xcell_ttest.txt", sep=""), append=T, sep="\t", row.names=F, col.names=F, quote=F)

test_data <- decon_temp
for (i in 1:(ncol(test_data)-1))
{
  # print(colnames(test_data)[i])
  ttest_result <- t.test(test_data[,i][test_data$group_status == "Tumor"], test_data[,i][test_data$group_status == "Dysplasia"], var.equal = TRUE, paired = FALSE)
  # print(ttest_result)
  ttest_list <- t(matrix(c(colnames(test_data)[i], ttest_result$estimate, ttest_result$p.value, ttest_result$p.value *65)))
  write.table(ttest_list, file=paste(date_time, grp_header, "_xcell_ttest.txt", sep=""), append=T, sep="\t", row.names=F, col.names=F, quote=F)
}


# cibersort
decon_temp <- as.data.frame(decon_cib[c(grp_A, grp_B),])
# View(decon_temp)
decon_temp$group_status <- NA
decon_temp[row.names(decon_temp) %in% grp_A,]$group_status <- "Tumor"
decon_temp[row.names(decon_temp) %in% grp_B,]$group_status <- "Dysplasia"

library(ggplot2)
pdf(file=paste(date_time, "_cib_TUM_vs_DYS.pdf", sep=""),height=6,width=6)
for (j in c(1:22))
{
  cell_type <- colnames(decon_temp[,j, drop=FALSE])
  print(
    ggplot(aes(y=decon_temp[,j], x=decon_temp$group_status), data=decon_temp) +
      geom_boxplot() +
      labs(title = paste("cib", cell_type, " "), x = "Region", y = cell_type) # +
    # theme(plot.title = element_text(hjust = 0.5))  # to center title
  )
}
dev.off()

# t-test
ttest_header1 <- t(matrix(c(grp_header, "CIBERSORT")))
ttest_header2 <- t(matrix(c("Cell_type", "Group1_mean", "Group2_mean", "p-value", "Bonferroni_corrected_p-value")))
write.table(ttest_header1, file=paste(date_time, grp_header, "cib_ttest.txt", sep="_"), append=F, sep="\t", row.names=F, col.names=F, quote=F)
write.table(ttest_header2, file=paste(date_time, grp_header, "cib_ttest.txt", sep="_"), append=T, sep="\t", row.names=F, col.names=F, quote=F)

test_data <- decon_temp
for (i in 1:22)
{
  # print(colnames(test_data)[i])
  ttest_result <- t.test(test_data[,i][test_data$group_status == "Tumor"], test_data[,i][test_data$group_status == "Dysplasia"], var.equal = TRUE, paired = FALSE)
  # print(ttest_result)
  ttest_list <- t(matrix(c(colnames(test_data)[i], ttest_result$estimate, ttest_result$p.value, ttest_result$p.value *22)))
  write.table(ttest_list, file=paste(date_time, grp_header, "cib_ttest.txt", sep="_"), append=T, sep="\t", row.names=F, col.names=F, quote=F)
}

# 9.3 Deconvoluted DYS vs PNS --------------------------

grp_header <- "decon_DYS_VS_PNS"
grp_A <- samples_DYS
grp_B <- samples_PNS

decon_temp <- as.data.frame(decon_xcell[c(grp_A, grp_B),])
# View(decon_temp)
decon_temp$group_status <- NA
decon_temp[row.names(decon_temp) %in% grp_A,]$group_status <- "Dysplasia"
decon_temp[row.names(decon_temp) %in% grp_B,]$group_status <- "Normal"

library(ggplot2)
pdf(file=paste(date_time, "_xcell_DYS_vs_PNS.pdf", sep=""),height=6,width=6)
for (j in c(1:67))
{
  cell_type <- colnames(decon_temp[,j, drop=FALSE])
  print(
    ggplot(aes(y=decon_temp[,j], x=decon_temp$group_status), data=decon_temp) +
      geom_boxplot() +
      labs(title = paste("xcell", cell_type, " "), x = "Region", y = cell_type) # +
    # theme(plot.title = element_text(hjust = 0.5))  # to center title
  )
}
dev.off()

# t-test
ttest_header1 <- t(matrix(c(grp_header, "xCell")))
ttest_header2 <- t(matrix(c("Cell_type", "Group1_mean", "Group2_mean", "p-value", "Bonferroni_corrected_p-value")))
write.table(ttest_header1, file=paste(date_time, grp_header, "_xcell_ttest.txt", sep=""), append=F, sep="\t", row.names=F, col.names=F, quote=F)
write.table(ttest_header2, file=paste(date_time, grp_header, "_xcell_ttest.txt", sep=""), append=T, sep="\t", row.names=F, col.names=F, quote=F)

test_data <- decon_temp
for (i in 1:(ncol(test_data)-1))
{
  # print(colnames(test_data)[i])
  ttest_result <- t.test(test_data[,i][test_data$group_status == "Dysplasia"], test_data[,i][test_data$group_status == "Normal"], var.equal = TRUE, paired = FALSE)
  # print(ttest_result)
  ttest_list <- t(matrix(c(colnames(test_data)[i], ttest_result$estimate, ttest_result$p.value, ttest_result$p.value *65)))
  write.table(ttest_list, file=paste(date_time, grp_header, "_xcell_ttest.txt", sep=""), append=T, sep="\t", row.names=F, col.names=F, quote=F)
}


# cibersort
decon_temp <- as.data.frame(decon_cib[c(grp_A, grp_B),])
# View(decon_temp)
decon_temp$group_status <- NA
decon_temp[row.names(decon_temp) %in% grp_A,]$group_status <- "Dysplasia"
decon_temp[row.names(decon_temp) %in% grp_B,]$group_status <- "Normal"

library(ggplot2)
pdf(file=paste(date_time, "_cib_DYS_vs_PNS.pdf", sep=""),height=6,width=6)
for (j in c(1:22))
{
  cell_type <- colnames(decon_temp[,j, drop=FALSE])
  print(
    ggplot(aes(y=decon_temp[,j], x=decon_temp$group_status), data=decon_temp) +
      geom_boxplot() +
      labs(title = paste("cib", cell_type, " "), x = "Region", y = cell_type) # +
    # theme(plot.title = element_text(hjust = 0.5))  # to center title
  )
}
dev.off()

# t-test
ttest_header1 <- t(matrix(c(grp_header, "CIBERSORT")))
ttest_header2 <- t(matrix(c("Cell_type", "Group1_mean", "Group2_mean", "p-value", "Bonferroni_corrected_p-value")))
write.table(ttest_header1, file=paste(date_time, grp_header, "cib_ttest.txt", sep="_"), append=F, sep="\t", row.names=F, col.names=F, quote=F)
write.table(ttest_header2, file=paste(date_time, grp_header, "cib_ttest.txt", sep="_"), append=T, sep="\t", row.names=F, col.names=F, quote=F)

test_data <- decon_temp
for (i in 1:22)
{
  # print(colnames(test_data)[i])
  ttest_result <- t.test(test_data[,i][test_data$group_status == "Dysplasia"], test_data[,i][test_data$group_status == "Normal"], var.equal = TRUE, paired = FALSE)
  # print(ttest_result)
  ttest_list <- t(matrix(c(colnames(test_data)[i], ttest_result$estimate, ttest_result$p.value, ttest_result$p.value *22)))
  write.table(ttest_list, file=paste(date_time, grp_header, "cib_ttest.txt", sep="_"), append=T, sep="\t", row.names=F, col.names=F, quote=F)
}






# 9.4 Deconvoluted TUM vs TME --------------------------

grp_header <- "decon_TUM_VS_TME"
grp_A <- samples_TUM_UNQ
grp_B <- samples_TME_UNQ

decon_temp <- as.data.frame(decon_xcell[c(grp_A, grp_B),])
# View(decon_temp)
decon_temp$group_status <- NA
decon_temp[row.names(decon_temp) %in% grp_A,]$group_status <- "Tumor"
decon_temp[row.names(decon_temp) %in% grp_B,]$group_status <- "Microenvironment"

library(ggplot2)
pdf(file=paste(date_time, "_xcell_TUM_vs_TME.pdf", sep=""),height=6,width=6)
for (j in c(1:67))
{
  cell_type <- colnames(decon_temp[,j, drop=FALSE])
  print(
    ggplot(aes(y=decon_temp[,j], x=decon_temp$group_status), data=decon_temp) +
      geom_boxplot() +
      labs(title = paste("xcell", cell_type, " "), x = "Region", y = cell_type) # +
    # theme(plot.title = element_text(hjust = 0.5))  # to center title
  )
}
dev.off()

# t-test
ttest_header1 <- t(matrix(c(grp_header, "xCell")))
ttest_header2 <- t(matrix(c("Cell_type", "Group1_mean", "Group2_mean", "p-value", "Bonferroni_corrected_p-value")))
write.table(ttest_header1, file=paste(date_time, grp_header, "_xcell_ttest.txt", sep=""), append=F, sep="\t", row.names=F, col.names=F, quote=F)
write.table(ttest_header2, file=paste(date_time, grp_header, "_xcell_ttest.txt", sep=""), append=T, sep="\t", row.names=F, col.names=F, quote=F)

test_data <- decon_temp
for (i in 1:(ncol(test_data)-1))
{
  # print(colnames(test_data)[i])
  ttest_result <- t.test(test_data[,i][test_data$group_status == "Tumor"], test_data[,i][test_data$group_status == "Microenvironment"], var.equal = TRUE, paired = FALSE)
  # print(ttest_result)
  ttest_list <- t(matrix(c(colnames(test_data)[i], ttest_result$estimate, ttest_result$p.value, ttest_result$p.value *65)))
  write.table(ttest_list, file=paste(date_time, grp_header, "_xcell_ttest.txt", sep=""), append=T, sep="\t", row.names=F, col.names=F, quote=F)
}


# cibersort
decon_temp <- as.data.frame(decon_cib[c(grp_A, grp_B),])
# View(decon_temp)
decon_temp$group_status <- NA
decon_temp[row.names(decon_temp) %in% grp_A,]$group_status <- "Tumor"
decon_temp[row.names(decon_temp) %in% grp_B,]$group_status <- "Microenvironment"

library(ggplot2)
pdf(file=paste(date_time, "_cib_TUM_vs_TME.pdf", sep=""),height=6,width=6)
for (j in c(1:22))
{
  cell_type <- colnames(decon_temp[,j, drop=FALSE])
  print(
    ggplot(aes(y=decon_temp[,j], x=decon_temp$group_status), data=decon_temp) +
      geom_boxplot() +
      labs(title = paste("cib", cell_type, " "), x = "Region", y = cell_type) # +
    # theme(plot.title = element_text(hjust = 0.5))  # to center title
  )
}
dev.off()

# t-test
ttest_header1 <- t(matrix(c(grp_header, "CIBERSORT")))
ttest_header2 <- t(matrix(c("Cell_type", "Group1_mean", "Group2_mean", "p-value", "Bonferroni_corrected_p-value")))
write.table(ttest_header1, file=paste(date_time, grp_header, "cib_ttest.txt", sep="_"), append=F, sep="\t", row.names=F, col.names=F, quote=F)
write.table(ttest_header2, file=paste(date_time, grp_header, "cib_ttest.txt", sep="_"), append=T, sep="\t", row.names=F, col.names=F, quote=F)

test_data <- decon_temp
for (i in 1:22)
{
  # print(colnames(test_data)[i])
  ttest_result <- t.test(test_data[,i][test_data$group_status == "Tumor"], test_data[,i][test_data$group_status == "Microenvironment"], var.equal = TRUE, paired = FALSE)
  # print(ttest_result)
  ttest_list <- t(matrix(c(colnames(test_data)[i], ttest_result$estimate, ttest_result$p.value, ttest_result$p.value *22)))
  write.table(ttest_list, file=paste(date_time, grp_header, "cib_ttest.txt", sep="_"), append=T, sep="\t", row.names=F, col.names=F, quote=F)
}




# 10. Deconvoluted progression --------------------------------------------

# xCell
cases_prog <- c(samples_PNS, samples_NAT, samples_DYS, samples_TUM, samples_TME)
decon_xcell_t <- t(decon_xcell)
decon_xcell_prog <- decon_xcell_t[,cases_prog]
decon_xcell_prog <- as.data.frame(t(decon_xcell[cases_prog,]))

pdf(file=paste(date_time, "_box_decon_progression.pdf", sep=""),height=2,width=3.5)
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
    theme_classic() +
    theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12))
  print(p)}
dev.off()

# cibersort
cases_prog <- c(samples_PNS, samples_NAT, samples_DYS, samples_TUM_UNQ, samples_TME)
decon_cib_t <- t(decon_cib)
decon_cib_prog <- decon_cib_t[,cases_prog]
decon_cib_prog <- as.data.frame(t(decon_cib[cases_prog,]))

pdf(file=paste(date_time, "_box_decon_progression_cibersort.pdf", sep=""),height=2,width=3.5)
res_barchart <- rownames(decon_cib_prog)
for (i in c(1:length(res_barchart))){
  library(tidyr)
  untidy <- as.data.frame(t(decon_cib_prog[res_barchart[i],]))
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
    theme_classic() +
    theme(axis.title.y=element_text(size=12), axis.text=element_text(size=12))
  print(p)}
dev.off()

# 11. TME correlates -----------------------------------------------

# remove any duplicates
tdata_paired_tme <- tdata_ns[samples_TME_UNQ,]
tdata_paired_tme_patient_index <- tdata_paired_tme$patient_index
unique(tdata_paired_tme_patient_index)

tdata_paired_tum <- tdata_tum_unq[tdata_tum_unq$cell_type == "TUM",] # from unique tumors
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



pdf(file=paste(date_time, "_corr_CXCLs_CCL20.pdf", sep=""),height=3.5,width=3.5)
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

 
pdf(file=paste(date_time, "_corr_CCL20_cib.pdf", sep=""),height=3.5,width=3.5)
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

pdf(file=paste(date_time, "_corr_CCL20_xcell.pdf", sep=""),height=3.5,width=3.5)
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


# 12. Deconvoluted Heatmap ----------------------------------------------------

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

tcdata_tme$T <- as.character(tcdata_tme$T)
tcdata_tme$N <- as.character(tcdata_tme$N)
tcdata_tme$M <- as.character(tcdata_tme$M) 
tcdata_tme$Stage <- recode(tcdata_tme$Stage, '2'="II", '3'="III", '4'="IV")

identical(colnames(decon_cib_tme_hm), rownames(tcdata_tme)) #TRUE

anno_cib_tme <- tcdata_tme[,c("Any recurrence", "Stage", "M", "N", "T", "Histology pattern", "LMP1 status"), drop=FALSE]
colorlist <- list('type' = c("tumor"="#EB6E80", "normal"="#008F95"), 
                  'T'=c("1"="#e6e6e6", "2"="#b3b3b3", "3"="#808080", "4"="#4d4d4d", "NA"="#ffffff"),
                  'N'=c("0"="#e6e6e6", "1"="#b3b3b3", "2"="#808080", "3"="#4d4d4d", "NA"="#ffffff"),
                  'M'=c("0"="#e6e6e6", "1"="#4d4d4d", "NA"="#ffffff"),
                  'Stage'=c("I"="#e6e6e6", "II"="#b3b3b3", "III"="#808080", "IV"="#4d4d4d", "NA"="#ffffff"),
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


# OO. Save --------------------------------------------------------------------
setwd(output_directory)
save.image("NPC_s3s_final.Rdata")

