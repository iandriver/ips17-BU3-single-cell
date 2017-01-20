#!/usr/bin/env Rscript
#R code for monocle2 ips17 BU3 single cell data Hawkins, et al.
#also refer to the monocle2 vignette http://www.bioconductor.org/packages/release/bioc/vignettes/monocle/inst/doc/monocle-vignette.pdf

library(monocle)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(BiocGenerics)
library(Hmisc)
library(igraph)

#if using modified monocle plotting options to produce figures from the paper
#provide path to modified plotting code https://github.com/iandriver/monocle-release/tree/master/R

#provide name add on for file saveing
fname="kmeans2_200_test"
new_dirname = file.path(getwd(), fname)
dir.create(new_dirname, showWarnings = TRUE)

#read in gene cell matrix file
fpkm_matrix <- read.delim("ips17_BU3_normalized_edgeR_counts_TMM.txt")

#assign gene symbol names as row.names
row.names(fpkm_matrix) <- fpkm_matrix$X

#rename the matrix as final
final_df <- fpkm_matrix
#filter out any blanks
final_df <- final_df[,!(names(final_df) %in% c('X'))]
fpkm_matrix <- final_df

#read in the Cell names and group assignments
sample_sheet <- read.delim("NKX_and_mito_and_origin_cells.txt")
#make the rownames the SampleIDs
row.names(sample_sheet) <- sample_sheet$SampleID
#create the pheno-data DataFrame for monocle2
pd <- new("AnnotatedDataFrame", data = sample_sheet)

#match the SampleIDs in the sample sheet to the matrix
fpkm2 <- fpkm_matrix[,match(row.names(sample_sheet), colnames(fpkm_matrix))]


#create the (gene) feature sheet by taking all of the genes in the matrix and making a dataframe (no extra data needed)
feature_sheet <- as.data.frame(row.names(fpkm2))
row.names(feature_sheet) <- feature_sheet$"row.names(fpkm2)"
feature_sheet$GeneID <- row.names(feature_sheet)
fd <- new("AnnotatedDataFrame", data = feature_sheet)

#create the monocle2 CellDataSet object
ips17_BU3_data <- newCellDataSet(as.matrix(fpkm2), phenoData = pd, featureData=fd)

# As per monocle2 vignette, use it to estimate RNA counts
rpc_matrix <- relative2abs(ips17_BU3_data)
# Now, make a new CellDataSet using the RNA counts
ips17_BU3_data <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
phenoData = pd,
featureData = fd,
lowerDetectionLimit=1,
expressionFamily=negbinomial.size())

#pre-calculate size factors and dispersions
ips17_BU3_data <- estimateSizeFactors(ips17_BU3_data)
ips17_BU3_data <- estimateDispersions(ips17_BU3_data)

#calulate the number of cells that each gene is expressed in
ips17_BU3_data <- detectGenes(ips17_BU3_data , min_expr = 0.1)
#create new feature data for the total sum of gene expression across all cells
fData(ips17_BU3_data)$GeneSum <- Matrix::rowSums(exprs(ips17_BU3_data))

#find expressed genes that are expressed in at least 6 cells and have a total expression >= 80
expressed_genes <- row.names(subset(fData(ips17_BU3_data), num_cells_expressed >= 6 & GeneSum >=80))

#filter cells with way to few or way to many genes expressed
pData(ips17_BU3_data)$Total_mRNAs <- Matrix::colSums(exprs(ips17_BU3_data))
upper_bound <- 10^(mean(log10(pData(ips17_BU3_data)$Total_mRNAs)) + 2*sd(log10(pData(ips17_BU3_data)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(ips17_BU3_data)$Total_mRNAs)) - 2*sd(log10(pData(ips17_BU3_data)$Total_mRNAs)))

#plot mRNA expression by Nkx_group or other
qplot(Total_mRNAs, data=pData(ips17_BU3_data), color=Nkx_group, geom="density") + geom_vline(xintercept=lower_bound) + geom_vline(xintercept=upper_bound)

ips17_BU3_data <- ips17_BU3_data[,pData(ips17_BU3_data)$Total_mRNAs > lower_bound & pData(ips17_BU3_data)$Total_mRNAs < upper_bound]
ips17_BU3_data <- detectGenes(ips17_BU3_data, min_expr = 0.1)

#create gene ids for classification
NKX_id <- row.names(subset(fData(ips17_BU3_data), GeneID == "NKX2-1"))
APOA2_id <- row.names(subset(fData(ips17_BU3_data), GeneID == "APOA2"))
CD47_id <- row.names(subset(fData(ips17_BU3_data), GeneID == "CD47"))
MSX1_id <- row.names(subset(fData(ips17_BU3_data), GeneID == "MSX1"))
SFTA3_id <- row.names(subset(fData(ips17_BU3_data), GeneID == "SFTA3"))
COL19A1_id <- row.names(subset(fData(ips17_BU3_data), GeneID == "COL19A1"))
#create empty newCellTypeHierarchy classifier
cth <- newCellTypeHierarchy()

#define "Lung" classifier expression rules
#expresses NKX2-1 and SFTA3 or
#expresses CD47 and SFTA3 or
#expresses NKX2-1 and CD47

cth <- addCellType(cth, "Lung", classify_func=function(x) {
              (x[NKX_id,] >= 1 & x[SFTA3_id,] >=1) |
              (x[CD47_id,] >= 1 & x[SFTA3_id,]>=1) |
              (x[NKX_id,] >= 1 & x[CD47_id,] >=1)
              })

#define "Liver" classifier expression rules
#doesnt express NKX2-1 and does express APOA2 or
#doesnt express NKX2-1 and does express MSX1
#doesnt express NKX2-1 and does express COL19A1 or
#doesnt express NKX2-1 or CD47
cth <- addCellType(cth, "Liver", classify_func=function(x) {
              (x[NKX_id,] < 1 & x[APOA2_id,] > 1) |
              (x[NKX_id,] < 1 & x[MSX1_id,] > 1) |
              (x[NKX_id,] < 1 &x[COL19A1_id,]>1) |
              (x[NKX_id,] < 1 & x[CD47_id,] < 1)
              })


#classify cells into Lung or Liver, with rest being ambiguous or unknown
ips17_BU3_data <- classifyCells(ips17_BU3_data, cth, 0.1)

#plot distribution of cell classifcations
pie <- ggplot(pData(ips17_BU3_data), aes(x = factor(1), fill = factor(CellType))) +
geom_bar(width = 1)
pie + coord_polar(theta = "y") +
theme(axis.title.x=element_blank(), axis.title.y=element_blank())

#significance test for genes that vary by the assigned cell type
diff_test_res <- differentialGeneTest(ips17_BU3_data[expressed_genes,], fullModelFormulaStr="~CellType", cores=3)
#take the top significant gene
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
#read in the file that contains cell cycle genes to exclude
cell_cycle_genes <- read.delim('cell_cycle_genes.txt')

#read in file with top significant genes bewtween the 2 kmeans cluster groups (scicast produced)
best_gene_list <- read.delim('kmeans_2_top200_allgene.txt')
best_genes <- best_gene_list$GeneID

#merge the two lists and remove duplicates
best_ordering <- as.factor(c(as.character(best_genes)))

#exclude cell cycle genes
best_ordering <- setdiff(best_ordering, cell_cycle_genes$GeneID)
best_ord <- best_ordering[!duplicated(best_ordering)]

best_ord <- best_ord[1:150]
#create cell trajectory based on the list of ordering genes
ips17_BU3_data <- setOrderingFilter(ips17_BU3_data, best_ord)
ips17_BU3_data <- reduceDimension(ips17_BU3_data, max_components=3)
ips17_BU3_data <- orderCells(ips17_BU3_data, reverse=TRUE)
source('plotting_larger_text.R')
#plot trajectories
plot_cell_trajectory(ips17_BU3_data, color_by="State", cell_size = 6, line_width = 1.8, show_branch_points=FALSE)
ggsave(file.path(new_dirname,paste(fname,"cell_traj_byState.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 8, units = c("in"), dpi = 300)
plot_cell_trajectory(ips17_BU3_data, color_by="CellType", cell_size = 6, line_width = 1.8, show_branch_points=FALSE)
ggsave(file.path(new_dirname,paste(fname,"cell_traj_byCellType.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 8, units = c("in"), dpi = 300)
plot_cell_trajectory(ips17_BU3_data, color_by="Nkx_group", cell_size = 6, line_width = 1.8, show_branch_points=FALSE)
ggsave(file.path(new_dirname,paste(fname,"cell_traj_byNkxGroup.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 8, units = c("in"), dpi = 300)
plot_cell_trajectory(ips17_BU3_data, color_by="Mitosis_Group", cell_size = 6, line_width = 1.8, show_branch_points=FALSE)
ggsave(file.path(new_dirname,paste(fname,"cell_traj_byMitosisGroup.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 8, units = c("in"), dpi = 300)
plot_cell_trajectory(ips17_BU3_data, color_by="CellOrigin", cell_size = 6, line_width = 1.8, show_branch_points=FALSE)
ggsave(file.path(new_dirname,paste(fname,"cell_traj_byCellOrigin.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 8, units = c("in"), dpi = 300)
traj <- plot_cell_trajectory(ips17_BU3_data, color_by="Pseudotime", cell_size = 6, line_width = 1.8, show_branch_points=FALSE)
traj + scale_colour_gradient(limits=c(0, 12), low="purple", high="green")
ggsave(file.path(new_dirname,paste(fname,"cell_traj_byPseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 8, units = c("in"), dpi = 300)

#plot genes in pseudotime
plot_genes_in_pseudotime(ips17_BU3_data['NKX2-1',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("NKX2-1",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['CD47',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("CD47",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['APOA2',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("APOA2",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['FGB',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("FGB",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['SOX9',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("SOX9",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['LAMA3',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("LAMA3",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['COL1A2',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("COL1A2",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['HAS2',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("HAS2",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['CPM',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("CPM",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['IGFBP5',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("IGFBP5",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['NFIB',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("NFIB",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['WNT5A',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("WNT5A",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['IRX2',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("IRX2",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)
plot_genes_in_pseudotime(ips17_BU3_data['SFTA3',], color_by = "State", cell_size = 6, relative_expr = F, trendline_width=1.8)
ggsave(file.path(new_dirname,paste("SFTA3",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height = 6, units = c("in"), dpi = 300)

#find top pseudtime genes
diff_test_res_ps <- differentialGeneTest(ips17_BU3_data[expressed_genes,], fullModelFormulaStr="~sm.ns(Pseudotime)", cores=3)
pstime_genes <- row.names (subset(diff_test_res_ps, qval < 0.001))
pst_df <- subset(fData(ips17_BU3_data[pstime_genes,]))
#save significant pseudtime genes
write.table(pst_df[order(-pst_df[,"GeneSum"]),], file=paste(new_dirname,'/',fname,'_',"significant_genes_by_pseudotime_ordered_by_expression.txt", sep=''), sep="\t", col.names=NA)
pst_df <- pst_df[order(-pst_df[,"GeneSum"]),]
diff_test_df_ps <- diff_test_res_ps[,c("GeneID", "pval", "qval")]
pst_df_qval <- diff_test_df_ps[match(row.names(pst_df), row.names(diff_test_df_ps)),]
pst_df_qval$GeneSum <- pst_df$GeneSum
pst_df_qval$num_cells_expressed <- pst_df$num_cells_expressed
top_30_expression_to_plot <- row.names(pst_df_qval[1:min(dim(pst_df_qval)[1],30),])
source('plotting_small_text.R')
plot_genes_in_pseudotime(ips17_BU3_data[top_30_expression_to_plot,], color_by = "State", nrow = 6, ncol = 5, cell_size = 2.4, relative_expr = F)
ggsave(file.path(new_dirname,paste("top30_byexpression_Pseudo_sig",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 14, height = 11, units = c("in"), dpi = 300)
pst_df_qval <- pst_df_qval[order(pst_df_qval[,"qval"]),]
top_30_qval_to_plot <- row.names(pst_df_qval[1:min(dim(pst_df_qval)[1],30),])
plot_genes_in_pseudotime(ips17_BU3_data[top_30_qval_to_plot,], color_by = "State", nrow = 6, ncol = 5, cell_size = 2.4, relative_expr = F)
ggsave(file.path(new_dirname,paste("top30_byqval_Pseudo_sig",fname,"Pseudotime.pdf", sep="_")), plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 14, height = 11, units = c("in"), dpi = 300)
#plot genes in jitter
plot_genes_jitter(ips17_BU3_data['CD47',], grouping="State", color_by="State", nrow=1, ncol=NULL, plot_trend = TRUE)
