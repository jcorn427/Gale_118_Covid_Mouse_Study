######### --Libraries--#########
setwd("/vol08/ngs/Gale_Lab_Projects/Gale118_CoV2MA10_C57BL6/Cornelius_analysis/de_work/DE_all_9_23")
library(stringr)
library(limma)
library(edgeR)
library(factoextra)
library(PCAtools)
library(biomaRt)
library(svglite)
library(CEMiTool)
library(qusage)
library(data.table)
library(WebGestaltR)
source("/vol08/ngs/Gale_Lab_Projects/Gale118_CoV2MA10_C57BL6/Cornelius_analysis/heatmap3LW_function.r")

# Set for heatmap top column color
colcolorlistgrodwn <- c(rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"), rep("white"))
colcolormatrix <- as.matrix(colcolorlistgrodwn)


# --Read in target files
message("STATUS: Load tables")
cm <- read.table("./count_matrix.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
target <- read.csv("./targetfile.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

# Rename samples
newsampleIDs <- c()
for (i in colnames(cm)) {
  i <- str_remove(i, "_RNA\\d+_Lib\\d+\\S*$")
  i <- str_replace_all(i, "-", "_")
  newsampleIDs <- c(newsampleIDs, i)
}
colnames(cm) <- newsampleIDs

# --generate figure of all counts
png("de_intensities_raw_counts.png", res = 100)
par(xpd = TRUE)
if (length(rownames(target)) > 10) {
  plotDensities(log2(cm + 0.1), legend = FALSE)
} else {
  plotDensities(log2(cm + 0.1),
                legend = "topright",
                inset = c(-0.2, 0), levels(rownames(target))
  )
}
dev.off()


# -- Normalize data
# order target and count matrix so they are the same (THIS IS IMPORTANT)
cm <- cm[, rownames(target)]

# CHECK IF ORDER IS THE SAME
if (all.equal(colnames(cm), rownames(target)) != TRUE) {
  print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
  print(rownames(target))
  print(colnames(cm))
}

# normalize
cm2 <- DGEList(counts = cm)
cm2 <- calcNormFactors(cm2, method = "TMM") # TMM normalization
png("mean_variance_norm.png")
Pi.CPM <- voom(counts = cm2, normalize.method = "none", plot = T, span = 0.1)
dev.off()
write.csv(Pi.CPM$E, "Pi.CPM$E_all.csv")

# box plot of unfiltered data
png("boxplot_vnorm_all.png", width = 10, height = 8, units = "in", res = 100)
# par(mar=c(1,1,1,1))
minvalue <- min(Pi.CPM$E)
maxvalue <- max(Pi.CPM$E)
boxplot(Pi.CPM$E,
        labels = target$GaleID, ylim = c(minvalue - 1, maxvalue + 1),
        ylab = "voom expression", main = "Count matrix", cex.axis = .6, las = 2,
        frame = FALSE
)
dev.off()

#-- filter out genes from each group that are below mean count of 10 across samples 
#Iteratively adjusted thresholds and decided that 10 was the best cutoff to get rid of the 
# mean-variance hook shown in the initial mean-variance plot (iterative results not saved just ran in R)
A <- rowMeans(cm)
isexpr <- A >= 10
cmfl_counts <- cm[isexpr, ]
write.csv(cmfl_counts, "count_matrix_renamed_fl.csv")


#Normalize again with filtering
cm2 <- DGEList(counts = cmfl_counts)
cm2 <- calcNormFactors(cm2, method = "TMM") # TMM normalization
png("mean_variance_norm_fl.png")
Pi.CPM <- voom(counts = cm2, normalize.method = "none", plot = T, span = 0.1)
dev.off()
write.csv(Pi.CPM$E, "Pi.CPM$E_all_fl.csv")

# Save voom normalized object to resume analysis
saveRDS(Pi.CPM, file = "Pi.CPM.rds")

# box plot of filtered data
png("boxplot_vnorm_all_fl.png", width = 10, height = 8, units = "in", res = 100)
# par(mar=c(1,1,1,1))
minvalue <- min(Pi.CPM$E)
maxvalue <- max(Pi.CPM$E)
boxplot(Pi.CPM$E,
        labels = target$GaleID, ylim = c(minvalue - 1, maxvalue + 1),
        ylab = "voom normalized expression", main = "Normalized count matrix", cex.axis = .6, las = 2,
        frame = FALSE
)
dev.off()

#Feature reduction
p <- PCAtools::pca(Pi.CPM$E, 
                   metadata = target)

PCAtools::biplot(p, x='PC1', y='PC2', 
                 lab=NULL, 
                 colby='Time_Point',
                 shape='Sample_Type',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 24,
                 legendLabSize = 20,
                 legendTitleSize = 24,
                 showLoadings = FALSE,
                 legendIconSize = 8,
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)
dev.off()

screeplot(p, xlab = "Principal component", title = "SCREE plot", ylim = c(0,100), components = getComponents(p, 1:10),
          colBar = "dodgerblue")

#Subsetting PCA by tissue
heart <- Pi.CPM$E[,grepl("_H", colnames(Pi.CPM$E))]
lung <- Pi.CPM$E[,grepl("_L", colnames(Pi.CPM$E))]
kidney <- Pi.CPM$E[,grepl("_K", colnames(Pi.CPM$E))]
spleen <- Pi.CPM$E[,grepl("_S", colnames(Pi.CPM$E))]

h_t <- target[grepl("_H", rownames(target)),]
l_t <- target[grepl("_L", rownames(target)),]
k_t <- target[grepl("_K", rownames(target)),]
s_t <- target[grepl("_S", rownames(target)),]

p_h <- PCAtools::pca(heart, 
                     metadata = h_t)
p_l <- PCAtools::pca(lung, 
                     metadata = l_t)
p_k <- PCAtools::pca(kidney, 
                     metadata = k_t)
p_s <- PCAtools::pca(spleen, 
                     metadata = s_t)

PCAtools::biplot(p_h, x='PC1', y='PC2', 
                 lab=NULL, 
                 colby='Time_Point',
                 shape='Sex',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 24,
                 legendLabSize = 20,
                 legendTitleSize = 24,
                 showLoadings = FALSE,
                 legendIconSize = 8,
                 title = 'Heart',
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)

PCAtools::biplot(p_l, x='PC1', y='PC2', 
                 lab=NULL, 
                 colby='Time_Point',
                 shape='Sex',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 24,
                 legendLabSize = 20,
                 legendTitleSize = 24,
                 showLoadings = FALSE,
                 legendIconSize = 8,
                 title = "lung",
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)

PCAtools::biplot(p_k, x='PC1', y='PC2', 
                 lab=NULL, 
                 colby='Time_Point',
                 shape='Sex',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 24,
                 legendLabSize = 20,
                 legendTitleSize = 24,
                 showLoadings = FALSE,
                 legendIconSize = 8,
                 title = 'kidney',
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)

PCAtools::biplot(p_s, x='PC1', y='PC2', 
                 lab=NULL, 
                 colby='Time_Point',
                 shape='Sex',
                 pointSize = 3,
                 legendPosition = 'right',
                 #hline=0,
                 #vline=0,
                 axisLabSize = 24,
                 legendLabSize = 20,
                 legendTitleSize = 24,
                 showLoadings = FALSE,
                 legendIconSize = 8,
                 title = 'Spleen',
                 gridlines.major = FALSE,
                 gridlines.minor = FALSE)

# Get mouse gene list
mouse <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
mouse_BM <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = mouse)
# mouse_BM <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "transcript_Hength"), mart = mouse)
## remove duplicate gene entries ##
mouse_BM <- mouse_BM[!duplicated(mouse_BM$ensembl_gene_id),]


all(row.names(Pi.CPM$E) %in% mouse_BM$ensembl_gene_id)
genomic_idx <- match(rownames(Pi.CPM$E), mouse_BM$ensembl_gene_id)
mouse_BM_ordered <- mouse_BM[genomic_idx,]
#mouse_BM_ordered <- as.matrix(mouse_BM_ordered)
rows <- nrow(mouse_BM_ordered)
rownames(mouse_BM_ordered) <- 1:rows
all(row.names(Pi.CPM$E) %in% mouse_BM_ordered$ensembl_gene_id)


## Generate model matrix and lmfit object for sex blind analysis ##

mn <- factor(target[, "Mouse_Number"])
age <- factor(target[, "Age"])
hd <- factor(target[, "Harvest_DPI"])
st <- factor(target[, "Sample_Type"])
mm <- model.matrix(~ 0 + age:hd:st)

rownames(mm) <- colnames(Pi.CPM$E)
colnames(mm) <- make.names(colnames(mm))
mm <- mm[, colnames(mm)[order(tolower(colnames(mm[, ])))]]
mm <- mm[, colSums(mm) > 0]

excludeAll <- nonEstimable(mm)
if (length(excludeAll) > 0) {
  message("WARNING: These samples are nonEstimatable, design matrix ", excludeAll)
}

if ("ti" %in% excludeAll) {
  return("interactions term non estimable")
}
mm <- mm[, !colnames(mm) %in% excludeAll]
if (!is.fullrank(mm)) {
  return("not full rank")
}


Pi.lmfit <- lmFit(Pi.CPM, design = mm, block = target$Animal_ID, correlation = dupcor$consensus)



######### --de analysis 10 wk lung, sex blind--#########

contrastsmatrix <- c(
  "age10.week.hd2.stL -age10.week.hdM.stL",
  "age10.week.hd4.stL -age10.week.hdM.stL",
  "age10.week.hd7.stL -age10.week.hdM.stL"
  # "age20.week.hd2.stL -age20.week.hdM.stL",
  # "age20.week.hd4.stL -age20.week.hdM.stL",
  # "age20.week.hd7.stL -age20.week.hdM.stL",
  # "age2.year.hd2.stL -age2.year.hdM.stL",
  # "age2.year.hd4.stL -age2.year.hdM.stL",
  # "age2.year.hd7.stL -age2.year.hdM.stL",
  # "age10.week.hd2.stK -age10.week.hdM.stK",
  # "age10.week.hd4.stK -age10.week.hdM.stK",
  # "age10.week.hd7.stK -age10.week.hdM.stK",
  # "age20.week.hd2.stK -age20.week.hdM.stK",
  # "age20.week.hd4.stK -age20.week.hdM.stK",
  # "age20.week.hd7.stK -age20.week.hdM.stK",
  # "age2.year.hd2.stK -age2.year.hdM.stK",
  # "age2.year.hd4.stK -age2.year.hdM.stK",
  # "age2.year.hd7.stK -age2.year.hdM.stK",
  # "age10.week.hd2.stS -age10.week.hdM.stS",
  # "age10.week.hd4.stS -age10.week.hdM.stS",
  # "age10.week.hd7.stS -age10.week.hdM.stS",
  # "age20.week.hd2.stS -age20.week.hdM.stS",
  # "age20.week.hd4.stS -age20.week.hdM.stS",
  # "age20.week.hd7.stS -age20.week.hdM.stS",
  # "age2.year.hd2.stS -age2.year.hdM.stS",
  # "age2.year.hd4.stS -age2.year.hdM.stS",
  # "age2.year.hd7.stS -age2.year.hdM.stS",
  # "age10.week.hd2.stH -age10.week.hdM.stH",
  # "age10.week.hd4.stH -age10.week.hdM.stH",
  # "age10.week.hd7.stH -age10.week.hdM.stH",
  # "age20.week.hd2.stH -age20.week.hdM.stH",
  # "age20.week.hd4.stH -age20.week.hdM.stH",
  # "age20.week.hd7.stH -age20.week.hdM.stH",
  # "age2.year.hd2.stH -age2.year.hdM.stH",
  # "age2.year.hd4.stH -age2.year.hdM.stH",
  # "age2.year.hd7.stH -age2.year.hdM.stH"
)
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- treat(fit)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.01)
summary(results)

topTreat_data <- topTreat(fit2, coef=1, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_lung_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=2, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_lung_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=3, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_lung_hd7.csv", quote = F)


dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

# filter for significant genes - up/down regulated
# sigMask <- subset(sigMask, rowSums(sigMask) != 0)
# fit2$genes <- data.frame(ID_REF=rownames(fit2))

write.csv(ExpressMatrixde, file = "expression_matrix_de.csv", quote = F)
write.csv(results, file = "results_de.csv", quote = F)
write.csv(dataMatrixde, file = "full_expression_matrix_de.csv", quote = F)
write.csv(fit2$t, file = "t_stats.csv", quote = F)
write.csv(fit2$p.value, file = "p_value.csv", quote = F)
for (i in 1:ncol(fit2$p.value)) fit2$p.value[, i] <- p.adjust(fit2$p.value[, i], method = "BH")
write.csv(fit2$p.value, file = "p_value_adj.csv", quote = F)

ExpressMatrixde_genes <- merge(mouse_BM_ordered, ExpressMatrixde,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)

ExpressMatrixde_genes <- ExpressMatrixde_genes[, !(names(ExpressMatrixde_genes) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes <- avereps(ExpressMatrixde_genes,
                                 ID = ExpressMatrixde_genes$external_gene_name
)
rownames(ExpressMatrixde_genes) <- ExpressMatrixde_genes[, "external_gene_name"]
ExpressMatrixde_genes <- ExpressMatrixde_genes[, !(colnames(ExpressMatrixde_genes) %in% c("external_gene_name"))]
ExpressMatrixde_genes <- as.matrix(data.frame(ExpressMatrixde_genes, check.names = FALSE))
class(ExpressMatrixde_genes) <- "numeric"
write.csv(ExpressMatrixde_genes, file = "expression_matrix_10wk_lung_extgenename.csv", quote = F)






dataMatrixde_genes <- merge(mouse_BM_ordered, dataMatrixde,
                            by.x = "ensembl_gene_id",
                            by.y = "row.names",
                            all.X = T, all.Y = T
)

dataMatrixde_genes <- dataMatrixde_genes[, !(names(dataMatrixde_genes) %in% c("ensembl_gene_id"))]
dataMatrixde_genes <- avereps(dataMatrixde_genes,
                              ID = dataMatrixde_genes$external_gene_name
)
rownames(dataMatrixde_genes) <- dataMatrixde_genes[, "external_gene_name"]
dataMatrixde_genes <- dataMatrixde_genes[, !(colnames(dataMatrixde_genes) %in% c("external_gene_name"))]
dataMatrixde_genes <- as.matrix(data.frame(dataMatrixde_genes, check.names = FALSE))
class(dataMatrixde_genes) <- "numeric"
write.csv(dataMatrixde_genes, file = "full_expression_matrix_10wk_lung_extgenename.csv", quote = F)

message(paste0("Dimensionality of DE genes ", dim(ExpressMatrixde)[1]))

new_colnames <- c(
  "10wk_lung_hd2", "10wk_lung_hd4", "10wk_lung_hd7"
)

colnames(results) <- new_colnames
colnames(dataMatrixde) <- new_colnames
colnames(dataMatrixde_genes) <- new_colnames
colnames(ExpressMatrixde) <- new_colnames
colnames(ExpressMatrixde_genes) <- new_colnames

# --heatmap
svglite("heatmap10w_lung.svg", width = 10, height = 10)
global_modules10wklung <- heatmap.L.4(ExpressMatrixde_genes,
                                      figmargins = c(20, 5),
                                      cutoff = 1, distmethod = "euclidean", cexcol = 2, colcolorlist = colcolormatrix,
                                      clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

# Barplot
results_t <- t(summary(results))
results_t <- results_t[, -2]

for (i in 1:(length(row.names(results_t)))) {
  results_t[i, 1] <- results_t[i, 1] * -1
}

DE <- as.data.frame(results_t)
DE <- setnames(DE,
               old = c("Var1", "Var2", "Freq"),
               new = c("Time_Point", "group", "DE_genes")
)

# Create plot
ggplot(DE, aes(
  x = Time_Point, y = DE_genes, fill = group,
  label = DE$DE_genes
)) +
  geom_bar(stat = "identity", position = "identity") +
  # geom_text(size = 5, position = position_stack(vjust = 0) )+
  # theme_light() +
  theme_minimal() +
  scale_fill_manual(values = c("#0808c4", "#da9618")) +
  # xlab("Time point")
  ylab("Number of Differentially Expressed Genes") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 15)
  )
ggsave(barplot10w_lung.svg, width = 6, height = 4, units="in", dpi = 300)

#######SAVE THESE FOR LATER ######
#Go back and re-run just 10wk and get sig de genes, then save ##
dataMatrixde_10wk <- dataMatrixde
dataMatrixde_genes_10wk <- dataMatrixde_genes
ExpressMatrixde_10wk <- ExpressMatrixde
ExpressMatrixde_10wk_genes <- ExpressMatrixde_genes

###WebGestaltR analysis####
for (cluster in unique(global_modules10wklung$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modules10wklung$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism="mmusculus",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(cm),
              referenceGeneType = "ensembl_gene_id",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(global_modules10wklung$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modules10wklung$modulesrows == cluster)
  write.csv(names(genes), paste0(cluster,"_genes_10wk.csv"))
  
  temp <- ExpressMatrixde_genes[row.names(ExpressMatrixde_genes)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_10wk.csv"))
}

## All lung, sex blind ##
contrastsmatrix <- c(
  "age10.week.hd2.stL -age10.week.hdM.stL",
  "age10.week.hd4.stL -age10.week.hdM.stL",
  "age10.week.hd7.stL -age10.week.hdM.stL",
  "age20.week.hd2.stL -age20.week.hdM.stL",
  "age20.week.hd4.stL -age20.week.hdM.stL",
  "age20.week.hd7.stL -age20.week.hdM.stL",
  "age2.year.hd2.stL -age2.year.hdM.stL",
  "age2.year.hd4.stL -age2.year.hdM.stL",
  "age2.year.hd7.stL -age2.year.hdM.stL"
  # "age10.week.hd2.stK -age10.week.hdM.stK",
  # "age10.week.hd4.stK -age10.week.hdM.stK",
  # "age10.week.hd7.stK -age10.week.hdM.stK",
  # "age20.week.hd2.stK -age20.week.hdM.stK",
  # "age20.week.hd4.stK -age20.week.hdM.stK",
  # "age20.week.hd7.stK -age20.week.hdM.stK",
  # "age2.year.hd2.stK -age2.year.hdM.stK",
  # "age2.year.hd4.stK -age2.year.hdM.stK",
  # "age2.year.hd7.stK -age2.year.hdM.stK",
  # "age10.week.hd2.stS -age10.week.hdM.stS",
  # "age10.week.hd4.stS -age10.week.hdM.stS",
  # "age10.week.hd7.stS -age10.week.hdM.stS",
  # "age20.week.hd2.stS -age20.week.hdM.stS",
  # "age20.week.hd4.stS -age20.week.hdM.stS",
  # "age20.week.hd7.stS -age20.week.hdM.stS",
  # "age2.year.hd2.stS -age2.year.hdM.stS",
  # "age2.year.hd4.stS -age2.year.hdM.stS",
  # "age2.year.hd7.stS -age2.year.hdM.stS",
  # "age10.week.hd2.stH -age10.week.hdM.stH",
  # "age10.week.hd4.stH -age10.week.hdM.stH",
  # "age10.week.hd7.stH -age10.week.hdM.stH",
  # "age20.week.hd2.stH -age20.week.hdM.stH",
  # "age20.week.hd4.stH -age20.week.hdM.stH",
  # "age20.week.hd7.stH -age20.week.hdM.stH",
  # "age2.year.hd2.stH -age2.year.hdM.stH",
  # "age2.year.hd4.stH -age2.year.hdM.stH",
  # "age2.year.hd7.stH -age2.year.hdM.stH"
)
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- treat(fit)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.01)
summary(results)

topTreat_data <- topTreat(fit2, coef=1, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_lung_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=2, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_lung_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=3, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_lung_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=4, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_20wk_lung_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=5, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_20wk_lung_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=6, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_20wk_lung_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=7, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_2yr_lung_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=8, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_2yr_lung_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=9, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_2yr_lung_hd7.csv", quote = F)


dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

# filter for significant genes - up/down regulated
# sigMask <- subset(sigMask, rowSums(sigMask) != 0)
# fit2$genes <- data.frame(ID_REF=rownames(fit2))

write.csv(ExpressMatrixde, file = "expression_matrix_de.csv", quote = F)
write.csv(results, file = "results_de.csv", quote = F)
write.csv(dataMatrixde, file = "full_expression_matrix_de.csv", quote = F)
write.csv(fit2$t, file = "t_stats.csv", quote = F)
write.csv(fit2$p.value, file = "p_value.csv", quote = F)
for (i in 1:ncol(fit2$p.value)) fit2$p.value[, i] <- p.adjust(fit2$p.value[, i], method = "BH")
write.csv(fit2$p.value, file = "p_value_adj.csv", quote = F)

ExpressMatrixde_genes <- merge(mouse_BM_ordered, ExpressMatrixde,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)

ExpressMatrixde_genes <- ExpressMatrixde_genes[, !(names(ExpressMatrixde_genes) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes <- avereps(ExpressMatrixde_genes,
                                 ID = ExpressMatrixde_genes$external_gene_name
)
rownames(ExpressMatrixde_genes) <- ExpressMatrixde_genes[, "external_gene_name"]
ExpressMatrixde_genes <- ExpressMatrixde_genes[, !(colnames(ExpressMatrixde_genes) %in% c("external_gene_name"))]
ExpressMatrixde_genes <- as.matrix(data.frame(ExpressMatrixde_genes, check.names = FALSE))
class(ExpressMatrixde_genes) <- "numeric"
write.csv(ExpressMatrixde_genes, file = "expression_matrix_all_lung_extgenename.csv", quote = F)






dataMatrixde_genes <- merge(mouse_BM_ordered, dataMatrixde,
                            by.x = "ensembl_gene_id",
                            by.y = "row.names",
                            all.X = T, all.Y = T
)

dataMatrixde_genes <- dataMatrixde_genes[, !(names(dataMatrixde_genes) %in% c("ensembl_gene_id"))]
dataMatrixde_genes <- avereps(dataMatrixde_genes,
                              ID = dataMatrixde_genes$external_gene_name
)
rownames(dataMatrixde_genes) <- dataMatrixde_genes[, "external_gene_name"]
dataMatrixde_genes <- dataMatrixde_genes[, !(colnames(dataMatrixde_genes) %in% c("external_gene_name"))]
dataMatrixde_genes <- as.matrix(data.frame(dataMatrixde_genes, check.names = FALSE))
class(dataMatrixde_genes) <- "numeric"
write.csv(dataMatrixde_genes, file = "full_expression_matrix_all_lung_extgenename.csv", quote = F)

message(paste0("Dimensionality of DE genes ", dim(ExpressMatrixde)[1]))

new_colnames <- c(
  "10wk_lung_hd2", "10wk_lung_hd4", "10wk_lung_hd7", "20wk_lung_hd2", "20wk_lung_hd4", "20wk_lung_hd7",
  "2yr_lung_hd2", "2yr_lung_hd4", "2yr_lung_hd7"
)

colnames(results) <- new_colnames
colnames(dataMatrixde) <- new_colnames
colnames(dataMatrixde_genes) <- new_colnames
colnames(ExpressMatrixde) <- new_colnames
colnames(ExpressMatrixde_genes) <- new_colnames

#######SAVE THESE FOR LATER ######
dataMatrixde_all_lung <- dataMatrixde
dataMatrixde_genes_all_lung <- dataMatrixde_genes
ExpressMatrixde_all_lung <- ExpressMatrixde
ExpressMatrixde_all_genes_lung <- ExpressMatrixde_genes

# --heatmap
svglite("heatmap_all_lung.svg", width = 10, height = 10)
global_modulesalllung <- heatmap.L.4(ExpressMatrixde_all_genes_lung,
                                      figmargins = c(20, 5),
                                      cutoff = 1, distmethod = "euclidean", cexcol = 2, colcolorlist = colcolormatrix,
                                      clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

# Barplot
results_t <- t(summary(results))
results_t <- results_t[, -2]

for (i in 1:(length(row.names(results_t)))) {
  results_t[i, 1] <- results_t[i, 1] * -1
}

DE <- as.data.frame(results_t)
DE <- setnames(DE,
               old = c("Var1", "Var2", "Freq"),
               new = c("Time_Point", "group", "DE_genes")
)

# Create plot
ggplot(DE, aes(
  x = Time_Point, y = DE_genes, fill = group,
  label = DE$DE_genes
)) +
  geom_bar(stat = "identity", position = "identity") +
  # geom_text(size = 5, position = position_stack(vjust = 0) )+
  # theme_light() +
  theme_minimal() +
  scale_fill_manual(values = c("#0808c4", "#da9618")) +
  # xlab("Time point")
  ylab("Number of Differentially Expressed Genes") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 15)
  )
ggsave("barplot_all_lung.svg", width = 6, height = 4, units="in", dpi = 300)


###WebGestaltR analysis####
for (cluster in unique(global_modulesalllung$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesalllung$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism="mmusculus",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(cm),
              referenceGeneType = "ensembl_gene_id",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(global_modulesalllung$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesalllung$modulesrows == cluster)
  write.csv(names(genes), paste0(cluster,"_genes_10wk.csv"))
  
  temp <- ExpressMatrixde_genes[row.names(ExpressMatrixde_genes)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_all_lung.csv"))
}







## Lung DDE ##

contrastsmatrix <- c(
  "(age20.week.hd2.stL -age20.week.hdM.stL) - (age10.week.hd2.stL -age10.week.hdM.stL)",
  "(age20.week.hd4.stL -age20.week.hdM.stL) - (age10.week.hd4.stL -age10.week.hdM.stL)",
  "(age20.week.hd7.stL -age20.week.hdM.stL) - (age10.week.hd7.stL -age10.week.hdM.stL)",
  "(age2.year.hd2.stL -age2.year.hdM.stL) - (age10.week.hd2.stL -age10.week.hdM.stL)",
  "(age2.year.hd4.stL -age2.year.hdM.stL) - (age10.week.hd4.stL -age10.week.hdM.stL)",
  "(age2.year.hd7.stL -age2.year.hdM.stL) - (age10.week.hd7.stL -age10.week.hdM.stL)",
  "(age2.year.hd2.stL -age2.year.hdM.stL) - (age20.week.hd2.stL -age20.week.hdM.stL)",
  "(age2.year.hd4.stL -age2.year.hdM.stL) - (age20.week.hd4.stL -age20.week.hdM.stL)",
  "(age2.year.hd7.stL -age2.year.hdM.stL) - (age20.week.hd7.stL -age20.week.hdM.stL)"
)
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- treat(fit)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.01)
summary(results)

#Generate topTreat tables for each comparison#
topTreat_data <- topTreat(fit2, coef=1, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_1 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_1, file = "lfc_p_dde_lung_20wk-10wk_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=2, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_2 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_2, file = "lfc_p_dde_lung_2y-10wk_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=3, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_3 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_3, file = "lfc_p_dde_lung_20wk-10wk_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=4, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_4 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_4, file = "lfc_p_dde_lung_2y-10wk_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=5, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_5 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_5, file = "lfc_p_dde_lung_20wk-10wk_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=6, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_6 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_6, file = "lfc_p_dde_lung_2y-10wk_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=7, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_7 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_7, file = "lfc_p_dde_lung_2y-20wk_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=8, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_8 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_8, file = "lfc_p_dde_lung_2y-20wk_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=9, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_9 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_9, file = "lfc_p_dde_lung_2y-20wk_hd7.csv", quote = F)




dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise
## Add p-values to dataMatrixde objects ##
for (i in 1:ncol(fit2$p.value)) fit2$p.value[, i] <- p.adjust(fit2$p.value[, i], method = "BH")
pvalsadjst <- fit2$p.value

merged_pval_matrix_all <- merge(dataMatrixde, pvalsadjst,
                                by.x = "row.names",
                                by.y = "row.names",
                                all.X = T, all.Y = T
)
## make first column the row names ##
merged_pval_matrix_all <- data.frame(merged_pval_matrix_all, row.names = 1)

## Use no pvalue expression matrix to filter sig genes in matrix with pvalue##
ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

ExpressMatrixde_merged_pval_all <- merged_pval_matrix_all[row.names(merged_pval_matrix_all)%in%row.names(ExpressMatrixde),]

# filter for significant genes - up/down regulated
sigMask <- subset(sigMask, rowSums(sigMask) != 0)
# fit2$genes <- data.frame(ID_REF=rownames(fit2))

write.csv(ExpressMatrixde_merged_pval_all, file = "expression_matrix_dde_pval_all.csv", quote = F)
write.csv(ExpressMatrixde, file = "expression_matrix_dde.csv", quote = F)
write.csv(results, file = "results_dde.csv", quote = F)
write.csv(dataMatrixde, file = "full_expression_matrix_dde.csv", quote = F)

## Add gene names to pval expression matrix ##
ExpressMatrixde_genes_merged_pval_all <- merge(mouse_BM_ordered, ExpressMatrixde_merged_pval_all,
                                               by.x = "ensembl_gene_id",
                                               by.y = "row.names",
                                               all.X = T, all.Y = T
)

## Average entries with same gene name then make gene names the row names ##
ExpressMatrixde_genes_merged_pval_all <- ExpressMatrixde_genes_merged_pval_all[, !(names(ExpressMatrixde_genes_merged_pval_all) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes_merged_pval_all <- avereps(ExpressMatrixde_genes_merged_pval_all,
                                                 ID = ExpressMatrixde_genes_merged_pval_all$external_gene_name
)
rownames(ExpressMatrixde_genes_merged_pval_all) <- ExpressMatrixde_genes_merged_pval_all[, "external_gene_name"]
ExpressMatrixde_genes_merged_pval_all <- ExpressMatrixde_genes_merged_pval_all[, !(colnames(ExpressMatrixde_genes_merged_pval_all) %in% c("external_gene_name"))]
ExpressMatrixde_genes_merged_pval_all <- as.matrix(data.frame(ExpressMatrixde_genes_merged_pval_all, check.names = FALSE))
class(ExpressMatrixde_genes_merged_pval_all) <- "numeric"
write.csv(ExpressMatrixde_genes_merged_pval_all, "expression_matrix_dde_lung_extgenename_pval.csv", quote = F)

## Add gene names to non-pval expression matrix  and remove duplicates##
ExpressMatrixde_genes_all <- merge(mouse_BM_ordered, ExpressMatrixde,
                                   by.x = "ensembl_gene_id",
                                   by.y = "row.names",
                                   all.X = T, all.Y = T
)

ExpressMatrixde_genes_all <- ExpressMatrixde_genes_all[, !(names(ExpressMatrixde_genes_all) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes_all <- avereps(ExpressMatrixde_genes_all,
                                     ID = ExpressMatrixde_genes_all$external_gene_name
)
rownames(ExpressMatrixde_genes_all) <- ExpressMatrixde_genes_all[, "external_gene_name"]
ExpressMatrixde_genes_all <- ExpressMatrixde_genes_all[, !(colnames(ExpressMatrixde_genes_all) %in% c("external_gene_name"))]
ExpressMatrixde_genes_all <- as.matrix(data.frame(ExpressMatrixde_genes_all, check.names = FALSE))
class(ExpressMatrixde_genes_all) <- "numeric"
write.csv(ExpressMatrixde_genes_all, "expression_matrix_genes_dde.csv", quote = F)


## Rename datamatrix and remove duplicates ##
dataMatrixde_genes <- merge(mouse_BM_ordered, dataMatrixde,
                            by.x = "ensembl_gene_id",
                            by.y = "row.names",
                            all.X = T, all.Y = T
)

dataMatrixde_genes <- dataMatrixde_genes[, !(names(dataMatrixde_genes) %in% c("ensembl_gene_id"))]
dataMatrixde_genes <- avereps(dataMatrixde_genes,
                              ID = dataMatrixde_genes$external_gene_name
)
rownames(dataMatrixde_genes) <- dataMatrixde_genes[, "external_gene_name"]
dataMatrixde_genes <- dataMatrixde_genes[, !(colnames(dataMatrixde_genes) %in% c("external_gene_name"))]
dataMatrixde_genes <- as.matrix(data.frame(dataMatrixde_genes, check.names = FALSE))
class(dataMatrixde_genes) <- "numeric"
write.csv(dataMatrixde_genes, "full_expression_matrix_dde_lung_extgenename.csv", quote = F)

message(paste0("Dimensionality of DDE genes ", dim(ExpressMatrixde)[1]))

new_colnames <- c("20wk_hd2-10wk_hd2_lung","20wk_hd4-10wk_hd4_lung", "20wk_hd7-10wk_hd7_lung", "2y_hd2-10wk_hd2_lung", "2y_hd4-10wk_hd4_lung", "2y_hd7-10wk_hd7_lung",
                  "2y_hd2-20wk_hd2_lung", "2y_hd4-20wk_hd4_lung", "2y_hd7-20wk_hd7_lung"
)

new_colnames_pval <- c("20wk_hd2-10wk_hd2_lung_lfc", "20wk_hd4-10wk_hd4_lung_lfc", "20wk_hd7-10wk_hd7_lung_lfc", "2y_hd2-10wk_hd2_lung_lfc", "2y_hd4-10wk_hd4_lung_lfc", 
                       "2y_hd7-10wk_hd7_lung_lfc", "2y_hd2-20wk_hd2_lung_lfc", "2y_hd4-20wk_hd4_lung_lfc", 
                       "2y_hd7-20wk_hd7_lung_lfc","20wk_hd2-10wk_hd2_lung_cpval", "20wk_hd4-10wk_hd4_lung_cpval", "20wk_hd7-10wk_hd7_lung_cpval", "2y_hd2-10wk_hd2_lung_cpval", 
                       "2y_hd4-10wk_hd4_lung_cpval", "2y_hd7-10wk_hd7_lung_cpval", "2y_hd2-20wk_hd2_lung_cpval", 
                       "2y_hd4-20wk_hd4_lung_cpval", "2y_hd7-20wk_hd7_lung_cpval"
)


colnames(results) <- new_colnames
colnames(dataMatrixde) <- new_colnames
colnames(dataMatrixde_genes) <- new_colnames
colnames(ExpressMatrixde) <- new_colnames
colnames(ExpressMatrixde_genes_all) <- new_colnames
colnames(ExpressMatrixde_genes_merged_pval_all) <- new_colnames_pval

# --heatmap
svglite("heatmap_all_lung_dde.svg", width = 10, height = 10)
global_modulesalllungdde <- heatmap.L.4(ExpressMatrixde_genes_all,
                                     figmargins = c(20, 5),
                                     cutoff = 1, distmethod = "euclidean", cexcol = 2, colcolorlist = colcolormatrix,
                                     clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

# Barplot
results_t <- t(summary(results))
results_t <- results_t[, -2]

for (i in 1:(length(row.names(results_t)))) {
  results_t[i, 1] <- results_t[i, 1] * -1
}

DE <- as.data.frame(results_t)
DE <- setnames(DE,
               old = c("Var1", "Var2", "Freq"),
               new = c("Time_Point", "group", "DE_genes")
)

# Create plot
ggplot(DE, aes(
  x = Time_Point, y = DE_genes, fill = group,
  label = DE$DE_genes
)) +
  geom_bar(stat = "identity", position = "identity") +
  # geom_text(size = 5, position = position_stack(vjust = 0) )+
  # theme_light() +
  theme_minimal() +
  scale_fill_manual(values = c("#0808c4", "#da9618")) +
  # xlab("Time point")
  ylab("Number of Differentially Expressed Genes") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 15)
  )
ggsave("barplot_all_lung_dde.svg", width = 6, height = 4, units="in", dpi = 300)


###WebGestaltR analysis####
for (cluster in unique(global_modulesalllungdde$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesalllungdde$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism="mmusculus",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(cm),
              referenceGeneType = "ensembl_gene_id",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(global_modulesalllungdde$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesalllungdde$modulesrows == cluster)
  write.csv(names(genes), paste0(cluster,"_genes_10wk.csv"))
  
  temp <- ExpressMatrixde_genes[row.names(ExpressMatrixde_genes)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_all_lung.csv"))
}

# Read in GMT file for CEMiTool #
gmt <- read_gmt("./m5.all.v2023.1.Mm.symbols.gmt")

# Read in sample annotation for CEMiTool
annot <- read.csv("CEMiTool_sample_annot_lung.csv")

## All lung CEMiTool ##
dataMatrixde_genes_all_lung <- as.data.frame(dataMatrixde_genes_all_lung)
cem <- cemitool(dataMatrixde_genes_all_lung, annot, gmt)
generate_report(cem)
diagnostic_report(cem)






# All kidney DE #
contrastsmatrix <- c(
  # "age10.week.hd2.stL -age10.week.hdM.stL",
  # "age10.week.hd4.stL -age10.week.hdM.stL",
  # "age10.week.hd7.stL -age10.week.hdM.stL",
  # "age20.week.hd2.stL -age20.week.hdM.stL",
  # "age20.week.hd4.stL -age20.week.hdM.stL",
  # "age20.week.hd7.stL -age20.week.hdM.stL",
  # "age2.year.hd2.stL -age2.year.hdM.stL",
  # "age2.year.hd4.stL -age2.year.hdM.stL",
  # "age2.year.hd7.stL -age2.year.hdM.stL"
  "age10.week.hd2.stK -age10.week.hdM.stK",
  "age10.week.hd4.stK -age10.week.hdM.stK",
  "age10.week.hd7.stK -age10.week.hdM.stK",
  "age20.week.hd2.stK -age20.week.hdM.stK",
  "age20.week.hd4.stK -age20.week.hdM.stK",
  "age20.week.hd7.stK -age20.week.hdM.stK",
  "age2.year.hd2.stK -age2.year.hdM.stK",
  "age2.year.hd4.stK -age2.year.hdM.stK",
  "age2.year.hd7.stK -age2.year.hdM.stK"
  # "age10.week.hd2.stS -age10.week.hdM.stS",
  # "age10.week.hd4.stS -age10.week.hdM.stS",
  # "age10.week.hd7.stS -age10.week.hdM.stS",
  # "age20.week.hd2.stS -age20.week.hdM.stS",
  # "age20.week.hd4.stS -age20.week.hdM.stS",
  # "age20.week.hd7.stS -age20.week.hdM.stS",
  # "age2.year.hd2.stS -age2.year.hdM.stS",
  # "age2.year.hd4.stS -age2.year.hdM.stS",
  # "age2.year.hd7.stS -age2.year.hdM.stS",
  # "age10.week.hd2.stH -age10.week.hdM.stH",
  # "age10.week.hd4.stH -age10.week.hdM.stH",
  # "age10.week.hd7.stH -age10.week.hdM.stH",
  # "age20.week.hd2.stH -age20.week.hdM.stH",
  # "age20.week.hd4.stH -age20.week.hdM.stH",
  # "age20.week.hd7.stH -age20.week.hdM.stH",
  # "age2.year.hd2.stH -age2.year.hdM.stH",
  # "age2.year.hd4.stH -age2.year.hdM.stH",
  # "age2.year.hd7.stH -age2.year.hdM.stH"
)
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- treat(fit)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.01)
summary(results)

topTreat_data <- topTreat(fit2, coef=1, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_kidney_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=2, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_kidney_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=3, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_kidney_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=4, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_20wk_kidney_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=5, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_20wk_kidney_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=6, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_20wk_kidney_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=7, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_2yr_kidney_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=8, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_2yr_kidney_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=9, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_2yr_kidney_hd7.csv", quote = F)


dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

# filter for significant genes - up/down regulated
# sigMask <- subset(sigMask, rowSums(sigMask) != 0)
# fit2$genes <- data.frame(ID_REF=rownames(fit2))

write.csv(ExpressMatrixde, file = "expression_matrix_de.csv", quote = F)
write.csv(results, file = "results_de.csv", quote = F)
write.csv(dataMatrixde, file = "full_expression_matrix_de.csv", quote = F)
write.csv(fit2$t, file = "t_stats.csv", quote = F)
write.csv(fit2$p.value, file = "p_value.csv", quote = F)
for (i in 1:ncol(fit2$p.value)) fit2$p.value[, i] <- p.adjust(fit2$p.value[, i], method = "BH")
write.csv(fit2$p.value, file = "p_value_adj.csv", quote = F)

ExpressMatrixde_genes <- merge(mouse_BM_ordered, ExpressMatrixde,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)

ExpressMatrixde_genes <- ExpressMatrixde_genes[, !(names(ExpressMatrixde_genes) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes <- avereps(ExpressMatrixde_genes,
                                 ID = ExpressMatrixde_genes$external_gene_name
)
rownames(ExpressMatrixde_genes) <- ExpressMatrixde_genes[, "external_gene_name"]
ExpressMatrixde_genes <- ExpressMatrixde_genes[, !(colnames(ExpressMatrixde_genes) %in% c("external_gene_name"))]
ExpressMatrixde_genes <- as.matrix(data.frame(ExpressMatrixde_genes, check.names = FALSE))
class(ExpressMatrixde_genes) <- "numeric"
write.csv(ExpressMatrixde_genes, file = "expression_matrix_all_kidney_extgenename.csv", quote = F)






dataMatrixde_genes <- merge(mouse_BM_ordered, dataMatrixde,
                            by.x = "ensembl_gene_id",
                            by.y = "row.names",
                            all.X = T, all.Y = T
)

dataMatrixde_genes <- dataMatrixde_genes[, !(names(dataMatrixde_genes) %in% c("ensembl_gene_id"))]
dataMatrixde_genes <- avereps(dataMatrixde_genes,
                              ID = dataMatrixde_genes$external_gene_name
)
rownames(dataMatrixde_genes) <- dataMatrixde_genes[, "external_gene_name"]
dataMatrixde_genes <- dataMatrixde_genes[, !(colnames(dataMatrixde_genes) %in% c("external_gene_name"))]
dataMatrixde_genes <- as.matrix(data.frame(dataMatrixde_genes, check.names = FALSE))
class(dataMatrixde_genes) <- "numeric"
write.csv(dataMatrixde_genes, file = "full_expression_matrix_all_kidney_extgenename.csv", quote = F)

message(paste0("Dimensionality of DE genes ", dim(ExpressMatrixde)[1]))

new_colnames <- c(
  "10wk_kidney_hd2", "10wk_kidney_hd4", "10wk_kidney_hd7", "20wk_kidney_hd2", "20wk_kidney_hd4", "20wk_kidney_hd7",
  "2yr_kidney_hd2", "2yr_kidney_hd4", "2yr_kidney_hd7"
)

colnames(results) <- new_colnames
colnames(dataMatrixde) <- new_colnames
colnames(dataMatrixde_genes) <- new_colnames
colnames(ExpressMatrixde) <- new_colnames
colnames(ExpressMatrixde_genes) <- new_colnames

#######SAVE THESE FOR LATER ######
dataMatrixde_all_kidney <- dataMatrixde
dataMatrixde_genes_all_kidney <- dataMatrixde_genes
ExpressMatrixde_all_kidney <- ExpressMatrixde
ExpressMatrixde_all_genes_kidney <- ExpressMatrixde_genes

# --heatmap
svglite("heatmap_all_kidney.svg", width = 10, height = 10)
global_modulesallkidney <- heatmap.L.4(ExpressMatrixde_genes,
                                     figmargins = c(20, 5),
                                     cutoff = 1, distmethod = "euclidean", cexcol = 2, colcolorlist = colcolormatrix,
                                     clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

# Barplot
results_t <- t(summary(results))
results_t <- results_t[, -2]

for (i in 1:(length(row.names(results_t)))) {
  results_t[i, 1] <- results_t[i, 1] * -1
}

DE <- as.data.frame(results_t)
DE <- setnames(DE,
               old = c("Var1", "Var2", "Freq"),
               new = c("Time_Point", "group", "DE_genes")
)

# Create plot
ggplot(DE, aes(
  x = Time_Point, y = DE_genes, fill = group,
  label = DE$DE_genes
)) +
  geom_bar(stat = "identity", position = "identity") +
  # geom_text(size = 5, position = position_stack(vjust = 0) )+
  # theme_light() +
  theme_minimal() +
  scale_fill_manual(values = c("#0808c4", "#da9618")) +
  # xlab("Time point")
  ylab("Number of Differentially Expressed Genes") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 15)
  )
ggsave("barplot_all_kidney.svg", width = 6, height = 4, units="in", dpi = 300)


###WebGestaltR analysis####
for (cluster in unique(global_modulesallkidney$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesallkidney$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism="mmusculus",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(cm),
              referenceGeneType = "ensembl_gene_id",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(global_modulesallkidney$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesallkidney$modulesrows == cluster)
  write.csv(names(genes), paste0(cluster,"_genes_all_kidney.csv"))
  
  temp <- ExpressMatrixde_genes[row.names(ExpressMatrixde_genes)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_all_kidney.csv"))
}




## Kidney DDE ##

contrastsmatrix <- c(
  "(age20.week.hd2.stK -age20.week.hdM.stK) - (age10.week.hd2.stK -age10.week.hdM.stK)",
  "(age20.week.hd4.stK -age20.week.hdM.stK) - (age10.week.hd4.stK -age10.week.hdM.stK)",
  "(age20.week.hd7.stK -age20.week.hdM.stK) - (age10.week.hd7.stK -age10.week.hdM.stK)",
  "(age2.year.hd2.stK -age2.year.hdM.stK) - (age10.week.hd2.stK -age10.week.hdM.stK)",
  "(age2.year.hd4.stK -age2.year.hdM.stK) - (age10.week.hd4.stK -age10.week.hdM.stK)",
  "(age2.year.hd7.stK -age2.year.hdM.stK) - (age10.week.hd7.stK -age10.week.hdM.stK)",
  "(age2.year.hd2.stK -age2.year.hdM.stK) - (age20.week.hd2.stK -age20.week.hdM.stK)",
  "(age2.year.hd4.stK -age2.year.hdM.stK) - (age20.week.hd4.stK -age20.week.hdM.stK)",
  "(age2.year.hd7.stK -age2.year.hdM.stK) - (age20.week.hd7.stK -age20.week.hdM.stK)"
)
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- treat(fit)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.01)
summary(results)

#Generate topTreat tables for each comparison#
topTreat_data <- topTreat(fit2, coef=1, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_1 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_1, file = "lfc_p_dde_kidney_20wk-10wk_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=2, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_2 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_2, file = "lfc_p_dde_kidney_2y-10wk_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=3, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_3 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_3, file = "lfc_p_dde_kidney_20wk-10wk_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=4, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_4 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_4, file = "lfc_p_dde_kidney_2y-10wk_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=5, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_5 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_5, file = "lfc_p_dde_kidney_20wk-10wk_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=6, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_6 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_6, file = "lfc_p_dde_kidney_2y-10wk_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=7, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_7 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_7, file = "lfc_p_dde_kidney_2y-20wk_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=8, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_8 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_8, file = "lfc_p_dde_kidney_2y-20wk_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=9, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_9 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_9, file = "lfc_p_dde_kidney_2y-20wk_hd7.csv", quote = F)




dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise
## Add p-values to dataMatrixde objects ##
for (i in 1:ncol(fit2$p.value)) fit2$p.value[, i] <- p.adjust(fit2$p.value[, i], method = "BH")
pvalsadjst <- fit2$p.value

merged_pval_matrix_all <- merge(dataMatrixde, pvalsadjst,
                                by.x = "row.names",
                                by.y = "row.names",
                                all.X = T, all.Y = T
)
## make first column the row names ##
merged_pval_matrix_all <- data.frame(merged_pval_matrix_all, row.names = 1)

## Use no pvalue expression matrix to filter sig genes in matrix with pvalue##
ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

ExpressMatrixde_merged_pval_all <- merged_pval_matrix_all[row.names(merged_pval_matrix_all)%in%row.names(ExpressMatrixde),]

# filter for significant genes - up/down regulated
sigMask <- subset(sigMask, rowSums(sigMask) != 0)
# fit2$genes <- data.frame(ID_REF=rownames(fit2))

write.csv(ExpressMatrixde_merged_pval_all, file = "expression_matrix_dde_pval_all.csv", quote = F)
write.csv(ExpressMatrixde, file = "expression_matrix_dde.csv", quote = F)
write.csv(results, file = "results_dde.csv", quote = F)
write.csv(dataMatrixde, file = "full_expression_matrix_dde.csv", quote = F)

## Add gene names to pval expression matrix ##
ExpressMatrixde_genes_merged_pval_all <- merge(mouse_BM_ordered, ExpressMatrixde_merged_pval_all,
                                               by.x = "ensembl_gene_id",
                                               by.y = "row.names",
                                               all.X = T, all.Y = T
)

## Average entries with same gene name then make gene names the row names ##
ExpressMatrixde_genes_merged_pval_all <- ExpressMatrixde_genes_merged_pval_all[, !(names(ExpressMatrixde_genes_merged_pval_all) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes_merged_pval_all <- avereps(ExpressMatrixde_genes_merged_pval_all,
                                                 ID = ExpressMatrixde_genes_merged_pval_all$external_gene_name
)
rownames(ExpressMatrixde_genes_merged_pval_all) <- ExpressMatrixde_genes_merged_pval_all[, "external_gene_name"]
ExpressMatrixde_genes_merged_pval_all <- ExpressMatrixde_genes_merged_pval_all[, !(colnames(ExpressMatrixde_genes_merged_pval_all) %in% c("external_gene_name"))]
ExpressMatrixde_genes_merged_pval_all <- as.matrix(data.frame(ExpressMatrixde_genes_merged_pval_all, check.names = FALSE))
class(ExpressMatrixde_genes_merged_pval_all) <- "numeric"
write.csv(ExpressMatrixde_genes_merged_pval_all, "expression_matrix_dde_kidney_extgenename_pval.csv", quote = F)

## Add gene names to non-pval expression matrix  and remove duplicates##
ExpressMatrixde_genes_all <- merge(mouse_BM_ordered, ExpressMatrixde,
                                   by.x = "ensembl_gene_id",
                                   by.y = "row.names",
                                   all.X = T, all.Y = T
)

ExpressMatrixde_genes_all <- ExpressMatrixde_genes_all[, !(names(ExpressMatrixde_genes_all) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes_all <- avereps(ExpressMatrixde_genes_all,
                                     ID = ExpressMatrixde_genes_all$external_gene_name
)
rownames(ExpressMatrixde_genes_all) <- ExpressMatrixde_genes_all[, "external_gene_name"]
ExpressMatrixde_genes_all <- ExpressMatrixde_genes_all[, !(colnames(ExpressMatrixde_genes_all) %in% c("external_gene_name"))]
ExpressMatrixde_genes_all <- as.matrix(data.frame(ExpressMatrixde_genes_all, check.names = FALSE))
class(ExpressMatrixde_genes_all) <- "numeric"
write.csv(ExpressMatrixde_genes_all, "expression_matrix_genes_dde.csv", quote = F)


## Rename datamatrix and remove duplicates ##
dataMatrixde_genes <- merge(mouse_BM_ordered, dataMatrixde,
                            by.x = "ensembl_gene_id",
                            by.y = "row.names",
                            all.X = T, all.Y = T
)

dataMatrixde_genes <- dataMatrixde_genes[, !(names(dataMatrixde_genes) %in% c("ensembl_gene_id"))]
dataMatrixde_genes <- avereps(dataMatrixde_genes,
                              ID = dataMatrixde_genes$external_gene_name
)
rownames(dataMatrixde_genes) <- dataMatrixde_genes[, "external_gene_name"]
dataMatrixde_genes <- dataMatrixde_genes[, !(colnames(dataMatrixde_genes) %in% c("external_gene_name"))]
dataMatrixde_genes <- as.matrix(data.frame(dataMatrixde_genes, check.names = FALSE))
class(dataMatrixde_genes) <- "numeric"
write.csv(dataMatrixde_genes, "full_expression_matrix_dde_kidney_extgenename.csv", quote = F)

message(paste0("Dimensionality of DDE genes ", dim(ExpressMatrixde)[1]))

new_colnames <- c("20wk_hd2-10wk_hd2_kidney","20wk_hd4-10wk_hd4_kidney", "20wk_hd7-10wk_hd7_kidney", "2y_hd2-10wk_hd2_kidney", "2y_hd4-10wk_hd4_kidney", "2y_hd7-10wk_hd7_kidney",
                  "2y_hd2-20wk_hd2_kidney", "2y_hd4-20wk_hd4_kidney", "2y_hd7-20wk_hd7_kidney"
)

new_colnames_pval <- c("20wk_hd2-10wk_hd2_kidney_lfc", "20wk_hd4-10wk_hd4_kidney_lfc", "20wk_hd7-10wk_hd7_kidney_lfc", "2y_hd2-10wk_hd2_kidney_lfc", "2y_hd4-10wk_hd4_kidney_lfc", 
                       "2y_hd7-10wk_hd7_kidney_lfc", "2y_hd2-20wk_hd2_kidney_lfc", "2y_hd4-20wk_hd4_kidney_lfc", 
                       "2y_hd7-20wk_hd7_kidney_lfc","20wk_hd2-10wk_hd2_kidney_cpval", "20wk_hd4-10wk_hd4_kidney_cpval", "20wk_hd7-10wk_hd7_kidney_cpval", "2y_hd2-10wk_hd2_kidney_cpval", 
                       "2y_hd4-10wk_hd4_kidney_cpval", "2y_hd7-10wk_hd7_kidney_cpval", "2y_hd2-20wk_hd2_kidney_cpval", 
                       "2y_hd4-20wk_hd4_kidney_cpval", "2y_hd7-20wk_hd7_kidney_cpval"
)


colnames(results) <- new_colnames
colnames(dataMatrixde) <- new_colnames
colnames(dataMatrixde_genes) <- new_colnames
colnames(ExpressMatrixde) <- new_colnames
colnames(ExpressMatrixde_genes_all) <- new_colnames
colnames(ExpressMatrixde_genes_merged_pval_all) <- new_colnames_pval

# --heatmap
svglite("heatmap_all_kidney_dde.svg", width = 10, height = 10)
global_modulesallkidneydde <- heatmap.L.4(ExpressMatrixde_genes_all,
                                        figmargins = c(20, 5),
                                        cutoff = 1, distmethod = "euclidean", cexcol = 2, colcolorlist = colcolormatrix,
                                        clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

# Barplot
results_t <- t(summary(results))
results_t <- results_t[, -2]

for (i in 1:(length(row.names(results_t)))) {
  results_t[i, 1] <- results_t[i, 1] * -1
}

DE <- as.data.frame(results_t)
DE <- setnames(DE,
               old = c("Var1", "Var2", "Freq"),
               new = c("Time_Point", "group", "DE_genes")
)

# Create plot
ggplot(DE, aes(
  x = Time_Point, y = DE_genes, fill = group,
  label = DE$DE_genes
)) +
  geom_bar(stat = "identity", position = "identity") +
  # geom_text(size = 5, position = position_stack(vjust = 0) )+
  # theme_light() +
  theme_minimal() +
  scale_fill_manual(values = c("#0808c4", "#da9618")) +
  # xlab("Time point")
  ylab("Number of Differentially Expressed Genes") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 15)
  )
ggsave("barplot_all_kidney_dde.svg", width = 6, height = 4, units="in", dpi = 300)


###WebGestaltR analysis####
for (cluster in unique(global_modulesallkidneydde$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesallkidneydde$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism="mmusculus",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(cm),
              referenceGeneType = "ensembl_gene_id",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(global_modulesallkidneydde$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesallkidneydde$modulesrows == cluster)
  write.csv(names(genes), paste0(cluster,"_genes_kidney_dde.csv"))
  
  temp <- ExpressMatrixde_genes[row.names(ExpressMatrixde_genes)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_all_kidney_dde.csv"))
}


# Read in GMT file for CEMiTool #
gmt <- read_gmt("./m5.all.v2023.1.Mm.symbols.gmt")

# Read in sample annotation for CEMiTool
annot <- read.csv("CEMiTool_sample_annot_kidney.csv")

## All kidney CEMiTool ##
dataMatrixde_genes_all_kidney <- as.data.frame(dataMatrixde_genes_all_kidney)
cem <- cemitool(dataMatrixde_genes_all_kidney, annot, gmt)
generate_report(cem)
diagnostic_report(cem)



# All spleen DE #
contrastsmatrix <- c(
  # "age10.week.hd2.stL -age10.week.hdM.stL",
  # "age10.week.hd4.stL -age10.week.hdM.stL",
  # "age10.week.hd7.stL -age10.week.hdM.stL",
  # "age20.week.hd2.stL -age20.week.hdM.stL",
  # "age20.week.hd4.stL -age20.week.hdM.stL",
  # "age20.week.hd7.stL -age20.week.hdM.stL",
  # "age2.year.hd2.stL -age2.year.hdM.stL",
  # "age2.year.hd4.stL -age2.year.hdM.stL",
  # "age2.year.hd7.stL -age2.year.hdM.stL"
  # "age10.week.hd2.stK -age10.week.hdM.stK",
  # "age10.week.hd4.stK -age10.week.hdM.stK",
  # "age10.week.hd7.stK -age10.week.hdM.stK",
  # "age20.week.hd2.stK -age20.week.hdM.stK",
  # "age20.week.hd4.stK -age20.week.hdM.stK",
  # "age20.week.hd7.stK -age20.week.hdM.stK",
  # "age2.year.hd2.stK -age2.year.hdM.stK",
  # "age2.year.hd4.stK -age2.year.hdM.stK",
  # "age2.year.hd7.stK -age2.year.hdM.stK",
  "age10.week.hd2.stS -age10.week.hdM.stS",
  "age10.week.hd4.stS -age10.week.hdM.stS",
  "age10.week.hd7.stS -age10.week.hdM.stS",
  "age20.week.hd2.stS -age20.week.hdM.stS",
  "age20.week.hd4.stS -age20.week.hdM.stS",
  "age20.week.hd7.stS -age20.week.hdM.stS",
  "age2.year.hd2.stS -age2.year.hdM.stS",
  "age2.year.hd4.stS -age2.year.hdM.stS",
  "age2.year.hd7.stS -age2.year.hdM.stS"
  # "age10.week.hd2.stH -age10.week.hdM.stH",
  # "age10.week.hd4.stH -age10.week.hdM.stH",
  # "age10.week.hd7.stH -age10.week.hdM.stH",
  # "age20.week.hd2.stH -age20.week.hdM.stH",
  # "age20.week.hd4.stH -age20.week.hdM.stH",
  # "age20.week.hd7.stH -age20.week.hdM.stH",
  # "age2.year.hd2.stH -age2.year.hdM.stH",
  # "age2.year.hd4.stH -age2.year.hdM.stH",
  # "age2.year.hd7.stH -age2.year.hdM.stH"
)
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- treat(fit)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.01)
summary(results)

topTreat_data <- topTreat(fit2, coef=1, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_spleen_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=2, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_spleen_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=3, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_spleen_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=4, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_20wk_spleen_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=5, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_20wk_spleen_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=6, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_20wk_spleen_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=7, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_2yr_spleen_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=8, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_2yr_spleen_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=9, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_2yr_spleen_hd7.csv", quote = F)


dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

# filter for significant genes - up/down regulated
# sigMask <- subset(sigMask, rowSums(sigMask) != 0)
# fit2$genes <- data.frame(ID_REF=rownames(fit2))

write.csv(ExpressMatrixde, file = "expression_matrix_de.csv", quote = F)
write.csv(results, file = "results_de.csv", quote = F)
write.csv(dataMatrixde, file = "full_expression_matrix_de.csv", quote = F)
write.csv(fit2$t, file = "t_stats.csv", quote = F)
write.csv(fit2$p.value, file = "p_value.csv", quote = F)
for (i in 1:ncol(fit2$p.value)) fit2$p.value[, i] <- p.adjust(fit2$p.value[, i], method = "BH")
write.csv(fit2$p.value, file = "p_value_adj.csv", quote = F)

ExpressMatrixde_genes <- merge(mouse_BM_ordered, ExpressMatrixde,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)

ExpressMatrixde_genes <- ExpressMatrixde_genes[, !(names(ExpressMatrixde_genes) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes <- avereps(ExpressMatrixde_genes,
                                 ID = ExpressMatrixde_genes$external_gene_name
)
rownames(ExpressMatrixde_genes) <- ExpressMatrixde_genes[, "external_gene_name"]
ExpressMatrixde_genes <- ExpressMatrixde_genes[, !(colnames(ExpressMatrixde_genes) %in% c("external_gene_name"))]
ExpressMatrixde_genes <- as.matrix(data.frame(ExpressMatrixde_genes, check.names = FALSE))
class(ExpressMatrixde_genes) <- "numeric"
write.csv(ExpressMatrixde_genes, file = "expression_matrix_all_spleen_extgenename.csv", quote = F)






dataMatrixde_genes <- merge(mouse_BM_ordered, dataMatrixde,
                            by.x = "ensembl_gene_id",
                            by.y = "row.names",
                            all.X = T, all.Y = T
)

dataMatrixde_genes <- dataMatrixde_genes[, !(names(dataMatrixde_genes) %in% c("ensembl_gene_id"))]
dataMatrixde_genes <- avereps(dataMatrixde_genes,
                              ID = dataMatrixde_genes$external_gene_name
)
rownames(dataMatrixde_genes) <- dataMatrixde_genes[, "external_gene_name"]
dataMatrixde_genes <- dataMatrixde_genes[, !(colnames(dataMatrixde_genes) %in% c("external_gene_name"))]
dataMatrixde_genes <- as.matrix(data.frame(dataMatrixde_genes, check.names = FALSE))
class(dataMatrixde_genes) <- "numeric"
write.csv(dataMatrixde_genes, file = "full_expression_matrix_all_spleen_extgenename.csv", quote = F)

message(paste0("Dimensionality of DE genes ", dim(ExpressMatrixde)[1]))

new_colnames <- c(
  "10wk_spleen_hd2", "10wk_spleen_hd4", "10wk_spleen_hd7", "20wk_spleen_hd2", "20wk_spleen_hd4", "20wk_spleen_hd7",
  "2yr_spleen_hd2", "2yr_spleen_hd4", "2yr_spleen_hd7"
)

colnames(results) <- new_colnames
colnames(dataMatrixde) <- new_colnames
colnames(dataMatrixde_genes) <- new_colnames
colnames(ExpressMatrixde) <- new_colnames
colnames(ExpressMatrixde_genes) <- new_colnames

#######SAVE THESE FOR LATER ######
dataMatrixde_all_spleen <- dataMatrixde
dataMatrixde_genes_all_spleen <- dataMatrixde_genes
ExpressMatrixde_all_spleen <- ExpressMatrixde
ExpressMatrixde_all_genes_spleen <- ExpressMatrixde_genes

# --heatmap
svglite("heatmap_all_spleen.svg", width = 10, height = 10)
global_modulesallspleen <- heatmap.L.4(ExpressMatrixde_genes,
                                       figmargins = c(20, 5),
                                       cutoff = 1, distmethod = "euclidean", cexcol = 2, colcolorlist = colcolormatrix,
                                       clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

# Barplot
results_t <- t(summary(results))
results_t <- results_t[, -2]

for (i in 1:(length(row.names(results_t)))) {
  results_t[i, 1] <- results_t[i, 1] * -1
}

DE <- as.data.frame(results_t)
DE <- setnames(DE,
               old = c("Var1", "Var2", "Freq"),
               new = c("Time_Point", "group", "DE_genes")
)

# Create plot
ggplot(DE, aes(
  x = Time_Point, y = DE_genes, fill = group,
  label = DE$DE_genes
)) +
  geom_bar(stat = "identity", position = "identity") +
  # geom_text(size = 5, position = position_stack(vjust = 0) )+
  # theme_light() +
  theme_minimal() +
  scale_fill_manual(values = c("#0808c4", "#da9618")) +
  # xlab("Time point")
  ylab("Number of Differentially Expressed Genes") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 15)
  )
ggsave("barplot_all_spleen.svg", width = 6, height = 4, units="in", dpi = 300)


###WebGestaltR analysis####
for (cluster in unique(global_modulesallspleen$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesallspleen$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism="mmusculus",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(cm),
              referenceGeneType = "ensembl_gene_id",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(global_modulesallspleen$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesallspleen$modulesrows == cluster)
  write.csv(names(genes), paste0(cluster,"_genes_all_spleen.csv"))
  
  temp <- ExpressMatrixde_genes[row.names(ExpressMatrixde_genes)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_all_spleen.csv"))
}


## spleen DDE ##

contrastsmatrix <- c(
  "(age20.week.hd2.stS -age20.week.hdM.stS) - (age10.week.hd2.stS -age10.week.hdM.stS)",
  "(age20.week.hd4.stS -age20.week.hdM.stS) - (age10.week.hd4.stS -age10.week.hdM.stS)",
  "(age20.week.hd7.stS -age20.week.hdM.stS) - (age10.week.hd7.stS -age10.week.hdM.stS)",
  "(age2.year.hd2.stS -age2.year.hdM.stS) - (age10.week.hd2.stS -age10.week.hdM.stS)",
  "(age2.year.hd4.stS -age2.year.hdM.stS) - (age10.week.hd4.stS -age10.week.hdM.stS)",
  "(age2.year.hd7.stS -age2.year.hdM.stS) - (age10.week.hd7.stS -age10.week.hdM.stS)",
  "(age2.year.hd2.stS -age2.year.hdM.stS) - (age20.week.hd2.stS -age20.week.hdM.stS)",
  "(age2.year.hd4.stS -age2.year.hdM.stS) - (age20.week.hd4.stS -age20.week.hdM.stS)",
  "(age2.year.hd7.stS -age2.year.hdM.stS) - (age20.week.hd7.stS -age20.week.hdM.stS)"
)
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- treat(fit)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.01)
summary(results)

#Generate topTreat tables for each comparison#
topTreat_data <- topTreat(fit2, coef=1, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_1 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_1, file = "lfc_p_dde_spleen_20wk-10wk_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=2, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_2 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_2, file = "lfc_p_dde_spleen_2y-10wk_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=3, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_3 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_3, file = "lfc_p_dde_spleen_20wk-10wk_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=4, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_4 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_4, file = "lfc_p_dde_spleen_2y-10wk_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=5, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_5 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_5, file = "lfc_p_dde_spleen_20wk-10wk_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=6, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_6 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_6, file = "lfc_p_dde_spleen_2y-10wk_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=7, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_7 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_7, file = "lfc_p_dde_spleen_2y-20wk_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=8, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_8 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_8, file = "lfc_p_dde_spleen_2y-20wk_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=9, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_9 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_9, file = "lfc_p_dde_spleen_2y-20wk_hd7.csv", quote = F)




dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise
## Add p-values to dataMatrixde objects ##
for (i in 1:ncol(fit2$p.value)) fit2$p.value[, i] <- p.adjust(fit2$p.value[, i], method = "BH")
pvalsadjst <- fit2$p.value

merged_pval_matrix_all <- merge(dataMatrixde, pvalsadjst,
                                by.x = "row.names",
                                by.y = "row.names",
                                all.X = T, all.Y = T
)
## make first column the row names ##
merged_pval_matrix_all <- data.frame(merged_pval_matrix_all, row.names = 1)

## Use no pvalue expression matrix to filter sig genes in matrix with pvalue##
ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

ExpressMatrixde_merged_pval_all <- merged_pval_matrix_all[row.names(merged_pval_matrix_all)%in%row.names(ExpressMatrixde),]

# filter for significant genes - up/down regulated
sigMask <- subset(sigMask, rowSums(sigMask) != 0)
# fit2$genes <- data.frame(ID_REF=rownames(fit2))

write.csv(ExpressMatrixde_merged_pval_all, file = "expression_matrix_dde_pval_all.csv", quote = F)
write.csv(ExpressMatrixde, file = "expression_matrix_dde.csv", quote = F)
write.csv(results, file = "results_dde.csv", quote = F)
write.csv(dataMatrixde, file = "full_expression_matrix_dde.csv", quote = F)

## Add gene names to pval expression matrix ##
ExpressMatrixde_genes_merged_pval_all <- merge(mouse_BM_ordered, ExpressMatrixde_merged_pval_all,
                                               by.x = "ensembl_gene_id",
                                               by.y = "row.names",
                                               all.X = T, all.Y = T
)

## Average entries with same gene name then make gene names the row names ##
ExpressMatrixde_genes_merged_pval_all <- ExpressMatrixde_genes_merged_pval_all[, !(names(ExpressMatrixde_genes_merged_pval_all) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes_merged_pval_all <- avereps(ExpressMatrixde_genes_merged_pval_all,
                                                 ID = ExpressMatrixde_genes_merged_pval_all$external_gene_name
)
rownames(ExpressMatrixde_genes_merged_pval_all) <- ExpressMatrixde_genes_merged_pval_all[, "external_gene_name"]
ExpressMatrixde_genes_merged_pval_all <- ExpressMatrixde_genes_merged_pval_all[, !(colnames(ExpressMatrixde_genes_merged_pval_all) %in% c("external_gene_name"))]
ExpressMatrixde_genes_merged_pval_all <- as.matrix(data.frame(ExpressMatrixde_genes_merged_pval_all, check.names = FALSE))
class(ExpressMatrixde_genes_merged_pval_all) <- "numeric"
write.csv(ExpressMatrixde_genes_merged_pval_all, "expression_matrix_dde_spleen_extgenename_pval.csv", quote = F)

## Add gene names to non-pval expression matrix  and remove duplicates##
ExpressMatrixde_genes_all <- merge(mouse_BM_ordered, ExpressMatrixde,
                                   by.x = "ensembl_gene_id",
                                   by.y = "row.names",
                                   all.X = T, all.Y = T
)

ExpressMatrixde_genes_all <- ExpressMatrixde_genes_all[, !(names(ExpressMatrixde_genes_all) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes_all <- avereps(ExpressMatrixde_genes_all,
                                     ID = ExpressMatrixde_genes_all$external_gene_name
)
rownames(ExpressMatrixde_genes_all) <- ExpressMatrixde_genes_all[, "external_gene_name"]
ExpressMatrixde_genes_all <- ExpressMatrixde_genes_all[, !(colnames(ExpressMatrixde_genes_all) %in% c("external_gene_name"))]
ExpressMatrixde_genes_all <- as.matrix(data.frame(ExpressMatrixde_genes_all, check.names = FALSE))
class(ExpressMatrixde_genes_all) <- "numeric"
write.csv(ExpressMatrixde_genes_all, "expression_matrix_genes_dde.csv", quote = F)


## Rename datamatrix and remove duplicates ##
dataMatrixde_genes <- merge(mouse_BM_ordered, dataMatrixde,
                            by.x = "ensembl_gene_id",
                            by.y = "row.names",
                            all.X = T, all.Y = T
)

dataMatrixde_genes <- dataMatrixde_genes[, !(names(dataMatrixde_genes) %in% c("ensembl_gene_id"))]
dataMatrixde_genes <- avereps(dataMatrixde_genes,
                              ID = dataMatrixde_genes$external_gene_name
)
rownames(dataMatrixde_genes) <- dataMatrixde_genes[, "external_gene_name"]
dataMatrixde_genes <- dataMatrixde_genes[, !(colnames(dataMatrixde_genes) %in% c("external_gene_name"))]
dataMatrixde_genes <- as.matrix(data.frame(dataMatrixde_genes, check.names = FALSE))
class(dataMatrixde_genes) <- "numeric"
write.csv(dataMatrixde_genes, "full_expression_matrix_dde_spleen_extgenename.csv", quote = F)

message(paste0("Dimensionality of DDE genes ", dim(ExpressMatrixde)[1]))

new_colnames <- c("20wk_hd2-10wk_hd2_spleen","20wk_hd4-10wk_hd4_spleen", "20wk_hd7-10wk_hd7_spleen", "2y_hd2-10wk_hd2_spleen", "2y_hd4-10wk_hd4_spleen", "2y_hd7-10wk_hd7_spleen",
                  "2y_hd2-20wk_hd2_spleen", "2y_hd4-20wk_hd4_spleen", "2y_hd7-20wk_hd7_spleen"
)

new_colnames_pval <- c("20wk_hd2-10wk_hd2_spleen_lfc", "20wk_hd4-10wk_hd4_spleen_lfc", "20wk_hd7-10wk_hd7_spleen_lfc", "2y_hd2-10wk_hd2_spleen_lfc", "2y_hd4-10wk_hd4_spleen_lfc", 
                       "2y_hd7-10wk_hd7_spleen_lfc", "2y_hd2-20wk_hd2_spleen_lfc", "2y_hd4-20wk_hd4_spleen_lfc", 
                       "2y_hd7-20wk_hd7_spleen_lfc","20wk_hd2-10wk_hd2_spleen_cpval", "20wk_hd4-10wk_hd4_spleen_cpval", "20wk_hd7-10wk_hd7_spleen_cpval", "2y_hd2-10wk_hd2_spleen_cpval", 
                       "2y_hd4-10wk_hd4_spleen_cpval", "2y_hd7-10wk_hd7_spleen_cpval", "2y_hd2-20wk_hd2_spleen_cpval", 
                       "2y_hd4-20wk_hd4_spleen_cpval", "2y_hd7-20wk_hd7_spleen_cpval"
)


colnames(results) <- new_colnames
colnames(dataMatrixde) <- new_colnames
colnames(dataMatrixde_genes) <- new_colnames
colnames(ExpressMatrixde) <- new_colnames
colnames(ExpressMatrixde_genes_all) <- new_colnames
colnames(ExpressMatrixde_genes_merged_pval_all) <- new_colnames_pval

# --heatmap
svglite("heatmap_all_spleen_dde.svg", width = 10, height = 10)
global_modulesallspleendde <- heatmap.L.4(ExpressMatrixde_genes_all,
                                          figmargins = c(20, 5),
                                          cutoff = 1, distmethod = "euclidean", cexcol = 2, colcolorlist = colcolormatrix,
                                          clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

# Barplot
results_t <- t(summary(results))
results_t <- results_t[, -2]

for (i in 1:(length(row.names(results_t)))) {
  results_t[i, 1] <- results_t[i, 1] * -1
}

DE <- as.data.frame(results_t)
DE <- setnames(DE,
               old = c("Var1", "Var2", "Freq"),
               new = c("Time_Point", "group", "DE_genes")
)

# Create plot
ggplot(DE, aes(
  x = Time_Point, y = DE_genes, fill = group,
  label = DE$DE_genes
)) +
  geom_bar(stat = "identity", position = "identity") +
  # geom_text(size = 5, position = position_stack(vjust = 0) )+
  # theme_light() +
  theme_minimal() +
  scale_fill_manual(values = c("#0808c4", "#da9618")) +
  # xlab("Time point")
  ylab("Number of Differentially Expressed Genes") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 15)
  )
ggsave("barplot_all_spleen_dde.svg", width = 6, height = 4, units="in", dpi = 300)


###WebGestaltR analysis####
for (cluster in unique(global_modulesallspleendde$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesallspleendde$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism="mmusculus",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(cm),
              referenceGeneType = "ensembl_gene_id",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(global_modulesallspleendde$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesallspleendde$modulesrows == cluster)
  write.csv(names(genes), paste0(cluster,"_genes_spleen_dde.csv"))
  
  temp <- ExpressMatrixde_genes[row.names(ExpressMatrixde_genes)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_all_spleen_dde.csv"))
}


# Read in GMT file for CEMiTool #
gmt <- read_gmt("./m5.all.v2023.1.Mm.symbols.gmt")

# Read in sample annotation for CEMiTool
annot <- read.csv("CEMiTool_sample_annot_spleen.csv")

## All spleen CEMiTool ##
dataMatrixde_genes_all_spleen <- as.data.frame(dataMatrixde_genes_all_spleen)
cem <- cemitool(dataMatrixde_genes_all_spleen, annot, gmt)
generate_report(cem)
diagnostic_report(cem)



# All heart DE #
contrastsmatrix <- c(
  # "age10.week.hd2.stL -age10.week.hdM.stL",
  # "age10.week.hd4.stL -age10.week.hdM.stL",
  # "age10.week.hd7.stL -age10.week.hdM.stL",
  # "age20.week.hd2.stL -age20.week.hdM.stL",
  # "age20.week.hd4.stL -age20.week.hdM.stL",
  # "age20.week.hd7.stL -age20.week.hdM.stL",
  # "age2.year.hd2.stL -age2.year.hdM.stL",
  # "age2.year.hd4.stL -age2.year.hdM.stL",
  # "age2.year.hd7.stL -age2.year.hdM.stL"
  # "age10.week.hd2.stK -age10.week.hdM.stK",
  # "age10.week.hd4.stK -age10.week.hdM.stK",
  # "age10.week.hd7.stK -age10.week.hdM.stK",
  # "age20.week.hd2.stK -age20.week.hdM.stK",
  # "age20.week.hd4.stK -age20.week.hdM.stK",
  # "age20.week.hd7.stK -age20.week.hdM.stK",
  # "age2.year.hd2.stK -age2.year.hdM.stK",
  # "age2.year.hd4.stK -age2.year.hdM.stK",
  # "age2.year.hd7.stK -age2.year.hdM.stK",
  # "age10.week.hd2.stS -age10.week.hdM.stS",
  # "age10.week.hd4.stS -age10.week.hdM.stS",
  # "age10.week.hd7.stS -age10.week.hdM.stS",
  # "age20.week.hd2.stS -age20.week.hdM.stS",
  # "age20.week.hd4.stS -age20.week.hdM.stS",
  # "age20.week.hd7.stS -age20.week.hdM.stS",
  # "age2.year.hd2.stS -age2.year.hdM.stS",
  # "age2.year.hd4.stS -age2.year.hdM.stS",
  # "age2.year.hd7.stS -age2.year.hdM.stS",
  "age10.week.hd2.stH -age10.week.hdM.stH",
  "age10.week.hd4.stH -age10.week.hdM.stH",
  "age10.week.hd7.stH -age10.week.hdM.stH",
  "age20.week.hd2.stH -age20.week.hdM.stH",
  "age20.week.hd4.stH -age20.week.hdM.stH",
  "age20.week.hd7.stH -age20.week.hdM.stH",
  "age2.year.hd2.stH -age2.year.hdM.stH",
  "age2.year.hd4.stH -age2.year.hdM.stH",
  "age2.year.hd7.stH -age2.year.hdM.stH"
)
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- treat(fit)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.01)
summary(results)

topTreat_data <- topTreat(fit2, coef=1, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_heart_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=2, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_heart_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=3, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_10wk_heart_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=4, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_20wk_heart_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=5, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_20wk_heart_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=6, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_20wk_heart_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=7, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_2yr_heart_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=8, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_2yr_heart_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=9, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes <- merge(mouse_BM_ordered, topTreat_data,
                             by.x = "ensembl_gene_id",
                             by.y = "row.names",
                             all.X = T, all.Y = T
)
write.csv(topTreat_data_genes, file = "lfc_p_2yr_heart_hd7.csv", quote = F)


dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise

ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

# filter for significant genes - up/down regulated
# sigMask <- subset(sigMask, rowSums(sigMask) != 0)
# fit2$genes <- data.frame(ID_REF=rownames(fit2))

write.csv(ExpressMatrixde, file = "expression_matrix_de.csv", quote = F)
write.csv(results, file = "results_de.csv", quote = F)
write.csv(dataMatrixde, file = "full_expression_matrix_de.csv", quote = F)
write.csv(fit2$t, file = "t_stats.csv", quote = F)
write.csv(fit2$p.value, file = "p_value.csv", quote = F)
for (i in 1:ncol(fit2$p.value)) fit2$p.value[, i] <- p.adjust(fit2$p.value[, i], method = "BH")
write.csv(fit2$p.value, file = "p_value_adj.csv", quote = F)

ExpressMatrixde_genes <- merge(mouse_BM_ordered, ExpressMatrixde,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)

ExpressMatrixde_genes <- ExpressMatrixde_genes[, !(names(ExpressMatrixde_genes) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes <- avereps(ExpressMatrixde_genes,
                                 ID = ExpressMatrixde_genes$external_gene_name
)
rownames(ExpressMatrixde_genes) <- ExpressMatrixde_genes[, "external_gene_name"]
ExpressMatrixde_genes <- ExpressMatrixde_genes[, !(colnames(ExpressMatrixde_genes) %in% c("external_gene_name"))]
ExpressMatrixde_genes <- as.matrix(data.frame(ExpressMatrixde_genes, check.names = FALSE))
class(ExpressMatrixde_genes) <- "numeric"
write.csv(ExpressMatrixde_genes, file = "expression_matrix_all_heart_extgenename.csv", quote = F)






dataMatrixde_genes <- merge(mouse_BM_ordered, dataMatrixde,
                            by.x = "ensembl_gene_id",
                            by.y = "row.names",
                            all.X = T, all.Y = T
)

dataMatrixde_genes <- dataMatrixde_genes[, !(names(dataMatrixde_genes) %in% c("ensembl_gene_id"))]
dataMatrixde_genes <- avereps(dataMatrixde_genes,
                              ID = dataMatrixde_genes$external_gene_name
)
rownames(dataMatrixde_genes) <- dataMatrixde_genes[, "external_gene_name"]
dataMatrixde_genes <- dataMatrixde_genes[, !(colnames(dataMatrixde_genes) %in% c("external_gene_name"))]
dataMatrixde_genes <- as.matrix(data.frame(dataMatrixde_genes, check.names = FALSE))
class(dataMatrixde_genes) <- "numeric"
write.csv(dataMatrixde_genes, file = "full_expression_matrix_all_heart_extgenename.csv", quote = F)

message(paste0("Dimensionality of DE genes ", dim(ExpressMatrixde)[1]))

new_colnames <- c(
  "10wk_heart_hd2", "10wk_heart_hd4", "10wk_heart_hd7", "20wk_heart_hd2", "20wk_heart_hd4", "20wk_heart_hd7",
  "2yr_heart_hd2", "2yr_heart_hd4", "2yr_heart_hd7"
)

colnames(results) <- new_colnames
colnames(dataMatrixde) <- new_colnames
colnames(dataMatrixde_genes) <- new_colnames
colnames(ExpressMatrixde) <- new_colnames
colnames(ExpressMatrixde_genes) <- new_colnames

#######SAVE THESE FOR LATER ######
dataMatrixde_all_heart <- dataMatrixde
dataMatrixde_genes_all_heart <- dataMatrixde_genes
ExpressMatrixde_all_heart <- ExpressMatrixde
ExpressMatrixde_all_genes_heart <- ExpressMatrixde_genes

# --heatmap
svglite("heatmap_all_heart.svg", width = 10, height = 10)
global_modulesallheart <- heatmap.L.4(ExpressMatrixde_all_genes_heart,
                                       figmargins = c(20, 5),
                                       cutoff = 1, distmethod = "euclidean", cexcol = 2, colcolorlist = colcolormatrix,
                                       clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

# Barplot
results_t <- t(summary(results))
results_t <- results_t[, -2]

for (i in 1:(length(row.names(results_t)))) {
  results_t[i, 1] <- results_t[i, 1] * -1
}

DE <- as.data.frame(results_t)
DE <- setnames(DE,
               old = c("Var1", "Var2", "Freq"),
               new = c("Time_Point", "group", "DE_genes")
)

# Create plot
ggplot(DE, aes(
  x = Time_Point, y = DE_genes, fill = group,
  label = DE$DE_genes
)) +
  geom_bar(stat = "identity", position = "identity") +
  # geom_text(size = 5, position = position_stack(vjust = 0) )+
  # theme_light() +
  theme_minimal() +
  scale_fill_manual(values = c("#0808c4", "#da9618")) +
  # xlab("Time point")
  ylab("Number of Differentially Expressed Genes") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 15)
  )
ggsave("barplot_all_heart.svg", width = 6, height = 4, units="in", dpi = 300)


###WebGestaltR analysis####
for (cluster in unique(global_modulesallheart$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesallheart$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism="mmusculus",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(cm),
              referenceGeneType = "ensembl_gene_id",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(global_modulesallheart$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesallheart$modulesrows == cluster)
  write.csv(names(genes), paste0(cluster,"_genes_all_heart.csv"))
  
  temp <- ExpressMatrixde_genes[row.names(ExpressMatrixde_genes)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_all_heart.csv"))
}


## heart DDE ##

contrastsmatrix <- c(
  "(age20.week.hd2.stH -age20.week.hdM.stH) - (age10.week.hd2.stH -age10.week.hdM.stH)",
  "(age20.week.hd4.stH -age20.week.hdM.stH) - (age10.week.hd4.stH -age10.week.hdM.stH)",
  "(age20.week.hd7.stH -age20.week.hdM.stH) - (age10.week.hd7.stH -age10.week.hdM.stH)",
  "(age2.year.hd2.stH -age2.year.hdM.stH) - (age10.week.hd2.stH -age10.week.hdM.stH)",
  "(age2.year.hd4.stH -age2.year.hdM.stH) - (age10.week.hd4.stH -age10.week.hdM.stH)",
  "(age2.year.hd7.stH -age2.year.hdM.stH) - (age10.week.hd7.stH -age10.week.hdM.stH)",
  "(age2.year.hd2.stH -age2.year.hdM.stH) - (age20.week.hd2.stH -age20.week.hdM.stH)",
  "(age2.year.hd4.stH -age2.year.hdM.stH) - (age20.week.hd4.stH -age20.week.hdM.stH)",
  "(age2.year.hd7.stH -age2.year.hdM.stH) - (age20.week.hd7.stH -age20.week.hdM.stH)"
)
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm)

fit <- contrasts.fit(Pi.lmfit, contr) # This results in a 3-dimensional array, which is not entirely analogous to fit2 from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
fit2 <- treat(fit)
results <- decideTests(fit2, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.01)
summary(results)

#Generate topTreat tables for each comparison#
topTreat_data <- topTreat(fit2, coef=1, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_1 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_1, file = "lfc_p_dde_heart_20wk-10wk_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=2, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_2 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_2, file = "lfc_p_dde_heart_2y-10wk_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=3, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_3 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_3, file = "lfc_p_dde_heart_20wk-10wk_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=4, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_4 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_4, file = "lfc_p_dde_heart_2y-10wk_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=5, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_5 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_5, file = "lfc_p_dde_heart_20wk-10wk_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=6, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_6 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_6, file = "lfc_p_dde_heart_2y-10wk_hd7.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=7, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_7 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_7, file = "lfc_p_dde_heart_2y-20wk_hd2.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=8, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_8 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_8, file = "lfc_p_dde_heart_2y-20wk_hd4.csv", quote = F)

topTreat_data <- topTreat(fit2, coef=9, sort.by="p", resort.by=NULL, number=length(fit2$coefficients))
topTreat_data_genes_9 <- merge(mouse_BM_ordered, topTreat_data,
                               by.x = "ensembl_gene_id",
                               by.y = "row.names",
                               all.X = T, all.Y = T
)
write.csv(topTreat_data_genes_9, file = "lfc_p_dde_heart_2y-20wk_hd7.csv", quote = F)




dataMatrixde <- fit2$coefficients
sigMask <- dataMatrixde * (results**2) # 1 if significant, 0 otherwise
## Add p-values to dataMatrixde objects ##
for (i in 1:ncol(fit2$p.value)) fit2$p.value[, i] <- p.adjust(fit2$p.value[, i], method = "BH")
pvalsadjst <- fit2$p.value

merged_pval_matrix_all <- merge(dataMatrixde, pvalsadjst,
                                by.x = "row.names",
                                by.y = "row.names",
                                all.X = T, all.Y = T
)
## make first column the row names ##
merged_pval_matrix_all <- data.frame(merged_pval_matrix_all, row.names = 1)

## Use no pvalue expression matrix to filter sig genes in matrix with pvalue##
ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

ExpressMatrixde_merged_pval_all <- merged_pval_matrix_all[row.names(merged_pval_matrix_all)%in%row.names(ExpressMatrixde),]

# filter for significant genes - up/down regulated
sigMask <- subset(sigMask, rowSums(sigMask) != 0)
# fit2$genes <- data.frame(ID_REF=rownames(fit2))

write.csv(ExpressMatrixde_merged_pval_all, file = "expression_matrix_dde_pval_all.csv", quote = F)
write.csv(ExpressMatrixde, file = "expression_matrix_dde.csv", quote = F)
write.csv(results, file = "results_dde.csv", quote = F)
write.csv(dataMatrixde, file = "full_expression_matrix_dde.csv", quote = F)

## Add gene names to pval expression matrix ##
ExpressMatrixde_genes_merged_pval_all <- merge(mouse_BM_ordered, ExpressMatrixde_merged_pval_all,
                                               by.x = "ensembl_gene_id",
                                               by.y = "row.names",
                                               all.X = T, all.Y = T
)

## Average entries with same gene name then make gene names the row names ##
ExpressMatrixde_genes_merged_pval_all <- ExpressMatrixde_genes_merged_pval_all[, !(names(ExpressMatrixde_genes_merged_pval_all) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes_merged_pval_all <- avereps(ExpressMatrixde_genes_merged_pval_all,
                                                 ID = ExpressMatrixde_genes_merged_pval_all$external_gene_name
)
rownames(ExpressMatrixde_genes_merged_pval_all) <- ExpressMatrixde_genes_merged_pval_all[, "external_gene_name"]
ExpressMatrixde_genes_merged_pval_all <- ExpressMatrixde_genes_merged_pval_all[, !(colnames(ExpressMatrixde_genes_merged_pval_all) %in% c("external_gene_name"))]
ExpressMatrixde_genes_merged_pval_all <- as.matrix(data.frame(ExpressMatrixde_genes_merged_pval_all, check.names = FALSE))
class(ExpressMatrixde_genes_merged_pval_all) <- "numeric"
write.csv(ExpressMatrixde_genes_merged_pval_all, "expression_matrix_dde_heart_extgenename_pval.csv", quote = F)

## Add gene names to non-pval expression matrix  and remove duplicates##
ExpressMatrixde_genes_all <- merge(mouse_BM_ordered, ExpressMatrixde,
                                   by.x = "ensembl_gene_id",
                                   by.y = "row.names",
                                   all.X = T, all.Y = T
)

ExpressMatrixde_genes_all <- ExpressMatrixde_genes_all[, !(names(ExpressMatrixde_genes_all) %in% c("ensembl_gene_id"))]
ExpressMatrixde_genes_all <- avereps(ExpressMatrixde_genes_all,
                                     ID = ExpressMatrixde_genes_all$external_gene_name
)
rownames(ExpressMatrixde_genes_all) <- ExpressMatrixde_genes_all[, "external_gene_name"]
ExpressMatrixde_genes_all <- ExpressMatrixde_genes_all[, !(colnames(ExpressMatrixde_genes_all) %in% c("external_gene_name"))]
ExpressMatrixde_genes_all <- as.matrix(data.frame(ExpressMatrixde_genes_all, check.names = FALSE))
class(ExpressMatrixde_genes_all) <- "numeric"
write.csv(ExpressMatrixde_genes_all, "expression_matrix_genes_dde.csv", quote = F)


## Rename datamatrix and remove duplicates ##
dataMatrixde_genes <- merge(mouse_BM_ordered, dataMatrixde,
                            by.x = "ensembl_gene_id",
                            by.y = "row.names",
                            all.X = T, all.Y = T
)

dataMatrixde_genes <- dataMatrixde_genes[, !(names(dataMatrixde_genes) %in% c("ensembl_gene_id"))]
dataMatrixde_genes <- avereps(dataMatrixde_genes,
                              ID = dataMatrixde_genes$external_gene_name
)
rownames(dataMatrixde_genes) <- dataMatrixde_genes[, "external_gene_name"]
dataMatrixde_genes <- dataMatrixde_genes[, !(colnames(dataMatrixde_genes) %in% c("external_gene_name"))]
dataMatrixde_genes <- as.matrix(data.frame(dataMatrixde_genes, check.names = FALSE))
class(dataMatrixde_genes) <- "numeric"
write.csv(dataMatrixde_genes, "full_expression_matrix_dde_heart_extgenename.csv", quote = F)

message(paste0("Dimensionality of DDE genes ", dim(ExpressMatrixde)[1]))

new_colnames <- c("20wk_hd2-10wk_hd2_heart","20wk_hd4-10wk_hd4_heart", "20wk_hd7-10wk_hd7_heart", "2y_hd2-10wk_hd2_heart", "2y_hd4-10wk_hd4_heart", "2y_hd7-10wk_hd7_heart",
                  "2y_hd2-20wk_hd2_heart", "2y_hd4-20wk_hd4_heart", "2y_hd7-20wk_hd7_heart"
)

new_colnames_pval <- c("20wk_hd2-10wk_hd2_heart_lfc", "20wk_hd4-10wk_hd4_heart_lfc", "20wk_hd7-10wk_hd7_heart_lfc", "2y_hd2-10wk_hd2_heart_lfc", "2y_hd4-10wk_hd4_heart_lfc", 
                       "2y_hd7-10wk_hd7_heart_lfc", "2y_hd2-20wk_hd2_heart_lfc", "2y_hd4-20wk_hd4_heart_lfc", 
                       "2y_hd7-20wk_hd7_heart_lfc","20wk_hd2-10wk_hd2_heart_cpval", "20wk_hd4-10wk_hd4_heart_cpval", "20wk_hd7-10wk_hd7_heart_cpval", "2y_hd2-10wk_hd2_heart_cpval", 
                       "2y_hd4-10wk_hd4_heart_cpval", "2y_hd7-10wk_hd7_heart_cpval", "2y_hd2-20wk_hd2_heart_cpval", 
                       "2y_hd4-20wk_hd4_heart_cpval", "2y_hd7-20wk_hd7_heart_cpval"
)


colnames(results) <- new_colnames
colnames(dataMatrixde) <- new_colnames
colnames(dataMatrixde_genes) <- new_colnames
colnames(ExpressMatrixde) <- new_colnames
colnames(ExpressMatrixde_genes_all) <- new_colnames
colnames(ExpressMatrixde_genes_merged_pval_all) <- new_colnames_pval

# --heatmap
svglite("heatmap_all_heart_dde.svg", width = 10, height = 10)
global_modulesallheartdde <- heatmap.L.4(ExpressMatrixde_genes_all,
                                          figmargins = c(20, 5),
                                          cutoff = 1, distmethod = "euclidean", cexcol = 2, colcolorlist = colcolormatrix,
                                          clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9
)
dev.off()

# Barplot
results_t <- t(summary(results))
results_t <- results_t[, -2]

for (i in 1:(length(row.names(results_t)))) {
  results_t[i, 1] <- results_t[i, 1] * -1
}

DE <- as.data.frame(results_t)
DE <- setnames(DE,
               old = c("Var1", "Var2", "Freq"),
               new = c("Time_Point", "group", "DE_genes")
)

# Create plot
ggplot(DE, aes(
  x = Time_Point, y = DE_genes, fill = group,
  label = DE$DE_genes
)) +
  geom_bar(stat = "identity", position = "identity") +
  # geom_text(size = 5, position = position_stack(vjust = 0) )+
  # theme_light() +
  theme_minimal() +
  scale_fill_manual(values = c("#0808c4", "#da9618")) +
  # xlab("Time point")
  ylab("Number of Differentially Expressed Genes") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 15)
  )
ggsave("barplot_all_heart_dde.svg", width = 6, height = 4, units="in", dpi = 300)


###WebGestaltR analysis####
for (cluster in unique(global_modulesallheartdde$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesallheartdde$modulesrows == cluster)
  WebGestaltR(enrichMethod="ORA",
              organism="mmusculus",
              interestGene = names(genes),
              interestGeneType="genesymbol",
              referenceGene = row.names(cm),
              referenceGeneType = "ensembl_gene_id",
              projectName = paste0(cluster,"_ORA"),
              minNum=5,
              enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                 "geneontology_Molecular_Function_noRedundant",
                                 "geneontology_Cellular_Component_noRedundant",
                                 "pathway_KEGG"
              ),
              sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05, nThreads=8)
}

## Generate list of de genes per cluster and LFCs ##

for (cluster in unique(global_modulesallheartdde$modulesrows)) {
  print(paste0("STATUS: gene enrichments for module ", cluster))
  genes <- which(global_modulesallheartdde$modulesrows == cluster)
  write.csv(names(genes), paste0(cluster,"_genes_heart_dde.csv"))
  
  temp <- ExpressMatrixde_genes[row.names(ExpressMatrixde_genes)%in%names(genes),]
  write.csv(temp, paste0(cluster,"_gene_lfcs_all_heart_dde.csv"))
}


# Read in GMT file for CEMiTool #
gmt <- read_gmt("./m5.all.v2023.1.Mm.symbols.gmt")

# Read in sample annotation for CEMiTool
annot <- read.csv("CEMiTool_sample_annot_heart.csv")

## All heart CEMiTool ##
dataMatrixde_genes_all_heart <- as.data.frame(dataMatrixde_genes_all_heart)
cem <- cemitool(dataMatrixde_genes_all_heart, annot, gmt)
generate_report(cem)
diagnostic_report(cem)