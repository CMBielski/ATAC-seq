#!/usr/bin/env Rscript

##########################################################################################
##########################################################################################

# Craig Bielski
# Pe'er Lab
# Detecting differential ATAC-seq peaks in CML DLI data

##########################################################################################
##########################################################################################

options(max.print = 100000)

'%!in%' <- function(x,y)!('%in%'(x,y))
'%!like%' <- function(x,y)!('%like%'(x,y))
luq <- function(x) {return(length(unique(x)))}

# auxiliary logging function
catverbose <- function(...){
  cat(format(Sys.time(), "%Y%m%d %H:%M:%S |"), ..., "\n")
}

# initialization
if('Rsubread' %!in% rownames(installed.packages())) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("Rsubread")
  require("Rsubread")
}

if('vsn' %!in% rownames(installed.packages())) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("vsn")
  require("vsn")
}

# required packages
pkgs = c('data.table',
         "readxl",
         'stringr',
         'plyr',
         'ggplot2',
         'gplots',
         'cowplot',
         'RColorBrewer',
         'viridis',
         'Rsubread',
         'DESeq2',
         'vsn',
         'fgsea',
         'biomaRt')

tmp <- lapply(pkgs, function (x) {
  suppressPackageStartupMessages(require(x, character.only = TRUE))
})
rm(tmp)

##########################################################################################
##########################################################################################

atac_samples <- as.data.table(readxl::read_excel("ATACseq_Samples.xlsx"))
bam_files <- as.data.table(list.files(pattern = ".bam$", recursive = T))
setnames(bam_files, "path")

# mapping
bam_files[, Sample_ID := sapply(str_split(path, "/"), "[[", 2)]
bam_files[, Cell_ID := sapply(str_split(path, "/"), "[[", 3)]
bam_files[, DFCI_ID := atac_samples$DFCI_ID[match(Sample_ID, atac_samples$Sample_ID)]]
bam_files[, Response := atac_samples$Response[match(Sample_ID, atac_samples$Sample_ID)]]
bam_files[, Condition := paste0(atac_samples$Timing[match(Sample_ID, atac_samples$Sample_ID)], "-Treatment")]
# bam_files$Cell <- "CD4"
bam_files$Cell <- "CD8"

hg38_genes <- fread("hg38_genes.gtf")
setnames(
  hg38_genes,
  c(
    "Chr",
    "Source",
    "Feature",
    "Start",
    "End",
    "Score",
    "Strand",
    "Phase",
    "Info"
  )
)
hg38_genes[, Gene := gsub("gene_name ", "", sapply(str_split(Info, ";"), "[[", 4))]
hg38_genes[, Gene := gsub('[\"]', "", Gene)]
hg38_genes[, Start := Start - 10000]
hg38_genes[, End := End + 10000]
hg38_genes <- hg38_genes[, .(Chr, Start, End, Gene)]
setkey(hg38_genes, Chr, Start, End)
  
ENCODE_blacklist <- fread("ENCODE_blacklist.bed")
DAC_blacklist <- fread("consensusBlacklist.bed")
DAC_blacklist[, c('V5', 'V6') := NULL]
setnames(DAC_blacklist, c("Chr", "Start", "End", "Info"))

setkey(DAC_blacklist, Chr, Start, End)

##########################################################################################
##########################################################################################

# merged bed file
merged_bed <- fread("all_peaks.merge.sort.bed")
setnames(merged_bed, c("chrom", "start", "end"))

# merged SAF
merged_saf <- fread("all_peaks.merge.sort.saf")
merged_saf[, Coord_ID := paste0(Chr, ":", Start, "-", End)]
setkey(merged_saf, Chr, Start, End)

fo <- foverlaps(merged_saf[, .(Chr, Start, End)],
                DAC_blacklist,
                type = "any")
fo[, coord_id := paste0(Chr, ":", i.Start, "-", i.End)]

# subset bam files
# cd4_bamfls <- bam_files[Cell %like% 'CD4']$path
cd8_bamfls <- bam_files[Cell %like% 'CD8']$path

# generate CD4 featurecounts
cd4_peak_fc = featureCounts(
  files = cd4_bamfls,
  annot.inbuilt = 'hg38',
  annot.ext = merged_saf,
  isPairedEnd = TRUE,
  allowMultiOverlap = TRUE
)

# generate CD8 featurecounts
cd8_peak_fc = featureCounts(
  files = cd8_bamfls,
  annot.inbuilt = 'hg38',
  annot.ext = merged_saf,
  isPairedEnd = TRUE,
  allowMultiOverlap = TRUE
)

cd4_counts <- cd4_peak_fc$counts
cd8_counts <- cd8_peak_fc$counts

##########################################################################################
##########################################################################################

# DESeq

# CD4 data
cd4_metadata <-
  data.frame(Condition = bam_files[Cell %like% 'CD4']$Condition,
             Response = bam_files[Cell %like% 'CD4']$Response,
             row.names = colnames(cd4_counts))

# CD8 data
cd8_metadata <-
  data.frame(Condition = bam_files[Cell %like% 'CD8']$Condition,
             Response = bam_files[Cell %like% 'CD8']$Response,
             row.names = colnames(cd8_counts))

##########################################################################################
##########################################################################################

# CD4

cd4_all_dds <-
  DESeqDataSetFromMatrix(cd4_counts[],
                         cd4_metadata[],
                         design =  ~ Response)

cd4_all_dds <-
  DESeq(
    cd4_all_dds,
    test = "LRT",
    full = ~ Response,
    reduced = ~ 1,
    minReplicatesForReplace = Inf
  ) # as per DESEeq2 documentation

cd4_pre_dds <-
  DESeqDataSetFromMatrix(cd4_counts[, which(cd4_metadata$Condition == "Pre-Treatment")], 
                         cd4_metadata[which(cd4_metadata$Condition == "Pre-Treatment"), ], 
                         design =  ~ Response)
cd4_pre_dds <-
  DESeq(
    cd4_pre_dds,
    test = "LRT",
    full = ~ Response,
    reduced = ~ 1,
    minReplicatesForReplace = Inf
  ) # as per DESEeq2 documentation

cd4_r_dds <-
  DESeqDataSetFromMatrix(cd4_counts[, which(cd4_metadata$Response == "R")], 
                         cd4_metadata[which(cd4_metadata$Response == "R"), ], 
                         design =  ~ Condition)
cd4_r_dds <-
  DESeq(
    cd4_r_dds,
    test = "LRT",
    full = ~ Condition,
    reduced = ~ 1,
    minReplicatesForReplace = Inf
  ) # as per DESEeq2 documentation

##########################################################################################
##########################################################################################

# custom sd vs. mean plot
my_msd_plot <- function(dds, dds_title = "") {
  # vsd::meanSdPlot
  plt_tmp = meanSdPlot(
    assay(dds),
    bins = 100,
    xlab = "Rank(mean)",
    ylab = "Standard deviation",
    plot = F
  )
  
  # customize
  plt_tmp <- plt_tmp$gg +
    # geom_smooth(color = "lightgrey") +
    scale_fill_viridis("Raw counts", limits = c(0, 150)) +
    ggtitle(dds_title) +
    # scale_fill_gradientn("Density",
    #                      limits = c(0, 150),
    #                      colours = c("black", "red", "orange", "yellow")) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      legend.position = 'bottom',
      title = element_text(size = 20),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 15),
      legend.text = element_text(size = 15)
    )
  
  return(plt_tmp)
  
}

##########################################################################################
##########################################################################################

# diagnostics
plotDispEsts(cd4_all_dds)
gg_pre <- suppressMessages(my_msd_plot(cd4_all_dds, dds_title = "Pre-normalization"))

# regularized log transform
# normalization with respect to library size, minimizes effects due to low counts, variance stabilization
cd4_rlog <- rlog(cd4_pre_dds)

gg_rlog <- suppressMessages(my_msd_plot(cd4_rlog, dds_title = "Regularized log transform"))

## log2(n + 1)
cd4_rtd <- normTransform(cd4_all_dds)
gg_rtd <- suppressMessages(my_msd_plot(cd4_rtd, dds_title = "log(count + 1)"))

## variance stabilizing transform
cd4_vsd <- varianceStabilizingTransformation(cd4_all_dds, blind = FALSE)
gg_vsf <- suppressMessages(my_msd_plot(cd4_vsd, dds_title = "Variance stabilization"))

## DESeq estimateSizeFactors
## estimates size factors via median ratio method (Anders & Huber 2010)
cd4_esf <- estimateSizeFactors(cd4_all_dds)
gg_esf <- suppressMessages(my_msd_plot(cd4_esf, dds_title = "Median ratio method"))

## peak normalization
## row-wise geometric means = 1
## http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# nz_idx = which(cd4_counts[, 1] > 0 & cd4_counts[, 2] > 0)
# normFactors <- cd4_counts[nz_idx, ] / exp(rowMeans(log(cd4_counts[nz_idx, ])))
# cd4_esfn <- estimateSizeFactors(cd4_all_dds[nz_idx, ], normMatrix = normFactors)
# gg_esfn <- suppressMessages(my_msd_plot(cd4_esfn, dds_title = "Median ratio method (normalized)"))

pdf("ATACseq_CD4_Normalization_Matrix.pdf", width = 16, height = 16, useDingbats = F)
plot_grid(gg_pre, gg_rlog, gg_rtd, gg_vsf, ncol = 2, align = 'vh')
dev.off()

##########################################################################################
##########################################################################################

cd4.pre_res <- results(cd4_pre_dds)
cd4.pre_res = cd4.pre_res[order(cd4.pre_res$pvalue),]

cd4.pre_resdata <-
  as.data.table(merge(
    as.data.frame(cd4.pre_res),
    as.data.frame(counts(cd4_pre_dds, normalized = TRUE)),
    by = "row.names",
    sort = FALSE
  ))
names(cd4.pre_resdata)[1] <- "Peak"
cd4.pre_resdata[, coord_id := merged_saf$Coord_ID[match(Peak, merged_saf$GeneID)]]
cd4.pre_resdata <- cd4.pre_resdata[coord_id %!in% fo[!is.na(Info)]$coord_id]
cd4.pre_resdata[, Chr := sapply(str_split(coord_id, "[:-]"), "[[", 1)]
cd4.pre_resdata[, Start := as.numeric(sapply(str_split(coord_id, "[:-]"), "[[", 2))]
cd4.pre_resdata[, End := as.numeric(sapply(str_split(coord_id, "[:-]"), "[[", 3))]
setkey(cd4.pre_resdata, Chr, Start, End)

cd4.pre_resdata <-
  foverlaps(cd4.pre_resdata, hg38_genes, type = "any")
cd4.pre_resdata[, c("Start", "End") := NULL]
setnames(
  cd4.pre_resdata,
  old = c("i.Start", "i.End"),
  new = c("Start", "End")
)
setkey(cd4.pre_resdata, Chr, Start, End)
cd4.pre_resdata = cd4.pre_resdata[order(padj)]

pdf("CD4_MA.pdf", width = 20, height = 14)
plotMA(cd4_pre_dds, alpha = 0.01, ylim = c(-3, 3))
dev.off()

cd4.r_res <- results(cd4_r_dds)
cd4.r_res = cd4.r_res[order(cd4.r_res$pvalue),]

cd4.r_resdata <-
  as.data.table(merge(
    as.data.frame(cd4.r_res),
    as.data.frame(counts(cd4_r_dds, normalized = TRUE)),
    by = "row.names",
    sort = FALSE
  ))
names(cd4.r_resdata)[1] <- "Peak"
cd4.r_resdata[, coord_id := merged_saf$Coord_ID[match(Peak, merged_saf$GeneID)]]
cd4.r_resdata[, Chr := sapply(str_split(coord_id, "[:-]"), "[[", 1)]
cd4.r_resdata[, Start := as.numeric(sapply(str_split(coord_id, "[:-]"), "[[", 2))]
cd4.r_resdata[, End := as.numeric(sapply(str_split(coord_id, "[:-]"), "[[", 3))]
setkey(cd4.r_resdata, Chr, Start, End)

cd4.r_resdata <-
  foverlaps(cd4.r_resdata, hg38_genes, type = "any")
cd4.r_resdata[, c("Start", "End") := NULL]
setnames(cd4.r_resdata,
         old = c("i.Start", "i.End"),
         new = c("Start", "End"))
setkey(cd4.r_resdata, Chr, Start, End)
cd4.r_resdata = cd4.r_resdata[order(padj)]

write.table(
  merged_saf[Coord_ID %in% cd4.r_resdata[pvalue < .05]$coord_id, .(Chr, Start, End)],
  "CD4_memory_PrePost.bed",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

cd4_r_bams = data.table(cd4_metadata[which(cd4_metadata$Response == "R"), ], keep.rownames = T)
setnames(cd4_r_bams, old = "rn", new = "path")
for (i in 1:3) {cd4_r_bams[, path := sub("\\.", "/", path)]}

write.table(
  cd4_r_bams,
  "CD4_memory_PrePost_bams.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

write.table(
  merged_saf[Coord_ID %in% cd4.pre_resdata[pvalue < .05]$coord_id, .(Chr, Start, End)],
  "CD4_memory_RnR.bed",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

cd4_pre_bams = data.table(cd4_metadata[which(cd4_metadata$Condition == "Pre-treatment"), ], keep.rownames = T)
setnames(cd4_pre_bams, old = "rn", new = "path")
for (i in 1:3) {cd4_pre_bams[, path := sub("\\.", "/", path)]}
  
write.table(
  cd4_pre_bams,
  "CD4_memory_RnR_bams.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

write.table(
  cd4.pre_resdata[, .(Gene, Peak, log2FoldChange, stat, pvalue, Peak_ID = coord_id)],
  "CD4_memory_RnR_diffPeaks.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

write.table(
  cd4.r_resdata[, .(Gene, Peak, log2FoldChange, stat, pvalue, Peak_ID = coord_id)],
  "CD4_memory_PrePost_diffPeaks.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)


##########################################################################################
##########################################################################################

# CD4

cd8_all_dds <-
  DESeqDataSetFromMatrix(cd8_counts[],
                         cd8_metadata[],
                         design =  ~ Response)

cd8_all_dds <-
  DESeq(
    cd8_all_dds,
    test = "LRT",
    full = ~ Response,
    reduced = ~ 1,
    minReplicatesForReplace = Inf
  ) # as per DESEeq2 documentation

cd8_pre_dds <-
  DESeqDataSetFromMatrix(cd8_counts[, which(cd8_metadata$Condition == "Pre-Treatment")], 
                         cd8_metadata[which(cd8_metadata$Condition == "Pre-Treatment"), ], 
                         design =  ~ Response)
cd8_pre_dds <-
  DESeq(
    cd8_pre_dds,
    test = "LRT",
    full = ~ Response,
    reduced = ~ 1,
    minReplicatesForReplace = Inf
  ) # as per DESEeq2 documentation

cd8_r_dds <-
  DESeqDataSetFromMatrix(cd8_counts[, which(cd8_metadata$Response == "R")], 
                         cd8_metadata[which(cd8_metadata$Response == "R"), ], 
                         design =  ~ Condition)
cd8_r_dds <-
  DESeq(
    cd8_r_dds,
    test = "LRT",
    full = ~ Condition,
    reduced = ~ 1,
    minReplicatesForReplace = Inf
  ) # as per DESEeq2 documentation

##########################################################################################
##########################################################################################

# diagnostics
plotDispEsts(cd8_all_dds)
gg_pre <- suppressMessages(my_msd_plot(cd8_all_dds, dds_title = "Pre-normalization"))

# regularized log transform
# normalization with respect to library size, minimizes effects due to low counts, variance stabilization
cd8_rlog <- rlog(cd8_pre_dds)

gg_rlog <- suppressMessages(my_msd_plot(cd8_rlog, dds_title = "Regularized log transform"))

## log2(n + 1)
cd8_rtd <- normTransform(cd8_all_dds)
gg_rtd <- suppressMessages(my_msd_plot(cd8_rtd, dds_title = "log(count + 1)"))

## variance stabilizing transform
cd8_vsd <- varianceStabilizingTransformation(cd8_all_dds, blind = FALSE)
gg_vsf <- suppressMessages(my_msd_plot(cd8_vsd, dds_title = "Variance stabilization"))

## DESeq estimateSizeFactors
## estimates size factors via median ratio method (Anders & Huber 2010)
cd8_esf <- estimateSizeFactors(cd8_all_dds)
gg_esf <- suppressMessages(my_msd_plot(cd8_esf, dds_title = "Median ratio method"))

## peak normalization
## row-wise geometric means = 1
## http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# nz_idx = which(cd8_counts[, 1] > 0 & cd8_counts[, 2] > 0)
# normFactors <- cd8_counts[nz_idx, ] / exp(rowMeans(log(cd8_counts[nz_idx, ])))
# cd8_esfn <- estimateSizeFactors(cd8_dds[nz_idx, ], normMatrix = normFactors)
# gg_esfn <- suppressMessages(my_msd_plot(cd8_esfn, dds_title = "Median ratio method (normalized)"))

pdf("ATACseq_CD8_Normalization_Matrix.pdf", width = 16, height = 16, useDingbats = F)
plot_grid(gg_pre, gg_rlog, gg_rtd, gg_vsf, ncol = 2, align = 'vh')
dev.off()

##########################################################################################
##########################################################################################

cd8.pre_res <- results(cd8_pre_dds)
cd8.pre_res = cd8.pre_res[order(cd8.pre_res$pvalue),]

cd8.pre_resdata <-
  as.data.table(merge(
    as.data.frame(cd8.pre_res),
    as.data.frame(counts(cd8_pre_dds, normalized = TRUE)),
    by = "row.names",
    sort = FALSE
  ))
names(cd8.pre_resdata)[1] <- "Peak"
cd8.pre_resdata[, coord_id := merged_saf$Coord_ID[match(Peak, merged_saf$GeneID)]]
cd8.pre_resdata <- cd8.pre_resdata[coord_id %!in% fo[!is.na(Info)]$coord_id]
cd8.pre_resdata[, Chr := sapply(str_split(coord_id, "[:-]"), "[[", 1)]
cd8.pre_resdata[, Start := as.numeric(sapply(str_split(coord_id, "[:-]"), "[[", 2))]
cd8.pre_resdata[, End := as.numeric(sapply(str_split(coord_id, "[:-]"), "[[", 3))]
setkey(cd8.pre_resdata, Chr, Start, End)

cd8.pre_resdata <-
  foverlaps(cd8.pre_resdata, hg38_genes, type = "any")
cd8.pre_resdata[, c("Start", "End") := NULL]
setnames(
  cd8.pre_resdata,
  old = c("i.Start", "i.End"),
  new = c("Start", "End")
)
setkey(cd8.pre_resdata, Chr, Start, End)
cd8.pre_resdata = cd8.pre_resdata[order(padj)]

pdf("CD8_MA.pdf", width = 20, height = 18)
plotMA(cd8_pre_dds, alpha = 0.01, ylim = c(-3, 3))
dev.off()

cd8.r_res <- results(cd8_r_dds)
cd8.r_res = cd8.r_res[order(cd8.r_res$pvalue),]

cd8.r_resdata <-
  as.data.table(merge(
    as.data.frame(cd8.r_res),
    as.data.frame(counts(cd8_r_dds, normalized = TRUE)),
    by = "row.names",
    sort = FALSE
  ))
names(cd8.r_resdata)[1] <- "Peak"
cd8.r_resdata[, coord_id := merged_saf$Coord_ID[match(Peak, merged_saf$GeneID)]]
cd8.r_resdata[, Chr := sapply(str_split(coord_id, "[:-]"), "[[", 1)]
cd8.r_resdata[, Start := as.numeric(sapply(str_split(coord_id, "[:-]"), "[[", 2))]
cd8.r_resdata[, End := as.numeric(sapply(str_split(coord_id, "[:-]"), "[[", 3))]
setkey(cd8.r_resdata, Chr, Start, End)

cd8.r_resdata <-
  foverlaps(cd8.r_resdata, hg38_genes, type = "any")
cd8.r_resdata[, c("Start", "End") := NULL]
setnames(cd8.r_resdata,
         old = c("i.Start", "i.End"),
         new = c("Start", "End"))
setkey(cd8.r_resdata, Chr, Start, End)
cd8.r_resdata = cd8.r_resdata[order(padj)]

write.table(
  merged_saf[Coord_ID %in% cd8.r_resdata[pvalue < .05]$coord_id, .(Chr, Start, End)],
  "CD8_memory_PrePost.bed",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

cd8_r_bams = data.table(cd8_metadata[which(cd8_metadata$Response == "R"), ], keep.rownames = T)
setnames(cd8_r_bams, old = "rn", new = "path")
for (i in 1:3) {cd8_r_bams[, path := sub("\\.", "/", path)]}

write.table(
  cd8_r_bams,
  "CD8_memory_PrePost_bams.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

write.table(
  merged_saf[Coord_ID %in% cd8.pre_resdata[pvalue < .05]$coord_id, .(Chr, Start, End)],
  "CD8_memory_RnR.bed",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

cd8_pre_bams = data.table(cd8_metadata[which(cd8_metadata$Condition == "Pre-treatment"), ], keep.rownames = T)
setnames(cd8_pre_bams, old = "rn", new = "path")
for (i in 1:3) {cd8_pre_bams[, path := sub("\\.", "/", path)]}

write.table(
  cd8_pre_bams,
  "CD8_memory_RnR_bams.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

write.table(
  cd8.pre_resdata[, .(Gene, Peak, log2FoldChange, stat, pvalue, Peak_ID = coord_id)],
  "CD8_memory_RnR_diffPeaks.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

write.table(
  cd8.r_resdata[, .(Gene, Peak, log2FoldChange, stat, pvalue, Peak_ID = coord_id)],
  "CD8_memory_PrePost_diffPeaks.txt",
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T
)

##########################################################################################
##########################################################################################

cd4_count_mx = as.data.frame(assay(cd4_rtd))
cd4_count_mx = cd4_count_mx[which(rownames(cd4_count_mx) %!in% merged_saf[Chr == "chrY"]$GeneID), ]

# cd4_tsne <- Rtsne(unique(scale(cd4_count_mx)), perplexity = 30, max_iter = 300)
# plot(cd4_tsne$Y[, 1], cd4_tsne$Y[, 2])

mycols <- viridis_pal(option = 'A')(9)[c(2, 7)]
pdf("CD4_PCA.pdf", width = 14, height = 14, useDingbats = F)
cd4_pca_plot <-
  plotPCA(cd4_rtd, intgroup = "Condition") + 
  geom_point(size = 5) +
  theme_classic() + 
  scale_color_manual(values = mycols) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 25))
cd4_pca_plot
dev.off()

pdf("CD4_Volcano.pdf", width = 14, height = 20)
ggplot(cd4.pre_resdata[baseMean > 10], aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(data = subset(cd4.pre_resdata[baseMean > 10], !(abs(log2FoldChange) > 1 &
                                                               -log10(pvalue) > 2)),
             color = "lightgrey",
             size = 3,
             alpha = 0.5) +
  geom_point(
    data = subset(cd4.pre_resdata[baseMean > 10], (abs(log2FoldChange) > 1 &
                                                     -log10(pvalue) > 2)),
    color = brewer.pal(n = 9, "OrRd")[7],
    size = 3,
    alpha = 0.8
  ) +
  xlab("log2(fold change)") +
  ylab("-log10(nominal pvalue)") +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 25)
  )
dev.off()

cd4.pca <- prcomp(t(cd4_count_mx), scale = T)
cd4.eigen.dat <- data.table(Variance = cd4.pca$sdev^2 / sum(cd4.pca$sdev^2))
cd4.eigen.dat[, PC := rank(order(Variance, decreasing = T))]

cd4.pca.dat <- as.data.table(cd4.pca$x[, 1:2], keep.rownames = T)
cd4.pca.dat[, PC1 := scale(PC1)]
cd4.pca.dat[, PC2 := scale(PC2)]
cd4.pca.dat[, Condition := cd4_metadata$Condition[match(rn, rownames(cd4_metadata))]]
cd4.pca.dat[, Response := cd4_metadata$Response[match(rn, rownames(cd4_metadata))]]

ggplot(cd4.pca.dat, aes(
  x = PC1,
  y = PC2,
  shape = Condition,
  color = Response
)) + 
  geom_point(size = 5) + 
  xlab(paste0("PC1: ", 100*round(cd4.eigen.dat$Variance[1], digits = 2), "% of Variance Explained")) + 
  ylab(paste0("PC2: ", 100*round(cd4.eigen.dat$Variance[2], digits = 2), "% of Variance Explained")) + 
  theme_classic() + 
  scale_color_viridis(discrete = T) +
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18), 
        legend.text = element_text(size = 18)) +
  guides(color = FALSE)

distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method="average")

zClust <-
  function(x,
           scale = "row",
           zlim = c(-3, 3),
           method = "average") {
    if (scale == "row")
      z <- t(scale(t(x)))
    if (scale == "col")
      z <- scale(x)
    z <- pmin(pmax(z, zlim[1]), zlim[2])
    hcl_row <- hclust(distCor(z), method = method)
    hcl_col <- hclust(distCor(t(z)), method = method)
    return(list(
      data = z,
      Rowv = as.dendrogram(hcl_row),
      Colv = as.dendrogram(hcl_col)
    ))
  }

cd4_rowmeans <- rowMeans(cd4_count_mx)
cd4_count_mx_filtered <- cd4_count_mx[which(cd4_rowmeans > 2), ]
rowVars <- apply(cd4_count_mx_filtered, 1, var)
cd4.mat <- zClust(cd4_count_mx_filtered[order(-rowVars)[1:2000], ], scale = "col")

side_cols <- cd4_metadata$Response
side_cols <- unlist(lapply(mapvalues(side_cols, from = c("NR", "R"), to = c("lightblue", "darkblue")), as.character))

pdf("CD4_memory_clust_new.pdf", width = 14, height = 14, useDingbats = F)
heatmap.2(
  cd4.mat$data,
  trace = 'none',
  col = viridis(n=299, option = 'A'),
  symm = F,
  symkey = F,
  symbreaks = F,
  Rowv = cd4.mat$Rowv,
  Colv = cd4.mat$Colv,
  dendrogram = 'both',
  labRow = "",
  ColSideColors = side_cols
)
dev.off()

##########################################################################################
##########################################################################################

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapTab <-
  getBM(
    attributes = c("hgnc_symbol", "entrezgene"),
    filters = "hgnc_symbol",
    values = resdata$Gene,
    mart = ensembl,
    uniqueRows = FALSE
  )
dupRows <-
  union(which(duplicated(mapTab[, 1])), which(duplicated(mapTab[, 2])))
mapTab <- as.data.table(mapTab[-dupRows,])

cd4.pre_resdata[, entrez_id := mapTab$entrezgene[match(Gene, mapTab$hgnc_symbol)]]
res.ranks <- cd4.pre_resdata[!is.na(stat)][order(stat)]$stat
names(res.ranks) <- cd4.pre_resdata[!is.na(stat)][order(stat)]$Gene

res.noise <- rnorm(sum(duplicated(res.ranks)), .01, .1)
res.ranks[duplicated(res.ranks)] <-
  res.ranks[duplicated(res.ranks)] + res.noise

# msigdb.all <- gmt2list("c2.all.v6.1.symbols.gmt")
msigdb.canon <- gmt2list("c2.cp.v6.1.symbols.gmt")
fgseaRes <- fgsea(
  pathways = msigdb.canon,
  stats = res.ranks,
  minSize = 15,
  maxSize = 500,
  nperm = 10000
)

fgseaRes = fgseaRes[order(pval)]

# plotEnrichment(msigdb.canon[[""]], res.ranks) + labs(title="")

##########################################################################################
##########################################################################################

cd8_count_mx = as.data.frame(assay(cd8_rtd))
cd8_count_mx = cd8_count_mx[which(rownames(cd8_count_mx) %!in% merged_saf[Chr == "chrY"]$GeneID), ]

# cd8_tsne <- Rtsne(unique(scale(cd8_count_mx)), perplexity = 30, max_iter = 300)
# plot(cd8_tsne$Y[, 1], cd8_tsne$Y[, 2])

mycols <- viridis_pal(option = 'A')(9)[c(2, 7)]
pdf("CD8_PCA.pdf", width = 14, height = 14, useDingbats = F)
cd8_pca_plot <-
  plotPCA(cd8_rtd, intgroup = "Condition") + 
  geom_point(size = 5) +
  theme_classic() + 
  scale_color_manual(values = mycols) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 25))
cd8_pca_plot
dev.off()

pdf("CD8_Volcano.pdf", width = 18, height = 20)
ggplot(cd8.pre_resdata[baseMean > 10], aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(data = subset(cd8.pre_resdata[baseMean > 10], !(abs(log2FoldChange) > 1 &
                                                               -log10(pvalue) > 2)),
             color = "lightgrey",
             size = 3,
             alpha = 0.5) +
  geom_point(
    data = subset(cd8.pre_resdata[baseMean > 10], (abs(log2FoldChange) > 1 &
                                                     -log10(pvalue) > 2)),
    color = brewer.pal(n = 9, "OrRd")[7],
    size = 3,
    alpha = 0.8
  ) +
  xlab("log2(fold change)") +
  ylab("-log10(nominal pvalue)") +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 25)
  )
dev.off()

cd8.pca <- prcomp(t(cd8_count_mx), scale = T)
cd8.eigen.dat <- data.table(Variance = cd8.pca$sdev^2 / sum(cd8.pca$sdev^2))
cd8.eigen.dat[, PC := rank(order(Variance, decreasing = T))]

cd8.pca.dat <- as.data.table(cd8.pca$x[, 1:2], keep.rownames = T)
cd8.pca.dat[, PC1 := scale(PC1)]
cd8.pca.dat[, PC2 := scale(PC2)]
cd8.pca.dat[, Condition := cd8_metadata$Condition[match(rn, rownames(cd8_metadata))]]
cd8.pca.dat[, Response := cd8_metadata$Response[match(rn, rownames(cd8_metadata))]]

ggplot(cd8.pca.dat, aes(
  x = PC1,
  y = PC2,
  shape = Condition,
  color = Response
)) + 
  geom_point(size = 5) + 
  xlab(paste0("PC1: ", 100*round(cd8.eigen.dat$Variance[1], digits = 2), "% of Variance Explained")) + 
  ylab(paste0("PC2: ", 100*round(cd8.eigen.dat$Variance[2], digits = 2), "% of Variance Explained")) + 
  theme_classic() + 
  scale_color_viridis(discrete = T) +
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18), 
        legend.text = element_text(size = 18)) +
  guides(color = FALSE)

distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method="average")

zClust <-
  function(x,
           scale = "row",
           zlim = c(-3, 3),
           method = "average") {
    if (scale == "row")
      z <- t(scale(t(x)))
    if (scale == "col")
      z <- scale(x)
    z <- pmin(pmax(z, zlim[1]), zlim[2])
    hcl_row <- hclust(distCor(z), method = method)
    hcl_col <- hclust(distCor(t(z)), method = method)
    return(list(
      data = z,
      Rowv = as.dendrogram(hcl_row),
      Colv = as.dendrogram(hcl_col)
    ))
  }

cd8_rowmeans <- rowMeans(cd8_count_mx)
cd8_count_mx_filtered <- cd8_count_mx[which(cd8_rowmeans > 2), ]
rowVars <- apply(cd8_count_mx_filtered, 1, var)
cd8.mat <- zClust(cd8_count_mx_filtered[order(-rowVars)[1:2000], ], scale = "col")

side_cols <- cd8_metadata$Response
side_cols <- unlist(lapply(mapvalues(side_cols, from = c("NR", "R"), to = c("lightblue", "darkblue")), as.character))

pdf("CD8_memory_clust_new.pdf", width = 14, height = 14, useDingbats = F)
heatmap.2(
  cd8.mat$data,
  trace = 'none',
  col = viridis(n=299, option = 'A'),
  symm = F,
  symkey = F,
  symbreaks = F,
  Rowv = cd8.mat$Rowv,
  Colv = cd8.mat$Colv,
  dendrogram = 'both',
  labRow = "",
  ColSideColors = side_cols
)
dev.off()


