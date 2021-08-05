setwd("/omics/groups/OE0246/internal/guz/cola_hc/examples/ReprocessedAllen")
library(cola)
library(scRNAseq)
data = readRDS('/omics/groups/OE0246/internal/guz/cola_hc/examples/ReprocessedAllen/ReprocessedAllen_data.rds')
mat = assays(data)$rsem_tpm
mat = log2(mat + 1)
mat = adjust_matrix(mat)

anno = colData(data)[, c("driver_1_s", "dissection_s", "Core.Type", "Primary.Type", "Secondary.Type")]
anno = as.data.frame(anno)

rh = hierarchical_partition(mat, subset = 500, cores = 4, anno = anno)
saveRDS(rh, file = "ReprocessedAllen_cola_rh.rds")

cola_report(rh, output = "ReprocessedAllen_cola_rh_report", title = "cola Report for Hierarchical Partitioning - 'ReprocessedAllen'")
