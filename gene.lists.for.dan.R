setwd("/")
# reading in Dan's deleterious variants
df <- readr::read_tsv("camp/lab/swantonc/working/blackj/Germline_project/tx421.germline.deleterious.tsv")
# Rahman (2015) list of cancer predisposition genes
load("camp/lab/swantonc/working/blackj/rahman_germline_gene_list.RData")
# How many of Dan's snps are within cancer predisposition genes
df$is_germline_cancer_gene[df$gene_refgene %in% germline_gene_list$`Gene  Symbol` | df$gene_refgene == "ERBB2"] <- "Predisposition_gene"
pattern <- "(.+,)(.+,)(.+,)(.+,)(.+)"
# Separating clinvar annotation - need to talk to Dan about this
df$clinvar_pathogenicity <- gsub(pattern = pattern, paste0("\\5"), df$clinvar_20190305)
# Variants that are pathogenic and within germline cancer genes
df$is_cancer_germline_variant[df$is_germline_cancer_gene == "Predisposition_gene" & df$clinvar_pathogenicity %in% 
                                c("Likely_pathogenic", "Pathogenic", "Pathogenic/Likely_pathogenic")] <- "Cancer_germline_variant"

# gnomad 'pLoF Metrics by Gene TSV' downloaded from https://gnomad.broadinstitute.org/downloads
gnomad_lof <- readr::read_tsv("camp/lab/swantonc/working/blackj/gnomad.v2.1.1.lof_metrics.by_gene.txt")
# I think this is LOEUF column
gnomad_to_merge <- gnomad_lof[,c("gene", "oe_lof_upper")]
df_merged <- merge(df, gnomad_to_merge, by.x = "gene_refgene", by.y = "gene")

# this now contains gnomad and Rahman data
save(df_merged, file = "camp/lab/swantonc/working/blackj/Germline_project/all.421.germline_jb_updated.RData")

load("camp/lab/swantonc/working/blackj/Germline_project/Ding_variants")
ding_variants$gene_start <- paste0(ding_variants$HUGO_Symbol, ding_variants$Start)

df_merged$gene_start <- paste0(df_merged$gene_refgene, df_merged$start)
df_merged$ding_variant[df_merged$gene_start %in% ding_variants$gene_start] <- "Ding Variant"
save(df_merged, file = "camp/lab/swantonc/working/blackj/Germline_project/all.421.germline_jb_updated.RData")
