library(SKAT)
library(fst)
# load SNP annotation table from Dan
annotation <- read_fst("~/Downloads/annotation.fst")
# load covariate data + phenotype (using test data for now)
example_pt_phenotype_binary <- c(1,0,0,0,0,1,0,1)
example_pt_phenotype_continuous <- c(0.1,1.1,0.5,0.5,0.2,0.8,1.2,1.0)
example_pt_covariate_1 <- c(0,0,0,0,0,0,1,1)
example_pt_covariate_2 <- c(-1.2,-3,0.1,5,1.2,3.1,1,0.8)

# load example data
data(SKAT.example)
# here, y.c is continuous phenotype and y.b is binary
# Z is genotype matrix (presumably, 2 is homozygous, 1 is heterozygous and 0 is absent)
# X is co-variate matrix (continuous)
attach(SKAT.example)

obj<-SKAT_Null_Model(y.c ~ X, out_type="C")
# there, out_type = C is for continuous (D for dichotomous)
SKAT(Z, obj)$p.value

# 1. Simple association testing
# 2. Assigning weights to SNPs
# 3. Combining burden testing and SKAT
### This is for optimal test: 
SKAT(Z, obj, method="SKATO")$p.value

# Step 1 - need to get SNP data into homo -, het, homo + (0,1,2)
annotation$SNP_code[annotation$vaf_gl <0.4] <- 0
annotation$SNP_code[annotation$vaf_gl >0.6] <- 2
annotation$SNP_code[!annotation$SNP_code %in% c(0,2)] <- 1

#Step 2 - create matrix with patients/regions as rows, and different SNPs as columns
annotation$chrpos <- paste0(annotation$chrom, annotation$pos)
annotation_matrix <- matrix(nrow = length(unique(annotation$sample_name_hash)), ncol = length(unique(annotation$chrpos)))
rownames(annotation_matrix) <- unique(annotation$sample_name_hash)
colnames(annotation_matrix) <- unique(annotation$chrpos)
for (i in 1:nrow(annotation_matrix)){
  for (j in 1:ncol(annotation_matrix)){
    annotation_matrix[i,j] <- annotation$SNP_code[annotation$sample_name_hash == rownames(annotation_matrix)[i] & annotation$chrpos == colnames(annotation_matrix)[j]]
  }
}

# Step 3 - make phenotype data matrix - can be binary or continuous
annotation_matrix_phenotype_binary <- matrix(nrow = length(unique(annotation$sample_name_hash)), ncol = 1)
rownames(annotation_matrix_phenotype_binary) <- unique(annotation$sample_name_hash)
annotation_matrix_phenotype_binary[1:8,1] <- example_pt_phenotype_binary
annotation_matrix_phenotype_continuous <- matrix(nrow = length(unique(annotation$sample_name_hash)), ncol = 1)
rownames(annotation_matrix_phenotype_continuous) <- unique(annotation$sample_name_hash)
annotation_matrix_phenotype_continuous[1:8,1] <- example_pt_phenotype_continuous

# Step 4 - make covariate matrix - here using 1 binary and 1 continuous
covariate_matrix <- matrix(nrow = length(unique(annotation$sample_name_hash)), ncol = 2)
covariate_matrix[,1] <- example_pt_covariate_1
covariate_matrix[,2] <- example_pt_covariate_2

# Step 5 - SKAT_Null_Model to estimate parameters under null model with no associations
# out_type = "C" for continuous, ("D" for dichotomous)
obj<-SKAT_Null_Model(annotation_matrix_phenotype_continuous ~ covariate_matrix, out_type="C")

# now run SKAT, or SKAT-O

# 1. Simple association testing
SKAT(annotation_matrix, obj)$p.value

# 2. Combining burden testing and SKAT - computes multiple different p values and uses minimum as test statistic
SKAT(annotation_matrix, obj, method = "SKATO")$p.value







