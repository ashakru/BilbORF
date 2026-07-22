library(BSgenome.Hsapiens.UCSC.hg38)
library(devtools)
library(rtracklayer)
library(dplyr)
library(ORFik)
library(GenomicFeatures)
library(tidyverse)

BSgenome <- BSgenome.Hsapiens.UCSC.hg38
gtf <- "inst/extdata/gencode.v35.annotation_chr10.gtf"
bed <- "inst/extdata/Ribo-seq_ORFs.bed"
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf)

# Prepare annotations
annotations <- prepare_annotations_fromGTF(gtf, BSgenome)

# Pull transcripts metadata
transcripts_meta <- import(gtf, format = "GTF") %>%
  as.data.frame() %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(gene_name, transcript_id, transcript_type, gene_id) %>%
  dplyr::distinct()

# Load ORFs
orfs <- import(bed, format = "BED")
names(orfs) <- orfs$name
orf_tab <- as.data.frame(orfs)

# Annotate ORFs
annotated_orfs <- annotate_orf_isoforms(annotations, orfs, BSgenome, transcripts_meta)

# Count translation in selected ORFs (not isoform strict currently)
offsets <- data.frame(fraction = c(25:30),
                      offsets_start = -12)
bam <- "inst/extdata/riboseq_chr10.bam"
orf_translation <- count_p_sites(bam, offsets, annotated_orfs)

# DRIMSeq quantification
library(DRIMSeq)
gtex_tx <- read_delim("inst/extdata/quantification_gencode.counts.txt")
gtex_meta <- read_csv("inst/extdata/GTEx_longReads_Seq_meta.csv") %>%
  dplyr::filter(sample_id %in% colnames(gtex_tx))
all(colnames(gtex_tx)[-1] %in% gtex_meta$sample_id)

transcripts_table <- read_csv("inst/extdata/gencode.v35.tx2gene.csv")
transcripts_table <- transcripts_table %>%
  dplyr::mutate(transcript_ensembl_id = sapply(strsplit(transcript_id, "[.]"), function(x){x[1]}))

heart_meta <- gtex_meta %>%
  dplyr::filter(grepl("Heart", tissue)) %>%
  mutate(group = factor(case_when(grepl("Ventricle", tissue) ~ "LV",
                                  TRUE ~ "AA")))  %>%
  as.data.frame()
heart_counts <- gtex_tx %>%
  dplyr::select(transcript, any_of(heart_meta$sample_id)) %>%
  rename(feature_id = transcript) %>%
  mutate(feature_id = sapply(strsplit(feature_id, "[.]"), function(x){x[1]})) %>%
  left_join(dplyr::select(transcripts_table, gene_id, transcript_ensembl_id),
            by = c("feature_id" = "transcript_ensembl_id")) %>%
  relocate(gene_id) %>%
  mutate(feature_id = make.unique(feature_id)) %>%
  as.data.frame() %>%
  dplyr::filter(!is.na(gene_id))


all(heart_meta$sample_id == colnames(heart_counts)[-1:-2])
levels(heart_meta$group)

d <- dmDSdata(counts = heart_counts, samples = heart_meta)
d <- dmFilter(d, min_samps_gene_expr = 7, min_samps_feature_expr = 3,
              min_gene_expr = 10, min_feature_expr = 10)
design_full <- model.matrix(~ group, data = samples(d))
set.seed(123)
d <- dmPrecision(d, design = design_full)
head(mean_expression(d), 3)
plotPrecision(d)
d <- dmFit(d, design = design_full, verbose = 1)
d <- dmTest(d, coef = "groupLV", verbose = 1)
drimseq_results <- results(d)
drimseq_results <- drimseq_results[order(drimseq_results$pvalue, decreasing = FALSE), ]
top_gene_id <- drimseq_results$gene_id[1]
drimseq_results_tx <- results(d, level = "feature")
drimseq_props_tx <- proportions(d)

annotated_orfs_tab <- annotated_orfs$table %>%
  mutate(feature_id = sapply(strsplit(transcript_id, "[.]"), function(x){x[1]}))

# Find differen
drimseq_results_tx <- drimseq_results_tx %>%
  mutate(signif = adj_pvalue < 0.05,
         hasORF = feature_id %in% annotated_orfs_tab$feature_id)

####### Important ORF examples #######
# c10riboseqorf127 ORF with carying frame translation between transcript isoforms








####### SANDBOX #######

# Debugging Kozak sequence function
seqName <- "Chromosome"
ORF1 <- GRanges(seqnames = seqName,
                ranges = IRanges(c(1007, 1096), width = 60),
                strand = c("+", "+"))
ORF2 <- GRanges(seqnames = seqName,
                ranges = IRanges(c(400, 100), width = 30),
                strand = c("-", "-"))
ORFs <- GRangesList(tx1 = ORF1, tx2 = ORF2)
ORFs <- makeORFNames(ORFs) # need ORF names
tx <- extendLeaders(ORFs, 100)

faFile <- FaFile(system.file("extdata/Danio_rerio_sample", "genome_dummy.fasta", package = "ORFik"))
kozakSequenceScore(ORFs, tx, faFile)

cds <- cdsBy(txdb, by = "tx", use.names = TRUE)[1:10]
tx <- exonsBy(txdb, by = "tx", use.names = TRUE)[names(cds)]
faFile <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens

kozakSequenceScore(cds, tx, faFile, species = "human")


x <- GRanges("chr1", IRanges(start = c(26, 29), end = c(27, 29)), "+")
names(x) <- rep("tx1_ORF1", length(x))
x <- groupGRangesBy(x)
# tx is the whole region
tx_gr <- GRanges("chr1", IRanges(c(5, 29), c(27, 30)), "+")
names(tx_gr) <- rep("tx1", length(tx_gr))
tx <- groupGRangesBy(tx_gr)
pmapToTranscriptF(x, tx)
