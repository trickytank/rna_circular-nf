# Prepare GTF data for only chromosome 1
grep -P '^chr1\t' data/reference/hg19_ens.gtf > data/reference/hg19_ens_chr1.gtf
