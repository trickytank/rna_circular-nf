/* Run circular RNA seq detection and characterization
 *
 * Rick Tankard
 */


//
nextflow.preview.dsl=2

// Define the default parameters

params.genome = "hg19"
params.fasta = "data/reference/chr1.fa"
params.gtf   = "data/reference/hg19_ens_chr1.gtf"
params.reads  = "../circRNA_detection_review/simu/pos_{1,2}.fastq.gz"

// Parse parameters
fasta_ref = file(params.fasta)
gtf_ref = file(params.gtf)


// Pipeline
/*
process fetch_ref_fasta {
  label 'circexplorer'

  input:
    val genome

  output:
    path("${genome}.fa")

  """
    # Download human reference genome sequence file
    fetch_ucsc.py ${genome} fa ${genome}.fa
  """
}

process fetch_ref_genes {
  label 'circexplorer'

  input:
    val genome

  output:
    path("${genome}_ens.gtf")
  """
    # Download human Ensembl gene annotation file
    fetch_ucsc.py ${genome} ens ${genome}_ens.txt

    # Convert gene annotation file to GTF format (require genePredToGtf)
    cut -f2-11 ${genome}_ens.txt|genePredToGtf file stdin ${genome}_ens.gtf

  """
}
*/


process create_star_index {
  label 'star'

  input:
      path genome // reference fasta
      path gtf // genes

  output:
      path "star_index_dir"

  script:
  """
  mkdir genome_dir

  STAR --runMode genomeGenerate \
       --genomeDir star_index_dir \
       --genomeFastaFiles ${genome} \
       --sjdbGTFfile ${gtf} \
       --runThreadN ${task.cpus}
  """

}

/*
process 'alignment_star' {

  script:
  """
    STAR \
        --chimSegmentMin 10 \
        --runThreadN ${task.cpus} \
        --genomeDir hg19_STAR_index \
        --readFilesIn read_1.fastq read_2.fastq
  """
}
*/

workflow {
  //reads_ch    = Channel.fromFilePairs(params.reads)

  // fetch_ref_genes(params.genome)
  // fetch_ref_fasta(params.genome)

  create_star_index(fasta_ref, gtf_ref)
  
}
