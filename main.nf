/* Run circular RNA seq detection and characterization
 *
 * Rick Tankard
 */


//
nextflow.preview.dsl=2

// Define the default parameters

params.genome = "hg19"
params.reads  = "../circRNA_detection_review/simu/pos_{1,2}.fastq.gz"



// Pipeline
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

/*
process 'create_star_index' {
  input:
      file genome from genome_file

  output:
      file "genome_dir" into star_index_ch

  script:
  """
  mkdir star_index_dir

  STAR --runMode genomeGenerate \
       --genomeDir genome_dir \
       --genomeFastaFiles ${genome} \
       --runThreadN ${task.cpus}
  """

}

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

  fetch_ref_genes(params.genome)
  fetch_ref_fasta(params.genome)
  
}
