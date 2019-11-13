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
       --runThreadN ${task.cpus} \
       --genomeSAindexNbases 12
  ## TODO: allow adjustment of --genomeSAindexNbases outside
  """

}

process 'alignment_star' {
  label 'star'
  tag 'TODO'

  input:
    path star_index_dir
    tuple sample_id, path(read1), path(read)
    
  output:
    tuple sample_id, path("out/$sample_id.Aligned.out.sam")
    tuple sample_id, path("out/$sample_id.Log.final.out")

  script:
  """
    mkdir out

    STAR \
        --chimSegmentMin 10 \
        --runThreadN ${task.cpus} \
        --genomeDir star_index_dir \
        --outSAMstrandField intronMotif \
        --readFilesIn $read1 $read2 \
        --outFileNamePrefix out/$sample_id.
  """
  // TODO: allow RNA strandedness to be specified.
  // The above is coded for unstranded RNA seq.
}
*/

workflow {
  //reads_ch    = Channel.fromFilePairs(params.reads)

  fetch_ref_genes(params.genome)
  fetch_ref_fasta(params.genome)
  
}
