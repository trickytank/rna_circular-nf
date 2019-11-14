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
params.annotation   = "data/reference/hg19_ens_chr1.txt"
params.reads  = "data/mini/pos_mini_{1,2}.fastq.gz"

// Parse parameters
fasta_ref = file(params.fasta)
gtf_ref = file(params.gtf)
annotation_ref = file(params.annotation)


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


process 'genome_fai' {
  label 'star'

  input:
      path genome

  output:
      path "${genome}.fai"

  script:
  """
  samtools faidx ${genome}
  """
}



process create_star_index {
  label 'star'

  input:
      path genome // reference fasta
      path gtf // genes

  output:
      path "star_index_dir"

  script:
  """
  mkdir star_index_dir

  STAR --runMode genomeGenerate \
       --genomeDir star_index_dir \
       --genomeFastaFiles ${genome} \
       --sjdbGTFfile ${gtf} \
       --runThreadN ${task.cpus} \
       --genomeSAindexNbases 12
  ## TODO: allow adjustment of --genomeSAindexNbases outside
  """

}

process alignment_star {
  label 'star'
  tag "$sample_id"

  input:
    path star_index_dir
    tuple sample_id, path(reads)
    
  output:
    tuple sample_id,
      path("out/${sample_id}.Aligned.out.sam"),
      path("out/${sample_id}.Chimeric.out.junction")
    path("out/${sample_id}.Log.final.out")

  script:
  """
    mkdir out

    STAR \
        --chimSegmentMin 10 \
        --runThreadN ${task.cpus} \
        --genomeDir star_index_dir \
        --outSAMstrandField intronMotif \
        --readFilesCommand gunzip -c \
        --readFilesIn $reads \
        --outFileNamePrefix out/${sample_id}.
  """
  // TODO: allow RNA strandedness to be specified.
  // TODO (low priority): allow uncompressed fastq to work
  // The above is coded for unstranded RNA seq.
}

process circ_parse {
  label 'circexplorer'
  tag "$sample_id"

  input:
    path genome // reference fasta
    path genome_fai
    path annotation // gene annotation txt file
    tuple sample_id,
      path("${sample_id}.Aligned.out.sam"),
      path("${sample_id}.Chimeric.out.junction")

  output:
    tuple sample_id,
      path("${sample_id}.back_spliced_junction.bed")
      // Also outputs circularRNA_known.txt, may not be important...

  script:
  """
    fast_circ.py parse \
      -r $annotation \
      -g $genome \
      -t STAR \
      ${sample_id}.Chimeric.out.junction

      mv back_spliced_junction.bed "${sample_id}.back_spliced_junction.bed"
  """

}

process circ_denovo {
  label 'circexplorer'
  tag "$sample_id"

  input:
    path genome // reference fasta
    path genome_fai
    path annotation // gene annotation txt file
    path gtf // genes
    tuple sample_id, path(reads)

  output:
    path '*' 

  script:
  """
    fast_circ.py denovo \
      -r $annotation \
      -g $genome \
      -G $gtf \
      -p ${task.cpus} \
      -f ${reads[0]}
  """
}

/*
process multiqc {
  label 'seq_qc'

}
*/

workflow {
  reads_ch    = Channel.fromFilePairs(params.reads)

  // fetch_ref_genes(params.genome)
  // fetch_ref_fasta(params.genome)

  genome_fai(fasta_ref)

  create_star_index(fasta_ref, gtf_ref)

  alignment_star(
    create_star_index.out,
    reads_ch
  )
  circ_parse(
    fasta_ref,
    genome_fai.out,
    annotation_ref,
    alignment_star.out[0]
  )

  circ_denovo(
    fasta_ref,
    genome_fai.out,
    annotation_ref,
    gtf_ref,
    reads_ch
  )

  publish:
    circ_parse.out to: 'results/circ_parse',
    mode: 'link'
}
