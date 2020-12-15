version 1.0

workflow debarcerWorkflow {
  input {
    String outdir = "./"    
    File bamFile
    File bedFile
    Int distance = 1
    Int position = 10
    String separator
    Int readCount = 0
    Boolean truncate = false
    Boolean ignoreOrphans = false
    Int maxDepth = 1000000
    String stepper = "nofilter"
    Int baseQuality = 25
    String familySize
    Float percentThreshold = 50 
    Int countThreshold = 1
    Float referenceThreshold = 95
    Float alternativeThreshold = 2
    Int filterThreshold = 10
    String extension = "png"
    Boolean report = true
    Int minCov = 1000
    Float minRatio = 0.1
    Int minUmis = 1000
    Int minChildren = 500
  }


  parameter_meta {





  }

  meta {
    author: "Richard Jovelin"
    email: "richard.jovelin@oicr.on.ca"
    description: "A package for de-barcoding and error correction of sequencing data containing molecular barcodes"
  }


  

  call assignSmmips {
    input:
      fastq1 = fastq1,
      fastq2 = fastq2,
      panel = panel,
      outdir = outdir,
      prefix = prefix,  
      maxSubs = maxSubs,
      upstreamNucleotides = upstreamNucleotides,
      umiLength = umiLength, 
      match = match,
      mismatch = mismatch,
      gapOpening = gapOpening,
      gapExtension = gapExtension,  
      alignmentOverlapThreshold = alignmentOverlapThreshold,
      matchesThreshold = matchesThreshold,
      remove = removeIntermediate
  }

  File assignedBam = assignSmmips.assignedBam 
  File assignedBamIndex = assignSmmips.assignedBamIndex 

  Boolean truncateColumn = if truncate then true else false
  Boolean ignoreOrphanReads = if ignoreOrphans then true else false

  call countVariants {
    input: 
      assignedBam = assignedBam,
      assignedBamIndex = assignedBamIndex,
      panel = panel,
      outdir = outdir,
      prefix = prefix,  
      truncate = truncateColumn,
      ignoreOrphans = ignoreOrphanReads,
      stepper = stepper,
      maxDepth = maxDepth,
      referenceName = referenceName,
      cosmicFile = cosmicFile
  }

  output {
    File outputExtractionMetrics = assignSmmips.extractionMetrics
    File outputReadCounts = assignSmmips.readCounts
    File outputSortedbam = assignSmmips.sortedbam
    File outputSortedbamIndex = assignSmmips.sortedbamIndex
    File outputAssignedBam = assignSmmips.assignedBam
    File outputAssignedBamIndex = assignSmmips.assignedBamIndex
    File outputUnassignedBam = assignSmmips.unassignedBam
    File outputUnassignedBamIndex = assignSmmips.unassignedBamIndex
    File outputEmptyBam = assignSmmips.emptyBam
    File outputEmptyBamIndex = assignSmmips.emptyBamIndex
    File outputCountTable = countVariants.countTable 
  }
}

task groupUmis {
  input {
    String modules = "debarcer/2.1.3"
    Int memory = 32
    Int timeout = 36
    String outdir = "./"    
    String region
    File bamFile
    Int distance
    Int position
    String separator
    Int readCount = 0
    Boolean truncate = false
    Boolean ignoreOrphans = false
  }

  command <<<
    set -euo pipefail
    debarcer group -o ~{outdir} -r ~{region} -b ~{bamFile} -d ~{distance} -p ~{position} -s ~{separator} -rc ~{readCount} -i ~{ignoreOrphans} -t ~{truncate}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
   File dataFile = "${outdir}/Datafiles/datafile_${region}.csv"
   File umiFamilies = "${outdir}/Umifiles/${region}.json"
   File umiRelationships = "${outdir}/Stats/UMI_relationships_${region}.txt"
   File mappedReadCounts = "${outdir}/Stats/Mapped_read_counts_${region}.json"
   File umis = "${outdir}/Stats/Umis_${region}_before_grouping.json"
  }
}


task collapseUmis {
  input {
    String modules = "debarcer/2.1.3  hg19/p13"
    Int memory = 32
    Int timeout = 36
    String outdir = "./"    
    File bamFile
    String region
    Int maxDepth = 1000000
    Boolean truncate = false
    Boolean ignoreOrphans = false
    String stepper = "nofilter"
    String separator
    Int baseQuality = 25
    String familySize
    File umiFile
    Float percentThreshold = 50 
    Int countThreshold = 1
    Int positionThreshold = 10
    File refFasta = $HG19_ROOT/hg19_random.fa
    File refDict = $HG19_ROOT/hg19_random.dict
    File refIndex = $HG19_ROOT/hg19_random.fa.fai
  }

  command <<<
    set -euo pipefail
    cp ~{refDict} .
    cp ~{refIndex} .
    debarcer collapse -o ~{outdir} -b ~{bamFile} -rf ~{refFasta} -r ~{region} -u ~{umiFile} -f ~{familySize} -ct ~{countThreshold} -pt ~{percentThreshold} -p ~{positionThreshold} -m ~{maxDepth} -t ~{truncate} -i ~{ignoreOrphans} -stp ~{stepper} -s ~{separator} -bq ~{baseQuality}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File coverage = "${outdir}/Stats/CoverageStats.yml"
    File consensus = "${outdir}/Consfiles/${region}.cons"
  }
}



task callVariants {
  input {
    String modules = "debarcer/2.1.3  hg19/p13"
    Int memory = 32
    Int timeout = 36
    String outdir = "./"    
    File refFasta = $HG19_ROOT/hg19_random.fa
    File refDict = $HG19_ROOT/hg19_random.dict
    File refIndex = $HG19_ROOT/hg19_random.fa.fai
    Float referenceThreshold = 95
    Float alternativeThreshold = 2
    Int filterThreshold = 10
    String familySize
  }

  command <<<
    set -euo pipefail
    cp ~{refDict} .
    cp ~{refIndex} .
    debarcer call -o ~{outdir} -rf ~{refFasta} -rt ~{referenceThreshold} -at ~{alternativeThreshold} -ft ~{filterThreshold} -f ~{familySize}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    Array[File] vcfFiles = glob("${outdir}/VCFfiles/*.vcf")
  }
}


task graph {
  input {
    String modules = "debarcer/2.1.3"
    Int memory = 32
    Int timeout = 36
    String outdir = "./"
    String extension = "png"
    String sample = sampleName
    Boolean report = true
    Int minCov = 1000
    Float minRatio = 0.1
    Int minUmis = 1000
    Int minChildren = 500
    Float refThreshold = 95
  }    

  command <<<
    set -euo pipefail
    debarcer plot -d ~{outdir} -e ~{extension} -s ~{sample} -r ~{report} -mv ~{minCov} -mr ~{minRatio} -mu ~{minUmis} -mc ~{minChildren} -rt ~{refThreshold}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    Array[File] figureFiles = glob("${outdir}/Figures/*.png")
    Array[File] imageFiles = glob("${outdir}/Figures/*.svg")
    File report = "${outdir}/Report/debarcer_report.html"
  }
}


