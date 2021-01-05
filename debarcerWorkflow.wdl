version 1.0

workflow debarcerWorkflow {
  input {
    String outdir = "./"    
    File bamFile
    File regionFile
    Int distance = 1
    Int position = 10
    String separator
    Int readCount = 0
    Boolean truncate = false
    Boolean ignoreOrphans = false
    Boolean ignore = false
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
    outdir: "Output directory with subdirectory structure"    
    bamFile: "Path to alignment file" 
    regionFile: "Path to file with genomic regions"
    distance: "Hamming distance threshold for connecting parent-children umis"
    position: "Umi position threshold for grouping umis together"
    separator: "String separating the UMI from the remaining of the read name"
    readCount: "Minimum number of reads in region required for grouping"
    truncate: "Only pileup columns in the exact region specificied are returned if True. (ie, discard reads overlapping the region)"
    ignoreOrphans: "Ignore orphans (paired reads that are not in a proper pair) if True"
    ignore: "Keep the most abundant family and ignore families at other positions within each group if True"
    maxDepth: "Maximum read depth"
    stepper: "Filter or include reads in the pileup. all: skip reads with BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP flags, nofilter: uses every single read turning off any filtering"
    baseQuality: "Base quality score threshold. Bases with quality scores below the threshold are not used in the consensus"
    familySize: " Comma-separated list of minimum umi family size to collapase on"
    percentThreshold: "Majority rule consensus threshold in pileup column"
    countThreshold: "Base count threshold in pileup column"
    referenceThreshold: "Positions with frequency below reference threshold are considered variable"
    alternativeThreshold: "Variants with frequency above alternative are considered alternative alleles"
    filterThreshold: "Minimum number of reads to pass alternative variants"
    extension: "Figure format"
    report: "Generate a report if true"
    minCov: "Minimum coverage value. Values below are plotted in red"
    minRatio: "Minimum children to parent umi ratio. Values below are plotted in red"
    minUmis: "Minimum umi count. Values below are plotted in red"
    minChildren: "Minimum children umi count. Values below are plotted in red"
  }

  meta {
    author: "Richard Jovelin"
    email: "richard.jovelin@oicr.on.ca"
    description: "A package for de-barcoding and error correction of sequencing data containing molecular barcodes"
    dependencies: [
      {
        name: "python/3.6",
        url: "https://www.python.org/downloads/"
      }
    ]
  }
  
  call regionFileIntoArray {
    input:
      regionFile = regionFile
  } 

   
  scatter(region in regionFileIntoArray.out) {
    call groupUmis {
      input:
        outdir = outdir,
        bamFile = bamFile,
        distance = distance,
        position = position,
        separator = seprator,   
        readCount = readCount,
        truncate = truncate,
        ignoreOrphans = ignoreOrphans,  
        region = region
    }
  }    
  
  #Array[File] dataFilesGroup = groupUmis.dataFile
  Array[File] umiFamiliesGroup = groupUmis.umiFamilies
  #Array[File] umiReleationshipsGroup = groupUmis.umiRelationships
  #Array[File] umiStats = groupUmis.umis
  #Array[File] readCountStats = groupUmis.mappedReadCounts
  
  Array[Pair[String, File]] pairedRegionsUmiFiles = zip(regionFileIntoArray.out, umiFamiliesGroup)
  
  scatter(i in pairedRegionsUmiFiles) {
    call collapseUmis {
      input:
        bamFile = bamFile,
        region = i.left,
        umiFile= i.right,
        maxDepth = maxDepth,
        truncate = truncate,
        ignoreOrphans = ignoreOrphans,
        stepper = stepper,
        separator = separator,
        baseQuality = baseQuality,
        familySize = familySize,
        percentThreshold = percentThreshold,
        countThreshold = countThreshold,
        positionThreshold = positionThreshold
    }
  }

  #File coverage = collapseUmis.coverage
  #Array[File] consensusFiles = collapseUmis.consensus


  call mergeConsensusFiles {
  input:
    outdir = outdir
  }
  
  #File mergedConsensus = mergeConsensusFiles.mergedConsensus

  call callVariants {
    input:
      outdir = outdir,
      referenceThreshold = referenceThreshold,
      alternativeThreshold = alternativeThreshold,
      filterThreshold = filterThreshold,
      familySize = familySize
  }


  #Array[File] vcfFiles = callVariants.vcfFiles

  
  call graph {
    input:
      outdir = outdir,
      extension = extension,
      report = report,
      minCov = minCov,
      minRatio = minRatio,
      minUmis = minUmis,
      minChildren = minChildren,
      refThreshold = refThreshold
  }    

  #File summaryReport = graph.summaryReport
  

  output {
    Array[File] consensusFiles = collapseUmis.consensus
    Array[File] mergedVcfs = "~{outdir}/VCFfiles/Merged_ConsensusFile_famsize*.vcf"
    File summaryReport = graph.summaryReport
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

  parameter_meta {
    modules: "Names and versions of modules to load"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    outdir: "Output directory with subdirectory structure"    
    region: "Genomic region in format chrA:xxx-xxx"
    bamFile: "Path to alignment file" 
    distance: "Hamming distance threshold for connecting parent-children umis"
    position: "Umi position threshold for grouping umis together"
    separator: "String separating the UMI from the remaining of the read name"
    readCount: "Minimum number of reads in region required for grouping"
    truncate: "Only pileup columns in the exact region specificied are returned if True. (ie, discard reads overlapping the region)"
    ignoreOrphans: "Ignore orphans (paired reads that are not in a proper pair) if True"
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

  parameter_meta {
    modules: "Names and versions of modules to load"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    outdir: "Output directory with subdirectory structure"    
    bamFile: "Path to alignment file" 
    region: "Genomic region in format chrA:xxx-xxx"
    maxDepth: "Maximum read depth"
    truncate: "Only pileup columns in the exact region specificied are returned if True. (ie, discard reads overlapping the region)"
    ignoreOrphans: "Ignore orphans (paired reads that are not in a proper pair) if True"
    stepper: "Filter or include reads in the pileup. all: skip reads with BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP flags, nofilter: uses every single read turning off any filtering"
    separator: "String separating the UMI from the remaining of the read name"
    baseQuality: "Base quality score threshold. Bases with quality scores below the threshold are not used in the consensus"
    familySize: " Comma-separated list of minimum umi family size to collapase on"
    umiFile: "File with UMI parent-children relationships"
    percentThreshold: "Majority rule consensus threshold in pileup column"
    countThreshold: "Base count threshold in pileup column"
    positionThreshold: "Umi position threshold for grouping umis together"
    refFasta: "Path to to the reference genome"
    refDict: "Path to the reference dictionary"
    refIndex: "Path to the reference index"
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


  parameter_meta {
    modules: "Names and versions of modules to load"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    outdir: "Output directory with subdirectory structure"    
    refFasta: "Path to to the reference genome"
    refDict: "Path to the reference dictionary"
    refIndex: "Path to the reference index"
    referenceThreshold: "Positions with frequency below reference threshold are considered variable"
    alternativeThreshold: "Variants with frequency above alternative are considered alternative alleles"
    filterThreshold: "Minimum number of reads to pass alternative variants"
    familySize: " Comma-separated list of minimum umi family size to collapase on"
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
    Boolean report = true
    Int minCov = 1000
    Float minRatio = 0.1
    Int minUmis = 1000
    Int minChildren = 500
    Float refThreshold = 95
  }    

 
  parameter_meta {
    modules: "Names and versions of modules to load"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    outdir: "Output directory with subdirectory structure"    
    extension: "Figure format"
    report: "Generate a report if true"
    minCov: "Minimum coverage value. Values below are plotted in red"
    minRatio: "Minimum children to parent umi ratio. Values below are plotted in red"
    minUmis: "Minimum umi count. Values below are plotted in red"
    minChildren: "Minimum children umi count. Values below are plotted in red"
    refThreshold: "Positions with frequency below reference threshold are considered variable"
  }

  command <<<
    set -euo pipefail
    debarcer plot -d ~{outdir} -e ~{extension} -r ~{report} -mv ~{minCov} -mr ~{minRatio} -mu ~{minUmis} -mc ~{minChildren} -rt ~{refThreshold}
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    Array[File] figureFiles = glob("${outdir}/Figures/*.png")
    Array[File] imageFiles = glob("${outdir}/Figures/*.svg")
    File summaryReport = "${outdir}/Report/debarcer_report.html"
  }
}


task mergeConsensusFiles {
  input {
    String modules = "debarcer/2.1.3"
    Int memory = 32
    Int timeout = 10
    String outdir = "./"
  }
  
  command <<<
    set -euo pipefail
    debarcer merge -d ~{outdir} -dt consensusfiles
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedConsensus = "${outdir}/Consfiles/Merged_ConsensusFile.cons"
  }
}



 
task regionFileIntoArray {
  input {
    File regionFile
    Int memory = 1
    Int timeout = 1
  }    

  command <<<
    set -euo pipefail
    arr=()
    for line in `cat "~{regionFile}" | while read line ; do echo $line | sed "s/ /:/" | sed "s/ /-/"; done;`; do arr+=("$line"); done;
  >>>

  runtime {
    memory:  "~{memory} GB"
    timeout: "~{timeout}"
  }

  output {
    Array[String] out = "${arr}"
  }
}




