version 1.0

workflow debarcerWorkflow {
  input {
    String outdir = "./"    
    File bamFile
    File bamIndex
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
  }


  parameter_meta {
    outdir: "Output directory with subdirectory structure"    
    bamFile: "Path to alignment file" 
    bamIndex: "Index file of the bam file" 
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
    output_meta: {
    dataFilesGroup: "List of data files with UMI parent-child counts from grouping",
    umiFamiliesGroup: "List of UMI files with UMI counts within families for each region",
    umiReleationshipsGroup: "list of summary files with parent-child umi relationships",
    umiStats: "List of files with UMI counts before grouping",
    readCountStats: "List of files with mapping read information for each region",
    coverage: "File with coverage information on each region",
    consensusFiles: "List of consensus files with nucleotide information after collapse",
    vcfFiles: "List of vcf files with SNP and indel information",
    mergedConsensus: "Merged consensus file across regions",
    mergedUmiFile: "Merged umi file with parent-child UMIs across regions",
    mergedDataFile: "Merged data file with UMI parent and child counts across regions"
    }
  }
  
  call regionFileIntoArray {
    input:
      regionFile = regionFile
  } 

  Array[String] genomic_regions = regionFileIntoArray.out
   
  scatter(interval in genomic_regions) {
    call groupUmis {
      input:
        outdir = outdir,
        bamFile = bamFile,
        bamIndex = bamIndex,
        distance = distance,
        position = position,
        separator = separator,   
        readCount = readCount,
        truncate = truncate,
        ignoreOrphans = ignoreOrphans,  
        region = interval
    }
  
    call collapseUmis {
      input:
        bamFile = bamFile,
        region = interval,
        umiFile= groupUmis.umiFamilies,
        maxDepth = maxDepth,
        truncate = truncate,
        ignoreOrphans = ignoreOrphans,
        stepper = stepper,
        separator = separator,
        baseQuality = baseQuality,
        familySize = familySize,
        percentThreshold = percentThreshold,
        countThreshold = countThreshold,
        position = position
    }
  }

  Array[File] dataFilesGroup = groupUmis.dataFile
  Array[File] umiFamiliesGroup = groupUmis.umiFamilies
  Array[File] umiReleationshipsGroup = groupUmis.umiRelationships
  Array[File] umiStats = groupUmis.umis
  Array[File] readCountStats = groupUmis.mappedReadCounts
  
  File coverage = collapseUmis.coverage[0]
  Array[File] consensusFiles = collapseUmis.consensus

  scatter(consensus in consensusFiles) {
    call callVariants {
      input:
        consensusFile = consensus,
        outdir = outdir,
        referenceThreshold = referenceThreshold,
        alternativeThreshold = alternativeThreshold,
        filterThreshold = filterThreshold,
        familySize = familySize
    }
  }

  Array[Array[File]] vcfFiles = callVariants.vcfFiles


  call mergeConsensusFiles {
    input:
      outdir = outdir,
      consensusFiles = consensusFiles
  }

  File mergedConsensus = mergeConsensusFiles.mergedConsensus
  
  call mergeUmiFiles {
    input: 
      outdir = outdir,
      umiFiles = umiFamiliesGroup
  }

  File mergedUmiFile = mergeUmiFiles.mergedUmiFile

  call mergeDataFiles {
    input:
      outdir = outdir,
      dataFiles = dataFilesGroup
  }

  File mergedDataFile = mergeDataFiles.mergedDataFile

  output {
    Array[File] outputdataFilesGroup =  dataFilesGroup
    Array[File] outputumiFamiliesGroup = umiFamiliesGroup
    Array[File] outputumiReleationshipsGroup = umiReleationshipsGroup
    Array[File] outputumiStats = umiStats
    Array[File] outputreadCountStats = readCountStats
    File outputcoverage = coverage
    Array[File] outputconsensusFiles = consensusFiles
    Array[Array[File]] outputvcfFiles = vcfFiles
    File outputmergedConsensus = mergedConsensus
    File outputmergedUmiFile = mergedUmiFile
    File outputmergedDataFile = mergedDataFile
  }
}


task groupUmis {
  input {
    String modules = "debarcer/2.1.4"
    Int memory = 32
    Int timeout = 36
    String outdir = "./"    
    String region
    File bamFile
    File bamIndex
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
    bamIndex: "Index file of the bam file" 
    distance: "Hamming distance threshold for connecting parent-children umis"
    position: "Umi position threshold for grouping umis together"
    separator: "String separating the UMI from the remaining of the read name"
    readCount: "Minimum number of reads in region required for grouping"
    truncate: "Only pileup columns in the exact region specificied are returned if True. (ie, discard reads overlapping the region)"
    ignoreOrphans: "Ignore orphans (paired reads that are not in a proper pair) if True"
  }


  command <<<
    set -euo pipefail
    echo ~{region}
    debarcer group -o ~{outdir} -b ~{bamFile} -d ~{distance} -p ~{position} -s "~{separator}" -rc ~{readCount} -i ~{ignoreOrphans} -t ~{truncate} -r ~{region}
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
    String modules = "debarcer/2.1.4  hg19/p13"
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
    Int position = 10
    File refFasta = "$HG19_ROOT/hg19_random.fa"
    File refDict = "$HG19_ROOT/hg19_random.dict"
    File refIndex = "$HG19_ROOT/hg19_random.fa.fai"
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
    position: "Umi position threshold for grouping umis together"
    refFasta: "Path to to the reference genome"
    refDict: "Path to the reference dictionary"
    refIndex: "Path to the reference index"
  }

  command <<<
    set -euo pipefail
    cp ~{refDict} .
    cp ~{refIndex} .
    debarcer collapse -o ~{outdir} -b ~{bamFile} -rf ~{refFasta} -r ~{region} -u ~{umiFile} -f ~{familySize} -ct ~{countThreshold} -pt ~{percentThreshold} -p ~{position} -m ~{maxDepth} -t ~{truncate} -i ~{ignoreOrphans} -stp ~{stepper} -s ~{separator} -bq ~{baseQuality}
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
    File consensusFile
    String modules = "debarcer/2.1.4  hg19/p13"
    Int memory = 32
    Int timeout = 36
    String outdir = "./"    
    File refFasta = "$HG19_ROOT/hg19_random.fa"
    File refDict = "$HG19_ROOT/hg19_random.dict"
    File refIndex = "$HG19_ROOT/hg19_random.fa.fai"
    Float referenceThreshold = 95
    Float alternativeThreshold = 2
    Int filterThreshold = 10
    String familySize
  }


  parameter_meta {
    consensusFile: "File with nucleotide counts"
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
    familySize: "Comma-separated list of minimum umi family size to collapase on"
  }

  command <<<
    set -euo pipefail
    cp ~{refDict} .
    cp ~{refIndex} .
    debarcer call -cf ~{consensusFile} -o ~{outdir} -rf ~{refFasta} -rt ~{referenceThreshold} -at ~{alternativeThreshold} -ft ~{filterThreshold} -f ~{familySize}
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


task mergeConsensusFiles {
  input {
    String modules = "debarcer/2.1.4"
    Int memory = 32
    Int timeout = 10
    String outdir = "./"
    Array[File] consensusFiles
  }

  parameter_meta {
    modules: "Names and versions of modules to load"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    outdir: "Output directory with subdirectory structure"
    consensusFiles: "List of consensus files to be merged"    
  }

  command <<<
    set -euo pipefail
    debarcer merge -o ~{outdir} -f ~{sep =" " consensusFiles} -dt consensusfiles
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


task mergeUmiFiles {
  input {
    String modules = "debarcer/2.1.4"
    Int memory = 32
    Int timeout = 10
    String outdir = "./"
    Array[File] umiFiles
  }

  parameter_meta {
    modules: "Names and versions of modules to load"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    outdir: "Output directory with subdirectory structure"
    umiFiles: "List of UMI files to be merged"    
  }

  command <<<
    set -euo pipefail
    debarcer merge -o ~{outdir} -f ~{sep =" " umiFiles} -dt umifiles
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedUmiFile = "${outdir}/Umifiles/Merged_UmiFile.json"
  }
}


task mergeDataFiles {
  input {
    String modules = "debarcer/2.1.4"
    Int memory = 32
    Int timeout = 10
    String outdir = "./"
    Array[File] dataFiles
  }

  parameter_meta {
    modules: "Names and versions of modules to load"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
    outdir: "Output directory with subdirectory structure"
    dataFiles: "List of data files to be merged"    
  }

  command <<<
    set -euo pipefail
    debarcer merge -o ~{outdir} -f ~{sep =" " dataFiles} -dt datafiles
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File mergedDataFile = "${outdir}/Datafiles/Merged_DataFile.csv"
  }
}


task regionFileIntoArray {
  input {
    File regionFile
    Int memory = 1
    Int timeout = 1
  }    


  parameter_meta {
    regionFile: "Path to file with genomic regions"
    memory: "Memory allocated for this job"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail
    cat ~{regionFile} | sed "s/\t/:/" | sed "s/\t/-/"
  >>>

  output {
    Array[String] out = read_lines(stdout())
  }

  runtime {
    memory:  "~{memory} GB"
    timeout: "~{timeout}"
  }
}


