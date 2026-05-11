#!/usr/bin/env nextflow

// Developer notes
//
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/ingress'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process alignReads {
    label "wfflu"
    cpus 2
    input:
        tuple val(meta), path("reads.fastq.gz")
        path "reference.fasta"
    output:
        tuple val(meta), path("align.bam"), path("align.bam.bai"), emit: alignments
        tuple val(meta), path("align.bamstats"), emit: bamstats

    shell:
    """
    # Check for input files
    if [[ ! -f "reads.fastq.gz" ]]; then
        echo "Error: reads.fastq.gz not found!" >&2
        exit 1
    fi

    if [[ ! -f "reference.fasta" ]]; then
        echo "Error: reference.fasta not found!" >&2
        exit 1
    fi

    mini_align -i reads.fastq.gz -r reference.fasta -p align_tmp -t 2 -m

    # keep only mapped reads
    samtools view --write-index -F 4 align_tmp.bam -o align.bam##idx##align.bam.bai

    # get stats from bam
    stats_from_bam -o align.bamstats -s align.bam.summary -t 2 align.bam
    """
}

process coverageCalc {
      label "wfflu"
      cpus 1
      input:
          tuple val(meta), path("align.bam"), path("align.bam.bai")
      output:
          tuple val(meta), path("depth.txt")
      """
      samtools depth -aa align.bam -Q 20 -q 1 > depth.txt
      """

}

process downSample {
    label 'wfflu'
    cpus 1
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai")
        path "reference.fasta"
    output:
        tuple val(meta), path("all_merged.sorted.bam"), path("all_merged.sorted.bam.bai"), emit: alignments
    """
    # get region info from fasta
    samtools faidx reference.fasta
    cut -f1-2 reference.fasta.fai > regions.txt
    
    # for every region we're going to downsample separatley
    while read -r region length;
    do

      # get upper and lower bounds of reference span

      upper=`echo \$((\${length}+(\${length}*10/100)))`
      lower=`echo \$((\${length}-(\${length}*10/100)))`

      # filter reads in region and covering region

      samtools view -bh align.bam -e "length(seq)>\${lower} && length(seq)<\${upper}" \${region} > \${region}.bam;

      # ignore regions with no reads
      count=`samtools view -c \${region}.bam`

      if [ "\${count}" -eq "0" ];
      then
        echo "no reads in \${region} so continuing"
        cp \${region}.bam \${region}_all.bam
        continue;
      fi

      lines=( ${params.downsample} / 2 )

      # get header 
      samtools view -H \${region}.bam > \${region}_all.sam

      samtools view -F16 \${region}.bam | shuf -n \${lines} >> \${region}_all.sam
      samtools view -f16 \${region}.bam | shuf -n \${lines} >> \${region}_all.sam
      samtools view -bh \${region}_all.sam > \${region}_all.bam

    done < regions.txt

    samtools merge all_merged.bam *_all.bam
    samtools sort all_merged.bam > all_merged.sorted.bam
    samtools index all_merged.sorted.bam
    echo "done"

    """
}

process medakaVariants {
    label "medaka"
    cpus 1
    input:
        tuple val(meta), path("downsample.bam"), path("downsample.bam.bai")
        path("reference.fasta")
    output:
        tuple val(meta), path("variants.annotated.filtered.vcf")
    script:
    // we use `params.override_basecaller_cfg` if present; otherwise use
    // `meta.basecall_models[0]` (there should only be one value in the list because
    // we're running ingress with `allow_multiple_basecall_models: false`; note that
    // `[0]` on an empty list returns `null`)
    String basecall_model = params.override_basecaller_cfg ?: meta.basecall_models[0]
    if (!basecall_model) {
        error "Found no basecall model information in the input data for " + \
            "sample '$meta.alias'. Please provide it with the " + \
            "`--override_basecaller_cfg` parameter."
    }
    """
    medaka consensus downsample.bam consensus.hdf --model "${basecall_model}:consensus"
    medaka variant --gvcf reference.fasta consensus.hdf variants.vcf --verbose
    medaka tools annotate --debug --pad 25 variants.vcf reference.fasta downsample.bam variants.annotated.vcf

    bcftools filter -e "ALT='.'" variants.annotated.vcf | bcftools filter -o variants.annotated.filtered.vcf -O v -e "INFO/DP<${params.min_coverage}" -
    """
}

process makeConsensus {
    label "wfflu"
    cpus 1
    input:
        tuple val(meta), path("variants.annotated.filtered.vcf"), path("depth.txt")
        path "reference.fasta"
    output:
        tuple val(meta), path("draft.consensus.fasta")
    """
    awk '{if (\$3<${params.min_coverage}) print \$1"\t"\$2+1}' depth.txt > mask.regions
    bgzip variants.annotated.filtered.vcf
    tabix variants.annotated.filtered.vcf.gz

    bcftools consensus --mask mask.regions  --mark-del '-' --mark-ins lc --fasta-ref reference.fasta -o draft.consensus.fasta variants.annotated.filtered.vcf.gz
    """
}

process typeFlu { 
    label "wfflutyping"
    cpus 1
    input:
        tuple val(meta), path("consensus.fasta")
        path("blast_db")
    output:
        tuple val(meta), path("rsv.typing.tsv"), emit: typing
        path "abricate.version", emit: version
    """
    abricate --version | sed 's/ /,/' > abricate.version
    abricate --datadir blast_db --db rsv --minid 80 --mincov 1 --debug consensus.fasta > rsv.typing.tsv
    """
}
// --mincov 1: this is set so that no hits are filtered out. The coverage here is dependent on gaps and length of the hit - not accurate. Will use GATK % coverage at 30X instead s

// original: abricate --datadir blast_db --db insaflu -minid 70 -mincov 60 --quiet consensus.fasta 1> insaflu.typing.txt


process processType {
    label "wfflu"
    cpus 1
    input:
        tuple val(meta), path("rsv.typing.tsv")
    output:
        tuple val(meta), path("processed_type.json"), emit: processed_type_json
        tuple val(meta), path("processed_type.tsv"), emit: processed_type_tsv
    script:
    """
    workflow-glue process_abricate --typing rsv.typing.tsv --output_json processed_type.json --output_tsv processed_type.tsv
    """
}


process prepNextclade {
    label "wfflu"
    cpus 1
    input:
        tuple val(meta), path("typing.json"), path("consensus.fasta")
    output:
        tuple val(meta), path("typed.consensus.fasta"), emit: consensus
    script:
    String alias = meta.alias
    """
    workflow-glue nextclade_helper \
        --typing typing.json \
        --consensus "consensus.fasta" \
        --output "typed.consensus.fasta"
    """
}



process nextclade {
    label "nextclade"
    cpus 1
    input:
        tuple val(meta), path("typed.consensus.fasta")
        path("nextclade_data")
    output:
        tuple val(meta), path("nextclade.tsv"), emit: nextclade_tsv
        tuple val(meta), path("nextclade.csv"), emit: nextclade_csv
        tuple val(meta), path("nextclade.json"), emit: nextclade_json
        tuple val(meta), path("nextclade.nwk"), emit: nextclade_nwk
        tuple val(meta), path("nextclade.aligned.fasta"), emit: nextclade_fa
    script:
    """
    nextclade --version

    HEADER=\$(head -n 1 typed.consensus.fasta)

    if [[ \$HEADER == *"typeA"* ]]; then
        echo "file is type A"
        SUBTYPE_DATA_DIR="${params.nextclade_data}/rsv_a"
    elif [[ \$HEADER == *"typeB"* ]]; then
        echo "file is type B"
        SUBTYPE_DATA_DIR="${params.nextclade_data}/rsv_b"
    elif [[ \$HEADER == "NA" ]]; then
        echo "No typing infromation"
        printf "NA\n" >> nextclade.tsv && 
        printf "NA\t%.0s" {1..69} >> nextclade.tsv && 
        printf "NA\n" >> nextclade.tsv
        printf "NA\n" >> nextclade.csv && 
        printf "NA; %.0s" {1..69} >> nextclade.csv && 
        printf "NA\n" >> nextclade.csv
        echo "NA" > nextclade.json
        echo "NA" > nextclade.nwk
        echo "NA" > nextclade.aligned.fasta
        exit 0
    elif [[ \$HEADER == "mixedAB" ]]; then
        echo "No typing infromation"
        printf "mixedAB\nNA\tmixedAB\t" > nextclade.tsv && 
        printf "NA\t%.0s" {1..67} >> nextclade.tsv && 
        printf "NA\n" >> nextclade.tsv
        printf "mixedAB\nNA; mixedAB; " > nextclade.csv && 
        printf "NA; %.0s" {1..69} >> nextclade.csv && 
        printf "NA\n" >> nextclade.csv
        echo "NA" > nextclade.json
        echo "NA" > nextclade.nwk
        echo "NA" > nextclade.aligned.fasta
        exit 0
    fi

    echo \$SUBTYPE_DATA_DIR

    nextclade run \
    --input-dataset \$SUBTYPE_DATA_DIR \
    --output-all=./ \
    typed.consensus.fasta 
    """
}


process getVersions {
    label "wfflu"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    bcftools --version | head -1 | sed 's/ /,/' >> versions.txt
    samtools --version | grep samtools | head -1 | sed 's/ /,/' >> versions.txt
    minimap2 --version | head -1 | sed 's/^/minimap2,/' >> versions.txt
    """
}


process getParams {
    label "wfflu"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process collectFilesInDir {
    label "wfflu"
    cpus 1
    input: tuple val(meta), path("staging_dir/*"), val(dirname)
    output: tuple val(meta), path(dirname)
    script:
    """
    mv staging_dir $dirname
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
    label "wfflu"
    cpus 1
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    echo $fname
    echo $dirname
    """
}



// workflow module
workflow pipeline {
    take:
        samples
        reference
        blastdb
        nextclade_data
    main:
        samples.multiMap{ meta, path, stats ->
            meta: meta
            stats: stats
        }.set { for_report }
    
        stats = for_report.stats.collect()

        alignment = alignReads(samples.map{ meta, reads, stats -> [ meta,reads ] }, reference)
        coverage = coverageCalc(alignment.alignments)

        // do crude downsampling
        if (params.rbk){
            println("RBK data - NOT Downsampling!!!")
            downsample = alignment
        } else if (params.downsample != null){
            println("Downsampling!!!")
            downsample = downSample(alignment.alignments, reference)
        } else {
            println("NOT Downsampling!!!")
            downsample = alignment
        }

        variants = medakaVariants(downsample.alignments, reference)

        for_draft = variants.join(coverage)

        draft = makeConsensus(for_draft, reference)
        flu_type = typeFlu(draft, blastdb)

        processed_type = processType(flu_type.typing)

        nextclade_prep = prepNextclade(processed_type.processed_type_json.join(draft, remainder: true))

        nextclade_result = nextclade(nextclade_prep, nextclade_data) // swu

        software_versions = getVersions()
        software_versions = software_versions.mix(flu_type.version.first()) | collectFile()
        workflow_params = getParams()

        // get all the per sample results together
        ch_per_sample_results = samples
        | join(coverage)
        | join(flu_type.typing)
        | join(processed_type.processed_type_json)
        | join(processed_type.processed_type_tsv)
        | join(nextclade_result.nextclade_tsv)
        | join(nextclade_result.nextclade_csv)

        // create channel with files to publish; the channel will have the shape `[file,
        // name of sub-dir to be published in]`.

        ch_to_publish = Channel.empty()
        | mix(
            software_versions | map { [it, null] },
            workflow_params | map { [it, null] },
            alignment.alignments
            | map { meta, bam, bai -> [bam, "$meta.alias/alignments"] },
            alignment.alignments
            | map { meta, bam, bai -> [bai, "$meta.alias/alignments"] }, 
            alignment.bamstats
            | map { meta, bamstats -> [bamstats, "$meta.alias/alignments"] },

            variants
            | map { meta, vcf -> [vcf, "$meta.alias/variants"]},
            coverage
            | map { meta, depth -> [depth, "$meta.alias/coverage"]},
            draft
            | map { meta, fa -> [fa, "$meta.alias/consensus"]},
            nextclade_prep //swu
            | map { meta, fasta -> [fasta, "$meta.alias/consensus"]},

            nextclade_result.nextclade_tsv //swu
            | map { meta, tsv -> [tsv, "$meta.alias/nextclade"]},
            nextclade_result.nextclade_csv //swu
            | map { meta, csv -> [csv, "$meta.alias/nextclade"]},
            nextclade_result.nextclade_json //swu
            | map { meta, json -> [json, "$meta.alias/nextclade"]},
            nextclade_result.nextclade_nwk //swu
            | map { meta, nwk -> [nwk, "$meta.alias/nextclade"]},
            nextclade_result.nextclade_fa //swu
            | map { meta, fasta -> [fasta, "$meta.alias/nextclade"]},

            processed_type.processed_type_json
            | map { meta, json -> [json, "$meta.alias/typing"]},
            processed_type.processed_type_tsv
            | map { meta, tsv -> [tsv, "$meta.alias/typing"]}
        )
        ch_to_publish.subscribe { println("Published item: ${it}") } // swu

    emit:
        results = ch_to_publish
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    // warn the user if overriding the basecall models found in the inputs
    if (params.override_basecaller_cfg) {
        log.warn \
            "Overriding basecall model with '${params.override_basecaller_cfg}'."
    }

    samples = fastq_ingress([
        "input": params.fastq,
        "stats": true,
        "sample_sheet": params.sample_sheet,
        "allow_multiple_basecall_models": false,
    ])


    //get reference
    if (params.reference == null){
      params.remove('reference')
      params._reference = projectDir.resolve("/hpf/largeprojects/pray/microbiology_testing/rsv/data/rsv_reference/rsv_ref.fasta").toString() // swu 
    } else {
      params._reference = file(params.reference, type: "file", checkIfExists:true).toString()
      params.remove('reference')
    }

    //get db
    if (params.blastdb == null){
      params.remove('blastdb')
      params._blastdb = projectDir.resolve("/hpf/largeprojects/pray/microbiology_testing/rsv/data/rsv_reference/blastdb").toString()
    } else {
      params._blastdb = file(params.blastdb, type: "directory", checkIfExists:true).toString() 
      params.remove('blastdb')
    }

    //get nextclade_data (swu)
    // params._nextclade_data = file(params.nextclade_data, type: "directory", checkIfExists:true).toString() 
    if (params.nextclade_data == null){
      params.remove('nextclade_data')
      params._nextclade_data = projectDir.resolve("/hpf/largeprojects/pray/microbiology_testing/rsv/data/nextclade").toString()
    } else {
      params._nextclade_data = file(params.nextclade_data, type: "directory", checkIfExists:true).toString() 
      params.remove('nextclade_data')
    }
    
    pipeline(samples, params._reference, params._blastdb, params._nextclade_data)
    pipeline.out.results
    | toList
    | flatMap
    | output
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
