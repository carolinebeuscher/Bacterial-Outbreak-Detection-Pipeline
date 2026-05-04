#!/usr/bin/env nextflow


/*
Nextflow pipeline for cleaning and assembling fastq files 
Usage: nextflow run nf_cmds.nf --sra_list sra.txt -c nf_cmds.config
*/

params.sra_list = './sra.txt'

process downloadSRA {
    debug true
    input:
    val sra_name
    
    output:
    path "${sra_name}/${sra_name}.sra"

    script:
    """
    echo "Fetching accession ${sra_name}..."
    #get accessions from ncbi 
    prefetch "${sra_name}"
    echo "Fetching accession ${sra_name} complete."

    """
}

process convertSRA {
    debug true
    input:
    path input_sra
    
    output:
    tuple val(input_sra.baseName), path("${input_sra.baseName}_1.fastq.gz"), path("${input_sra.baseName}_2.fastq.gz")

    script:
    """
    echo "Converting to ${input_sra} fastq..."
    # convert to all sra files to fastq files
    fasterq-dump \
    "${input_sra}" \
    --outdir . \
    --split-files \
    --skip-technical

    # compress
    pigz -9 *.fastq
    echo "Converting ${input_sra} complete."
    """
}

process cleanFastq {
    debug true
    input:
    tuple val(sample_id), path(r1), path(r2)
    
    output:
    tuple val(sample_id), path("${sample_id}.R1.fq.gz"), path("${sample_id}.R2.fq.gz")
    path "${sample_id}.json"
    path "${sample_id}.html"

    script:
    """
    echo "Cleaning ${sample_id}..."
    # Clean with fastp
    fastp \
    -i "${r1}" \
    -I "${r2}" \
    -o "${sample_id}.R1.fq.gz" \
    -O "${sample_id}.R2.fq.gz" \
    --json "${sample_id}.json" \
    --html "${sample_id}.html"
    echo "Cleaning ${sample_id} complete."
    """
}


process assembleGenome {
    debug true
    input:
    tuple val(sample_id), path(r1), path(r2)
    
    output:
    path "${sample_id}.fna"

    script:
    """
    echo "Assembling ${sample_id}..."
    # Assemble with skesa
    skesa \
    --reads "${r1}","${r2}" \
    --cores 4 \
    --min_contig 1000 \
    --contigs_out "${sample_id}".fna
    echo "Assembling ${sample_id} complete."
    """
}


process fastqMetrics {
    debug true
    input:
    tuple val(sample_id), path(r1), path(r2)
    
    output:
    path "${sample_id}_stats.tsv"

    script:
    """
    echo "Computing fastq stats for ${sample_id}..."
    #-a All Statistics (including N50, Q30)
    #-b grabs basename
    seqkit stats -a -b "${r1}" "${r2}" > "${sample_id}_stats.tsv"
    echo "Computing fastq stats for ${sample_id} complete."
    """
}



workflow {

    main:
    // parse SRA input list
    // .splitText() means each line as a separate item
    // .trim() removes newlines
    sra_ch = channel.fromPath(params.sra_list)
        .splitText() { it.trim() }

    // download SRA files
    downloadSRA(sra_ch)

    // convert sra files to fastq.gz
    convertSRA(downloadSRA.out)

    // clean fastq files
    cleanFastq(convertSRA.out)
    // [0] selects only the tuple (r1, r2), not the json/html outputs
    read_pairs_ch = cleanFastq.out[0]

    // assembly fastq files into genome
    assembleGenome(read_pairs_ch)

    // get fastq metrics
    fastqMetrics(read_pairs_ch)


    publish:
    first_output = downloadSRA.out
    second_output = convertSRA.out
    third_output = cleanFastq.out[0]
    fourth_output = assembleGenome.out
    fifth_output = fastqMetrics.out
}

output {
    first_output {
        path 'raw_fastq'
        mode 'copy'
    }

    second_output {
        path 'raw_fastq'
        mode 'copy'
    }

    third_output {
        path 'clean_fastq'
        mode 'copy'
    }

    fourth_output {
        path 'Assemblies'
        mode 'copy'
    }

    fifth_output {
        path 'clean_fastq'
        mode 'copy'
    }

}