#!/usr/bin/env nextflow
/*
Nextflow pipeline for outbreak detection of bacterial isolates from Illumina short-read data
Usage: nextflow run nf_cmds.nf --sra_list sra.txt --input-dir ./input --config nf_cmds.config
*/

params.sra_list = './sra.txt'
params.input_dir = './input'
params.mash_db_dir = "${projectDir}/data/mash"
// 16S rRNA identification is fragile for fragmented short-read assemblies
// so it's kept as an optional, informational side-branch; 
// Mash against the whole assembly is the primary genus/species identification method.
params.run_16S_id = false

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
    tuple val(sample_id), path("${sample_id}.fna")

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
    tuple val(sample_id), path("${sample_id}_stats.tsv")

    script:
    """
    echo "Computing fastq stats for ${sample_id}..."
    #-a All Statistics (including N50, Q30)
    #-b grabs basename
    seqkit stats -a -b "${r1}" "${r2}" > "${sample_id}_stats.tsv"
    echo "Computing fastq stats for ${sample_id} complete."
    """
}

process filterContigs {
    debug true
    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.fna")

    script:
    """
    #filter contigs to remove those < 1000 bp and with < 10X coverage
    filter.contigs.py \
    	--infile "${assembly}" \
		--outfile "${sample_id}.fna" \
		--discarded "${sample_id}_discarded_contigs.fna" \
		--cov 10 \
		--len 1000 \
		1> "${sample_id}_contig-filtering.stdout.log" \
		2> "${sample_id}_contig-filtering.stderr.log"
    """

}

process genePrediction {
    debug true
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    path "${sample_id}.gff"

    script:
    """
    #perform ab initio coding sequence prediction using prodigal
    prodigal \
    -i "${assembly}" \
    -c \
    -m \
    -f gff \
    -o "${sample_id}.gff" \
    2>&1 | tee "./${sample_id}.log"
    """
}

process extract16S {
    debug true
    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_16S.gff"), emit: gff
    tuple val(sample_id), path("${sample_id}_16S.fa"), emit: fasta

    script:
    """
    #identify the 16S rRNA gene sequence coordinates (GFF format)
    #short-read assemblies often fragment multi-copy rRNA operons, so a
    #default-threshold barrnap run can legitimately find zero 16S hits;
    #`|| true` stops grep's no-match exit (1) from killing the script here
    barrnap \
    "${assembly}" \
    | grep "Name=16S_rRNA;product=16S ribosomal RNA" \
    > 16S.gff || true

    if [[ ! -s 16S.gff ]]; then
        echo "No 16S rRNA hit at default barrnap thresholds; retrying with relaxed --reject/--evalue" >&2
        barrnap --reject 0.1 --evalue 1e-3 \
        "${assembly}" \
        | grep "Name=16S_rRNA;product=16S ribosomal RNA" \
        > 16S.gff || true
    fi

    if [[ ! -s 16S.gff ]]; then
        echo "ERROR: no 16S rRNA gene detected in ${assembly} (sample ${sample_id}), even with relaxed barrnap thresholds (--reject 0.1 --evalue 1e-3)." >&2
        exit 1
    fi

    #extract the nucleotide sequence (as FastA format)
    bedtools getfasta \
    -fi "${assembly}" \
    -bed 16S.gff -s \
    -fo 16S.fa

    mv 16S.gff "${sample_id}_16S.gff"
    mv 16S.fa "${sample_id}_16S.fa"
    """
}


process genusIdentification {
    debug true
    input:
    tuple val(sample_id), path(sixteenS_fasta)

    output:
    tuple val(sample_id), path("${sample_id}_16S_rRNA_blastn.tsv")

    script:
    """

    # Obtain the 16S rRNA database from BLAST server for use globally
        update_blastdb.pl --decompress 16S_ribosomal_RNA
		rm -f tax*

    # perform blastn search against 16S ribosomal RNA database
    # genus: minimum 95% identity is cutoff for 16S rRNA
	{
		echo -e "Organism Name\tQuery ID\tSubject ID\tPercent Identity\tAlignment Length\tMismatches\tGap Openings\tQuery Start\tQuery End\tSubject Start\tSubject End\tE Value\tBit Score\tQuery Coverage"
		blastn -query "${sixteenS_fasta}" \
			-db 16S_ribosomal_RNA \
			-perc_identity 95 \
			-max_target_seqs 10 \
			-outfmt "6 stitle qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" | \
		sed \$'s/ 16S ribosomal RNA[^\t]*/\t/'
	} > "${sample_id}_16S_rRNA_blastn.tsv"
    """
}

process readGenusIdentification {
    debug true
    input:
    tuple val(sample_id), path(blastn_results)
    
    output:
    tuple val(sample_id), stdout

    script:
    """
    #!/usr/bin/env python3
    ## Read the BLASTn results and identify the genus based on the top hit
    import pandas as pd

    df = pd.read_csv("${blastn_results}", sep='\t')
    top_hit = df.sort_values('Bit Score', ascending=False).iloc[0]
    genus = top_hit['Organism Name'].split()[0]

    print(genus, end='')
    """
}

process downloadMashDB {
    debug true
    storeDir params.mash_db_dir

    output:
    path "refseq.genomes.k21s1000.msh", emit: sketch
    path "assembly_summary.txt", emit: summary

    script:
    """
    echo "Downloading pre-sketched NCBI RefSeq genome database (~750MB)..."
    curl -LO https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh

    echo "Downloading RefSeq bacterial assembly summary (~200MB)..."
    curl -LO https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
    """
}

process mashIdentification {
    debug true
    input:
    tuple val(sample_id), path(assembly), path(mash_sketch), path(assembly_summary)

    output:
    tuple val(sample_id), path("${sample_id}_mash_taxonomy_results.tsv")

    script:
    """
    # Mash compares the whole assembly against a local RefSeq sketch database,
    # so it's robust to genomes where a specific gene (eg. 16S) didn't assemble
    mash sketch -m 2 "${assembly}"
    mash dist "${mash_sketch}" "${assembly}.msh" | sort -gk3 | head -n 10 > "${sample_id}_microbial_distances.tsv"
    mash screen "${mash_sketch}" "${assembly}" | sort -gr > "${sample_id}_microbial_screen.tsv"

    {
        echo -e "Organism Name\tSample ID\tPercent Identity\tShared Hashes\tP value"

        awk -F'\t' 'BEGIN {OFS="\t"}

            # Track which file we are currently reading
            FNR == 1 { file_count++ }

            # Step A: Read the assembly_summary.txt taxonomy file
            file_count == 1 && /^GCF/ {
                full_name = \$8;
                if (\$9 != "") {
                    s_val = \$9;
                    sub(/^(strain|strain[:=])/, "", s_val);
                    sub(/^[= ]+/, "", s_val);
                    full_name = full_name " (strain: " s_val ")";
                }
                tax[\$1] = full_name;
                next;
            }

            # Step B: Process the Jaccard distances TSV file
            file_count == 2 && /GCF/ {
                match(\$1, /GCF_[0-9]+\\.[0-9]+/);
                gcf = substr(\$1, RSTART, RLENGTH);
                split(\$2, p, "/");
                sample_map[gcf] = p[length(p)-1];
                next;
            }

            # Step C: Process the microbial screening TSV file produced by Mash
            file_count == 3 && (\$5 ~ /GCF/) {
                match(\$5, /GCF_[0-9]+\\.[0-9]+/);
                gcf = substr(\$5, RSTART, RLENGTH);
                if (gcf in tax && gcf in sample_map) {
                    print tax[gcf], sample_map[gcf], \$1*100, \$2, \$4
                }
            }' "${assembly_summary}" "${sample_id}_microbial_distances.tsv" "${sample_id}_microbial_screen.tsv"
    } > "${sample_id}_mash_taxonomy_results.tsv"
    """
}

process readMashIdentification {
    debug true
    input:
    tuple val(sample_id), path(mash_results)

    output:
    tuple val(sample_id), stdout

    script:
    """
    #!/usr/bin/env python3
    ## Read the Mash taxonomy results and identify genus/species from the top hit
    ## (mash_taxonomy_results.tsv is already sorted best-hit-first upstream)
    import pandas as pd

    df = pd.read_csv("${mash_results}", sep='\t')
    top_hit = df.iloc[0]
    words = top_hit['Organism Name'].split()
    genus = words[0]
    species = f"{words[0]} {words[1]}" if len(words) > 1 else genus

    print(f"{genus}\t{species}", end='')
    """
}

process genotyping {
    debug true
    input:
    tuple val(sample_id), path(assembly), val(genus), val(species)

    output:
    tuple val(sample_id), path("${sample_id}-genotyping.tsv")

    script:
    """
    # mlst's own BLAST-based scheme detection is reliable on its own; genus/species
    # (from Mash) are recorded here for traceability but not forced as --scheme
    echo "Mash-identified organism for ${sample_id}: ${species} (genus: ${genus})"
    mlst \
        "${assembly}" \
        > "${sample_id}-genotyping.tsv"
    """
}

process genomeQA{
        debug true
    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}-busco_summary.txt")

    script:
    """
    ## download BUSCO database if not already present in this work dir
    if [[ ! -f "data/busco/lineages/bacteria_odb10/dataset.cfg" ]]; then
	busco --download_path data/busco --download bacteria_odb10 --quiet
    fi

    ## run BUSCO on the assembly
    busco \
        -i "${assembly}" \
        -l bacteria_odb10 \
        -o "${sample_id}-busco" \
        --out_path . \
        --mode genome

    cp "${sample_id}-busco/short_summary.specific.bacteria_odb10.${sample_id}-busco.txt" "${sample_id}-busco_summary.txt"
    """
}

process contigQA{
        debug true
    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.tsv")

    script:
    """
    # DEBUG RUN: using GUNC's small 'test_data' database instead of the real
    # ~7GB progenomes_2.1 database, to validate this step without a huge
    # download. For production use, switch back to the full database:
    # `gunc download_db data/gunc` (no -db flag) and --db
    # data/gunc/gunc_db_progenomes2.1.dmnd below.
    if [[ ! -f "data/gunc/ci_test.dmnd" ]]; then
	mkdir -p data/gunc
	gunc download_db data/gunc -db test_data
    fi

    # run GUNC on the assembly; output is a directory with a tsv per sample
    mkdir -p "${sample_id}_gunc_out"
    gunc run \
		--input_fasta "${assembly}" \
		--db "data/gunc/ci_test.dmnd" \
		--out_dir "${sample_id}_gunc_out" \
		--detailed_output \
		--threads 8 \
	> "${sample_id}"_gunc.log 2>&1

    find "${sample_id}_gunc_out" -name "*maxCSS_level.tsv" -exec cp {} "${sample_id}.tsv" \\;
"""
}

process phylogeneticAnalysis {
    debug true
    // parsnp has no native osx-arm64 conda build; runs via Docker instead
    // (amd64 image, transparently emulated by Docker Desktop on Apple Silicon)
    container 'staphb/parsnp:1.5.6'
    containerOptions '--platform linux/amd64'

    input:
    // grouped by genus (see workflow) so a tree only ever compares taxonomically
    // related samples - mixing distant genera makes parsnp silently drop the
    // poorly-aligning genome from the tree instead of erroring
    tuple val(genus), path(assemblies)

    output:
    tuple val(genus), path("${genus}-parsnp.tree")

    script:
    """
    mkdir -p parsnp_input_assemblies
    cp ${assemblies} parsnp_input_assemblies/

    # RAxML (parsnp's ML tree builder, used whenever there are enough real SNPs)
    # hard-fails below 3 taxa ("TOO FEW SPECIES") - a 2-leaf tree has no real
    # topology to infer anyway, so treat < 3 genomes like the singleton case
    n_genomes=\$(ls parsnp_input_assemblies | wc -l)
    if [ "\$n_genomes" -lt 3 ]; then
        echo "Only \$n_genomes genome(s) identified as ${genus}; Parsnp/RAxML need at least 3 to build a comparative tree." >&2
        echo "(no tree: only \$n_genomes sample(s) identified as ${genus})" > "${genus}-parsnp.tree"
    else
        # run Parsnp to generate a core SNP phylogenetic tree from all samples of this genus
        parsnp \
        -d parsnp_input_assemblies \
        -r ! \
        -o parsnp_outdir \
        -p 4

        cp parsnp_outdir/parsnp.tree "${genus}-parsnp.tree"
    fi
    """
}

process renderPhylogeneticTree {
    debug true
    input:
    tuple val(genus), path(tree)

    output:
    tuple val(genus), path("${genus}-parsnp.svg")

    script:
    """
    if grep -q "^(no tree" "${tree}"; then
        echo '<svg xmlns="http://www.w3.org/2000/svg" width="420" height="40">' > "${genus}-parsnp.svg"
        echo '<text x="10" y="20" font-family="sans-serif" font-size="12">'"\$(cat "${tree}")"'</text></svg>' >> "${genus}-parsnp.svg"
    else
        figtree -graphic SVG "${tree}" "${genus}-parsnp.svg"
    fi
    """
}

process summarizeResults {
    debug true
    input:
    tuple val(sample_id), path(genotyping_results), path(fastq_metrics), path(genomeQA_results), path(contigQA_results), val(genus), val(species), path(phylogeneticAnalysis_results)

    output:
    path("${sample_id}-summary_report.tsv")

    script:
    """
    # summarize results into a single report
    summarize_results.py \
        --sample_id "${sample_id}" \
        --genotyping "${genotyping_results}" \
        --fastq_metrics "${fastq_metrics}" \
        --genomeQA "${genomeQA_results}" \
        --contigQA "${contigQA_results}" \
        --phylogenetic "${phylogeneticAnalysis_results}" \
        --species "${species}" \
        --genus "${genus}" \
        --output "${sample_id}-summary_report.tsv"
    """
}

process combineSummaryReports {
    debug true
    input:
    path(reports)

    output:
    path("summary_report.tsv")

    script:
    """
    # merge every per-sample report into one file, keeping a single header
    first=\$(ls *-summary_report.tsv | sort | head -n1)
    head -n1 "\$first" > summary_report.tsv
    for f in \$(ls *-summary_report.tsv | sort); do
        tail -n +2 "\$f" >> summary_report.tsv
    done
    """
}


workflow {

    main:
    // Reads can come from an sra.txt of accessions to download, a local
    // input/ directory of fastq.gz read pairs, or both at once - whichever
    // of the two sources is actually present gets used.
    def sra_list_file = file(params.sra_list)
    def input_dir = file(params.input_dir)

    if (!sra_list_file.exists() && !input_dir.isDirectory()) {
        error "No input found: expected an SRA accession list at '${params.sra_list}' and/or a directory of fastq.gz read pairs at '${params.input_dir}'."
    }

    if (sra_list_file.exists()) {
        // parse SRA input list
        // .splitText() means each line as a separate item
        // .trim() removes newlines
        sra_ch = channel.fromPath(params.sra_list)
            .splitText() { it.trim() }
            .filter { it }

        // download SRA files
        downloadSRA(sra_ch)

        // convert sra files to fastq.gz
        convertSRA(downloadSRA.out)
        sra_read_pairs_ch = convertSRA.out
        downloadSRA_out = downloadSRA.out
        convertSRA_out = convertSRA.out
    } else {
        sra_read_pairs_ch = channel.empty()
        downloadSRA_out = channel.empty()
        convertSRA_out = channel.empty()
    }

    if (input_dir.isDirectory()) {
        // pick up local Illumina-style read pairs (eg. Sample_S01_L001_R1_001.fastq.gz);
        // fromFilePairs keys on everything before the R1/R2 marker, so strip the
        // trailing _S##_L### to recover a clean sample_id
        local_read_pairs_ch = channel.fromFilePairs("${params.input_dir}/*_R{1,2}_*.fastq.gz", checkIfExists: false)
            .map { key, reads -> [key.replaceAll(/_S\d+_L\d+$/, ''), reads[0], reads[1]] }
    } else {
        local_read_pairs_ch = channel.empty()
    }

    // clean fastq files (both sources go through the same cleaning step)
    cleanFastq(sra_read_pairs_ch.mix(local_read_pairs_ch))
    // [0] selects only the tuple (r1, r2), not the json/html outputs
    read_pairs_ch = cleanFastq.out[0]

    // assembly fastq files into genome
    assembleGenome(read_pairs_ch)

    // get fastq metrics
    fastqMetrics(read_pairs_ch)

    // filter contigs
    filterContigs(assembleGenome.out)

    // gene annotation
    genePrediction(filterContigs.out)

    // optional 16S rRNA identification (informational only, off by default -
    // see params.run_16S_id). Doesn't feed genotyping or the final summary;
    // Mash below is the primary genus/species identification method.
    if (params.run_16S_id) {
        extract16S(filterContigs.out)
        genusIdentification(extract16S.out.fasta)
        readGenusIdentification(genusIdentification.out)
        sixteenS_gff_out = extract16S.out.gff
        sixteenS_blastn_out = genusIdentification.out
    } else {
        sixteenS_gff_out = channel.empty()
        sixteenS_blastn_out = channel.empty()
    }

    // download/reuse the local Mash RefSeq sketch database (once for the whole run)
    downloadMashDB()

    // identify genus+species per sample via Mash against the whole assembly -
    // far more robust to fragmented assemblies than a single-gene BLAST, and
    // avoids querying NCBI live for every sample
    mash_input_ch = filterContigs.out
        .combine(downloadMashDB.out.sketch)
        .combine(downloadMashDB.out.summary)
    mashIdentification(mash_input_ch)
    readMashIdentification(mashIdentification.out)

    // split the combined "genus\tspecies" stdout into two per-sample value channels
    genus_ch = readMashIdentification.out.map { id, gs -> [id, gs.tokenize('\t')[0]] }
    species_ch = readMashIdentification.out.map { id, gs -> [id, gs.tokenize('\t')[1]] }

    // join the filtered assembly with its Mash-derived genus and species by sample_id
    genotyping_ch = filterContigs.out
        .join(genus_ch)
        .join(species_ch)

    // perform genotyping
    genotyping(genotyping_ch)

    // perform genome quality assessment with BUSCO
    genomeQA(filterContigs.out)

    // perform contig quality assessment with GUNC
    contigQA(filterContigs.out)

    // group filtered assemblies by their Mash-identified genus, so each
    // phylogenetic tree only ever compares taxonomically related samples
    genus_assembly_ch = filterContigs.out
        .join(genus_ch)
        .map { id, assembly, genus -> [genus, assembly] }
        .groupTuple()

    phylogeneticAnalysis(genus_assembly_ch)

    // render each genus's Newick tree to SVG (native arm64 conda; parsnp itself is Docker-only)
    renderPhylogeneticTree(phylogeneticAnalysis.out)

    // map each sample back to its own genus's tree, keyed by sample_id again
    // .join() only reliably pairs unique keys; genus isn't unique when multiple
    // samples share it, so cross every sample against the (small) set of trees
    // and filter to its own genus instead of joining on genus directly
    sample_tree_ch = genus_ch
        .combine(renderPhylogeneticTree.out)
        .filter { id, sample_genus, tree_genus, svg -> sample_genus == tree_genus }
        .map { id, sample_genus, tree_genus, svg -> [id, svg] }

    // join all per-sample results by sample_id, including each sample's own tree
    summary_ch = genotyping.out
        .join(fastqMetrics.out)
        .join(genomeQA.out)
        .join(contigQA.out)
        .join(genus_ch)
        .join(species_ch)
        .join(sample_tree_ch)

    // summarize all results, then merge every sample's row into one report
    summarizeResults(summary_ch)
    combineSummaryReports(summarizeResults.out.collect())

    publish:
    first_output = downloadSRA_out
    second_output = convertSRA_out
    third_output = cleanFastq.out[0]
    fourth_output = assembleGenome.out
    fifth_output = fastqMetrics.out
    sixth_output = filterContigs.out
    seventh_output = genePrediction.out
    eighth_output = sixteenS_gff_out.mix(sixteenS_blastn_out)
    ninth_output = mashIdentification.out
    tenth_output = genotyping.out
    eleventh_output = genomeQA.out
    twelfth_output = contigQA.out
    thirteenth_output = phylogeneticAnalysis.out.mix(renderPhylogeneticTree.out)
    fourteenth_output = combineSummaryReports.out
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

    sixth_output {
        path 'filtered_assemblies'
        mode 'copy'
    }

    seventh_output {
        path 'gene_prediction'
        mode 'copy'
    }

    eighth_output {
        path '16S_rRNA'
        mode 'copy'
    }

    ninth_output {
        path 'mash_identification'
        mode 'copy'
    }

    tenth_output {
        path 'genotyping_results'
        mode 'copy'
    }

    eleventh_output {
        path 'genomeQA'
        mode 'copy'
    }

    twelfth_output {
        path 'contigQA'
        mode 'copy'
    }

    thirteenth_output {
        path 'phylogenetic_analysis'
        mode 'copy'
    }

    fourteenth_output {
        path 'summary_reports'
        mode 'copy'
    }

}