
SAMPLESHEET="samplesheet.csv"
RESULTS_DIR="run1"
GENOME_FA="/dados01/jorge/rnaseq_diatraea/reference_genomes/sugarcane/assembly/SofficinarumxspontaneumR570_771_v2.0.fa.gz"
GENOME_GTF="/dados01/jorge/rnaseq_diatraea/reference_genomes/sugarcane/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf"
MAX_MEMORY="500GB"
MAX_CPUS="300"
WORK_DIR="/dados02/jorge/comparative_saccharum/"


nextflow run nf-core/rnaseq \
        --input "$SAMPLESHEET" \
        --outdir "$RESULTS_DIR" \
        --fasta "$GENOME_FA" \
        --gtf "$GENOME_GTF" \
        --skip_alignment \
        --pseudo_aligner salmon \
        --extra_salmon_quant_args '--seqBias --gcBias --numGibbsSamples 30 --validateMappings' \
        --skip_qc false \
        --skip_fastqc false \
        --skip_rseqc false \
        --skip_multiqc false \
        --max_memory "$MAX_MEMORY" \
        --max_cpus "$MAX_CPUS" \
        $INDEX_FLAGS \
        -w "$WORK_DIR/rnaseq_work" \
        -profile docker \
        -resume
