# Repository for comparative gene co-expression netowrks in sugarcane reutilizing public RNAseq data

Our main goal is to show the potential of reusing RNAseq data in sugarcane genomics.

We have selected a dataset composed of two studys Muñoz-Perez et al. 2025 (https://onlinelibrary.wiley.com/doi/10.1111/ppl.70612) and Ta Quang Kiet et al. 2025 (https://onlinelibrary.wiley.com/doi/10.1111/ppl.70612)
The first one have two genotypes with contrasting NUE (nitrogen responsive and non-responsive), in two contrasting nitrogen avaiability conditions (low-high), samples are taken across leaf segments.
The second one have samples of Saccharum officinarum  and Saccharum robustum under low and high nitrogen (gradient).

Both studies perform RNAseq experiments confronting contransting NUE genotypes with contrasting nitrogen conditions giving some degree of compatibility of both experimental designs.
Both studies claim to have found important players in nitrogen response in their respective experiments. We want to compare the gene co-expression networks resulting from this genotypes and to explore the common and unique processes in nitrogen responose for this genotypes underlying this gene co-expression networks. For this we are going to explore the conservation of the topological features of the networks, the conservation of the neighbors of important players..

(Thinking about making only two networks)

Table 1. Data summary.
| Network | Study                    | Genotype | Trait           | Samples | Soil | Organism                | Notes                           | Reference                   |
|---------|--------------------------|---------|-----------------|---------|------|-------------------------|----------------------------------|-----------------------------| 
| 1       | Muñoz-Perez et al. 2025  | RB937570 | non-responsive   | 24      | sand | Sugarcane               | Include high and low nitrogen conditions | R570              |
| 1       | Muñoz-Perez et al. 2025  | RB975375 | responsive       | 24      | sand | Sugarcane               | Include high and low nitrogen conditions | R570              |
| 2       | Ta Quang Kiet et al. 2025| 51NG3    | responsive       | 9       | sand | Saccharum robustum      | Include high and low nitrogen conditions | LA purple         |
| 2       | Ta Quang Kiet et al. 2025| TAGZ     | non-responsive   | 9       | sand | Saccharum officinarum   | Include high and low nitrogen conditions | LA purple         |


1) Data retrieval and processing

To get Muñoz-Perez et al. 2025 samples we used Sweet Recycler precomputed count matrices (...) 
To get Ta Quang Kiet et al. 2025 we download the data from the The SugarCane multi-Omics Database (https://ngdc.cncb.ac.cn/scod/browse/genome)
we used the same methodology to get the count matrices from that used Sweet Recycler in Ta Quang Kiet et al. 2025, ie:

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

2) Gene co-expression inference
We get rid of transcripts with a coefficient of variation (cv) less than 15% in raw counts.

## Gene filtering

| Network | before filtering (genes) | after filtering cv > 15% | after filtering pcor(cor > 0.8 and FDR < 0.05) | 
|---------|---------------------------|-------------------------|------------------------------------------------|
| 1       |       190,973             |     170,790             |      102,020                                   |
| 2       |       215,183             |     170,736             |      170,103                                   |

## Correlation calculation

We filter out correlations |p| < 0.8, we calculated p-values and correct them using FDR < 0.05
We set the max_perason to 0.9999 this aims to remove suspicious perfect correlations.

(plots of correlations distributions)

| Network | Edges      | Nodes   | Transitivity | Connected Components | Giant Component Size |
|---------|------------|---------|--------------|----------------------|----------------------|
| 1       |75,380,961  | 102,020 | 0.690058     | 958                  | 99,881               |
| 2       |681,090,855 | 170,103 | 0.73957      | 44                   | 170,012              |





