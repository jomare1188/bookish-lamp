samples="samplesheet.csv"
outdir="run1"
annotation="reference/Mper_Falcon.PsGc.gff"
genome="reference/Mper_Falcon.PsGc.fasta"

#nextflow run nf-core/rnaseq --input $samples --outdir $outdir --gff $annotation --fasta $genome --aligner star_salmon --skip_qc false -resume -profile conda --min_mapped_reads 0 --seq_platform ILLUMINA


nextflow run nf-core/rnaseq --input $samples --outdir $outdir --gff $annotation --fasta $genome --aligner star_salmon --skip_qc false -resume -profile docker --min_mapped_reads 0 --seq_platform ILLUMINA -c custom.config --skip_dupradar




