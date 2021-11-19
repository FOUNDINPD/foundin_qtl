# foundin_qtl
FOUNDIN-PD multi-omics QTL analysis

## This project contains the per-omic quantitative trait loci (QTL) analyses
- rnab-qtl: bulk RNAseq eQTL analysis
    - includes interaction eQTL analysis for bulk RNA with [SCADEN](https://github.com/KevinMenden/scaden) estimated cell-type fractions
    - includes quantitative trait score (QTS) analysis based on PD risk
    - includes colocalization with common PD risk loci
- atac-qtl: bulk ATACseq QTL (caQTL) analysis
    - includes quantitative trait score (QTS) analysis based on PD risk
    - includes colocalization with common PD risk loci
- scrn-qtl: single-cell RNAseq eQTL analysis
    - includes quantitative trait score (QTS) analysis based on PD risk
    - includes colocalization with common PD risk loci
- rnab3a-qtl: alternative polyadenylation (APA) quantative trait loci (3'aQTL) analysis from bulk RNAseq
    - includes quantitative trait score (QTS) analysis based on PD risk
    - includes colocalization with common PD risk loci