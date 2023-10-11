## Genome installation

signeR accepts genome from any species installed with [BSgenome](https://bioconductor.org/packages/release/BiocViews.html#___BSgenome).

For exemple, to install the hg19 human genome, you can use the command:
```R
BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
```
