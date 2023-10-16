## SNV matrix

SNV matrix is a text file with a (tab-delimited) matrix of SNV counts found on analyzed genomes. It must contain one row for each genome sample and 97 columns, the first one with sample IDs and, after that, one column for each mutation type. Mutations should be specified in the column names (headers), by both the base change and the trinucleotide context where it occurs (for example: C>A:ACA). The table below shows an example of the SNV matrix structure.

<table class="table table-striped table-hover" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> </th>
   <th style="text-align:left;"> C>A:ACA </th>
   <th style="text-align:left;"> C>A:ACC </th>
   <th style="text-align:left;"> C>A:ACG </th>
   <th style="text-align:left;"> C>A:ACT </th>
   <th style="text-align:left;"> C>A:CCA </th>
   <th style="text-align:left;"> ... </th>
   <th style="text-align:left;"> T>G:TTT </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PD3851a </td>
   <td style="text-align:left;"> 31 </td>
   <td style="text-align:left;"> 34 </td>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:left;"> 21 </td>
   <td style="text-align:left;"> 24 </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PD3904a </td>
   <td style="text-align:left;"> 110 </td>
   <td style="text-align:left;"> 91 </td>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:left;"> 87 </td>
   <td style="text-align:left;"> 108 </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> 77 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PD3890a </td>
   <td style="text-align:left;"> 122 </td>
   <td style="text-align:left;"> 112 </td>
   <td style="text-align:left;"> 13 </td>
   <td style="text-align:left;"> 107 </td>
   <td style="text-align:left;"> 99 </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> 50 </td>
  </tr>
</tbody>
</table>

You can create an SNV matrix from a VCF file using the method 
```R
genCountMatrixFromVcf
```
or from a MAF file using the method
```R
genCountMatrixFromMAF
```
from signeR package. See the [documentation](https://bioconductor.org/packages/release/bioc/vignettes/signeR/inst/doc/signeR-vignette.html#toc3) for more details. 

## Columns
The first column needs to contain the sample ID and the other columns contain the 96 trinucleotide contexts.

## Rows
Each row contains the sample ID and the counts for each trinucleotide context.

## Example file

[21 breast cancer](https://raw.githubusercontent.com/TojalLab/signeR/devel/inst/extdata/21_breast_cancers.mutations.txt)

[VCF example](https://raw.githubusercontent.com/TojalLab/signeR/devel/inst/extdata/example.vcf)

[MAF example](https://raw.githubusercontent.com/TojalLab/signeR/devel/inst/extdata/example.maf)
