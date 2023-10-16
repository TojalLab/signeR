# Opportunity matrix

The opportunity matrix is a tab-delimited text file matrix of counts of trinucleotide contexts found in the studied genomes. It must structured as the SNV matrix, with mutations specified on the columns (for each SNV count, the Opportunity matrix shows the total number of genomic loci where the referred mutation could have occurred). The Opportunity matrix can be generated from a BSgenome, see signeR documentation. The table below shows an example of the opportunity matrix structure.

<table class="table table-striped table-hover" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> 366199887 </th>
   <th style="text-align:left;"> 211452373 </th>
   <th style="text-align:left;"> 45626142 </th>
   <th style="text-align:left;"> 292410567 </th>
   <th style="text-align:left;"> 335391892 </th>
   <th style="text-align:left;"> 239339768 </th>
   <th style="text-align:left;"> ... </th>
   <th style="text-align:left;"> 50233875 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 202227618 </td>
   <td style="text-align:left;"> 116207171 </td>
   <td style="text-align:left;"> 25138239 </td>
   <td style="text-align:left;"> 161279580 </td>
   <td style="text-align:left;"> 184193767 </td>
   <td style="text-align:left;"> 131051208 </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> 177385805 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 225505378 </td>
   <td style="text-align:left;"> 130255706 </td>
   <td style="text-align:left;"> 28152934 </td>
   <td style="text-align:left;"> 179996700 </td>
   <td style="text-align:left;"> 206678032 </td>
   <td style="text-align:left;"> 147634427 </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> 199062504 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 425545790 </td>
   <td style="text-align:left;"> 245523433 </td>
   <td style="text-align:left;"> 53437284 </td>
   <td style="text-align:left;"> 339065644 </td>
   <td style="text-align:left;"> 389386002 </td>
   <td style="text-align:left;"> 278770926 </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> 375075216 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 452332390 </td>
   <td style="text-align:left;"> 259934779 </td>
   <td style="text-align:left;"> 55862550 </td>
   <td style="text-align:left;"> 361010972 </td>
   <td style="text-align:left;"> 412168035 </td>
   <td style="text-align:left;"> 292805460 </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> 396657807 </td>
  </tr>
</tbody>
</table>

You can create an opportunity matrix from the reference genome using the method 
```R
genOpportunityFromGenome
```
from signeR package. See the [documentation](https://bioconductor.org/packages/release/bioc/vignettes/signeR/inst/doc/signeR-vignette.html#toc3) for more details. 

## Columns
There is no header in this file and each column represents a trinucleotide context.

## Rows
Each row contains the count frequency of the trinucleotides in the whole analyzed region for each sample.

## Example file

[21 breast cancer](https://raw.githubusercontent.com/TojalLab/signeR/devel/inst/extdata/21_breast_cancers.opportunity.txt)
