# Previous signatures matrix

Previous signatures is a tab-delimited text file with a matrix of previously known signatures. It must contain one column for each signature and one row for each of the 96 SNV types (considering trinucleotide contexts). Mutation types should be contained on the first column, in the same form as the column names of the SNV matrix. The table below shows a example of the previous signatures matrix structure.

<table class="table table-striped table-hover" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> </th>
   <th style="text-align:left;"> Signature 2 </th>
   <th style="text-align:left;"> Signature 3 </th>
   <th style="text-align:left;"> Signature 5 </th>
   <th style="text-align:left;"> Signature 6 </th>
   <th style="text-align:left;"> ... </th>
   <th style="text-align:left;"> Signature 8 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> C>A:ACA </td>
   <td style="text-align:left;"> 0.01110 </td>
   <td style="text-align:left;"> 0.00067 </td>
   <td style="text-align:left;"> 0.02218 </td>
   <td style="text-align:left;"> 0.01494 </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> 0.03672 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C>A:ACC </td>
   <td style="text-align:left;"> 0.00915 </td>
   <td style="text-align:left;"> 0.00062 </td>
   <td style="text-align:left;"> 0.01788 </td>
   <td style="text-align:left;"> 0.00896 </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> 0.03324 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C>A:ACG </td>
   <td style="text-align:left;"> 0.00150 </td>
   <td style="text-align:left;"> 0.00010 </td>
   <td style="text-align:left;"> 0.00213 </td>
   <td style="text-align:left;"> 0.00221 </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> 0.00252 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> ... </td>
  </tr>
  <tr>
   <td style="text-align:left;"> T>G:TTT </td>
   <td style="text-align:left;"> 0.00403 </td>
   <td style="text-align:left;"> 2.359E-05 </td>
   <td style="text-align:left;"> 0.0130 </td>
   <td style="text-align:left;"> 0.01337 </td>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> 0.00722 </td>
  </tr>
</tbody>
</table>


## Columns
The first columnn needs to contain the trinucleotide contexts and other columns contain the known signatures.

## Rows
Each row contains the expected frequency of the given mutation in the apponted trinucleotide context.

## Example file

[21 breast cancer](https://raw.githubusercontent.com/TojalLab/signeR/devel/inst/extdata/Cosmic_signatures_BRC.txt)
