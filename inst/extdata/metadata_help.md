# Metadata matrix

Clinical data is a tab-delimited text file with a matrix of available metadata (clinical and/or survival) for each sample. It must have a first column of sample IDs, named “SampleID”, whose entries match the row names of the **SNV matrix**. The number and title of the remaining columns are optional, however if *survival* data is included it must be organized in a column named **time** (in months) and another named **status** (which contains 1 for death events and 0 for censored samples). The table below shows an example of the clinical data matrix structure.

<table class="table table-striped table-hover" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> SampleID </th>
   <th style="text-align:left;"> gender </th>
   <th style="text-align:left;"> ajcc_pathologic_stage </th>
   <th style="text-align:left;"> ethnicity </th>
   <th style="text-align:left;"> race </th>
   <th style="text-align:left;"> status </th>
   <th style="text-align:left;"> time </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PD3851a </td>
   <td style="text-align:left;"> male </td>
   <td style="text-align:left;"> Stage I </td>
   <td style="text-align:left;"> not hispanic or latino </td>
   <td style="text-align:left;"> white </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 236 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PD3890a </td>
   <td style="text-align:left;"> male </td>
   <td style="text-align:left;"> Stage II </td>
   <td style="text-align:left;"> not hispanic or latino </td>
   <td style="text-align:left;"> black or african american </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 199 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PD3904a </td>
   <td style="text-align:left;"> female </td>
   <td style="text-align:left;"> Stage II </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 745 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PD3905a </td>
   <td style="text-align:left;"> female </td>
   <td style="text-align:left;"> Stage IV </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> white </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 299 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PD3945a </td>
   <td style="text-align:left;"> male </td>
   <td style="text-align:left;"> Stage IV </td>
   <td style="text-align:left;"> not hispanic or latino </td>
   <td style="text-align:left;"> asian </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 799 </td>
  </tr>
</tbody>
</table>


## Columns
The first column must contain the sample ID. Other columns may contain sample groupings or other features that you would like to co-analyze with exposure data.

## Rows
Each row contains clinical information for one sample: its ID and all other data of interest.

## Example file

[21 breast cancer](https://raw.githubusercontent.com/TojalLab/signeR/devel/inst/extdata/clinical-test-signerflow.tsv)
