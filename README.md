# vcf-annotator

vcf-annotator is a method for annotating vcf files using the [Ensembl Rest API](http://useast.ensembl.org/index.html).

Inside this Github repository, the annotated variants are included as the file "output_annotated_vcf_data.csv"

I have also included the original, unannotated vcf file as "test_vcf_data.txt" for the reviewer's convenience.

# Requirements
* Install `git` from [https://git-scm.com](https://git-scm.com/).
* Install `python3` from [https://www.python.org/](https://www.python.org/) (I used version 3.9.16).
* Install `pip3` from [https://pip.pypa.io/en/stable/installing/](https://pip.pypa.io/en/stable/installing/).

# Installation

Install `vcf-annotator` by cloning this github repository and then using pip to install the program locally:

```sh
git clone git@github.com:aps120/variant-annotator.git
cd variant-annotator
pip install ./
```
# Basic usage

Currently, this program only supports vcf files for which genomic coordinates are assigned based
on the human assembly GRCh37 (hg19). So genomic coordinates must be given in reference to GRCh37 
(if you want to use GRCh38, edit line 61 in `functions.py` to point to the GRCh38 Ensembl Server). 
You will also need to be connected to the internet to work with the Ensembl Rest API.

To annotate your vcf, type the following into your terminal:


```sh
vcf-annotator <input vcf> <output prefix>
```

This will produce a csv file named `<output prefix>_annotated_vcf_data.csv` containing the following columns:
* Chrom: chromosome number
* Pos: genomic position
* Ref: reference sequence
* Alt: sequence of variant(s)
* hgvs of variant: variant notated in [standard HGVS nomenclature](https://varnomen.hgvs.org/bg-material/simple/)
* Depth of coverage: Depth of sequence coverage at the site of variation
* No. reads for variant: number of reads supporting the variant
* % reads supporting variant: Percentage of reads supporting variant (or variant frequency as a percentage)
* gene: gene of variant in [HGNC](https://www.genenames.org/) nomenclature (intergenic regions labeled "intergenic")
* variation: the type of variation for variant
* effect: effect of variant
* min allele freq (%): The minor allele frequency as a percent ('n/a' if there's no minor allele)
* min allele hgsv: minor allele notated in [standard HGVS nomenclature](https://varnomen.hgvs.org/bg-material/simple/) ('n/a' if there's no minor allele)
* min allele effect: effect of the minor allele (n/a if there's no minor allele)  

To reproduce my annotated variant csv file with the test data, type the following into your terminal:
 
 ```sh
vcf-annotator test_vcf_data.txt test
```
This should produce the same annotated variant csv file with the name  `test_annotated_vcf_data.csv`.
It took about ~2 h to run this on my Macbook Pro (The step which slows this down is downloading data).