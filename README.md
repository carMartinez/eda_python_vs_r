A reference guide for exploratory data analysis (EDA) in Python and R. Genetic
data from the breast cancer project of The Cancer Genome Atlas (TCGA) is
used to walk through some routine tasks when you first get your hands on a
new (semi-clean) dataset. See the
[blog post](http://cmartinez.io/eda-python-vs-r/)
for explanations and side-by-side Python vs R comparisons.

To download the data, first download the
[GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)
and the `MANIFEST.txt` file in this repo, then simply run

`./gdc-client download -m MANIFEST.txt`

from a terminal. Any dataset from TCGA could be used instead of breast
cancer as well.
