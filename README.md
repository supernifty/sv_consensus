## VCF consensus

Given multiple VCF files containing SNV or SV calls, build a consensus for evaluation.

## Install
```
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```
## Usage

The following functionality is available:
* snv_concordance.py: count agreement for each SNP across VCFs
* snv_concordance_vcf.py: generate a VCF of variants with a minimum level of concordance
* snv_intersect.py: find how many variants are in common across every combination of VCFs
