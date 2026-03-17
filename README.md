# Tumor Mutational Burden

**Institut Curie - TMB analysis**

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Version](https://img.shields.io/badge/version-1.6.0-green.svg)](CHANGELOG)
[![PyPI Package](https://img.shields.io/badge/PyPI-pytmb--curie-blue.svg)](https://pypi.org/project/pytmb-curie/)

This tool was designed to calculate a **Tumor Mutational Burden (TMB)** score from a VCF file.

The TMB is usually defined as the total number of non-synonymous mutations per coding area of a tumor genome. This metric is mainly used as a biomarker in clinical practice to determine whether to use or not immunomodulatory anticancer drugs (immune checkpoint inhibitors such as Nivolumab). Whole Exome Sequencing (WES) allows comprehensive measurement of TMB and is considered the gold standard. In practice, as TMB is mainly used in the routine of the clinic, and due to the high cost of WES, TMB calculation based on gene panels is preferred.

Currently, the main limitation of TMB calculation is the lack of standard for its calculation. Therefore, we decided to propose a very **versatile** tool allowing the user to define exactly which type of variants to use or filter.

---

## Tool summary

### Installation

#### Option 1 — pip from PyPI (recommended)

```bash
# Install the package
pip install pytmb-curie

# To also enable pyEffGenomeSize (requires pybedtools + pandas)
pip install pytmb-curie[effgenomesize]
```

After installation two CLI commands are available directly in your `$PATH`:

```
pyTMB
pyEffGenomeSize
```

#### Option 2 — conda environment + pip (development)

```bash
# 1. Create and activate the conda environment
conda env create -f environment.yml
conda activate pyTMB_new

# 2. Install the package (regular)
pip install .

# 3. Or install in editable / development mode
pip install -e .

# 4. To also enable pyEffGenomeSize (requires pybedtools + pandas)
pip install -e ".[effgenomesize]"
```

After installation two CLI commands are available directly in your `$PATH`:

```
pyTMB
pyEffGenomeSize
```

#### Option 3 — conda only (pre-built bioconda package)

```bash
conda env create -n pytmb
conda activate pytmb
conda install -c bioconda -c conda-forge tmb=1.6.0
```

#### Option 4 — run scripts directly (backward-compatible shims)

The `bin/` directory still contains thin wrapper scripts that delegate to the
installed package.  After installing with `pip install .` you can still call:

```bash
python bin/pyTMB.py [options]
python bin/pyEffGenomeSize.py [options]
```

### Recommendations

In order to have homogenous VCF entry files and to avoid VCF ambiguities, we recommend to normalize the VCF files before calculating the TMB. This is especially useful if the VCF file contains Multi Nucleotide Variants (MNVs) or multiallelic variants.

```bash
bcftools norm -f FASTA -m -o file_norm.vcf file
```

### Implementation

The idea behind this tool is quite simple. All variants are scanned and filtered according to the criteria provided by the user. If a variant passes all the filters, it is therefore used for the TMB calculation. In other words, if no filters are provided, the tool will simply count the number of variants.

The TMB is defined as the number of variants over the size of the genomic region (in Mb).  
To calculate the effective genome size, the user can provide a BED file (`--bed`) with the design of the assay.  
This BED file should be ordered, 0-based and with no header.  
Another alternative is to specify the size directly using `--effGenomeSize`.  
Importantly, **it is the user's responsibility to provide the BED corresponding to the VCF input file.**

We also provide the `pyEffGenomeSize` command to calculate the effective size from a BAM file using annotations, coverage and mapping quality thresholds defined by the user.

---

## Package structure

Since version 1.6.0, the code is organized as an installable Python package:

```
pytmb/
├── __init__.py           # version + public API
├── config.py             # loadConfig()
├── vcf_utils.py          # getTag(), getMultiAlleleHeader()
├── genome_size.py        # getEffGenomeSizeFromBed(), getEffGenomeSizeFromMosdepth()
├── filters.py            # isAnnotatedAs(), isPolym(), isCancerHotspot(), …
├── tmb.py                # calculate_tmb() — importable library function
└── cli/
    ├── run_tmb.py        # pyTMB entry point
    └── run_effgenomesize.py  # pyEffGenomeSize entry point
```

The core `calculate_tmb()` function can also be used programmatically:

```python
from pytmb import calculate_tmb, loadConfig

db_flags     = loadConfig("config/snpeff.yml")
caller_flags = loadConfig("config/mutect2.yml")

results = calculate_tmb(
    vcf_path="sample.vcf.gz",
    db_flags=db_flags,
    caller_flags=caller_flags,
    eff_genome_size=33_280_000,
    vaf=0.05,
    maf=0.001,
    min_depth=20,
    min_alt_depth=2,
    filter_low_qual=True,
    filter_non_coding=True,
    filter_syn=True,
    filter_polym=True,
    polym_db="1k,gnomad",
)
print("TMB =", results["tmb"])
```

---

## Quick help

```bash
pyTMB -h

usage: pyTMB [-h] -i VCF --dbConfig DBCONFIG --varConfig VARCONFIG
             [--sample SAMPLE] [--effGenomeSize EFFGENOMESIZE] [--bed BED]
             [--vaf VAF] [--maf MAF] [--minDepth MINDEPTH]
             [--minAltDepth MINALTDEPTH] [--filterLowQual] [--filterIndels]
             [--filterCoding] [--filterSplice] [--filterNonCoding]
             [--filterSyn] [--filterNonSyn] [--filterCancerHotspot]
             [--filterPolym] [--filterRecurrence] [--polymDb POLYMDB]
             [--cancerDb CANCERDB] [--verbose] [--debug] [--export EXPORT]
             [--version]

Calculate a Tumour Mutational Burden (TMB) score from a VCF file.

options:
  -h, --help            show this help message and exit
  -i VCF, --vcf VCF     Input file (.vcf, .vcf.gz, .bcf, .bcf.gz) (default: None)
  --dbConfig DBCONFIG   Databases config file (YAML) (default: None)
  --varConfig VARCONFIG Variant calling config file (YAML) (default: None)
  --sample SAMPLE       Specify the sample ID to focus on (default: None)
  --effGenomeSize EFFGENOMESIZE
                        Effective genome size (bp) (default: None)
  --bed BED             Capture design BED file (default: None)
  --vaf VAF             Filter variants with Allelic Ratio < vaf (default: 0)
  --maf MAF             Filter variants with MAF > maf (default: 1)
  --minDepth MINDEPTH   Filter variants with depth < minDepth (default: 1)
  --minAltDepth MINALTDEPTH
                        Filter variants with alt-allele depth < minAltDepth (default: 1)
  --filterLowQual       Filter low quality (not PASS) variants (default: False)
  --filterIndels        Filter insertions/deletions (default: False)
  --filterCoding        Filter coding variants (default: False)
  --filterSplice        Filter splice variants (default: False)
  --filterNonCoding     Filter non-coding variants (default: False)
  --filterSyn           Filter synonymous variants (default: False)
  --filterNonSyn        Filter non-synonymous variants (default: False)
  --filterCancerHotspot Filter variants annotated as cancer hotspots (default: False)
  --filterPolym         Filter polymorphism variants (see --maf) (default: False)
  --filterRecurrence    Filter on run-level recurrence values (default: False)
  --polymDb POLYMDB     Databases for polymorphism detection, comma-separated (default: gnomad)
  --cancerDb CANCERDB   Databases for cancer hotspot annotation, comma-separated (default: cosmic)
  --verbose             Activate verbose mode (default: False)
  --debug               Export original VCF with TMB_FILTER tag (default: False)
  --export EXPORT       Export a VCF with passing variants to this path (default: None)
  --version             Version number
```

---

## Configs

Working with VCF files is usually not straightforward, and mainly depends on the variant caller and annotation tools/databases used.
In order to make this tool as flexible as possible, we set up **two configuration files** to define which fields have to be checked and in which case.

The `--dbConfig` file describes all details about annotation. We provide configurations for:
- **Annovar** — `config/annovar.yml`
- **snpEff** — `config/snpeff.yml`
- **VEP** — `config/vep.yml`

These files can be customized by the user.

The `--varConfig` file contains all variant-caller-specific parameters. Config files for:
- **Varscan2** — `config/varscan2.yml`
- **Mutect2** — `config/mutect2.yml`
- **Strelka** — `config/strelka.yml`

are provided as examples.

The `yaml` config files list the different **key:values** for each function.
For example, to assess whether a variant is coding (for Annovar):

```yaml
isCoding:
  Func.refGene:
    - exonic
```

Regarding databases, the polymorphism fields for Annovar are:

```yaml
polymDb:
  1k:
    - 1000g2015aug_all
  gnomad:
    - gnomAD_exome_ALL
  esp:
    - esp6500siv2_all
  exac:
    - ExAC_ALL
```

The user can then choose databases with `--polymDb 1k,gnomad,esp,exac`.  
The same logic applies for `--cancerDb`.

---

## Usage

### `pyTMB`

#### General parameters

##### `-i`
Input file (`.vcf`, `.vcf.gz`, `.bcf`, `.bcf.gz`)

##### `--sample`
Specify the sample ID to focus on. Required when dealing with multi-sample VCFs.

##### `--bed` and `--effGenomeSize`
Specify either a sorted BED file with no header, or the size of the effective genome directly.

### Filters

#### `--vaf MINVAF`
Filter variants with Allelic Ratio < minVAF. The field name is defined in `config/caller.yml`.
The tool first checks the **FORMAT** field and then falls back to the **INFO** field.

#### `--maf MAXMAF`
Filter variants with MAF > maxMAF. The databases to use are set with `--polymDb`
and the `config/databases.yml` file.

#### `--minDepth MINDEPTH`
Filter variants with depth < minDepth. The field name is defined in `config/caller.yml`.
The tool first checks the **FORMAT** field and then falls back to the **INFO** field.

#### `--minAltDepth MINALTDEPTH`
Filter variants with alternative allele depth < minAltDepth. Checked in the **FORMAT** field.

#### `--filterLowQual`
Filter variants for which the **FILTER** field is not **PASS** or for which the **QUAL** value is not null.

#### `--filterIndels`
Filter insertion/deletion variants.

#### `--filterCoding`
Filter coding variants as defined in the `config/databases.yml` file.

#### `--filterSplice`
Filter splice variants as defined in the `config/databases.yml` file.

#### `--filterNonCoding`
Filter non-coding variants as defined in the `config/databases.yml` file.

#### `--filterSyn`
Filter synonymous variants as defined in the `config/databases.yml` file.

#### `--filterNonSyn`
Filter non-synonymous variants as defined in the `config/databases.yml` file.

#### `--filterCancerHotspot`
Filter variants annotated as cancer hotspots as defined in the `config/databases.yml` file.
All variants with a cancer annotation (e.g. a COSMIC ID) will be removed.

#### `--filterPolym`
Filter polymorphism variants from genome databases. The databases can be listed with `--polymDb`.
The fields to scan for each database are defined in `config/databases.yml` and the population
frequency is compared against `--maf`.

#### `--filterRecurrence`
Filter on run-level recurrence values. The VCF must already contain recurrence information
as defined in the `config/databases.yml` file.

### Outputs

By default, the tool prints a summary with the calculated TMB value.

#### `--export PATH`
Export a VCF file containing only the variants used for TMB calculation.

#### `--debug`
Export a VCF file with the tag **TMB_FILTERS** in the **INFO** field. This tag contains the
reason why each variant would be filtered.

---

### `pyEffGenomeSize`

This tool calculates the effective genome size from a BAM file. This parameter has a strong
impact on the TMB result. For instance, if only coding variants are used it makes sense to
restrict the denominator to coding regions only.

```bash
pyEffGenomeSize -h

usage: pyEffGenomeSize [-h] --bed BED [--gtf GTF] [--bam BAM]
                       [--minCoverage MINCOVERAGE] [--minMapq MINMAPQ]
                       [--filterNonCoding] [--filterCoding]
                       [--featureTypes FEATURETYPES [FEATURETYPES ...]]
                       [--saveIntermediates] [-t THREAD]
                       [--oprefix OPREFIX] [--verbose] [--version]
```

#### General parameters

##### `--bed`
The input BED file to filter. Should be 0-based, sorted, and with no header. **Required.**

##### `--gtf`
A sorted GTF file for genome annotation (e.g. `gencode.v19.annotation.gtf` or `.gtf.gz`).

##### `--bam`
A BAM file from your experiment to extract mapping quality and coverage information.
When provided, mosdepth is automatically run to filter regions based on coverage and mapping quality.

#### Filters

##### `--minCoverage`
Minimum coverage per region of the BED file. Requires `--bam`.

##### `--minMapq`
Mapping quality threshold. Reads below this value are ignored. Requires `--bam`.

##### `--filterNonCoding`
Remove regions considered non-coding from the GTF/BED intersection to keep only exonic regions.
Requires `--gtf`.

##### `--filterCoding`
Remove regions considered coding based on the `transcript_type` field in the GTF.
Requires `--gtf` and `--featureTypes`.

##### `--featureTypes`
Choose one or more feature types from `exon`, `gene`, `transcript`, `UTR`, `CDS` to
retain in the final BED file. Default: `exon`. Required with `--filterCoding`.

#### Output / Misc

##### `--saveIntermediates`
Keep intermediate files (mosdepth output, filtered GTF, intersect BED) instead of deleting them.

##### `-t`, `--thread`
Number of threads for mosdepth. Default: 1.

##### `--oprefix`
Output file prefix. Default: `pyeffg`.

##### `--verbose`
Activate verbose mode.

##### `--version`
Show version number.

---

## Usage examples and recommendations

### Gene Panel

Calculate the TMB on a gene panel VCF (coding size = 1.59 Mb, caller = Varscan2, annotation = Annovar)
with the following criteria:
- minDepth at 100×
- non-synonymous, coding and splice variants only
- polymorphisms filtered using 1K, gnomAD databases at MAF 0.001 (0.1%)

```bash
pyTMB -i ${VCF} \
    --effGenomeSize 1590000 \
    --dbConfig config/annovar.yml \
    --varConfig config/varscan2.yml \
    --maf 0.001 --minDepth 100 --minAltDepth 2 \
    --filterLowQual \
    --filterNonCoding \
    --filterSplice \
    --filterSyn \
    --filterPolym --polymDb 1k,gnomad \
    --export panel_tmb_variants.vcf.gz \
    > TMB_results.log
```

### Whole Exome Sequencing

For WES, filter low-quality, non-coding, synonymous and polymorphic variants.
Indels and splicing variants are kept. An effective genome size of 33 Mb is used.

For Mutect2 + snpEff:

```bash
pyTMB -i ${VCF} \
    --effGenomeSize 33280000 \
    --dbConfig config/snpeff.yml \
    --varConfig config/mutect2.yml \
    --vaf 0.05 --maf 0.001 --minDepth 20 --minAltDepth 2 \
    --filterLowQual \
    --filterNonCoding \
    --filterSyn \
    --filterPolym --polymDb 1k,gnomad \
    > TMB_results.log
```

---

## Credits

This pipeline has been written by the bioinformatics core facility in close collaboration with the Clinical Bioinformatics and the Genetics Service of the Institut Curie. Many thanks to the seqOIA-IT team for their help in the development and for the extensive testing of the tool!

If you are using this tool for your own research, please cite:  
Dupain, C., Gutman, T., Girard, E. et al. *Tumor mutational burden assessment and standardized bioinformatics approach using custom NGS panels in clinical routine.* BMC Biol 22, 43 (2024). https://doi.org/10.1186/s12915-024-01839-8

## AI Disclosure: Augmented

This project is **AI-augmented** and utilized AI (e.g., Claude) to:
* **Generate** boilerplate code and specific utility functions.
* **Refactor** existing code for better performance and readability.
* **Draft** unit tests and technical documentation.

**Verification:** Every AI-generated contribution was manually reviewed, debugged, and integrated into the final codebase.

## Contacts

For any question, bug or suggestion, please use the issues system or contact the bioinformatics core facility.
