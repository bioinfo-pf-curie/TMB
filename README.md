# Tumor Mutational Burden

**Institut Curie - TMB analysis**

[![Install with](https://anaconda.org/anaconda/conda-build/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda)


This tool was designed to calculate a **Tumor Mutational Burden (TMB)** score from a VCF file.

The TMB is usually defined as the total number of non-synonymous mutations per coding area of a tumor genome. This metric is mainly used as a biomarker in clinical practice to determine whether to use or not immunomodulatory anticancer drugs (immune checkpoint inhibitors such as Nivolumab). Whole Exome Sequencing (WES) allows comprehensive measurement of TMB and is considered the gold standard. In practice, as TMB is mainly used in the routine of the clinic, and due to the high cost of WES, TMB calculation based on gene panels is preferred.

Currently, the main limitation of TMB calculation is the lack of standard for its calculation. Therefore, we decided to propose a very **versatile** tool allowing the user to define exactly which type of variants to use or filter.

## Tool summary

### Requirements

The tool was implemented in `python3`, and require the librairies `cyvcf2` and `yaml`.  
We provide a `conda` file to build a simple python environment.  
To do so, simply use:

```
conda env create -f environment.yml -p PATH_TO_INSTALL
```

### Recommendations

In order to have homogenous VCF entry files, avoid VCF ambiguities and have reproducible results we recommend to normalize your VCF before running this script, especially if the VCF file contains Multi Nucleotide Variants (MNVs) or multiallelic variants.

For that we suggest to use `bcftools norm -f FASTA -m- -o file_norm.vcf file` command.

### Implementation

The idea behind this script is quite simple. All variants are scanned and filtered according to the criteria provided by the user. If a variant passes all the filters, it is therefore used for the TMB calculation. In other words, if no filters are provided, the script will simply count the number of variants.

The TMB is defined as the number of variants over the size of the genomic region (in Mb). In order to calculate the size of the genome (ie. the `effectiveGenomeSize`), the user can provide a BED file (`--bed`) with the design of the assay. This tool also provides a feature to create a bed file with specific coverage and mapping quality thresholds defined by the user (script in `bin/pyEffGenomeSize.py`).

However, it is usually recommended to adapt the genome size to the filters applied. For instance, if only coding variants are used, it would make sense to use only the genomic size of coding region for the TMB calculation. So far, **this is the user responsability to provide an intial bed file with corresponding genomic features.** and to specify it to the `pyEffGenomeSize.py` script of directly to the `--bed` parameter. The user can also provide the size of the bed with the `--effGenomeSize` parameter.


## Quick help
```bash
python3 ~/Documents/Tom/Pipelines/tmb/bin/pyEffGenomeSize.py -h

usage: pyEffGenomeSize.py [-h] [--bed BED] [--gtf GTF] [--bam BAM]
                          [--minCoverage MINCOVERAGE] [--minMapq MINMAPQ]
                          [--filterNonCoding] [--filterCoding]
                          [--featureTypes FEATURETYPES [FEATURETYPES ...]]
                          [--saveIntermediates] [-t THREAD]
                          [--oprefix OPREFIX] [--verbose] [--version]

optional arguments:
  -h, --help            show this help message and exit
  --bed BED             BED file (.bed) (default: None)
  --gtf GTF             GTF file for genome annotation (.gtf) (default: None)
  --bam BAM             BAM file for mapping statistics (.bam) (default: None)
  --minCoverage MINCOVERAGE
                        minimum coverage of the region (default: 0)
  --minMapq MINMAPQ     minimum coverage of the region (default: 0)
  --filterNonCoding     Filter regions associated with non coding annotations
                        (default: False)
  --filterCoding        Filter regions associated with coding annotations
                        (default: False)
  --featureTypes FEATURETYPES [FEATURETYPES ...]
                        List of Features (exon, gene, transcript, UTR, CDS) to
                        keep (3rd column from gtf/gff). Required with
                        --filterCoding argument (default: [])
  --saveIntermediates   Save mosdepth intermediate files (default: False)
  -t THREAD, --thread THREAD
                        Number of threads for mosdepth run (default: 1)
  --oprefix OPREFIX     Suffix for filtered bed (default: pyeffg)
  --verbose             Active verbose mode (default: False)
  --version             Version number


```

```bash
python bin/pyTMB.py -h

usage: pyTMB.py [-h] [-i VCF] [--dbConfig DBCONFIG] [--varConfig VARCONFIG]
                [--sample SAMPLE] [--effGenomeSize EFFGENOMESIZE] [--bed BED]
                [--vaf VAF] [--maf MAF] [--minDepth MINDEPTH]
                [--minAltDepth MINALTDEPTH] [--filterLowQual] [--filterIndels]
                [--filterCoding] [--filterSplice] [--filterNonCoding]
                [--filterSyn] [--filterNonSyn] [--filterCancerHotspot]
                [--filterPolym] [--filterRecurrence] [--polymDb POLYMDB]
                [--cancerDb CANCERDB] [--verbose] [--debug] [--export]
                [--version]

optional arguments:
  -h, --help            show this help message and exit
  -i VCF, --vcf VCF     Input file (.vcf, .vcf.gz, .bcf) (default: None)
  --dbConfig DBCONFIG   Databases config file (default:
                        ./config/databases.yml)
  --varConfig VARCONFIG
                        Variant calling config file (default:
                        ./config/calling.yml)

  --effGenomeSize EFFGENOMESIZE
                        Effective genome size (default: None)
  --bed BED             Capture design to use if effGenomeSize is not defined
                        (BED file) (default: None)
  --vaf VAF             Filter variants with Allelic Ratio <= vaf (default:
                        0.05)
  --maf MAF             Filter variants with MAF > maf (default: 0.001)
  --minDepth MINDEPTH   Filter variants with depth < minDepth (default: 5)
  --minAltDepth MINALTDEPTH
                        Filter variants with alternative allele depth <=
                        minAltDepth (default: 2)
  --filterLowQual       Filter low quality (i.e not PASS) variant (default:
                        False)
  --filterIndels        Filter insertions/deletions (default: False)
  --filterCoding        Filter Coding variants (default: False)
  --filterSplice        Filter Splice variants (default: False)
  --filterNonCoding     Filter Non-coding variants (default: False)
  --filterSyn           Filter Synonymous variants (default: False)
  --filterNonSyn        Filter Non-Synonymous variants (default: False)
  --filterCancerHotspot
                        Filter variants annotated as cancer hotspots (default:
                        False)
  --filterPolym         Filter polymorphism variants in genome databases. See
                        --maf (default: False)
  --filterRecurrence    Filter on recurrence values (default: False)

  --polymDb POLYMDB     Databases used for polymorphisms detection (comma
                        separated) (default: gnomad)
  --cancerDb CANCERDB   Databases used for cancer hotspot annotation (comma
                        separated) (default: cosmic)

  --sample SAMPLE       Specify the sample ID to focus on (default: None)
  --debug               Export original VCF with TMB_FILTER tag (default:
                        False)
  --export              Export a VCF with the considered variants (default:
                        False)
  --verbose             Active verbose mode (default: False)
  --version             Version number
```

## Configs

Working with vcf files is usually not straighforward, and mainly depends on the variant caller and annotation tools/databases used.
In order to make this tool as flexible as possible, we decided to set up **two configurations files** to defined which fields have to be checked and in which case.

The `--dbConfig` file described all details about annotation. As an exemple, we provide some configurations for **Annovar** (*config/annovar.yml*)
and **snpEff** (*config/snpeff.yaml*) tool.  
These files can be customized by the user.

In the same way, all parameters which are variant caller specific can be set up in another config file using the `--varConfig` parameter.
Config files for **Varscan2** (*config/varscan2.yml*) and **Mutect2** (*config/mutect2.yml*) are provided as examples.

The `yaml` config files must list the different **keys:values** for each function.
As an exemple, to assess whether a variant is coding, the programm will used (for Annovar);

```
isCoding:
  Func.refGene:
    - exonic			
```

It will therefore search in the `INFO` field the key `Func.refGene` and the value `exonic`.

Regarding the databases, this is the same idea. Here is the list of databases, and fields used to check the `MAF` value by default (for Annovar) :

```
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

The user can thus choose to scan all databases with the `--polymDb 1k,gnomad,esp,exac` parameter.  
The same is true for the `--cancerDb` parameter.

## Usage

## `pyTMB.py`:
### Filters

#### `--vaf MINVAF`
Filter variants with Allelic Ratio < minVAF. Note the field used to get the Allelic Ratio field is defined in the *config/caller.yml* file.
In this case, the programm will first look for this information in the **FORMAT** field, and then in the **INFO** field.

#### `--maf MAXMAF`
Filter variants with MAF < maf. Note the databases used to check the Min Allele Frequency are set using the `--polymDb`
parameters and the *config/databases.yml* file.

#### `--minDepth MINDEPTH`
Filter variants with depth < minDepth. Note the field used to get the depth is defined in the *config/caller.yml* file.
In this case, the programm will first look for this information in the **FORMAT** field, and then in the **INFO** field.

#### `--minAltDepth MINALTDEPTH`
FIlter alternative allele with depth < minAltDepth. The programm will look in the **FORMAT** field exclusively.

#### `--filterLowQual`
Filter variants for which is the **FILTER** field is not **PASS** or for which the **QUAL** value is not null.

#### ` --filterIndels`
Filter insertions/deletions variants.

#### `--filterCoding`
Filter Coding variants as defined in the *config/databases.yml* field.

#### `--filterSplice`
Filter Splice variants as defined in the *config/databases.yml* field.

#### `--filterNonCoding`
Filter Non-coding variants as defined in the *config/databases.yml* field.

#### `--filterSyn`
Filter Synonymous variants as defined in the *config/databases.yml* field.

#### `--filterNonSyn`
Filter Non-Synonymous variants as defined in the *config/databases.yml* field.

#### `--filterCancerHotspot`
Filter variants annotated as cancer hotspots as defined in the *config/databases.yml* field.
So far, all variants with a 'cancer' annotation (for instance with a COSMIC Id) will be removed.

#### `--filterPolym`
Filter polymorphism variants from genome databases. The databases to considered can be listed with the `--polymDb` parameter.
The fields to scan for each database are defined in the *config/databases.yml* file and the population frequency is compared with the `--minMAF` field.

#### `--filterRecurrence`
Filter on recurrence values (for instance, intra-run occurence). In this case, the vcf file must contains the recurrence information
which can be defined the *config/databases.yml* file.

## Outputs

By default, the script outputs a few information with the calculated TMB value.

#### `--export`

This option allows to export a vcf file which only contains the variants used for TMB calculation.

#### `--debug`

The option allows to export a vcf file with the tag **TMB_FILTERS** in the **INFO** field. This tag therefore contains the reason for which a variant would be filtered.

## `pyTMB.py`:

### Filters

#### `--minCoverage`

Define the minimum coverage accepted for each region of the bed file

#### `--minMapq`

Mapping quality threshold. reads with a mapping quality less than this are ignored

#### `--filterNonCoding`

This filter removes regions considered as non coding from the gtf and bed files to only keep exonic regions.

#### `--filterCoding`

This filter removes regions considered as coding based on the transcript_type field in the gtf.
This filter **requires** the parameter `featureTypes`

#### `--featureTypes`

This parameter offers the possibility to choose one or multiple features to select from the following ("exon", "gene", "transcript", "UTR", "CDS") to keep in the final bed file.

## Usage and recommendations

Here is a list of recommended parameters for different user cases.

### Gene Panel

Let's calculated the TMB on a gene panel vcf file (coding size = 1.59Mb, caller = varscan, annotation = Annovar) with the following criteria:
- minDepth at 100X
- non-synonymous
- coding and splice
- non polymorphism variants using 1K, gnomAD databases and a MAF of 0.001

In this case, a typical usage would be :

```
python pyTMB.py -i ${VCF} --effGenomeSize 1590000 \
--dbConfig config/annovar.yml \
--varConfig config/varscan.yml \
--minMAF 0.001 --minDepth 100 --minAltDepth 2\
--filterLowQual \
--filterNonCoding \
--filterSplice \
--filterSyn \
--filterPolym --polymDb 1k,gnomad \
--export > TMB_results.log
```

### Exome / Whole Genome Sequencing

For WES, we recommend filtering low quality, non coding, synonymous, polymorphic variants. Here, indels and splicing variants are kept. For WES, an effective Genome size of 33Mb is used but a tailored size depending on the variants and regions is preferred.


In the case of a WES variant calling using Mutect2 as variant caller and Snpeff as annotation tool, a typical usage would be :

```
python pyTMB.py -i ${VCF} --effGenomeSize 33280000 \
--dbConfig config/snpeff.yml \
--varConfig config/mutect2.yml \
--vaf 0.05 --maf 0.001 --minDepth 20 --minAltDepth 2 \
--filterLowQual \
--filterNonCoding \
--filterSyn \
--filterPolym --polymDb 1k,gnomad  > TMB_results.log
```

### Credits

This pipeline has been written by the bioinformatics core facility in close collaboration with the Clinical Bioinformatics and the Genetics Service of the Institut Curie.

### Contacts

For any question, bug or suggestion, please use the issues system or contact the bioinformatics core facility.
