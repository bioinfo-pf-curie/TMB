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

### Implementation

The idea behind this script is quite simple. All variants are scanned and filtered according to the criteria provided by the user. If a variant passes all the filters, it is therefore used for the TMB calculation. In other words, if no filters are provided, the script will simply count the number of variants.

The TMB is defined as the number of variants over the size of the genomic region (in Mb). In order to calculate the size of the genome (ie. the `effectiveGenomeSize`), the user can provide a BED file (`--bed`) with the design of the assay. However, it is usually recommanded to adapt the genome size to the filters applied. For instance, if only coding variants are used, it would make sense to use only the genomic size of coding region for the TMB calculation. So far, **this is the user responsability to calculate the appropriate genome size** and to specify it with the `--effGenomeSize` parameter.


## Quick help

```bash
python bin/pyTMB.py --help
usage: pyTMB.py [-h] [-i VCF] [--dbConfig DBCONFIG] [--varConfig VARCONFIG] 
                [--effGenomeSize EFFGENOMESIZE] [--bed BED]
                [--minVAF MINVAF] [--minMAF MINMAF] [--minDepth MINDEPTH] [--minAltDepth MINALTDEPTH]
                [--filterLowQual] [--filterIndels] [--filterCoding]
                [--filterSplice] [--filterNonCoding] [--filterSyn]
                [--filterNonSyn] [--filterCancerHotspot] [--filterPolym]
                [--filterRecurrence] [--polymDb POLYMDB] [--cancerDb CANCERDB]
                [--verbose] [--debug] [--export] [--version]

optional arguments:
  -h, --help                          Show this help message and exit
  -i VCF, --vcf VCF                   Input file (.vcf, .vcf.gz, .bcf)
  --dbConfig DBCONFIG                 Databases config file
  --varConfig VARCONFIG               Variant calling config file

  --effGenomeSize EFFGENOMESIZE       Effective genome size
  --bed BED                           Capture design to use if effGenomeSize is not defined (BED file)

  --minVAF MINVAF                     Filter variants with Allelic Ratio < minVAF
  --minMAF MINMAF                     Filter variants with MAF < minMAF
  --minDepth MINDEPTH                 Filter variants with depth < minDepth
  --minAltDepth MINALTDEPTH           FIlter alternative allele with depth < minAltDepth
  --filterLowQual                     Filter low quality (i.e not PASS) variant
  --filterIndels                      Filter insertions/deletions
  --filterCoding                      Filter Coding variants
  --filterSplice                      Filter Splice variants
  --filterNonCoding                   Filter Non-coding variants
  --filterSyn                         Filter Synonymous variants
  --filterNonSyn                      Filter Non-Synonymous variants
  --filterCancerHotspot               Filter variants annotated as cancer hotspots
  --filterPolym                       Filter polymorphism variants in genome databases. See --minMAF
  --filterRecurrence                  Filter on recurrence values
  
  --polymDb POLYMDB                   Databases used for polymorphisms detection (comma separated)
  --cancerDb CANCERDB                 Databases used for cancer hotspot annotation (comma separated)
  
  --verbose
  --debug
  --export
  --version
		  
```

## Configs

Working with vcf files is usually not straighforward, and mainly depends on the variant caller and annotation tools/databases used.
In order to make this tool as flexible as possible, we decided to set up **two configurations files** to defined with fields as to be checked and in which case.

The `--dbConfig` file described all details about annotation. As an exemple, we provide some configurations for **Annovar** (*conf/annovar.yml*) 
and **snpEff** (*conf/snpeff.yaml*) tool.  
These files can be customized by the user.

In the same way, all parameters which are variant caller specific can be set up in another config file using the `--varConfig` parameter.
Config files for **Varscan2** (*conf/varscan2.yml*) and **Mutect2** (*conf/mutect2.yml*) are provided as examples.

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

### Filters

#### `--minVAF MINVAF`
Filter variants with Allelic Ratio < minVAF. Note the field used to get the Allelic Ratio field is defined in the *conf/caller.yml* file.
In this case, the programm will first look for this information in the **FORMAT** field, and then in the **INFO** field.

#### `--minMAF MINMAF`
Filter variants with MAF < minMAF. Note the databases used to check the Min Allele Frequency are set using the `--polymDb` 
parameters and the *conf/databases.yml* file.

#### `--minDepth MINDEPTH`
Filter variants with depth < minDepth. Note the field used to get the depth is defined in the *conf/caller.yml* file. 
In this case, the programm will first look for this information in the **FORMAT** field, and then in the **INFO** field. 

#### `--minAltDepth MINALTDEPTH`
FIlter alternative allele with depth < minAltDepth. The programm will look in the **FORMAT** field exclusively.

#### `--filterLowQual`
Filter variants for which is the **FILTER** field is not **PASS** or for which the **QUAL** value is not null.

#### ` --filterIndels`
Filter insertions/deletions variants.

#### `--filterCoding`
Filter Coding variants as defined in the *conf/databases.yml* field.

#### `--filterSplice`
Filter Splice variants as defined in the *conf/databases.yml* field.

#### `--filterNonCoding`
Filter Non-coding variants as defined in the *conf/databases.yml* field.

#### `--filterSyn`
Filter Synonymous variants as defined in the *conf/databases.yml* field.

#### `--filterNonSyn`
Filter Non-Synonymous variants as defined in the *conf/databases.yml* field.

#### `--filterCancerHotspot`
Filter variants annotated as cancer hotspots as defined in the *conf/databases.yml* field.
So far, all variants with a 'cancer' annotation (for instance with a COSMIC Id) will be removed.

#### `--filterPolym`
Filter polymorphism variants from genome databases. The databases to considered can be listed with the `--polymDb` parameter.
The fields to scan for each database are defined in the *conf/databases.yml* file and the population frequency is compared with the `--minMAF` field.

#### `--filterRecurrence`
Filter on recurrence values (for instance, intra-run occurence). In this case, the vcf file must contains the recurrence information 
which can be defined the *conf/databases.yml* file.

## Outputs

By default, the script outputs a few information with the calculated TMB value.

#### `--export`

This option allows to export a vcf file which only contains the variants used for TMB calculation.

#### `--debug`

The option allows to export a vcf file with the tag **TMB_FILTERS** in the **INFO** field. This tag therefore contains the reason for which a variant would be filtered.


## Usage and recommandations

Here is a list of recommanded parameters for different user cases.

### Gene Panel

Let's calculated the TMB on a gene panel vcf file (coding size = 1.59Mb, caller = varscan, annotation = Annovar) with the following criteria: 
- minDepth at 100X
- non-synonymous
- coding and splice
- non polymorphism variants using 1K, gnomAD databases and a MAF of 0.001

In this case, a typical usage would be :

```
python pyTMB.py -i ${VCF} \
--dbConfig conf/annovar.yml \
--varConfig conf/varscan.yml \
--minDepth 100 \
--filterLowQual \
--filterNonCoding \
--filterSplice \
--filterSyn \
--filterPolym --minMAF 0.001 --polymDb 1k,gnomad \
--effGenomeSize 1590000 \
--export > TMB_results.log
```

### Exome / Whole Genome Sequencing

For larger panel, we recommand the following parameters ;

```
python pyTMB.py -i ${VCF} \
--dbConfig conf/snpeff.yml \
--varConfig conf/mutect2.yml \
--minVAF 0.05 \
--filterLowQual \
--minDepth 100 \
--filterNonCoding \
--filterIndel \
--filterSyn \
--filterPolym --minMAF 0.001 --polymDb 1k,gnomad
```

### Credits

This pipeline has been written by the bioinformatics core facility in close collaboration with the Clinical Bioinformatics and the Genetics Service of the Institut Curie.

### Contacts

For any question, bug or suggestion, please use the issues system or contact the bioinformatics core facility.

