# Aliases

<h1>Genome aliases</h1>

[TOC]

TLDR; **The genome alias system in refgenie allows users to refer to assets with arbitrary strings managed with `refgenie alias` command.**

## Motivation

Many systems rely on human-readable identifiers of genomes, such as "hg38". However, two users may refer to different things with the same identifier, such as the many slight variations of the *hg38* genome assembly. Such identifier mismatches lead to compatibility issues that incur the wrath of bioinformaticians everywhere. A step toward solving this problem is to use unique identifiers that unambiguously identify a particular assembly, such as those provided by the NCBI Assembly database; however, this approach relies on a central authority, and therefore can not apply to custom genomes or custom assets. Besides, human-readable identifiers persist because there's something simple and satisfying about referring to a genome or piece of data with a simple string that makes some sense and is easy to remember, like *hg38*. 

## Solutions

### Sequence-derived identifiers

Refgenie’s approach extends the [refget](http://samtools.github.io/hts-specs/refget.html) algorithm by GA4GH, introduced in 2019 to *collections of annotated sequences*. This means that the unique sequence-derived genome identifier calculated by refgenie captures not only sequence content, but also related metadata like sequence names and length. So, instead of referring to human genome as, e.g. "hg38" refgenie unambiguously identifies it as `58de7f33a36ccd9d6e3b1b3afe6b9f37cd5b2867bbfb929a`. 

#### Genome namespace initialization

The genome digest is calculated based on a FASTA file once the genome namespace is first created. This can happen when the `fasta` asset is pulled or built.

To start, initialize an empty refgenie configuration file from the shell and subscribe to the desired asset server:



```bash
export REFGENIE=$(pwd)/refgenie.yaml
refgenie init -c $REFGENIE -s http://rg.databio.org
```

    Initialized genome configuration file: /Users/mstolarczyk/code/refgenie/docs_jupyter/refgenie.yaml
    Created directories:
     - /Users/mstolarczyk/code/refgenie/docs_jupyter/data
     - /Users/mstolarczyk/code/refgenie/docs_jupyter/alias


Now, let's pull a `fasta` asset, which is one way to initialize a genome:


```bash
refgenie pull rCRSd/fasta --force
```

    Compatible refgenieserver instances: ['http://rg.databio.org']
    No local digest for genome alias: rCRSd
    Setting 'rCRSd' identity with server: http://rg.databio.org/v3/genomes/genome_digest/rCRSd
    Determined server digest for local genome alias (rCRSd): 94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4
    Set genome alias (94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4: rCRSd)
    Created alias directories: 
     - /Users/mstolarczyk/code/refgenie/docs_jupyter/alias/rCRSd
    Downloading URL: http://rg.databio.org/v3/assets/archive/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta
    [2K94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta:default  [35m100.0%[0m •  •  • …[0m •  • ? • … …
    [?25hDownload complete: /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/fasta__default.tgz
    Extracting asset tarball: /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/fasta__default.tgz
    Default tag for '94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta' set to: default
    Initializing genome: rCRSd
    Loaded AnnotatedSequenceDigestList (1 sequences)
    Set genome alias (94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4: rCRSd)
    Created alias directories: 
     - /Users/mstolarczyk/code/refgenie/docs_jupyter/alias/rCRSd/fasta/default


Following the `refgenie pull` command logs we notice that multiple steps happened:
1. refgenie used the human-readable genome name from the `refgenie pull` call (`rCRSd`) to query the server for any digest associated with it
2. refgenie set the digest it got back from the server as the genome identifier and set the human-readable genome name as an alias
3. refgenie used the genome idenfitier (not the user-specified name) to query the server for the `fasta` asset

From now on, the unique sequence-derived genome identifier will be used to query asset servers


```bash
refgenie pull rCRSd/bowtie2_index --force
```

    Compatible refgenieserver instances: ['http://rg.databio.org']
    Downloading URL: http://rg.databio.org/v3/assets/archive/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index
    [2K94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index:defau…  [35m100.…[0m   • •
    [?25hDownload complete: /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/bowtie2_index__default.tgz
    Extracting asset tarball: /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/bowtie2_index__default.tgz
    Default tag for '94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index' set to: default
    Created alias directories: 
     - /Users/mstolarczyk/code/refgenie/docs_jupyter/alias/rCRSd/bowtie2_index/default


### Genome aliases

To make the user's life easier genome aliases system in refgenie allows to set arbitrary genome aliases that can be then used to refer to a genome. Users can interact with genome aliases using `refgenie alias` command:




```bash
refgenie alias --help
```

    usage: refgenie alias [-h] {remove,set,get} ...
    
    Interact with aliases.
    
    positional arguments:
      {remove,set,get}
        remove          Remove aliases.
        set             Set aliases.
        get             Get aliases.
    
    optional arguments:
      -h, --help        show this help message and exit


#### Set aliases

To set an alias "mito" for genome identified by digest `94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4` one needs to issue the command below:


```bash
refgenie alias set --aliases mito --digest 94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4
```

    Set genome alias (94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4: mito)
    Created alias directories: 
     - /Users/mstolarczyk/code/refgenie/docs_jupyter/alias/mito


#### Get aliases

To get see the entire aliases collection managed by refgenie one needs to issue the command below:



```bash
refgenie alias get
```

    [3m                          Genome aliases                          [0m
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┓
    ┃[1m [0m[1mgenome                                          [0m[1m [0m┃[1m [0m[1malias      [0m[1m [0m┃
    ┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━┩
    │ 94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4 │ rCRSd, mito │
    └──────────────────────────────────────────────────┴─────────────┘


### `alias` and `data` directories

Refgenie stores asset data in two directories: `alias` and `data`. The `data` directory consists of the actual asset files, which are built or pulled from asset servers. The files in this directory are named using the digests, which helps refgenie to unambigously identify genomes. The `alias` holds symbolic links to asset data in `data` directory. **This way users do not need to be aware of the digest-named files at all and there is no waste of disk space due to symbolic links**. 

Here's a general view of the contents of both directories:


```bash
tree alias -d 
tree data -d 
```

    [01;34malias[00m
    ├── [01;34mmito[00m
    │   ├── [01;34mbowtie2_index[00m
    │   │   └── [01;34mdefault[00m
    │   └── [01;34mfasta[00m
    │       └── [01;34mdefault[00m
    └── [01;34mrCRSd[00m
        ├── [01;34mbowtie2_index[00m
        │   └── [01;34mdefault[00m
        └── [01;34mfasta[00m
            └── [01;34mdefault[00m
    
    10 directories
    [01;34mdata[00m
    └── [01;34m94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4[00m
        ├── [01;34mbowtie2_index[00m
        │   └── [01;34mdefault[00m
        └── [01;34mfasta[00m
            └── [01;34mdefault[00m
    
    5 directories


As you can see, the `alias` directory holds *both* of the defined aliases. Let's take a closer look at one of them


```bash
tree alias/rCRSd/fasta
```

    [01;34malias/rCRSd/fasta[00m
    └── [01;34mdefault[00m
        ├── [01;36mrCRSd.chrom.sizes[00m -> ../../../../data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.chrom.sizes
        ├── [01;36mrCRSd.fa[00m -> ../../../../data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.fa
        └── [01;36mrCRSd.fa.fai[00m -> ../../../../data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.fa.fai
    
    1 directory, 3 files


This explicitly shows that the files inside `alias/rCRSd/fasta/default` are in fact symbolic links that point to the actual asset files in `data` directory.
