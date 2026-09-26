# Tutorial


I assume you've already installed refgenie. In this tutorial I'll show you a few ways to use refgenie from the command line (commands that start with a `!`), and also some Python commands.

To start, initialize an empty refgenie configuration file from the shell and subscribe to the desired asset server:


```python
!refgenie init -c refgenie.yaml -s http://rg.databio.org
```

    Initialized genome configuration file: /Users/mstolarczyk/code/refgenie/docs_jupyter/refgenie.yaml
    Created directories:
     - /Users/mstolarczyk/code/refgenie/docs_jupyter/data
     - /Users/mstolarczyk/code/refgenie/docs_jupyter/alias


Here's what it looks like:


```python
!cat refgenie.yaml
```

    config_version: 0.4
    genome_folder: /Users/mstolarczyk/code/refgenie/docs_jupyter
    genome_servers: 
     - http://rg.databio.org
    genomes: null



```python
!refgenie listr -c refgenie.yaml
```

    [3m                             Remote refgenie assets                             [0m
    [3m                       Server URL: http://rg.databio.org                        [0m
    ┏━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    ┃[1m [0m[1mgenome          [0m[1m [0m┃[1m [0m[1massets                                                   [0m[1m [0m┃
    ┡━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
    │ rCRSd            │ fasta, bowtie2_index, bwa_index, hisat2_index,            │
    │                  │ star_index, bismark_bt2_index                             │
    │ hg18_cdna        │ fasta, kallisto_index                                     │
    │ hs38d1           │ fasta, suffixerator_index, bowtie2_index, bwa_index,      │
    │                  │ tallymer_index, hisat2_index, star_index,                 │
    │                  │ bismark_bt2_index                                         │
    │ hg38_cdna        │ fasta, kallisto_index, salmon_index                       │
    │ human_repeats    │ fasta, suffixerator_index, bowtie2_index, bwa_index,      │
    │                  │ tallymer_index, hisat2_index, star_index,                 │
    │                  │ bismark_bt2_index                                         │
    │ rn6_cdna         │ fasta, kallisto_index, salmon_index                       │
    │ mm10_cdna        │ fasta, kallisto_index, salmon_index                       │
    │ hg38_chr22       │ fasta, suffixerator_index, bowtie2_index, bwa_index,      │
    │                  │ tallymer_index, hisat2_index, star_index,                 │
    │                  │ bismark_bt2_index                                         │
    │ hg38             │ fasta, gencode_gtf, ensembl_gtf, refgene_anno,            │
    │                  │ fasta_txome, ensembl_rb, feat_annotation,                 │
    │                  │ suffixerator_index, cellranger_reference, bowtie2_index,  │
    │                  │ bwa_index, tallymer_index, hisat2_index, star_index,      │
    │                  │ bismark_bt2_index, salmon_partial_sa_index                │
    │ hg19_cdna        │ fasta, kallisto_index, salmon_index                       │
    │ human_rDNA       │ fasta, suffixerator_index, bowtie2_index, bwa_index,      │
    │                  │ tallymer_index, hisat2_index, star_index,                 │
    │                  │ bismark_bt2_index                                         │
    │ human_alu        │ fasta, suffixerator_index, bowtie2_index, bwa_index,      │
    │                  │ tallymer_index, hisat2_index, bismark_bt2_index           │
    │ human_alphasat   │ fasta, suffixerator_index, bowtie2_index, bwa_index,      │
    │                  │ tallymer_index, hisat2_index, star_index,                 │
    │                  │ bismark_bt2_index                                         │
    │ mouse_chrM2x     │ fasta, suffixerator_index, bowtie2_index, bwa_index,      │
    │                  │ tallymer_index, hisat2_index, star_index,                 │
    │                  │ bismark_bt2_index                                         │
    │ t7               │ fasta, bowtie2_index                                      │
    │ mm10             │ fasta, gencode_gtf, ensembl_gtf, refgene_anno,            │
    │                  │ fasta_txome, ensembl_rb, feat_annotation,                 │
    │                  │ suffixerator_index, cellranger_reference, bwa_index,      │
    │                  │ bowtie2_index, hisat2_index, tallymer_index, star_index,  │
    │                  │ bismark_bt2_index, salmon_partial_sa_index                │
    │ dm6              │ fasta, gencode_gtf, ensembl_gtf, refgene_anno,            │
    │                  │ bowtie2_index                                             │
    │ hg18             │ fasta, gencode_gtf, fasta_txome, suffixerator_index,      │
    │                  │ cellranger_reference, bwa_index, bowtie2_index,           │
    │                  │ tallymer_index, hisat2_index, star_index,                 │
    │                  │ bismark_bt2_index                                         │
    │ hg19             │ fasta, gencode_gtf, ensembl_gtf, refgene_anno,            │
    │                  │ fasta_txome, ensembl_rb, feat_annotation,                 │
    │                  │ suffixerator_index, cellranger_reference, bwa_index,      │
    │                  │ bowtie2_index, tallymer_index, hisat2_index, star_index,  │
    │                  │ salmon_partial_sa_index, bismark_bt2_index                │
    │ rn6              │ fasta, ensembl_gtf, refgene_anno, fasta_txome,            │
    │                  │ suffixerator_index, bwa_index, bowtie2_index,             │
    │                  │ tallymer_index, hisat2_index, star_index,                 │
    │                  │ bismark_bt2_index, salmon_partial_sa_index                │
    │ hg38_noalt_decoy │ fasta, suffixerator_index, bowtie2_index, bwa_index,      │
    │                  │ tallymer_index, hisat2_index, bismark_bt2_index           │
    │ mm10_primary     │ fasta, bowtie2_index, bwa_index                           │
    │ hg38_primary     │ fasta, bowtie2_index, bwa_index                           │
    │ hg38_mm10        │ fasta, bwa_index                                          │
    └──────────────────┴───────────────────────────────────────────────────────────┘
    [2;3m             use refgenie listr -g <genome> for more detailed view              [0m


Now let's enter python and do some stuff.


```python
import refgenconf
rgc = refgenconf.RefGenConf(filepath="refgenie.yaml")
```

Use `pull` to download one of the assets:


```python
rgc.pull("mouse_chrM2x", "fasta", "default")
```


    Output()





    (['43f14ba8beed34d52edb244e26f193df6edbb467bd55d37a', 'fasta', 'default'],
     {'asset_path': 'fasta',
      'asset_digest': '8dfe402f7d29d5b036dd8937119e4404',
      'archive_digest': 'bfb7877ee114c61a17a50bd471de47a2',
      'asset_size': '39.4KB',
      'archive_size': '9.1KB',
      'seek_keys': {'fasta': '43f14ba8beed34d52edb244e26f193df6edbb467bd55d37a.fa',
       'fai': '43f14ba8beed34d52edb244e26f193df6edbb467bd55d37a.fa.fai',
       'chrom_sizes': '43f14ba8beed34d52edb244e26f193df6edbb467bd55d37a.chrom.sizes'},
      'asset_parents': [],
      'asset_children': ['43f14ba8beed34d52edb244e26f193df6edbb467bd55d37a/suffixerator_index:default',
       '43f14ba8beed34d52edb244e26f193df6edbb467bd55d37a/bowtie2_index:default',
       '43f14ba8beed34d52edb244e26f193df6edbb467bd55d37a/bwa_index:default',
       '43f14ba8beed34d52edb244e26f193df6edbb467bd55d37a/tallymer_index:default',
       '43f14ba8beed34d52edb244e26f193df6edbb467bd55d37a/hisat2_index:default',
       '43f14ba8beed34d52edb244e26f193df6edbb467bd55d37a/star_index:default',
       '43f14ba8beed34d52edb244e26f193df6edbb467bd55d37a/bismark_bt2_index:default']},
     'http://rg.databio.org')



Once it's downloaded, use `seek` to retrieve a path to it.


```python
rgc.seek("mouse_chrM2x", "fasta")
```




    '/Users/mstolarczyk/code/refgenie/docs_jupyter/alias/mouse_chrM2x/fasta/default/mouse_chrM2x.fa'



You can get the unique asset identifier with `id()`


```python
rgc.id("mouse_chrM2x", "fasta")
```




    '8dfe402f7d29d5b036dd8937119e4404'



## Building and pulling from the command line

Here, we can build a fasta asset instead of pulling one. Back to the shell, we'll grab the Revised Cambridge Reference Sequence (human mitochondrial genome, because it's small):


```python
!wget -O rCRSd.fa.gz http://big.databio.org/refgenie_raw/files.rCRSd.fasta.fasta
```

    --2021-03-09 12:22:40--  http://big.databio.org/refgenie_raw/files.rCRSd.fasta.fasta
    Resolving big.databio.org (big.databio.org)... 128.143.245.181, 128.143.245.182
    Connecting to big.databio.org (big.databio.org)|128.143.245.181|:80... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 8399 (8.2K) [application/octet-stream]
    Saving to: ‘rCRSd.fa.gz’
    
    rCRSd.fa.gz         100%[===================>]   8.20K  --.-KB/s    in 0.006s  
    
    2021-03-09 12:22:40 (1.35 MB/s) - ‘rCRSd.fa.gz’ saved [8399/8399]
    



```python
!refgenie build rCRSd/fasta -c refgenie.yaml  --files fasta=rCRSd.fa.gz -R
```

    Using 'default' as the default tag for 'rCRSd/fasta'
    Recipe validated successfully against a schema: /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/refgenie/schemas/recipe_schema.yaml
    Building 'rCRSd/fasta:default' using 'fasta' recipe
    Initializing genome: rCRSd
    Loaded AnnotatedSequenceDigestList (1 sequences)
    Set genome alias (94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4: rCRSd)
    Created alias directories: 
     - /Users/mstolarczyk/code/refgenie/docs_jupyter/alias/rCRSd
    Saving outputs to:
    - content: /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4
    - logs: /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/_refgenie_build
    ### Pipeline run code and environment:
    
    *              Command:  `/Library/Frameworks/Python.framework/Versions/3.6/bin/refgenie build rCRSd/fasta -c refgenie.yaml --files fasta=rCRSd.fa.gz -R`
    *         Compute host:  MichalsMBP
    *          Working dir:  /Users/mstolarczyk/code/refgenie/docs_jupyter
    *            Outfolder:  /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/_refgenie_build/
    *  Pipeline started at:   (03-09 12:22:41) elapsed: 0.0 _TIME_
    
    ### Version log:
    
    *       Python version:  3.6.5
    *          Pypiper dir:  `/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/pypiper`
    *      Pypiper version:  0.12.1
    *         Pipeline dir:  `/Library/Frameworks/Python.framework/Versions/3.6/bin`
    *     Pipeline version:  None
    
    ### Arguments passed to pipeline:
    
    * `asset_registry_paths`:  `['rCRSd/fasta']`
    *             `assets`:  `None`
    *            `command`:  `build`
    *        `config_file`:  `refgenie.yaml`
    *             `docker`:  `False`
    *              `files`:  `[['fasta=rCRSd.fa.gz']]`
    *             `genome`:  `None`
    *      `genome_config`:  `refgenie.yaml`
    * `genome_description`:  `None`
    *             `logdev`:  `False`
    *          `new_start`:  `False`
    *          `outfolder`:  `/Users/mstolarczyk/code/refgenie/docs_jupyter/data`
    *             `params`:  `None`
    *             `recipe`:  `None`
    *            `recover`:  `True`
    *       `requirements`:  `False`
    *             `silent`:  `False`
    *     `skip_read_lock`:  `False`
    *    `tag_description`:  `None`
    *          `verbosity`:  `None`
    *            `volumes`:  `None`
    
    ----------------------------------------
    
    Target to produce: `/Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/_refgenie_build/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4_fasta__default.flag`  
    
    > `cp rCRSd.fa.gz /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.fa.gz` (63575)
    <pre>
    psutil.ZombieProcess process still exists but it's a zombie (pid=63575)
    Warning: couldn't add memory use for process: 63575
    </pre>
    Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
      PID: 63575;	Command: cp;	Return code: 0;	Memory used: 0GB
    
    
    > `gzip -df /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.fa.gz` (63576)
    <pre>
    psutil.ZombieProcess process still exists but it's a zombie (pid=63576)
    Warning: couldn't add memory use for process: 63576
    </pre>
    Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
      PID: 63576;	Command: gzip;	Return code: 0;	Memory used: 0GB
    
    
    > `samtools faidx /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.fa` (63577)
    <pre>
    </pre>
    Command completed. Elapsed time: 0:00:00. Running peak memory: 0.001GB.  
      PID: 63577;	Command: samtools;	Return code: 0;	Memory used: 0.001GB
    
    
    > `cut -f 1,2 /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.fa.fai > /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.chrom.sizes` (63578)
    <pre>
    psutil.ZombieProcess process still exists but it's a zombie (pid=63578)
    Warning: couldn't add memory use for process: 63578
    </pre>
    Command completed. Elapsed time: 0:00:00. Running peak memory: 0.001GB.  
      PID: 63578;	Command: cut;	Return code: 0;	Memory used: 0GB
    
    
    > `touch /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/_refgenie_build/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4_fasta__default.flag` (63580)
    <pre>
    psutil.ZombieProcess process still exists but it's a zombie (pid=63580)
    Warning: couldn't add memory use for process: 63580
    </pre>
    Command completed. Elapsed time: 0:00:00. Running peak memory: 0.001GB.  
      PID: 63580;	Command: touch;	Return code: 0;	Memory used: 0GB
    
    Asset digest: 4eb430296bc02ed7e4006624f1d5ac53
    Default tag for '94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta' set to: default
    
    ### Pipeline completed. Epilogue
    *        Elapsed time (this run):  0:00:00
    *  Total elapsed time (all runs):  0:00:00
    *         Peak memory (this run):  0.0015 GB
    *        Pipeline completed time: 2021-03-09 12:22:41
    Finished building 'fasta' asset
    Created alias directories: 
     - /Users/mstolarczyk/code/refgenie/docs_jupyter/alias/rCRSd/fasta/default


The asset should be available for local use, let's call `refgenie list` to check it:


```python
!refgenie list -c refgenie.yaml --genome rCRSd
```

    [3m                        Local refgenie assets                         [0m
    [3m             Server subscriptions: http://rg.databio.org              [0m
    ┏━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━┓
    ┃[1m [0m[1mgenome   [0m[1m [0m┃[1m [0m[1masset ([0m[1;3mseek_keys[0m[1m)                         [0m[1m [0m┃[1m [0m[1mtags     [0m[1m [0m┃
    ┡━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━┩
    │ rCRSd     │ fasta ([3mfasta, fai, chrom_sizes[0m)            │ default   │
    └───────────┴────────────────────────────────────────────┴───────────┘


We can retrieve the path to this asset with:


```python
!refgenie seek rCRSd/fasta -c refgenie.yaml
```

    /Users/mstolarczyk/code/refgenie/docs_jupyter/alias/rCRSd/fasta/default/rCRSd.fa


Naturally, we can do the same thing from within Python:


```python
rgc = refgenconf.RefGenConf("refgenie.yaml")
rgc.seek("rCRSd", "fasta")
```




    '/Users/mstolarczyk/code/refgenie/docs_jupyter/alias/rCRSd/fasta/default/rCRSd.fa'



Now, if we have bowtie2-build in our `$PATH` we can build the `bowtie2_index` asset with no further requirements.

Let's check the requirements with `refgenie build --requirements`:



```python
!refgenie build rCRSd/bowtie2_index -c refgenie.yaml --requirements
```

    'bowtie2_index' recipe requirements: 
    - assets:
    	fasta (fasta asset for genome); default: fasta


Since I already have the fasta asset, that means I don't need anything else to build the bowtie2_index.


```python
!refgenie build rCRSd/bowtie2_index -c refgenie.yaml
```

    Using 'default' as the default tag for 'rCRSd/bowtie2_index'
    Recipe validated successfully against a schema: /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/refgenie/schemas/recipe_schema.yaml
    Building 'rCRSd/bowtie2_index:default' using 'bowtie2_index' recipe
    Saving outputs to:
    - content: /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4
    - logs: /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/default/_refgenie_build
    ### Pipeline run code and environment:
    
    *              Command:  `/Library/Frameworks/Python.framework/Versions/3.6/bin/refgenie build rCRSd/bowtie2_index -c refgenie.yaml`
    *         Compute host:  MichalsMBP
    *          Working dir:  /Users/mstolarczyk/code/refgenie/docs_jupyter
    *            Outfolder:  /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/default/_refgenie_build/
    *  Pipeline started at:   (03-09 12:22:45) elapsed: 0.0 _TIME_
    
    ### Version log:
    
    *       Python version:  3.6.5
    *          Pypiper dir:  `/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/pypiper`
    *      Pypiper version:  0.12.1
    *         Pipeline dir:  `/Library/Frameworks/Python.framework/Versions/3.6/bin`
    *     Pipeline version:  None
    
    ### Arguments passed to pipeline:
    
    * `asset_registry_paths`:  `['rCRSd/bowtie2_index']`
    *             `assets`:  `None`
    *            `command`:  `build`
    *        `config_file`:  `refgenie.yaml`
    *             `docker`:  `False`
    *              `files`:  `None`
    *             `genome`:  `None`
    *      `genome_config`:  `refgenie.yaml`
    * `genome_description`:  `None`
    *             `logdev`:  `False`
    *          `new_start`:  `False`
    *          `outfolder`:  `/Users/mstolarczyk/code/refgenie/docs_jupyter/data`
    *             `params`:  `None`
    *             `recipe`:  `None`
    *            `recover`:  `False`
    *       `requirements`:  `False`
    *             `silent`:  `False`
    *     `skip_read_lock`:  `False`
    *    `tag_description`:  `None`
    *          `verbosity`:  `None`
    *            `volumes`:  `None`
    
    ----------------------------------------
    
    Target to produce: `/Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/default/_refgenie_build/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4_bowtie2_index__default.flag`  
    
    > `bowtie2-build /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.fa /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4` (63609)
    <pre>
    Settings:
      Output files: "/Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.*.bt2"
      Line rate: 6 (line is 64 bytes)
      Lines per side: 1 (side is 64 bytes)
      Offset rate: 4 (one in 16)
      FTable chars: 10
      Strings: unpacked
      Max bucket size: default
      Max bucket size, sqrt multiplier: default
      Max bucket size, len divisor: 4
      Difference-cover sample period: 1024
      Endianness: little
      Actual local endianness: little
      Sanity checking: disabled
      Assertions: disabled
      Random seed: 0
      Sizeofs: void*:8, int:4, long:8, size_t:8
    Input files DNA, FASTA:
      /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.fa
    Building a SMALL index
    Reading reference sizes
      Time reading reference sizes: 00:00:00
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:00
    bmax according to bmaxDivN setting: 8284
    Using parameters --bmax 6213 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 6213 --dcv 1024
    Constructing suffix-array element generator
    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:00:00
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:00
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:00:00
      Sanity-checking and returning
    Building samples
    Reserving space for 12 sample suffixes
    Generating random suffixes
    QSorting 12 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 12 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 33136 (target: 6212)
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering Ebwt loop
    Getting block 1 of 1
      No samples; assembling all-inclusive block
      Sorting block of length 33136 for bucket 1
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 33137 for bucket 1
    Exited Ebwt loop
    fchr[A]: 0
    fchr[C]: 10248
    fchr[G]: 20610
    fchr[T]: 24948
    fchr[$]: 33136
    Exiting Ebwt::buildToDisk()
    Returning from initFromVector
    Wrote 4205567 bytes to primary EBWT file: /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.1.bt2
    Wrote 8292 bytes to secondary EBWT file: /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.2.bt2
    Re-opening _in1 and _in2 as input streams
    Returning from Ebwt constructor
    Headers:
        len: 33136
        bwtLen: 33137
        sz: 8284
        bwtSz: 8285
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 20
        eftabSz: 80
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 2072
        offsSz: 8288
        lineSz: 64
        sideSz: 64
        sideBwtSz: 48
        sideBwtLen: 192
        numSides: 173
        numLines: 173
        ebwtTotLen: 11072
        ebwtTotSz: 11072
        color: 0
        reverse: 0
    Total time for call to driver() for forward index: 00:00:00
    Reading reference sizes
      Time reading reference sizes: 00:00:00
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:00
      Time to reverse reference sequence: 00:00:00
    bmax according to bmaxDivN setting: 8284
    Using parameters --bmax 6213 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 6213 --dcv 1024
    Constructing suffix-array element generator
    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:00:00
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:00
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:00:00
      Sanity-checking and returning
    Building samples
    Reserving space for 12 sample suffixes
    Generating random suffixes
    QSorting 12 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 12 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 33136 (target: 6212)
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering Ebwt loop
    Getting block 1 of 1
      No samples; assembling all-inclusive block
      Sorting block of length 33136 for bucket 1
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 33137 for bucket 1
    Exited Ebwt loop
    fchr[A]: 0
    fchr[C]: 10248
    fchr[G]: 20610
    fchr[T]: 24948
    fchr[$]: 33136
    Exiting Ebwt::buildToDisk()
    Returning from initFromVector
    Wrote 4205567 bytes to primary EBWT file: /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.rev.1.bt2
    Wrote 8292 bytes to secondary EBWT file: /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.rev.2.bt2
    Re-opening _in1 and _in2 as input streams
    Returning from Ebwt constructor
    Headers:
        len: 33136
        bwtLen: 33137
        sz: 8284
        bwtSz: 8285
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 20
        eftabSz: 80
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 2072
        offsSz: 8288
        lineSz: 64
        sideSz: 64
        sideBwtSz: 48
        sideBwtLen: 192
        numSides: 173
        numLines: 173
        ebwtTotLen: 11072
        ebwtTotSz: 11072
        color: 0
        reverse: 1
    Total time for backward call to driver() for mirror index: 00:00:00
    </pre>


    Command completed. Elapsed time: 0:00:00. Running peak memory: 0.003GB.  
      PID: 63609;	Command: bowtie2-build;	Return code: 0;	Memory used: 0.003GB
    
    
    > `touch /Users/mstolarczyk/code/refgenie/docs_jupyter/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/default/_refgenie_build/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4_bowtie2_index__default.flag` (63611)
    <pre>
    psutil.ZombieProcess process still exists but it's a zombie (pid=63611)
    Warning: couldn't add memory use for process: 63611
    </pre>
    Command completed. Elapsed time: 0:00:00. Running peak memory: 0.003GB.  
      PID: 63611;	Command: touch;	Return code: 0;	Memory used: 0GB
    
    Asset digest: 1262e30d4a87db9365d501de8559b3b4
    Default tag for '94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index' set to: default
    
    ### Pipeline completed. Epilogue
    *        Elapsed time (this run):  0:00:01
    *  Total elapsed time (all runs):  0:00:00
    *         Peak memory (this run):  0.0028 GB
    *        Pipeline completed time: 2021-03-09 12:22:46
    Finished building 'bowtie2_index' asset
    Created alias directories: 
     - /Users/mstolarczyk/code/refgenie/docs_jupyter/alias/rCRSd/bowtie2_index/default


We can see a list of available recipes like this:


```python
!refgenie list -c refgenie.yaml --recipes
```

    bismark_bt1_index, bismark_bt2_index, blacklist, bowtie2_index, bwa_index, cellranger_reference, dbnsfp, dbsnp, ensembl_gtf, ensembl_rb, epilog_index, fasta, fasta_txome, feat_annotation, gencode_gtf, hisat2_index, kallisto_index, refgene_anno, salmon_index, salmon_partial_sa_index, salmon_sa_index, star_index, suffixerator_index, tallymer_index, tgMap


We can get the unique digest for any asset with `refgenie id`:


```python
!refgenie id rCRSd/fasta -c refgenie.yaml
```

    4eb430296bc02ed7e4006624f1d5ac53


## Versions


```python
from platform import python_version 
python_version()
```




    '3.6.5'




```python
!refgenie --version
```

    refgenie 0.10.0-dev | refgenconf 0.10.0-dev

