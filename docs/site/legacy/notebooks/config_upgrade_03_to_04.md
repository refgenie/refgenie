# Configuration file upgrade demonstration


In the following tutorial we will present the process of upgrading the refgenie configuration file and asset files from version **0.3** to version **0.4**.

First, let's install the refgenie and refgenconf Python packages that support version 0.3 of refgenie configuration file

## Working environment setup

Let's install the legacy refgenconf and refgenie Python packages


```bash
pip install refgenconf==0.9.3
pip install refgenie==0.9.3
```

    Collecting refgenconf==0.9.3
      Using cached https://files.pythonhosted.org/packages/52/c3/6aed361205272e30cd3570ca1c33feae6ad977ad32ddff8e509752046272/refgenconf-0.9.3-py3-none-any.whl
    Requirement already satisfied: requests in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from refgenconf==0.9.3) (2.21.0)
    Requirement already satisfied: attmap>=0.12.5 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenconf==0.9.3) (0.12.12.dev0)
    Requirement already satisfied: tqdm>=4.38.0 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenconf==0.9.3) (4.47.0)
    Requirement already satisfied: pyyaml in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from refgenconf==0.9.3) (5.1)
    Requirement already satisfied: yacman>=0.6.9 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenconf==0.9.3) (0.7.0)
    Requirement already satisfied: pyfaidx in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenconf==0.9.3) (0.5.9.1)
    Requirement already satisfied: urllib3<1.25,>=1.21.1 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from requests->refgenconf==0.9.3) (1.24.1)
    Requirement already satisfied: idna<2.9,>=2.5 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from requests->refgenconf==0.9.3) (2.8)
    Requirement already satisfied: certifi>=2017.4.17 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from requests->refgenconf==0.9.3) (2019.3.9)
    Requirement already satisfied: chardet<3.1.0,>=3.0.2 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from requests->refgenconf==0.9.3) (3.0.4)
    Requirement already satisfied: ubiquerg>=0.2.1 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from attmap>=0.12.5->refgenconf==0.9.3) (0.6.1)
    Requirement already satisfied: oyaml in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from yacman>=0.6.9->refgenconf==0.9.3) (0.9)
    Requirement already satisfied: setuptools>=0.7 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from pyfaidx->refgenconf==0.9.3) (41.0.1)
    Requirement already satisfied: six in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from pyfaidx->refgenconf==0.9.3) (1.12.0)
    Installing collected packages: refgenconf
      Found existing installation: refgenconf 0.10.0.dev0
        Uninstalling refgenconf-0.10.0.dev0:
          Successfully uninstalled refgenconf-0.10.0.dev0
    Successfully installed refgenconf-0.9.3
    [33mWARNING: You are using pip version 19.2.3, however version 20.2.3 is available.
    You should consider upgrading via the 'pip install --upgrade pip' command.[0m
    Collecting refgenie==0.9.3
      Using cached https://files.pythonhosted.org/packages/af/52/c1e1bc63b3543f591ebdf44caccfaab3c730708256d926b9f4b1c34d1865/refgenie-0.9.3-py3-none-any.whl
    Requirement already satisfied: pyfaidx>=0.5.5.2 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenie==0.9.3) (0.5.9.1)
    Requirement already satisfied: refgenconf>=0.9.1 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenie==0.9.3) (0.9.3)
    Requirement already satisfied: logmuse>=0.2.6 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenie==0.9.3) (0.2.6)
    Requirement already satisfied: piper>=0.12.1 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenie==0.9.3) (0.12.1)
    Requirement already satisfied: six in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from pyfaidx>=0.5.5.2->refgenie==0.9.3) (1.12.0)
    Requirement already satisfied: setuptools>=0.7 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from pyfaidx>=0.5.5.2->refgenie==0.9.3) (41.0.1)
    Requirement already satisfied: attmap>=0.12.5 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenconf>=0.9.1->refgenie==0.9.3) (0.12.12.dev0)
    Requirement already satisfied: requests in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from refgenconf>=0.9.1->refgenie==0.9.3) (2.21.0)
    Requirement already satisfied: tqdm>=4.38.0 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenconf>=0.9.1->refgenie==0.9.3) (4.47.0)
    Requirement already satisfied: pyyaml in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from refgenconf>=0.9.1->refgenie==0.9.3) (5.1)
    Requirement already satisfied: yacman>=0.6.9 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenconf>=0.9.1->refgenie==0.9.3) (0.7.0)
    Requirement already satisfied: ubiquerg>=0.4.5 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from piper>=0.12.1->refgenie==0.9.3) (0.6.1)
    Requirement already satisfied: psutil in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from piper>=0.12.1->refgenie==0.9.3) (5.6.1)
    Requirement already satisfied: pandas in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from piper>=0.12.1->refgenie==0.9.3) (1.0.3)
    Requirement already satisfied: urllib3<1.25,>=1.21.1 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from requests->refgenconf>=0.9.1->refgenie==0.9.3) (1.24.1)
    Requirement already satisfied: chardet<3.1.0,>=3.0.2 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from requests->refgenconf>=0.9.1->refgenie==0.9.3) (3.0.4)
    Requirement already satisfied: certifi>=2017.4.17 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from requests->refgenconf>=0.9.1->refgenie==0.9.3) (2019.3.9)
    Requirement already satisfied: idna<2.9,>=2.5 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from requests->refgenconf>=0.9.1->refgenie==0.9.3) (2.8)
    Requirement already satisfied: oyaml in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from yacman>=0.6.9->refgenconf>=0.9.1->refgenie==0.9.3) (0.9)
    Requirement already satisfied: pytz>=2017.2 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from pandas->piper>=0.12.1->refgenie==0.9.3) (2018.9)
    Requirement already satisfied: python-dateutil>=2.6.1 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from pandas->piper>=0.12.1->refgenie==0.9.3) (2.8.0)
    Requirement already satisfied: numpy>=1.13.3 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from pandas->piper>=0.12.1->refgenie==0.9.3) (1.17.3)
    Installing collected packages: refgenie
      Found existing installation: refgenie 0.10.0.dev0
        Uninstalling refgenie-0.10.0.dev0:
          Successfully uninstalled refgenie-0.10.0.dev0
    Successfully installed refgenie-0.9.3
    [33mWARNING: You are using pip version 19.2.3, however version 20.2.3 is available.
    You should consider upgrading via the 'pip install --upgrade pip' command.[0m


Now let's set up a directory that we will use for the config file and refgenie assets


```bash
export WORKDIR=~/Desktop/testing/refgenie/upgrade_test
rm -r $WORKDIR # remove first just to make sure the directory does not exist
mkdir -p $WORKDIR
cd $WORKDIR
```

Let's set `$REFGENIE` environment variable to point refgenie to the configuration file location and initialize it


```bash
export REFGENIE=$WORKDIR/g.yml
refgenie init -c $REFGENIE -s http://rg.databio.org:82/
```

    Initialized genome configuration file: /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/g.yml


Note that we subscribe to a test instance of refgenieserver, that supports both the old and new refgenie clients. This is because it exposes different API versions, that these clients use: `v2` (refgenie v0.9.3) and `v3` (refgenie v0.10.0-dev)

## Pull/build test assets

Next, let's retrieve couple of assets. As mentioned above, `v2` API is used to retrieve the asset.


```bash
refgenie pull rCRSd/fasta human_repeats/fasta rCRSd/bowtie2_index human_repeats/bwa_index
```

    Downloading URL: http://rg.databio.org:82/v2/asset/rCRSd/fasta/archive
    Download complete: /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/rCRSd/fasta__default.tgz
    Extracting asset tarball and saving to: /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/rCRSd/fasta/default
    Default tag for 'rCRSd/fasta' set to: default
    Downloading URL: http://rg.databio.org:82/v2/asset/human_repeats/fasta/archive
    Download complete: /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_repeats/fasta__default.tgz
    Extracting asset tarball and saving to: /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_repeats/fasta/default
    Default tag for 'human_repeats/fasta' set to: default
    Downloading URL: http://rg.databio.org:82/v2/asset/rCRSd/bowtie2_index/archive
    Download complete: /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/rCRSd/bowtie2_index__default.tgz
    Extracting asset tarball and saving to: /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/rCRSd/bowtie2_index/default
    Default tag for 'rCRSd/bowtie2_index' set to: default
    Downloading URL: http://rg.databio.org:82/v2/asset/human_repeats/bwa_index/archive
    Download complete: /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_repeats/bwa_index__default.tgz
    Extracting asset tarball and saving to: /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_repeats/bwa_index/default
    Default tag for 'human_repeats/bwa_index' set to: default


Now, let's download a small FASTA file and build a fasta asset for an arbitrary genome, which is not available at `http://rg.databio.org:82/`


```bash
wget -O human_alu.fa.gz http://big.databio.org/refgenie_raw/files.human_alu.fasta.fasta

```

    --2020-10-12 17:39:25--  http://big.databio.org/refgenie_raw/files.human_alu.fasta.fasta
    Resolving big.databio.org (big.databio.org)... 128.143.245.182, 128.143.245.181
    Connecting to big.databio.org (big.databio.org)|128.143.245.182|:80... connected.
    HTTP request sent, awaiting response... 200 OK
    Length: 501 [application/octet-stream]
    Saving to: ‘human_alu.fa.gz’
    
    human_alu.fa.gz     100%[===================>]     501  --.-KB/s    in 0s      
    
    2020-10-12 17:39:25 (1.19 MB/s) - ‘human_alu.fa.gz’ saved [501/501]
    



```bash
refgenie build human_alu/fasta --files fasta=human_alu.fa.gz
```

    Using 'default' as the default tag for 'human_alu/fasta'
    Building 'human_alu/fasta:default' using 'fasta' recipe
    Saving outputs to:
    - content: /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_alu
    - logs: /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_alu/fasta/default/_refgenie_build
    ### Pipeline run code and environment:
    
    *              Command:  `/Library/Frameworks/Python.framework/Versions/3.6/bin/refgenie build human_alu/fasta --files fasta=human_alu.fa.gz`
    *         Compute host:  MichalsMBP
    *          Working dir:  /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test
    *            Outfolder:  /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_alu/fasta/default/_refgenie_build/
    *  Pipeline started at:   (10-12 17:39:27) elapsed: 0.0 _TIME_
    
    ### Version log:
    
    *       Python version:  3.6.5
    *          Pypiper dir:  `/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/pypiper`
    *      Pypiper version:  0.12.1
    *         Pipeline dir:  `/Library/Frameworks/Python.framework/Versions/3.6/bin`
    *     Pipeline version:  None
    
    ### Arguments passed to pipeline:
    
    * `asset_registry_paths`:  `['human_alu/fasta']`
    *             `assets`:  `None`
    *            `command`:  `build`
    *        `config_file`:  `/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/refgenie/refgenie.yaml`
    *             `docker`:  `False`
    *              `files`:  `[['fasta=human_alu.fa.gz']]`
    *             `genome`:  `None`
    *      `genome_config`:  `None`
    * `genome_description`:  `None`
    *             `logdev`:  `False`
    *          `new_start`:  `False`
    *          `outfolder`:  `/Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test`
    *             `params`:  `None`
    *             `recipe`:  `None`
    *            `recover`:  `False`
    *       `requirements`:  `False`
    *             `silent`:  `False`
    *    `tag_description`:  `None`
    *          `verbosity`:  `None`
    *            `volumes`:  `None`
    
    ----------------------------------------
    
    Target to produce: `/Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_alu/fasta/default/_refgenie_build/human_alu_fasta__default.flag`  
    
    > `cp human_alu.fa.gz /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_alu/fasta/default/human_alu.fa.gz` (70063)
    <pre>
    psutil.ZombieProcess process still exists but it's a zombie (pid=70063)
    Warning: couldn't add memory use for process: 70063
    </pre>
    Command completed. Elapsed time: 0:00:00. Running peak memory: 0.001GB.  
      PID: 70063;	Command: cp;	Return code: 0;	Memory used: 0.001GB
    
    
    > `gzip -df /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_alu/fasta/default/human_alu.fa.gz` (70064)
    <pre>
    psutil.ZombieProcess process still exists but it's a zombie (pid=70064)
    Warning: couldn't add memory use for process: 70064
    </pre>
    Command completed. Elapsed time: 0:00:00. Running peak memory: 0.001GB.  
      PID: 70064;	Command: gzip;	Return code: 0;	Memory used: 0.0GB
    
    
    > `samtools faidx /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_alu/fasta/default/human_alu.fa` (70065)
    <pre>
    psutil.ZombieProcess process still exists but it's a zombie (pid=70065)
    Warning: couldn't add memory use for process: 70065
    </pre>
    Command completed. Elapsed time: 0:00:00. Running peak memory: 0.001GB.  
      PID: 70065;	Command: samtools;	Return code: 0;	Memory used: 0.0GB
    
    
    > `cut -f 1,2 /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_alu/fasta/default/human_alu.fa.fai > /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_alu/fasta/default/human_alu.chrom.sizes` (70066)
    <pre>
    psutil.ZombieProcess process still exists but it's a zombie (pid=70066)
    Warning: couldn't add memory use for process: 70066
    </pre>
    Command completed. Elapsed time: 0:00:00. Running peak memory: 0.001GB.  
      PID: 70066;	Command: cut;	Return code: 0;	Memory used: 0GB
    
    
    > `touch /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_alu/fasta/default/_refgenie_build/human_alu_fasta__default.flag` (70068)
    <pre>
    psutil.ZombieProcess process still exists but it's a zombie (pid=70068)
    Warning: couldn't add memory use for process: 70068
    </pre>
    Command completed. Elapsed time: 0:00:00. Running peak memory: 0.001GB.  
      PID: 70068;	Command: touch;	Return code: 0;	Memory used: 0GB
    
    
    > `cd /Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test/human_alu/fasta/default; find . -type f -not -path './_refgenie_build*' -exec md5sum {} \; | sort -k 2 | awk '{print $1}' | md5sum`
    Asset digest: 9e8fa06e6125f89be4fb974879cb91a6
    Default tag for 'human_alu/fasta' set to: default
    
    ### Pipeline completed. Epilogue
    *        Elapsed time (this run):  0:00:00
    *  Total elapsed time (all runs):  0:00:00
    *         Peak memory (this run):  0.001 GB
    *        Pipeline completed time: 2020-10-12 17:39:27
    Computing initial genome digest...
    Initializing genome...
    Finished building 'fasta' asset


Let's inspect the asset inventory


```bash
refgenie list
```

    Server subscriptions: http://rg.databio.org:82
    Local genomes: human_alu, human_repeats, rCRSd
    Local recipes: bismark_bt1_index, bismark_bt2_index, blacklist, bowtie2_index, bwa_index, cellranger_reference, dbnsfp, dbsnp, ensembl_gtf, ensembl_rb, epilog_index, fasta, fasta_txome, feat_annotation, gencode_gtf, hisat2_index, kallisto_index, refgene_anno, salmon_index, salmon_partial_sa_index, salmon_sa_index, star_index, suffixerator_index, tallymer_index
    Local assets:
               human_alu/   fasta.chrom_sizes:default, fasta.fai:default, fasta:default
           human_repeats/   bwa_index:default, fasta.chrom_sizes:default, fasta.fai:default, fasta:default
                   rCRSd/   bowtie2_index:default, fasta.chrom_sizes:default, fasta.fai:default, fasta:default


As you can see, assets for all three genomes are available.

## Upgrade refgenie software

Now, let's upgrade to refgenie==0.10.0-dev, which introduces the concept of sequence-derived genome identifiers to uniqly identify genomes.


```bash
pip install git+https://github.com/refgenie/refgenconf.git@dev_config_upgrade
pip install git+https://github.com/refgenie/refgenie.git@dev_config_upgrade
```

    Collecting git+https://github.com/refgenie/refgenconf.git@dev_config_upgrade
      Cloning https://github.com/refgenie/refgenconf.git (to revision dev_config_upgrade) to /private/var/folders/3f/0wj7rs2144l9zsgxd3jn5nxc0000gn/T/pip-req-build-kxmw8i6n
      Running command git clone -q https://github.com/refgenie/refgenconf.git /private/var/folders/3f/0wj7rs2144l9zsgxd3jn5nxc0000gn/T/pip-req-build-kxmw8i6n
      Running command git checkout -b dev_config_upgrade --track origin/dev_config_upgrade
      Switched to a new branch 'dev_config_upgrade'
      Branch 'dev_config_upgrade' set up to track remote branch 'dev_config_upgrade' from 'origin'.
    Requirement already satisfied: attmap>=0.12.5 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenconf==0.10.0.dev0) (0.12.12.dev0)
    Requirement already satisfied: pyyaml in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from refgenconf==0.10.0.dev0) (5.1)
    Requirement already satisfied: requests in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from refgenconf==0.10.0.dev0) (2.21.0)
    Requirement already satisfied: yacman>=0.7.0 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenconf==0.10.0.dev0) (0.7.0)
    Requirement already satisfied: future in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenconf==0.10.0.dev0) (0.17.1)
    Requirement already satisfied: jsonschema in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from refgenconf==0.10.0.dev0) (3.0.1)
    Requirement already satisfied: rich in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenconf==0.10.0.dev0) (3.3.0)
    Requirement already satisfied: ubiquerg>=0.2.1 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from attmap>=0.12.5->refgenconf==0.10.0.dev0) (0.6.1)
    Requirement already satisfied: chardet<3.1.0,>=3.0.2 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from requests->refgenconf==0.10.0.dev0) (3.0.4)
    Requirement already satisfied: certifi>=2017.4.17 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from requests->refgenconf==0.10.0.dev0) (2019.3.9)
    Requirement already satisfied: idna<2.9,>=2.5 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from requests->refgenconf==0.10.0.dev0) (2.8)
    Requirement already satisfied: urllib3<1.25,>=1.21.1 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from requests->refgenconf==0.10.0.dev0) (1.24.1)
    Requirement already satisfied: oyaml in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from yacman>=0.7.0->refgenconf==0.10.0.dev0) (0.9)
    Requirement already satisfied: attrs>=17.4.0 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from jsonschema->refgenconf==0.10.0.dev0) (19.1.0)
    Requirement already satisfied: pyrsistent>=0.14.0 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from jsonschema->refgenconf==0.10.0.dev0) (0.14.11)
    Requirement already satisfied: six>=1.11.0 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from jsonschema->refgenconf==0.10.0.dev0) (1.12.0)
    Requirement already satisfied: setuptools in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from jsonschema->refgenconf==0.10.0.dev0) (41.0.1)
    Requirement already satisfied: colorama<0.5.0,>=0.4.0 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from rich->refgenconf==0.10.0.dev0) (0.4.1)
    Requirement already satisfied: pprintpp<0.5.0,>=0.4.0 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from rich->refgenconf==0.10.0.dev0) (0.4.0)
    Requirement already satisfied: typing-extensions<4.0.0,>=3.7.4 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from rich->refgenconf==0.10.0.dev0) (3.7.4.2)
    Requirement already satisfied: pygments<3.0.0,>=2.6.0 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from rich->refgenconf==0.10.0.dev0) (2.6.1)
    Requirement already satisfied: dataclasses<0.8,>=0.7; python_version >= "3.6" and python_version < "3.7" in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from rich->refgenconf==0.10.0.dev0) (0.7)
    Requirement already satisfied: commonmark<0.10.0,>=0.9.0 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from rich->refgenconf==0.10.0.dev0) (0.9.1)
    Building wheels for collected packages: refgenconf
      Building wheel for refgenconf (setup.py) ... [?25ldone
    [?25h  Created wheel for refgenconf: filename=refgenconf-0.10.0.dev0-cp36-none-any.whl size=64959 sha256=37191046ce6136b2bd777b1aa274a2d6a5ffb508af7e4969ac0ae97c1682b1f5
      Stored in directory: /private/var/folders/3f/0wj7rs2144l9zsgxd3jn5nxc0000gn/T/pip-ephem-wheel-cache-516dw93w/wheels/a8/b1/82/f79eaabaad4cf5c64fb4914e06dd04726c5c226785974aee4e
    Successfully built refgenconf
    Installing collected packages: refgenconf
      Found existing installation: refgenconf 0.9.3
        Uninstalling refgenconf-0.9.3:
          Successfully uninstalled refgenconf-0.9.3
    Successfully installed refgenconf-0.10.0.dev0
    [33mWARNING: You are using pip version 19.2.3, however version 20.2.3 is available.
    You should consider upgrading via the 'pip install --upgrade pip' command.[0m
    Collecting git+https://github.com/refgenie/refgenie.git@dev_config_upgrade
      Cloning https://github.com/refgenie/refgenie.git (to revision dev_config_upgrade) to /private/var/folders/3f/0wj7rs2144l9zsgxd3jn5nxc0000gn/T/pip-req-build-3i4zdr4w
      Running command git clone -q https://github.com/refgenie/refgenie.git /private/var/folders/3f/0wj7rs2144l9zsgxd3jn5nxc0000gn/T/pip-req-build-3i4zdr4w
      Running command git checkout -b dev_config_upgrade --track origin/dev_config_upgrade
      Switched to a new branch 'dev_config_upgrade'
      Branch 'dev_config_upgrade' set up to track remote branch 'dev_config_upgrade' from 'origin'.
    Requirement already satisfied: logmuse>=0.2.6 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenie==0.10.0.dev0) (0.2.6)
    Requirement already satisfied: piper>=0.12.1 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenie==0.10.0.dev0) (0.12.1)
    Requirement already satisfied: pyfaidx>=0.5.5.2 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from refgenie==0.10.0.dev0) (0.5.9.1)
    Requirement already satisfied: ubiquerg>=0.4.5 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from piper>=0.12.1->refgenie==0.10.0.dev0) (0.6.1)
    Requirement already satisfied: yacman in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from piper>=0.12.1->refgenie==0.10.0.dev0) (0.7.0)
    Requirement already satisfied: attmap>=0.12.5 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from piper>=0.12.1->refgenie==0.10.0.dev0) (0.12.12.dev0)
    Requirement already satisfied: psutil in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from piper>=0.12.1->refgenie==0.10.0.dev0) (5.6.1)
    Requirement already satisfied: pandas in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from piper>=0.12.1->refgenie==0.10.0.dev0) (1.0.3)
    Requirement already satisfied: setuptools>=0.7 in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from pyfaidx>=0.5.5.2->refgenie==0.10.0.dev0) (41.0.1)
    Requirement already satisfied: six in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from pyfaidx>=0.5.5.2->refgenie==0.10.0.dev0) (1.12.0)
    Requirement already satisfied: oyaml in /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages (from yacman->piper>=0.12.1->refgenie==0.10.0.dev0) (0.9)
    Requirement already satisfied: pyyaml>=3.13 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from yacman->piper>=0.12.1->refgenie==0.10.0.dev0) (5.1)
    Requirement already satisfied: pytz>=2017.2 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from pandas->piper>=0.12.1->refgenie==0.10.0.dev0) (2018.9)


    Requirement already satisfied: python-dateutil>=2.6.1 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from pandas->piper>=0.12.1->refgenie==0.10.0.dev0) (2.8.0)
    Requirement already satisfied: numpy>=1.13.3 in /Users/mstolarczyk/Library/Python/3.6/lib/python/site-packages (from pandas->piper>=0.12.1->refgenie==0.10.0.dev0) (1.17.3)
    Building wheels for collected packages: refgenie
      Building wheel for refgenie (setup.py) ... [?25ldone
    [?25h  Created wheel for refgenie: filename=refgenie-0.10.0.dev0-cp36-none-any.whl size=29266 sha256=d78485a0207036ddd91c36eb66b1973bdb3588aaff925d165d5e5aed483f968c
      Stored in directory: /private/var/folders/3f/0wj7rs2144l9zsgxd3jn5nxc0000gn/T/pip-ephem-wheel-cache-wmsjgl78/wheels/07/12/55/f50538357799dd2938a702a2f9e8b84a849975e61b0c59e7a0
    Successfully built refgenie
    Installing collected packages: refgenie
      Found existing installation: refgenie 0.9.3
        Uninstalling refgenie-0.9.3:
          Successfully uninstalled refgenie-0.9.3
    Successfully installed refgenie-0.10.0.dev0
    [33mWARNING: You are using pip version 19.2.3, however version 20.2.3 is available.
    You should consider upgrading via the 'pip install --upgrade pip' command.[0m



```bash
refgenie --version
```

    refgenie 0.10.0-dev | refgenconf 0.10.0-dev


Execution of refgenie commands fails since the config is incompatible:


```bash
refgenie list 
```

    Traceback (most recent call last):
      File "/Library/Frameworks/Python.framework/Versions/3.6/bin/refgenie", line 10, in <module>
        sys.exit(main())
      File "/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/refgenie/refgenie.py", line 821, in main
        skip_read_lock=skip_read_lock)
      File "/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/refgenconf/refgenconf.py", line 110, in __init__
        raise ConfigNotCompliantError(msg)
    refgenconf.exceptions.ConfigNotCompliantError: This genome config (v0.3) is not compliant with v0.4 standards. 
    To use current refgenconf, please use upgrade_config function to upgrade, ordowngrade refgenconf: 'pip install "refgenconf>=0.7.0,<0.10.0"'. 
    If refgenie is installed, you can use 'refgenie upgrade --target-version 0.4'; For config format documentation please see http://refgenie.databio.org/en/latest/genome_config/




## Upgrade refgenie configuration file

Let's upgrade the config to v0.4, just as the error message suggests. We will use `--force` option to run the command in batch mode.


```bash
refgenie upgrade --force --target-version 0.4
```

    Upgrading v0.3 config file format to v0.4
    Retrieved rCRSd digest from the server (511fb1178275e7d529560d53b949dba40815f195623bce8e)
    Retrieved human_repeats digest from the server (ebf26d2f064462bea7029e6b4d2298967d7435bff82ed224)
    Genome digest for human_alu is not available on any of the servers. Generating the digest from a local fasta file
    Loaded AnnotatedSequenceDigestList (8 sequences)
    Creating 'data' and 'alias' directories in '/Users/mstolarczyk/Desktop/testing/refgenie/upgrade_test'.
    Copying assets to 'data' and creating alias symlinks in 'alias'. Genomes that the digest could not be determined for 'will be ignored.
    Removing genome assets that have been copied to 'data' directory.


The upgrade succeded for all the assets that were previously managed by refgenie, regardless of the fact if the sequence-derived genome identifiers were avialable on the server. For ones that were not (`human_alu` genome) refgenie calculated the digest from the locally available FASTA file using the same algorithm that has been used to generate digests for the genomes on the server.


```bash
refgenie list
```

    [3m         Local refgenie assets          [0m
    [3m         Server subscriptions:          [0m
    [3m        http://rg.databio.org:82        [0m
    ┏━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┓
    ┃[1m [0m[1mgenome       [0m[1m [0m┃[1m [0m[1massets              [0m[1m [0m┃
    ┡━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━┩
    │ rCRSd         │ fasta, bowtie2_index │
    │ human_repeats │ fasta, bwa_index     │
    │ human_alu     │ fasta                │
    └───────────────┴──────────────────────┘



```bash

```
