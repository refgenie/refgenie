# Asset build output example

Command:
```
$ refgenie build -c genomes.yaml hg38/fasta --fasta hg38.fa.gz
```
Output:
```
Output to: hg38 /Users/mstolarczyk/Desktop/testing/test_genomes /Users/mstolarczyk/Desktop/testing/test_genomes/hg38
Removed existing flag: '/Users/mstolarczyk/Desktop/testing/test_genomes/hg38/refgenie_failed.flag'
### Pipeline run code and environment:

*              Command:  `/Library/Frameworks/Python.framework/Versions/3.6/bin/refgenie build -c genomes.yaml hg38/fasta --fasta hg38.fa.gz`
*         Compute host:  MichalsMBP
*          Working dir:  /Users/mstolarczyk/Desktop/testing/test_genomes
*            Outfolder:  /Users/mstolarczyk/Desktop/testing/test_genomes/hg38/
*  Pipeline started at:   (09-17 08:42:19) elapsed: 0.0 _TIME_

### Version log:

*       Python version:  3.6.5
*          Pypiper dir:  `/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.0dev
*         Pipeline dir:  `/Library/Frameworks/Python.framework/Versions/3.6/bin`
*     Pipeline version:  None

### Arguments passed to pipeline:

*            `command`:  `build`
*             `silent`:  `False`
*          `verbosity`:  `None`
*             `logdev`:  `False`
*      `genome_config`:  `genomes.yaml`
*            `recover`:  `False`
*        `config_file`:  `/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/refgenie/refgenie.yaml`
*          `new_start`:  `False`
*             `docker`:  `False`
*               `tags`:  `None`
*            `volumes`:  `None`
*          `outfolder`:  `/Users/mstolarczyk/Desktop/testing/test_genomes`
*       `requirements`:  `False`
*             `genome`:  `None`
* `asset_registry_paths`:  `['hg38/fasta']`
*              `fasta`:  `hg38.fa.gz`
*        `ensembl_gtf`:  `None`
*        `gencode_gtf`:  `None`
*                `gff`:  `None`
*            `context`:  `None`
*            `refgene`:  `None`

----------------------------------------

MissingAssetError: using 'default' as the default tag
Inputs required to build 'fasta': fasta
Building asset 'fasta'
Target to produce: `/Users/mstolarczyk/Desktop/testing/test_genomes/hg38/fasta/default/build_complete.flag`

> `cp hg38.fa.gz /Users/mstolarczyk/Desktop/testing/test_genomes/hg38/fasta/default/hg38.fa.gz` (38283)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.002GB.
  PID: 38283;	Command: cp;	Return code: 0;	Memory used: 0.002GB


> `gzip -d /Users/mstolarczyk/Desktop/testing/test_genomes/hg38/fasta/default/hg38.fa.gz` (38284)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 0.002GB.
  PID: 38284;	Command: gzip;	Return code: 0;	Memory used: 0.001GB


> `samtools faidx /Users/mstolarczyk/Desktop/testing/test_genomes/hg38/fasta/default/hg38.fa` (38285)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 0.005GB.
  PID: 38285;	Command: samtools;	Return code: 0;	Memory used: 0.005GB


> `cut -f 1,2 /Users/mstolarczyk/Desktop/testing/test_genomes/hg38/fasta/default/hg38.fa.fai > /Users/mstolarczyk/Desktop/testing/test_genomes/hg38/fasta/default/hg38.chrom.sizes` (38286)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.005GB.
  PID: 38286;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `touch /Users/mstolarczyk/Desktop/testing/test_genomes/hg38/fasta/default/build_complete.flag` (38288)
<pre>
psutil.ZombieProcess process still exists but it's a zombie (pid=38288)
Warning: couldn't add memory use for process: 38288
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.005GB.
  PID: 38288;	Command: touch;	Return code: 0;	Memory used: 0GB


> `cd /Users/mstolarczyk/Desktop/testing/test_genomes/hg38/fasta/default; find . -type f -exec md5sum {} \; | sort -k 2 | awk '{print $1}' | md5sum`
Default tag for 'hg38/fasta' set to: default
Computing initial genome digest...
Initializing genome...
Finished building asset 'fasta'

### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:10:23
*  Total elapsed time (all runs):  0:16:17
*         Peak memory (this run):  0.01 GB
*        Pipeline completed time: 2019-09-17 08:52:42

```
