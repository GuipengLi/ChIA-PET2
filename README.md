ChIA-PET2
=========

ChIA-PET2 is a versatile and flexible pipeline for analysing different variants of ChIA-PET data from raw sequencing reads to chromatin loops.

ChIA-PET2 was named not only because it is a tool for ChIA-PET data analysis, but also because it supports at least 2 different ChIA-PET protocols (bridge linker protocol or half-linkers protocol) data, 2 modes of read alignments (short or long read alignment) and 2 dimensional contact map output.

ChIA-PET2 integrates all steps required for ChIA-PET data analysis, including linkers trimming, reads mapping, duplicates removing, peaks calling and chromatin loops calling. It supports different kinds of ChIA-PET data generated from different ChIA-PET protocols. It also provides quality controls for different steps of ChIA-PET analysis. In addition, ChIA-PET2 can use phased genotype data to call allele-specific chromatin interactions. We applied ChIA-PET2 to different ChIA-PET datasets, demonstrating its significant improved performance as well as its ability to easily process ChIA-PET raw data.

Install
-------

ChIA-PET2 could be installed in a linux-like system. The ChIA-PET2 pipeline requires the following dependencies, which are usually already installed in a bioinformatics cluster.

- [MACS2](https://github.com/taoliu/MACS) v2.1.0+ recommended: for peaks calling
- [samtools](https://github.com/samtools/samtools) v1.3+: for sam file manipulation
- [bedtools](https://github.com/arq5x/bedtools2) v2.25.0+: for bed/bedpe file manipulation
- [R](https://www.r-project.org/) with the ggplot2 and VGAM packages : for calling MICC and plotting Quality Control figure

To install ChIA-PET2:

    tar -zxvf ChIA-PET2.tar.gz
    cd ChIA-PET2
    make

ChIA-PET2 will be installed in the bin directory in user home (~/bin) by default. It's recommended to set the bin/ directory to the PATH.


Usage
-----

Just type in **' ChIA-PET2 -h '** for detailed usage.

    $ ChIA-PET2 -h
    usage : ChIA-PET2 -g genomeindex -b bedtoolsgenome -f fq1 -r fq2 -A linkerA -B linkerB -o OUTdir -n prefixname
    Use option -h|--help for more information

    ChIA-PET2 0.9.2   2016.04.19
    ----------------------------
    OPTIONS

      -s|--start:     start from which step(1:8): 1:Trim Linkers; 2:Map Reads; 3:Build PETs; 4:Call Peaks; 5:Find Interactions
                      6:Plot QC; 7:Estimate statistical confidence; 8:Phase PETs(optional), default=1
      -g|--genome:    genome index for bwa
      -b|--bedtoolsgenome: chromsomes size file for bedtools
      -f|--forward:   one fastq(.gz) file
      -r|--reverse:   the other fastq(.gz) file
      -A|--linkerA:   one linker sequence, default=GTTGGATAAG
      -B|--linkerB:   the other linker sequence, default=GTTGGAATGT
      -o|--output:    output folder, default=output
      -n|--name:      output prefix name, default=out
      -m|--mode:      0,1,2;  0: A/B linkers; 1: bridge liker; 2: Enzyme site, default=0
      -e|--err:       Maximum mismatches allowed in linker sequence, default=0
      -k|--keepempty: 0,1,2; 0:No linker-empty reads; 1:keep 1 linker-empty read; 2:keep 2 linker-empty reads. default=0
      -t|--thread:    threads to run, default=1
      -d|--short:     short reads (0 or 1), default=0
      -M|--macs2 parameters, default="-q 0.05"
      -Q|--mapq:      mapq cutoff, default=30
      -C|--cutoffPET: PET count cutoff before running MICC, default=1
      -S|--slop:      slop length, default=100
      -E|--extend:    extend length on both sides, default=500
      -l|--length:    min length of reads after linker trimming. default=15
      -P|--phased:    optional phased genotype file: 'chr1\tstart\tend\tA\tC'
      [-h|--help]:    help
      [-v|--version]: version


Output
------

Important result files:

- **prefixname.interactions.intra.bedpe**: inra-chromosomal loops (11 columns)
- **prefixname.interactions.inter.bedpe**: iner-chromosomal loops (11 columns)
- **prefixname.interactions.MICC**: significance estimation for each loops (13 columns)
- **prefixname.QCplot.pdf**: Quality control figure for different steps of analysis.

prefixname.interactions.MICC file has 11+2 columns. The last 2 columns ( **-log10(1-PostProbability)** and **FDR** ) are estimated by MICC based on a Bayesian mixture model. See the toy example below. This means the chromatin loop between peak_1(chr:9118-10409) with peak depth 3330 and peak_3(chr:89064-90360) with peak depth 3814 has 39 supportive pair-end tags(PETs). **-log10(1-PostProbability)=4.5** and **FDR=0.03** .

|chr |start|end  |chr |start |end  |peak1 |peak2 |depth1|depth2|#PET|PP | FDR|
|----|-----|-----|----|------|-----|------|------|------|------|----|---|----|
|chr1|9118 |10409|chr1|89064 |90360|peak_1|peak_3|3330  |3814  |39  |4.5|0.03|


Other intermediated files include:

- prefixname_1.valid.sam
- prefixname_2.valid.sam
- prefixname.bedpe
- prefixname.rmdup.bedpe
- prefixname_peaks.slopPeak
- prefixname.trim.stat
- prefixname.bedpe.stat


Toolkits
--------
**trimLinker:** trim linker sequences in the pair-end fastq(.gz) files. There will be also a summary file named output.trim.stat

    $ Usage:   trimLinker [-option] [argument]
    option:  -h  show help information
         -t  threads
         -e  number of mismatch allowed in linker, default=0
         -k  keep empty: 0, 1, 2
         -l  min length of timmed reads: default=15
         -o  output directory
         -m  mode: 0 or 1, A/B linkers(0) or bridge linker(1)
         -A  linkerA
         -B  linkerB
         -n  output name prefix
    example: trimLinker -t 4 -o outdir -m 0 -n name -k 0 -A AAAAAAAT -B TTTTACGG fq1.fq.gz fq2.fq.gz


**buildBedpe:** build bedpe file from two sam files. The read names should be in the same order. keepseq_flag indicates whether to keep sequences in the bedpe file (default=0, no need to keep sequences unless we want to do allele-specific analysis).

    $ buldBedpe file1.sam file2.sam output  MAPQ_cutoff thread keepseq_flag


**removeDup:** remove duplicate PETs

    $ removeDup in.bedpe out.bedpe


**buildTagAlign:** build tag file from bedpe file for MACS2 input.

    $ buildTagAlign in.bedpe out.tag.bed


**bedpe2Interaction:** Detect chromatin interactions from bedpe file.

    $ pairToBed -a in.bedpe -b peaks.bed -type both > tmp.bedpe
    $ bedpe2Interaction tmp.bedpe interaction.out.bedpe out.bedpe.stat


**QCplots.R:** Plot the Quality Control figure.

    $ Rscript QCplots.R directory name


**MICC2.R:** Call MICC package to estimate the statistically significance of chromatin loops.

    $ Rscript MICC2.R intra.bedpe inter.bedpe miccOUT CUTOFFPET


**bedpe2Matrix:** Generate the Hi-C style matrix.

    $ bedpe2Matrix --binsize 10000 --chrsizes bedtoolsgenome.txt --ifile in.rmDup.bedpe --oprefix PREFIX --progress


**bedpe2Phased:** Generate the allele-specific PETs. Each line of the phased genotype file phased.bed has 5 columns separated with tabs '\t': "chr1 start end A C"

    $ grep -v "*" in.bedpe > bothmapped.bedpe
    $ pairToBed -a bothmapped.bedpe -b phased.bed -type either > bedpe.flt.bed
    $ bedpe2Phased bedpe.flt.bed outputprefix


Citation
--------

Please cite the following article if you use ChIA-PET2 in your research:

- ChIA-PET2: an versatile and flexible pipeline for ChIA-PET data analysis (submitted)

If you use MICC to call significant chromatin loops (step 7 in the pipeline), please cite the following:

- He C, Zhang MQ, Wang X: MICC: an R package for identifying chromatin interactions from ChIA-PET data. Bioinformatics 2015, 31:3832-3834.


Contact
-------

Author: Guipeng Li

Email:  guipeng.lee(AT)gmail.com
