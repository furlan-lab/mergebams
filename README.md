<p align="center"><img src="mergebams.png" alt="" width="250"></a></p>
<hr>

version 0.2

Merge sam/bam files with intelligent cell barcode preservation.  This has been tested on bam file and tsv output from the 10X Genomics Cellranger program.  The implementation of mergebams was motivated by and primarily designed for working with Cellranger output.

## Updates

version 0.2, mergebams will check the headers of all bams to ensure they contain the same sequence names.  mergebams will not merge bams that do not have the same reference names.  mergebams will use the header from the first bam in the input argument. mergebams will now add a header comment that includes the filenames of the bams that were merged.

version 0.1 - first version in rust - implemented using the python version as a reference. Rust implementation is faster and supports multi-threading.

## Requirements

1. Rust >  (mergebams uses the bam crate as a dependency)
2. Samtools is required to run the vignette below

## Installation of rust (only if needed and only tested on Mac/Linux)

To install rust visit https://www.rust-lang.org/tools/install.


## Compile mergebams

#### Clone repo 

```bash
git clone https://github.com/furlan-lab/mergebams.git
cd mergebams
```

#### Compile

mergebams will compile creating an executable file in thge target/release directory called mergebams.  This binary is a self-contained file and can be moved to a location that is suitable.  In this vignette we will use mergebames in the target/release directory

```bash
cargo build --release
```

## Test installation of mergeBams

You should then be able to test installation on the provided toy bam and barcode files.  After running the following, you should see the help screen displayed.

```bash
target/release/mergebams -h
```

## Help

```bash
Merge sam/bam files with intelligent cell barcode preservation. This has been tested on bam file and tsv output from the
10X Genomics Cellranger program. The implementation of mergebams was motivated by and primarily designed for working
with Cellranger output.

USAGE:
    mergebams [FLAGS] [OPTIONS] --inputs <inputs> --labels <labels>

FLAGS:
    -c, --celltag    set if cell barcode tag should not be CB
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -b, --bcs <bcs>          barcodes files, comma-separated gzipped or not
    -i, --inputs <inputs>    input bams, comma-separated
    -l, --labels <labels>    strings for prepending cell barcode (i.e. sample name), comma-separated
    -o, --out <out>          folder for output (merged_bam.bam, merged_bcs.tsv.gz)
    -t, --out <threads>      threads
```

## Usage

The following is an example of merging two bam files and two barcodes.tsv files that were derived from them.  They are provided in the test directory

```bash
cd test
../target/release/mergebams -i bam1.bam,bam2.bam \
          -l t1_,t2_ \
          -b barcodes1.tsv.gz,barcodes2.tsv.gz \
          -o .
```

## Expected output

**In the above example mergebams will take input bams bam1.bam and bam2.bam which have the following data...**

```bash
samtools view bam1.bam | head -n 3 -
```

```bash
VH00738:4:AAAW2TWHV:2:2504:18932:30136  16  1 633374  3 90M *0  0 TGCCCATCGTCCTAGAATTAATTCCCCTAAAAATCTTTGAAATAGGGCCCGTATTTACCCTATAGCACCCCCTCTACCCCCTCTAGAGCC  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  NH:i:2  HI:i:1  AS:i:88 nM:i:0  RG:Z:ITS_D383_3:0:1:AAAW2TWHV:2 RE:A:I  xf:i:0  CR:Z:ATTGGACAGTCATGCT CY:Z:CCCCCCCCCCCCCCCC CB:Z:ATTGGACAGTCATGCT-1 UR:Z:CGGATCTGGT UY:Z:CCCCCCCCCC UB:Z:CGGATCTGGT
VH00738:4:AAAW2TWHV:1:2502:17909:53799  16  1 1014103 255 90M *0  0 CCAGCAGCGTCTGGCTGTCCACCCGAGCGGTGTGGCGCTGCAGGACAGGGTCCCCCTTGCCAGCCAGGGCCTGGGCCCCGGCAGCACGGT  CCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  NH:i:1  HI:i:1  AS:i:88 nM:i:0  RG:Z:ITS_D383_3:0:1:AAAW2TWHV:1 TX:Z:ENST00000379389,+273,90M;ENST00000624652,+324,90M;ENST00000624697,+349,90M GX:Z:ENSG00000187608  GN:Z:ISG15  fx:Z:ENSG00000187608  RE:A:E  xf:i:25 CR:Z:TTTACTGAGTCGATAA CY:Z:CCCCCCCCCCCCCCCC CB:Z:TTTACTGAGTCGATAA-1 UR:Z:GCCTCTTCCG UY:Z:CCCCCCCCCC UB:Z:GCCTCTTCCG
VH00738:4:AAAW2TWHV:1:2214:43700:49596  0 1 1374390 255 56M259N34M  * 0 0 GGGCCCGCAGACCCGGCTGCCCAGCACTCCAGAGACGGGCCAAGGCGGGCGGCCGCCTGCCCAAGGAACGGCCCTCAACAGCTGGGAAGT  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCC  NH:i:1  HI:i:1  AS:i:90 nM:i:0  RG:Z:ITS_D383_3:0:1:AAAW2TWHV:1 TX:Z:ENST00000321751,+147,90M;ENST00000338338,+394,90M;ENST00000338370,+419,90M;ENST00000378853,+229,90M  GX:Z:ENSG00000175756  GN:Z:AURKAIP1 fx:Z:ENSG00000175756  RE:A:E  xf:i:25CR:Z:ATCATGGCAGACGCTC  CY:Z:CCCCCCCCCCCCCCCC CB:Z:ATCATGGCAGACGCTC-1 UR:Z:GCATTATAGC UY:Z:CCCCCCCCCC UB:Z:GCATTATAGC
```

**AND**

```bash
samtools view bam2.bam | head -n 3 -
```

```bash
VH00738:4:AAAW2TWHV:2:1212:32395:39166  0 1 631094  3 90M *0  0 ATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTAT  CCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  NH:i:2  HI:i:1  AS:i:88 nM:i:0  RG:Z:ITS_D544_3:0:1:AAAW2TWHV:2 RE:A:I  xf:i:0  CR:Z:CCCAATCTCCTAAGTG CY:Z:CCCCCCCCCCCCCCCC CB:Z:CCCAATCTCCTAAGTG-1 UR:Z:TATATGTTTG UY:Z:CCCCCCCC;C UB:Z:TATATGTTTG
VH00738:4:AAAW2TWHV:2:1402:46048:13060  16  1 633977  255 90M *0  0 AACCACCCAACTATCTATAAACCTAGCCATGGCCATCCCCTTATGAGCGGGCGCAGTGATTATAGGCTTTCGCTCTAAGATTAAAAATGC  ;CCC-CC;CCC;;-CCC;C;C-CCCCCC;CCCC;C-CC;-;-C;;;CCC--C;CCCCCCCCCCC;C-CCC-C-CCCCCCCCCCCCCCCCC  NH:i:1  HI:i:1  AS:i:88 nM:i:0  RG:Z:ITS_D544_3:0:1:AAAW2TWHV:2 RE:A:I  xf:i:0  CR:Z:CATTCGCTCCTGCTTG CY:Z:CCCCCCCCC;CCCCCC CB:Z:CATTCGCTCCTGCTTG-1 UR:Z:TGTCATCAGA UY:Z:C--;CCCCCC UB:Z:TGTCATCAGA
VH00738:4:AAAW2TWHV:2:1605:40556:12681  0 1 633982  255 90M *0  0 CCCAACTATCTATAAACCTAGCCATGGCCATCCCCTTATGAGCGGGCGCAGTGATTATAGGCTTTCGCTCTAAGATTAAAAATGCCCTAG  -C;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCC;CCCCCCCCCCCC-CC  NH:i:1  HI:i:1  AS:i:88 nM:i:0  RG:Z:ITS_D544_3:0:1:AAAW2TWHV:2 RE:A:I  xf:i:0  CR:Z:CACATTTTCTTTAGGG CY:Z:CCCCCCCCCC;CCCCC CB:Z:CACATTTTCTTTAGGG-1 UR:Z:CTTAAACGGT UY:Z:CCC;CCCC-C UB:Z:CTTAAACGGT
```

**These bam files will be concatenated but will prepend the cell barcode (CB tag) with the label supplied in the program call using the -l flag**

```bash
(samtools view out_bam.bam | head -n 3 -; samtools view out_bam.bam | tail -n 3 -) > topandbottom.txt
cat topandbottom.txt
```

```bash
VH00738:4:AAAW2TWHV:2:2504:18932:30136  16  1 633374  3 90M *0  0 TGCCCATCGTCCTAGAATTAATTCCCCTAAAAATCTTTGAAATAGGGCCCGTATTTACCCTATAGCACCCCCTCTACCCCCTCTAGAGCC  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  NH:i:2  HI:i:1  AS:i:88 nM:i:0  RG:Z:ITS_D383_3:0:1:AAAW2TWHV:2 RE:A:I  xf:i:0  CR:Z:ATTGGACAGTCATGCT CY:Z:CCCCCCCCCCCCCCCC UR:Z:CGGATCTGGT UY:Z:CCCCCCCCCC UB:Z:CGGATCTGGT CB:Z:t1_ATTGGACAGTCATGCT-1
VH00738:4:AAAW2TWHV:1:2502:17909:53799  16  1 1014103 255 90M *0  0 CCAGCAGCGTCTGGCTGTCCACCCGAGCGGTGTGGCGCTGCAGGACAGGGTCCCCCTTGCCAGCCAGGGCCTGGGCCCCGGCAGCACGGT  CCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  NH:i:1  HI:i:1  AS:i:88 nM:i:0  RG:Z:ITS_D383_3:0:1:AAAW2TWHV:1 TX:Z:ENST00000379389,+273,90M;ENST00000624652,+324,90M;ENST00000624697,+349,90M GX:Z:ENSG00000187608  GN:Z:ISG15  fx:Z:ENSG00000187608  RE:A:E  xf:i:25 CR:Z:TTTACTGAGTCGATAA CY:Z:CCCCCCCCCCCCCCCC UR:Z:GCCTCTTCCG UY:Z:CCCCCCCCCC UB:Z:GCCTCTTCCG CB:Z:t1_TTTACTGAGTCGATAA-1
VH00738:4:AAAW2TWHV:1:2214:43700:49596  0 1 1374390 255 56M259N34M  * 0 0 GGGCCCGCAGACCCGGCTGCCCAGCACTCCAGAGACGGGCCAAGGCGGGCGGCCGCCTGCCCAAGGAACGGCCCTCAACAGCTGGGAAGT  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCC  NH:i:1  HI:i:1  AS:i:90 nM:i:0  RG:Z:ITS_D383_3:0:1:AAAW2TWHV:1 TX:Z:ENST00000321751,+147,90M;ENST00000338338,+394,90M;ENST00000338370,+419,90M;ENST00000378853,+229,90M  GX:Z:ENSG00000175756  GN:Z:AURKAIP1 fx:Z:ENSG00000175756  RE:A:E  xf:i:25CR:Z:ATCATGGCAGACGCTC  CY:Z:CCCCCCCCCCCCCCCC UR:Z:GCATTATAGC UY:Z:CCCCCCCCCCUB:Z:GCATTATAGC  CB:Z:t1_ATCATGGCAGACGCTC-1
VH00738:4:AAAW2TWHV:2:2502:67066:55086  4 * 0 0 * *0  0 TCGTACGCTTCCGTGTTCCTCATTAAAGGCCGAGCCCATATAAGAAATATTAGACAGACGTTGTGAGATGTTCAGATCGGAAGAGCGTCG  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  RG:Z:ITS_D544_3:1:1:AAAW2TWHV:2fr:Z:CCGTGTTCCTCATTA fq:Z:CCCCCCCCCCCCCCC  fb:Z:CCGTGTTCCTCATTA  fx:Z:CD71 xf:i:24 CR:Z:GAACATCTCACAACGT CY:Z:CCCCCCCCCCCCCCCC UR:Z:CTGTCTAATGUY:Z:CCCCCCCCCC  UB:Z:CTGTCTAATG CB:Z:t2_GAACATCTCACAACGT-1
VH00738:4:AAAW2TWHV:2:2309:76742:23037  4 * 0 0 * *0  0 ATTTTCTTTTGCGATGGTAGATTATGGTAAGTGGCCCATATAAGAAAATTCGGACATGAAGTTAGACGCAAAGAGATCGGAAGAGCGTCG  CCCCCCCC;C;-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CC;CCCCCCCCCCCC;CCCCC;CCCCCCCCCCCCCCCCCC  RG:Z:ITS_D544_3:1:1:AAAW2TWHV:2fr:Z:GCGATGGTAGATTAT fq:Z:;-CCCCCCCCCCCCC  fb:Z:GCGATGGTAGATTAT  fx:Z:CXCR3  xf:i:24 CR:Z:CTTTGCGTCTAACTTC CY:Z:CCCCC-CC;CCC;CCC UR:Z:ATGTCCGAATUY:Z:CC;C-CCCC-  UB:Z:ATGTCCGAAT CB:Z:t2_CTTTGCGTCTAACTTC-1
VH00738:4:AAAW2TWHV:2:1605:38966:50732  4 * 0 0 * *0  0 CGTCAATGTGCCGTGTTCCTCATTATAAGGGCCCCCCATATAAGAAATATTTGGTGGACACTCAGACGAATGGAGATCGGAAGAGCGTCG  CCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  RG:Z:ITS_D544_3:1:1:AAAW2TWHV:2fr:Z:CCGTGTTCCTCATTA fq:Z:CCCCCCCCCCCCCCC  fb:Z:CCGTGTTCCTCATTA  fx:Z:CD71 xf:i:24 CR:Z:CCATTCGTCTGAGTGT CY:Z:CCCCCCCCCCCCC;CC UR:Z:CCACCAAATAUY:Z:CCCCCCCCCC  UB:Z:CCACCAAATA CB:Z:t2_CCATTCGTCTGAGTGT-1
```

Similarly and if desired, mergebams will concatenate and add labels to barcodes.tsv files (for compressed barcodes.tsv.gz see below for an explanation of how compression of barcodes files are handled).  For example, in the above case...


```bash
zcat barcodes1.tsv.gz | head -n 3
```

```bash
AAACCTGAGCACCGTC-1
AAACCTGAGCGATATA-1
AAACCTGAGCTCTCGG-1
```

**AND**

```bash
zcat barcodes2.tsv.gz | tail -n 3
```

```bash
TTTGTCATCCAAACAC-1
TTTGTCATCCGCGGTA-1
TTTGTCATCTCGCTTG-1
```

Will be joined and given labels.

```bash
(zcat out_barcodes.tsv.gz | head -n 3; zcat out_barcodes.tsv.gz | tail -n 3 ) > topandbottombc.txt
cat topandbottombc.txt
```

```bash
t1_AAACCTGAGCACCGTC-1
t1_AAACCTGAGCGATATA-1
t1_AAACCTGAGCTCTCGG-1
t2_TTTGTCATCCAAACAC-1
t2_TTTGTCATCCGCGGTA-1
t2_TTTGTCATCTCGCTTG-1
```


Note that this program is compression aware and will detect compression of the barcodes file.  



## Acknowledgements

Written by Scott Furlan with help from cfcooldood and rcguy
