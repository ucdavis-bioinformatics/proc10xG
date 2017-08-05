# proc10G set of python scripts

A series of python scripts to process data generated using the 10x genomics DNA system. scripts are designed to extract and trim the reads of gem barocode informaton and primer sequence respectively. Compare the gem barcodes to a whitelist (allowing for 1 mismatch) and annotating the read as having a barcode which MATCHES a barcode int the whitelist, contains MISMATCH1, is AMBIGUOUS (with 1 edit distance, matches multiple whitelisted barcodes), or UNKONWN does not match any whitelisted barcodes at edit distance 1.

Scripts ready for use
* proc10xGenomics.py - process fastq files generated from bcl2fastq, longranger mkfastq, or supernova mkfastq
* samConcat2Tag.py - extract the FASTA/FASTQ comment appended to SAM output from bwa mem and generates sam tags

Scripts in progress, not ready for use
* profile_mapping.py - profile the gem barcode alignments
* process_mapping.py - remap ambiguous alignment using gem barocode to identify correct placement 

## proc10xGenomics.py

	Usage: usage process_10xReads.py -o [output file prefix (path + name)] -(aipbtg) --quiet -1 [read1] -2 [read2]
	process_10xReads.py will process read file produced by 10x genomics and do some stuff

	Options:
	  --version             show program's version number and exit
	  -h, --help            show this help message and exit
	  -o OUTPUT_DIR, --output=OUTPUT_DIR
	                        Directory + prefix to output reads, if stdout
	  -a, --all             output all reads, not just those with valid gem
	                        barcode, STATUS will be UNKNOWN, or AMBIGUOUS
	  -i                    output in interleaved format, if -o stdout,
	                        interleaved will be chosen automatically
	  -p                    profile the reads and barcodes, FUTURE
	  -b BCTRIM, --bctrim=BCTRIM
	                        trim gem barcode
	  -t TRIM, --trim=TRIM  trim addional bases after the gem barcode
	  -g, --nogzip          do not gzip the output, ignored if output is stdout
	  --quiet               turn off verbose output

	  Inputs:
	    10x fastq files to input

	    -1 read1, --read1=read1
	                        read1 of a pair, multiple files can be specified
	                        separated by comma
	    -2 read2, --read2=read2
	                        read2 of a pair, multiple files can be specified
	                        separated by comma

## Examples

python proc10xGenomics.py -o testing -1 testdata/CaCon-sm_R1_001.fastq.gz -2 testdata/CaCon-sm_R2_001.fastq.gz

bwa mem -C testdata/polished_p_ctg.fa testing_R1_001.fastq testing_R2_001.fastq | python samConcat2Tag.py | samtools sort -n - | samtools view -h -o mapping.sam -

bwa mem -t 1 -p -C testdata/polished_p_ctg.fa testing_R1_001.fastq testing_R2_001.fastq | python samConcat2Tag.py | samtools sort -n -o mapping.bam -

python samConcat2Tag.py saved.sam | samtools sort - | samtools view

python proc10xGenomics.py -a -1 data/CaCon-sm_R1_001.fastq.gz -2 data/CaCon-sm_R2_001.fastq.gz | bwa mem -t 1 -p -C data/polished_p_ctg.fa - | python samConcat2Tag.py | samtools sort -m 768M --threads 0 -n -o mapping.bam -

#################################################################################################################
## first process reads with process_10xReads.py which extracts the GEM barcode and primer sequence,
##   then compares the barcode to a white list, marking reads with status
##	   MATCH - a perfect match of bc to the whitelist
##     MISMATCH1 - a single mismatch to only one bc on the whitelist
##     AMBIGUOUS - a single mismatch to more than one bc on the whitelist
##     UNKNOWN - greater than 1 mismatch from any bc on the whitelist
##   then appends the status, library barcocde, GEM barcode, primer sequences and cooresponding
##		quality scores to the comment of the read ID and the whitelisted barcode to the
##      beginning of the read, in interleaved format
##  Then map to the genome using bwa mem with flags -p (interleaved) and -C (appends comment to the sam file)
##  Next process with samContcat2Tag.py which extracts the appended commend and add the following tags
##     ST:Z - Read status
##     BX:Z - GEM Barcode ID (with appended sample '-1'), whitelisted ID
##     BC:Z - Library Barcode
##     QT:Z - Library Barcode Quality (if Index read not provided then all '!')
##     RX:Z - GEM Barcode Sequence
##     QX:Z - GEM Barcode Quality
##     TR:Z - Primer Sequence
##     TQ:Z - Primer Quality
##  Finally, sort using samtools sort, sorting on reads ID (GEM Barcode)
#################################################################################################################

python process_10xReads.py -a -1 data/CaCon-sm_R1_001.fastq.gz \
  -2 data/CaCon-sm_R2_001.fastq.gz | \
  bwa mem -t 1 -p -C data/polished_p_ctg.fa - | python samConcat2Tag.py | samtools sort -n -o mapping.bcmapped.bam - 2> stderr.out > stdout.out
