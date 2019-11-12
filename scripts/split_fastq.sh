awk \
	-v OUTFASTQ="$1_OUTFOLDER/split_fastq" \
	'
	BEGIN {
		CUROUT = OUTFASTQ"/test.fastq";
	}
	((FNR) % 4 == 1) {
		CUROUT = OUTFASTQ"/"substr($1,2)".fastq";
	}
	{
		print $0 > CUROUT
	}
	' $2_INPUT
