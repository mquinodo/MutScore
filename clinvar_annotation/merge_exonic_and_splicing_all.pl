#!/usr/bin/perl
use strict;                # Terminate when error occurs
use warnings;              # Display all warning messages

###################
### EXPLANATION ###
###################
### script that merge the exonic variants + the splicing variants into 1 file

##################
# get "input file" and "$fishtank" and the "output file" from command line
my ( $fishtank, $outputFile, $bait) = @ARGV;

#####################################
### Optional but recommended
#####################################

if (!(-f $fishtank)){	die "Invalid input file";	}
if ( !defined($outputFile) || $outputFile eq ''){	die "Missing output file";	}
if (!(-f $bait)){	die "Invalid bait file";	}

#####################################
### END OF OPTIONAL PART
#####################################


open( IN, $fishtank ) or die "Unable to open file $fishtank :$!"; # open the fishtank file (in reading mode)
open( OUT, ">".$outputFile ) or die "Unable to open file $outputFile :$!"; # open (or create) the output file in writing mode
open( PIPPO, $bait ) or die "Unable to open file $bait :$!"; # open the fishtank file (in reading mode)

## Header
print OUT "##filtering\texonic and splicing variants\n";
print OUT "#line\tVar_type\tvariant_info\tchrom\tbegin\tend\tref\talt\tzyg\t.\t.\tchrom\tbegin\trs\tref\talt\t.\t.\tvcf\theader_Qual\tQual\n";


#### EXTRACT THE SPLICING VARIANTS
my $number_splicing=0;
my $number_exonic=0;
while (<PIPPO>){
	# if(/^splicing|^exonic;splicing/){ ## do not take ncRNA_splicing into account
	# 	print OUT "line\t",$_;
	# 	$number_splicing++;
	# }
	# elsif(/^exonic/){
	# 	$number_exonic++;
	# }
	if(/^inensdldldl/){
		#do nothing
	}
	else{
		#if(/^exonic|^intergenic|^ncRNA|^downstream|^upstream|^UTR/){
		if(/^exonic/){
			$number_exonic++;	
		}
		else{
			print OUT "line\t",$_;
			$number_splicing++;
		}
	}
}
close(PIPPO);


## print the exonic variants
my $variants=0;
while (<IN>){
	$variants++;
	print OUT $_;
}
close(IN);
close(OUT);


print "INFO: $number_splicing splicing variants and $variants exonic variants in the file ($number_exonic from variant_function_file)\n";	
print "test\n";


exit;
