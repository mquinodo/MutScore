#!/usr/bin/perl
use strict;                # Terminate when error occurs
use warnings;              # Display all warning messages

###################
### EXPLANATION ###
###################
### script that adds VEST3 and PROVEAN SCORES of deleteriousness

##################
# get "input file" and "$fishtank" and the "output file" from command line
my ( $fishtank, $outputFile, $bait) = @ARGV;

#####################################
### Optional but recommended
#####################################

if (!(-f $fishtank)){	die "Invalid input file";	}
if ( !defined($outputFile) || $outputFile eq ''){	die "Missing output file";	}
if (!(-f $bait)){	die "Missing bait file";	}

#####################################
### END OF OPTIONAL PART
#####################################

open( IN, $fishtank ) or die "Unable to open file $fishtank :$!"; # open the fishtank file (in reading mode)
open( OUT, ">".$outputFile ) or die "Unable to open file $outputFile :$!"; # open (or create) the output file in writing mode
open( PIPPO, $bait ) or die "Unable to open file $bait :$!"; 

print OUT "##Annotation\tadd annovar information\n";

### CREATE THE DICT OF VARIANTS
my $header="";
my $number_columns=0;
my %list_of_position; # dict containing the variants with specific position
my %list_of_rs; #in case it is recognazed by the SNP ID only
while (<PIPPO>){
	chomp;
	if (/^Chr/){ #select the header for the new columns
		my @array=split(/\t/);
		$header=join("\t",@array[5..$#array]);
		$number_columns=$#array+1-5;
	}
	else{ #extract the variants information
		my @array=split(/\t/);
		my ($chr, $begin, $end, $ref, $alt)=split(/\t/);
		my $info=join("\t",@array[5..$#array]);
		# $begin=$begin+1; #from 0-based to 1-based
		$list_of_position{"$chr\t$begin\t$ref\t$alt"}=$info;	
	}
}
close(PIPPO);

# print $number_columns,"\n";

my $variants=0; #number of variants in the input file
my $number=0; #number of variants found in the bait
my $rs=0;
# BROWSE THE VARIANTS AND ADD THE INFORMATION IF PRESENT IN THE BAIT FILE
while (<IN>){
	if (/^#line|^>q|^##Var_type/){ # HEADER
		chomp;
		print OUT $_,"\t$header\n";# add the column in the header
	}
	elsif(/^##/){ #Annotation header
		print OUT; 
	}
	if (/^\w/){ ## VARIANTS
	# else{
		chomp;
		$variants++;
		my @array=split("\t");
		my ($line,$type,$var_info,$chr,$begin,$end,$ref,$alt) = split(/\t/); #list of elements in the record


		my @record = split(/\t/); #list of elements in the record
		my $annotation="NA\t"x($number_columns-1)."NA"; #empty: list of NA (the two 0 are for # patients in Control and Nexome)

		#if variant exists in the bait based of the position and ref / alt
		if( exists $list_of_position{"$chr\t$begin\t$ref\t$alt"}) { 
			$number++;
			$annotation = $list_of_position{"$chr\t$begin\t$ref\t$alt"};
		}
		$annotation =~ s/\t\t/\t.\t/g;
		print OUT $_,"\t$annotation\n";
	}
}
close(IN);
close(OUT);

print "INFO: $number variants on $variants are annotated for their variants.\n"; # ($rs based only on rs number).\n";
exit;
