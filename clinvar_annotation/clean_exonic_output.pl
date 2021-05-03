#!/usr/bin/perl
#!/bin/bash
use strict;                # Terminate when error occurs
use warnings;              # Display all warning messages

###################
### EXPLANATION ###
###################
### script that cleans the input file into a more readable file (arrange the format, delete some doublet column, ...)


##################
# get "input file" and "$fishtank" and the "output file" from command line
my ( $fishtank, $outputFile) = @ARGV;

#####################################
### Optional but recommended
#####################################

if (!(-f $fishtank)){	die "Invalid input file";	}
if ( !defined($outputFile) || $outputFile eq ''){	die "Missing output file";	}

#####################################
### END OF OPTIONAL PART
#####################################

open( IN, $fishtank ) or die "Unable to open file $fishtank :$!"; # open the fishtank file (in reading mode)
open( OUT, ">".$outputFile ) or die "Unable to open file $outputFile :$!"; # open (or create) the output file in writing mode

print OUT "##Cleaning\tarrange the output file\n";
# header of the output files

## DICT OF GENES
while (<IN>){
	if(/^#line/){
		chomp;
		my @record = split(/\t/); #list of elements in the record
		print OUT "##$record[1]\tgene\t$record[2]\tchr\tbegin\tend\t$record[14]\t$record[15]\t$record[8]\t$record[13]\n";
		#print OUT "\tread_depth\tread_ref_alt\t%_alt\tgenotype_quality\tQD\tMQ\tBQRS\tSOR\tFS\tMQRS\tRPRS\n";
	}
	if(/^##/){
		print OUT;
	}
	if (/^\w/){
		chomp;
		my @record = split(/\t/); #list of elements in the record
		my $annotationFinal=join("\t",@record[21..$#record]);
		my $variant=$record[1];
		my $gene_list="";

		if(/exonic;splicing/){
			$gene_list=(split(/;/,$record[2]))[0];
			my $gene_splicing=(split(/\(/,(split(/;/,$record[2]))[1]))[0];
			unless($gene_list =~ /\Q$gene_splicing\E/){
				$gene_list.=";$gene_splicing";
			}
		}
		elsif(/splicing/){
			$gene_list=(split(/\(/,$record[2]))[0];
		}
		else{
			my $temp = $record[2];
			$temp =~ s/\([^())]*\)//g;
			for my $gene_cut (split(/,/,$temp)){
				my $gene_name=(split(/:/,$gene_cut))[0];
				unless($gene_list =~ /\Q$gene_name\E/){
					if($gene_list eq ""){
						$gene_list="$gene_name";	
					}
					else{
						$gene_list.=";$gene_name";	
					}
				}
			}
		}

		my $iso=$record[2];
		my $chr=$record[3];
		my $zyg=$record[8];
		my $start=$record[12];
		my $rs=$record[9];
		my $ref=$record[14];
		my $alt=$record[15];
		my $end=$start+length($alt);
		

		print OUT "$variant\t$gene_list\t$iso\t$chr\t$start\t$end\t$ref\t$alt\t$zyg\t$rs\n";
		#print OUT "$read_depth\t$read_ref_alt\t$perc_alt\t$genotype_quality\t$QD\t$MQ\t$BQRS\t$SOR\t$FS\t$MQRS\t$RPRS\n";
	}
}
close(IN);
close(OUT);

exit;
