#!/usr/bin/env perl
# $Id: get_parent_terms.pl 2015-10-29 erick.antezana $
#
# Script  : get_parent_terms.pl
# Purpose : Collects the parent terms (not all the ancestors) of a given term in the given OBO ontology
# Usage   : get_parent_terms.pl my_ontology.obo term_id > parent_terms.txt
# License : Copyright (c) 2006-2015 by Erick Antezana. All rights reserved.
#           This program is free software; you can redistribute it and/or
#           modify it under the same terms as Perl itself.
# Contact : Erick Antezana <erick.antezana -@- gmail.com>
#
################################################################################

use Carp;
use strict;
use warnings;
use Data::Dumper;

use OBO::Parser::OBOParser;
use OBO::Core::Relationship;

use Getopt::Long;

my %opts = ();
GetOptions (\%opts,
	'i=s{1,1}',
	'f=s{1,1}',
	'o=s{1,1}',
	'help|h')
or die("Error in command line arguments, ask for help: -h\n");

my $file    = $opts{f};
my $initFile = $opts{i};
my $outputFile = $opts{o};

unless ($file and $initFile and $outputFile) {print_help()};

sub print_help {
	print "\n";
	print "\tdescription: Get all ancesters defined in init file.\n";
	print "\tusage      : getOntologyFile.pl [options]\n";
	print "\toptions    :\n";
	print "\t\t-f  	 OBO input file\n";
	print "\t\t-i 	 init file\n";
	print "\t\t-o 	 output file\n";
	print "\texample:\n";
	#print "\t\tperl get_parent_terms.pl -f go.obo -t GO:0000234\n";
	exit;
}

sub parseInitFile
{
	my($pathFile) = @_;
	open my $FILE, "<", $pathFile or die $!;
	my @termArray;
	while(my $line = <$FILE>)
	{
		my @lineTerms;
		if($line =~ /^\s*([^\s]+)/g)
		{
			push(@lineTerms, $1);
			my %terms;
			while($line =~ /(\w{4}:\d+)/g)
			{
				$terms{$1}++;
			}
			push(@lineTerms, keys %terms);
		}
		push(@termArray, \@lineTerms);
	}
	close($FILE);
	return @termArray;
}

my @allInitTerms = &parseInitFile($initFile);
my $my_parser = OBO::Parser::OBOParser->new();
my $ontology = $my_parser->work($file);

open my $OUTPUT, ">", $outputFile or die $!;
foreach my $termsLine (@allInitTerms)
{
	my @termArray;
	my $originalTerms = " ";
	for(my $i = 1; $i < scalar @$termsLine; ++$i)
	{
		push(@termArray, @{$ontology->get_ancestor_terms($ontology->get_term_by_id($$termsLine[$i]))});
		$originalTerms .= $$termsLine[$i] ." ";
	}
	print "getting " . scalar  @termArray ." from " .$$termsLine[0] ." ... \n";
	print $OUTPUT $$termsLine[0];
	my $tmpLine = $originalTerms;
	foreach my $ontologyTerm (@termArray)
	{
		$tmpLine .= $ontologyTerm->id() . " "  if (defined $ontologyTerm->name());
	}
	chomp $tmpLine;
	print $OUTPUT $tmpLine . "\n";
}

close($OUTPUT);

#foreach my $term (@{$ontology->get_parent_terms($ontology->get_term_by_id($term_id))}) {
#	print $term->id();
#	print "\t", $term->name() if (defined $term->name());
#	print "\n";
#	my @rel = @{$ontology->get_relationships_by_target_term($term)};
#	
#	#print Dumper(@rel);
#	foreach my $x (@rel)
#	{
#		#print $x->[2]->{'TYPE'} . " \n";
#		#print OBO::Core::Relationship::type($x->[0]);
#		print "\tid:" . $x->id() . ",";
#		print "\ttype:" . $x->type(). ",";
#		#print "name:" . $x->name();
#		
#		#print Dumper($x);
#		print " ----------------\n";
#		
#	}
#	exit;
#	print "\n";
#	#print "\t", $ontology->get_relationships_by_target_term($term) if (defined $term->subnamespace());
#	#get_relationships_by_target_term
#	#get_relationships_by_source_term
#	#get_descendent_terms_by_relationship_type
#	print "\n";
#
#
exit 0;
