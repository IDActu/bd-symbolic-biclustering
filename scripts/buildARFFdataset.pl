#!/usr/bin/perl
########################################################################
# Script transforms data from an original binary matrix into ARFF file.
# Frantisek Malinka
########################################################################

#load modules
use strict;
use warnings;
use OBO::Parser::OBOParser;
use OBO::Core::Relationship;

use IO::Handle;
STDOUT->autoflush(1);

#arguments
my $pathTrain_input = "";
my $pathTest_input = "";
my $min_cover = "";
my $pathTrain = "train.arff";
my $pathTest = "test.arff";
my $LOCFLAG = 0;
my $locPath = "";
my $geneFLAG = 0;
my $genePath = "";
my $rowOBOPath = "";
my $rowOBOFLAG = 0;
my $fbgn2GOFLAG = 0;
my $fbgn2GOPath = "";

#Check program parameters
for(my $iarg = 0; $iarg <= $#ARGV+1; ++$iarg)
{
	if(defined $ARGV[$iarg] and defined $ARGV[$iarg+1] and $ARGV[$iarg] eq "--pathTrain")
	{
		$pathTrain_input = $ARGV[$iarg+1];
	}
	if(defined $ARGV[$iarg] and defined $ARGV[$iarg+1] and $ARGV[$iarg] eq "--pathTest")
	{
		$pathTest_input = $ARGV[$iarg+1];
	}
	if(defined $ARGV[$iarg] and defined $ARGV[$iarg+1] and $ARGV[$iarg] eq "--min_cover")
	{
		$min_cover = $ARGV[$iarg+1];
	}
	if(defined $ARGV[$iarg] and defined $ARGV[$iarg+1] and $ARGV[$iarg] eq "--outputTrain")
	{
		$pathTrain = $ARGV[$iarg+1];
	}
	if(defined $ARGV[$iarg] and defined $ARGV[$iarg+1] and $ARGV[$iarg] eq "--outputTest")
	{
		$pathTest = $ARGV[$iarg+1];
	}
	if(defined $ARGV[$iarg] and defined $ARGV[$iarg+1] and $ARGV[$iarg] eq "--location")
	{
		$LOCFLAG = 1;
		$locPath = $ARGV[$iarg+1];		
	}
	if(defined $ARGV[$iarg] and defined $ARGV[$iarg+1] and $ARGV[$iarg] eq "--gene")
	{
		$geneFLAG = 1;
		$genePath = $ARGV[$iarg+1];		
	}
	if(defined $ARGV[$iarg] and defined $ARGV[$iarg+1] and $ARGV[$iarg] eq "--rowOBO")
	{
		$rowOBOFLAG = 1;
		$rowOBOPath = $ARGV[$iarg+1];		
	}
	if(defined $ARGV[$iarg] and defined $ARGV[$iarg+1] and $ARGV[$iarg] eq "--fbgn2GO")
	{
		$fbgn2GOFLAG = 1;
		$fbgn2GOPath= $ARGV[$iarg+1];
	}
}

#print Help
if($pathTrain_input eq "" or $min_cover eq "" or $pathTest_input eq "" or $LOCFLAG == 0 or $geneFLAG ==  0 or $rowOBOFLAG == 0)
{
	print STDERR "All parameters have to be filled!\n";
	print "HELP - buildARFFdataset.pl\n
	Script transforms data from an original binary matrix into ARFF file.\n
	--pathTrain [file.csv]           train matrix
	--pathTest [file.csv]            test matrix
	--min_cover [int]                minimum coverage of all GO terms
	--outputTrain file [train.arff]  output path for the train dataset
	--outputTest file [test.arff]    output path for the test dataset
	--location [file]                location ontology file
	--gene [file]                    gene ontology file
	--rowOBO [file]                  row ontology file in OBO format
	--fbgn2GO [file]                 mapping FBgn to GO
	\n"; 
	exit;
}


#***********************************************************************
#parse Location file
sub parseLocFile
{
	my($FILE, $locTerms_ref, $locIDs_ref) = @_;
	my %locOntology;
	my %allTerms;
	my %isTrainLoc;
	foreach my $item (@$locIDs_ref)
	{
		$isTrainLoc{$item}++;
	}
	while(my $line = <$FILE>)
	{
		#parse line
		if($line =~ /([^\s]+)/g)
		{
			my $baseid = $1;
			my @terms;
			while($line =~ /([^\s]+)/g)
			{
				push(@terms, $1);
				if(exists $isTrainLoc{$baseid})
				{
					$allTerms{$1}++;
				}
			}
			$locOntology{$baseid} = \@terms;
		}
	}
	@$locTerms_ref = keys %allTerms;
	return %locOntology;
}


#***********************************************************************
#get IDs from a file
sub getIDs
{
	my($FILE, $fbgnIDs_ref, $locIDs_ref) = @_;
	my @matrix;
	my $locIDs = <$FILE>;
	my @loc_tmp = split(",", $locIDs);
	push(@matrix, \@loc_tmp);
	
	foreach my $lc (@loc_tmp)
	{
		if($lc =~ /([^\,"']+)/)
		{
			push(@$locIDs_ref, $1);
		}
	}
	
	while(my $fbgnID = <$FILE>)
	{
		my @tmp = split(",", $fbgnID);
		if($tmp[0] =~ /"?([^,"]+)"?/)
		{
			$tmp[0] = $1;
			push(@$fbgnIDs_ref, $1);
		}
		push(@matrix, \@tmp);
	}
	return @matrix;
}

#***********************************************************************
#get FBgn identifiers
sub getFBgn
{
	my($fbgnTerms_ref, $terms_ref, $mappingFbgn2GO_ref) = @_;
	my $iid = 1;
	if($fbgn2GOFLAG)
	{
		foreach my $term (@$fbgnTerms_ref)
		{
			foreach my $item (@{$$mappingFbgn2GO_ref{$term}})
			{
				$$terms_ref{$item}++;
			}
		}
	}
}

sub getKeysFromArray
{
	my($terms_ref) = @_;
	my @fbgnTerms;
	while(my($key,$value) = each %$terms_ref)
	{
		push(@fbgnTerms, $key);
	}
	return @fbgnTerms;
}

#***********************************************************************
#get all FBgn IDs
sub getAllFbgnIDs
{
	my($fbgnIDs_ref, $ontology_ref, $mappingFbgn2GO_ref) = @_;
	my %terms;
	my $iid = 1;
	
	my $previous = 0;
	my @fbgnTerms = @$fbgnIDs_ref;
	&getFBgn(\@fbgnTerms, \%terms, $mappingFbgn2GO_ref);
	@fbgnTerms = keys %terms;
	
	for my $id (@fbgnTerms)
	{
		for my $termid (@{$$ontology_ref->get_ancestor_terms($$ontology_ref->get_term_by_id($id))})
		{
			$terms{$termid->id()}++;
		}
	}
	@fbgnTerms = keys %terms;
	print "(Total number of terms: " . scalar @fbgnTerms . ")";

	return @fbgnTerms;
}

#***********************************************************************
#load train matrix
sub loadMatrix
{
	my($FILE) = @_;
	my @matrix;
	while(my $line = <$FILE>)
	{
		my @splitline = split(",",$line);
		push(@matrix, \@splitline);
	}
	return @matrix;
}


#***********************************************************************
#build location data in ARFF format
sub createLocDataArff
{
	my($dataMatrix_ref, $locIDs_ref, $locTerms_ref, $locationOntology_ref) = @_;
	my @locdata;
	my %ID_posInArray;
	for(my $i = 0; $i < scalar (@$locTerms_ref); ++$i)
	{
		my $loc = $$locTerms_ref[$i];
		$loc =~ /([^\s]+)/;
		my $locID = $1;
		$ID_posInArray{$locID} = $i;
	}

	foreach my $myLoc (@$locIDs_ref)
	{
		$myLoc =~ /([^\s]+)/;
		my $myLocID = $1;		
		my @tmp_loc;
		for(my $i = 0; $i < scalar (@$locTerms_ref); ++$i)
		{
			push(@tmp_loc, "'-'");
		}
		
		#LOCATION FILE ONTOLOGY
		my @list = @{$$locationOntology_ref{$myLocID}};
		foreach my $locTerm (@list)
		{
			if(exists $ID_posInArray{$locTerm})
			{
				$tmp_loc[$ID_posInArray{$locTerm}] = "'+'";
			}
		}
		push(@locdata, join(",", @tmp_loc));	
	}
	return @locdata;
}

#***********************************************************************
#build location data in ARFF format for TEST data
sub createLocDataArff_test
{
	my($dataMatrix_ref, $locIDs_ref, $locIDs_test_ref, $locationOntology_ref, $locTerms_ref) = @_;
	my @locdata;
	my %ID_posInArray;
	for(my $i = 0; $i < scalar (@$locTerms_ref); ++$i)
	{
		my $loc = $$locTerms_ref[$i];
		$loc =~ /([^\s]+)/;
		my $locID = $1;
		$ID_posInArray{$locID} = $i;
	}

	foreach my $myLoc (@$locIDs_test_ref)
	{
		$myLoc =~ /([^\s]+)/;
		my $myLocID = $1;		
		my @tmp_loc;
		for(my $i = 0; $i < scalar (@$locTerms_ref); ++$i)
		{
			push(@tmp_loc, "'-'");
		}

		my @list = @{$$locationOntology_ref{$myLocID}};
		foreach my $locTerm (@list)
		{
			if(exists $ID_posInArray{$locTerm})
			{
				$tmp_loc[$ID_posInArray{$locTerm}] = "'+'";
			}
		}
		push(@locdata, join(",", @tmp_loc));	
	}
	return @locdata;
}

#***********************************************************************
#build FBgn (row) data in ARFF format
sub createFBgnDataArff
{
	my($dataMatrix_ref, $fbgnTerms_ref, $geneOntology_ref, $ontology_ref, $mappingFbgn2GO_ref) = @_;
	my $row_number = (scalar @$dataMatrix_ref);
	
	my @fbgndata;
	for(my $irow = 1; $irow < $row_number; ++$irow)
	{
	
		my $pos = 0; my $neg = 0;
		my %termOccurrences;

		foreach my $ref (@{$$mappingFbgn2GO_ref{${$dataMatrix_ref->[$irow]}[0]}})
		{
			for my $termid (@{$$ontology_ref->get_ancestor_terms($$ontology_ref->get_term_by_id($ref))})
			{
				$termOccurrences{$termid->id()}++;
			}
		}
			
		my $tmp_line;
		foreach my $fbgnTerm (@$fbgnTerms_ref)
		{
			$fbgnTerm =~ /(GO:\d+)/;
			$fbgnTerm = $1;
			
			if(exists $termOccurrences{$fbgnTerm})
			{
				$tmp_line .= "'+',";
				++$pos;
			}
			else
			{
				$tmp_line .= "'-',";
				++$neg;
			}
		}
		push(@fbgndata, $tmp_line);
		print "[$irow\/". eval ($row_number - 1) ."] ".${$dataMatrix_ref->[$irow]}[0]." positive: $pos, negative: $neg\n";
	}
	return @fbgndata;
}

#***********************************************************************
#filter uncovered FBgn terms
sub filterfbgnterms
{
	my($fbgndata_ref, $fbgnTerms_ref, $bannedID_ref) = @_;
	my @freq;
	foreach my $line (@$fbgndata_ref)
	{
		my @items = split(",", $line);
		for(my $i = 0; $i < scalar @items; ++$i)
		{
			if($items[$i] =~ /\+/)
			{
				$freq[$i]++;
			}
		}
	}
	my $i = 0;
	for(my $count = scalar(@freq)-1; $count >= 0; --$count)
	{
		if(!defined $freq[$count] or $freq[$count] <= $min_cover)
		{
			$$bannedID_ref{$count}++;
			splice(@$fbgnTerms_ref, $count, 1);
		}
	}
}

#***********************************************************************
#parse Gene file (row)
sub parseGeneFile
{
	my($FILE, $locTerms_ref, $locIDs_ref) = @_;
	my %locOntology;
	my %allTerms;
	my %isTrainLoc;
	foreach my $item (@$locIDs_ref)
	{
		$isTrainLoc{$item}++;
	}
	while(my $line = <$FILE>)
	{
		#parse line
		if($line =~ /([^\s]+)/g)
		{
			my $baseid = $1;
			my @terms;
			while($line =~ /([^\s]+)/g)
			{
				push(@terms, $1);
				if(exists $isTrainLoc{$baseid})
				{
					$allTerms{$1}++;
				}
			}
			$locOntology{$baseid} = \@terms;
		}
	}
	@$locTerms_ref = keys %allTerms;
	return %locOntology;
}

#***********************************************************************
#create FBgn + KEGG data (row) in ARFF format
sub createFBgnKeggDataArff
{
	my($dataMatrix_ref, $keggTerms_ref,$geneOntology_ref, $fbgndata_ref) = @_;
	my $row_number = (scalar @$dataMatrix_ref);
	for(my $irow = 1; $irow < $row_number; ++$irow)
	{
		my $pos = 0; my $neg = 0;
		my %occur;
		my $tmp_line;
		my $identifier = ${$dataMatrix_ref->[$irow]}[0];
		$identifier =~ s/\"|\s//g;
		my @array = $$geneOntology_ref{$identifier};
		$occur{$_}++ for (@{$array[0]});
		
		foreach my $keggTerm (@$keggTerms_ref)
		{
			if(exists $occur{$keggTerm})
			{					
				$tmp_line .= "'+',";
				$pos++;
			}
			else
			{
				$tmp_line .= "'-',";
				$neg++;
			}				
		}
		$$fbgndata_ref[$irow-1] .= $tmp_line;
		print "KEGG[$irow\/". eval ($row_number - 1) ."] ".${$dataMatrix_ref->[$irow]}[0]." positive: $pos, negative: $neg\n";
	}
}
#***********************************************************************
#return hash of fbgn IDs and Go IDs
sub parseFbgn2GO
{
	my %mapping;
	open my $FILE, "<", $fbgn2GOPath or die $!;
	while( my $line = <$FILE>)
	{
		if($line =~ /FB\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+/)
		{
			my $fbgn = $1;
			my $go = $3;
			push(@{$mapping{$fbgn}}, $go);
		}
	}
	return %mapping;
}

#***********************************************************************
############################# MAIN ####################################

#prepare file descriptors
open my $FILE_TRAIN, "<", $pathTrain_input or die $!;
open my $FILEARFF, ">", $pathTrain or die $!;
my $headTEST = "";
open my $FILE_LOCATION, "<", $locPath or die $!;
open my $FILE_GENE, "<", $genePath or die $!;

my %mappingFbgn2GO;
if($fbgn2GOFLAG)
{
	print "Mapping FBgn to GO terms...";
	%mappingFbgn2GO = &parseFbgn2GO();
	print "DONE\n";
}

######################### TRAIN PART ######
my @locIDs, my @fbgnIDs;
my @dataMatrix = &getIDs($FILE_TRAIN, \@fbgnIDs, \@locIDs);
my %geneOntology; my @geneTerms;
my @fbgnTerms;
my @keggTerms;
print "Parsing Gene file...";
%geneOntology = &parseGeneFile($FILE_GENE, \@geneTerms, \@fbgnIDs);
print "DONE\n";
print "Parsing OBO file...";
@keggTerms = @geneTerms;
#load OBO row ontology
my $my_parser = OBO::Parser::OBOParser->new();
my $ontology = $my_parser->work($rowOBOPath);
print "DONE\n";
print "Getting universum..."; 
@fbgnTerms = &getAllFbgnIDs(\@fbgnIDs, \$ontology, \%mappingFbgn2GO);
print "DONE\n";

print "Creating row TRAIN ARFF data...\n";
my @fbgndata = &createFBgnDataArff(\@dataMatrix, \@fbgnTerms,\%geneOntology, \$ontology, \%mappingFbgn2GO);
print "DONE\n";

my %bannedID;
#filter uncovered terms
print "Filtering GO terms...";
&filterfbgnterms(\@fbgndata, \@fbgnTerms, \%bannedID);
print "DONE\n";
print "Creating row TRAIN ARFF data with KEGG...\n";
&createFBgnKeggDataArff(\@dataMatrix, \@keggTerms,\%geneOntology, \@fbgndata);
print "DONE\n";

my %locationOntology; my @locTerms;
print "Parsing location file...";
%locationOntology = &parseLocFile($FILE_LOCATION, \@locTerms, \@locIDs);
print "DONE\n";
print "Creating column TRAIN ARFF data...";
my @locdata = &createLocDataArff(\@dataMatrix, \@locIDs, \@locTerms, \%locationOntology);
print "DONE\n";

print "Creating $pathTrain...";
#print ARFF header
print $FILEARFF "\@relation 'drosophila'\n";
$headTEST .= "\@relation 'drosophila'\n";
foreach my $fbgnTerm (@fbgnTerms)
{
	print $FILEARFF "\@attribute '$fbgnTerm' {+,-}\n";
	$headTEST .= "\@attribute '$fbgnTerm' {+,-}\n";
}

foreach my $keggTerm (@keggTerms)
{
	print $FILEARFF "\@attribute '$keggTerm' {+,-}\n";
	$headTEST .= "\@attribute '$keggTerm' {+,-}\n";
}

foreach my $mylocid (@locTerms)
{
	print $FILEARFF "\@attribute '$mylocid' {+,-}\n";
	$headTEST .= "\@attribute '$mylocid' {+,-}\n";
}
print $FILEARFF "\@attribute 'classification' {+,-}\n\@data\n";
$headTEST .= "\@attribute 'classification' {+,-}\n\@data\n";

#PRINT  TRAIN
my $dataARFF;
my $datorow = 1;
foreach my $myfbgn (@fbgndata)
{
	my @tmp_fbgn = split(",", $myfbgn);
	my $new_fbgnvector = "";
	for(my $iban = 0; $iban < scalar @tmp_fbgn; ++$iban)
	{
		if(!exists $bannedID{$iban})
		{
			$new_fbgnvector .= $tmp_fbgn[$iban].",";
		}
	}
	my $datocol = 1;
	foreach my $myloc (@locdata)
	{
		print $FILEARFF "$new_fbgnvector$myloc,";
		if($dataMatrix[$datorow][$datocol] > 0)
		{
			print $FILEARFF  "'+'\n";
		}
		else
		{
			print $FILEARFF  "'-'\n";
		}
		++$datocol;
	}
	++$datorow;
}
print "DONE\n";

######################### TEST PART ######
#prepare test file descriptors
open my $FILE_TEST, "<", $pathTest_input or die $!;
open my $FILEARFF_TEST, ">", $pathTest or die $!;

print $FILEARFF_TEST $headTEST;

my @locIDs_test, my @fbgnIDs_test;
my @dataMatrix_test = &getIDs($FILE_TEST, \@fbgnIDs_test, \@locIDs_test);
print "Creating row TEST ARFF data...\n";
my @fbgndata_test = &createFBgnDataArff(\@dataMatrix_test, \@fbgnTerms, \%geneOntology, \$ontology, \%mappingFbgn2GO);
print "DONE\n";

print "Creating row TEST ARFF data with KEGG...\n";
&createFBgnKeggDataArff(\@dataMatrix_test, \@keggTerms,\%geneOntology, \@fbgndata_test);
print "DONE\n";

print "Creating column TEST ARFF data...";
my @locdata_test = &createLocDataArff_test(\@dataMatrix_test, \@locIDs, \@locIDs_test, \%locationOntology, \@locTerms);
print "DONE\n";

print "Creating $pathTest...";
$datorow = 1;
foreach my $myfbgn (@fbgndata_test)
{
	my $datocol = 1;
	foreach my $myloc (@locdata_test)
	{
		if(!($dataMatrix_test[$datorow][$datocol] =~ /NA/))
		{
			print $FILEARFF_TEST "$myfbgn$myloc,";
			if($dataMatrix_test[$datorow][$datocol] > 0)
			{
				print $FILEARFF_TEST  "'+'\n";
			}
			else
			{
				print $FILEARFF_TEST  "'-'\n";
			}
		}
		++$datocol;
	}
	++$datorow;
}
print "DONE\n";
close($FILE_TEST);
close($FILEARFF_TEST);

print "Elapsed time: ";
print time() - $^T;
print "sec\n";

close($FILE_TRAIN);
close($FILEARFF);
