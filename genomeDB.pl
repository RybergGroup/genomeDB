#! /usr/bin/perl -w

use strict;
use FindBin ();
use lib "$FindBin::Bin";
use DBI;
#use Bio::DB::Sam;
use LWP::Simple;
use Text::Wrap;
$Text::Wrap::columns = 80;
use sequences;

####################################################################
### Global variables that may or may not be changed by arguments ###
####################################################################

my $database_name;
my $taxon;
my $note;
my %things_to_do;
my $score_cut_off = 35;
my $start_time = time;
my $last_time = $start_time;
my $extraSQLcondition='';
my $SEQOUT = *STDOUT;
my $genetic_code = 'standard';
my $include_mate_in_other_group = 'N';
my $BLASTtype = 'taxonomy';
my $protein_query = 'F';
my %categories = ( 'fungi' => 'Eukaryota; Fungi; Dikarya; Basidiomycota; Agaricomycotina', 'bact' => 'Bacteria' ); # Give your organisms annotations as keys
my $ambig_critearia = 3; # The criteria to decide if the classification is ambiguous
                         # The function that use the criteria is evaluate_ambiguity
my $overwrite = 'CF';
my $delimiter = '\|';
my $position = 1;
my $BLAST_GBdatabase = 'protein';
my $e_cut_off;
my $bam_for_annotations = 'F';
my $SNP_for_annotations = 'F';
###############################################
### Parsing arguments                       ###
###############################################

for (my $i=0; $i < scalar @ARGV; ++$i) {
    if ($ARGV[$i] eq '-t' || $ARGV[$i] eq '--taxon') {
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
	    $taxon = $ARGV[++$i];
	}
	else { die "-t/--taxon requere a taxon name as next argument.\n"; }
    }
    elsif ($ARGV[$i] eq '-N' || $ARGV[$i] eq '--note') {
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
	    $note = $ARGV[++$i];
	}
	else { die "-N/--note requere a note annotatione as next argument.\n"; }
    }
    elsif ($ARGV[$i] eq '-db' || $ARGV[$i] eq '--database_name') {
 	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
            $database_name = $ARGV[++$i];
	}
        else { die "-db/--database_name requere a database name as next argument.\n"; }
    }
    elsif ($ARGV[$i] eq '-C' || $ARGV[$i] eq '--create_database') {
	$things_to_do{"Create database"}='T';
    }
    elsif ($ARGV[$i] eq '-s' || $ARGV[$i] eq '--read_scaffolds') {
        $things_to_do{"read_scaffolds"}='T';
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
	    $things_to_do{"read_scaffolds"}=$ARGV[++$i];
	}
    }
    elsif ($ARGV[$i] eq '-a' || $ARGV[$i] eq '--read_annotations') {
	$things_to_do{"read_annotations"}='T';
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
            $things_to_do{"read_annotations"}=$ARGV[++$i];
        }
    }
    elsif ($ARGV[$i] eq '-T' || $ARGV[$i] eq '--get_reads_of_taxon') {
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
            $things_to_do{"get_reads_of_taxon"}=$ARGV[++$i];
	}
	else { die "-T/--get_reads_of_taxon require a taxon name and bam file as second respectively third argument.\n"; }
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
            $things_to_do{"get_reads_of_taxon"}.= " ;###; " . $ARGV[++$i];
        }
    }
    elsif ($ARGV[$i] eq '-R' || $ARGV[$i] eq '--get_reads_of_scaffold') {
        if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
            $things_to_do{"get_reads_of_scaffold"}=$ARGV[++$i];
        }
	else { die "-R/--get_reads_of_scaffold require a scaffold/contig name and bam file as second respectively third argument.\n"; }
        if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
            $things_to_do{"get_reads_of_scaffold"}.= " ;###; " . $ARGV[++$i];
        }
    }
    elsif ($ARGV[$i] eq '-b' || $ARGV[$i] eq '--read_bam') {
        $things_to_do{"read_bam"}='T';
        if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
            $things_to_do{"read_bam"}=$ARGV[++$i];
        }
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
	    ++$i;
	    if ($ARGV[$i] =~ /anno/i) { $bam_for_annotations = 'T'; }
	}
    }
    elsif ($ARGV[$i] eq '-v' || $ARGV[$i] eq '--read_vcf') {
        $things_to_do{"readVCF"}='T';
        if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
            $things_to_do{"readVCF"}=$ARGV[++$i];
        }
        if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
            ++$i;
            if ($ARGV[$i] =~ /anno/i) { $SNP_for_annotations = 'T'; }
        }
    }
    elsif ($ARGV[$i] eq '--N50') {
	$things_to_do{"N50"} = 'all';
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
	    $things_to_do{"N50"} = $ARGV[++$i];
	}
    }
    elsif ($ARGV[$i] eq '--set_genetic_code') {
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
	    $genetic_code = $ARGV[++$i];
	    if (!&sequences::gene_codes_available($genetic_code)) { die "'$genetic_code' is not an available genetic code. Available options are: ", &sequences::available_gene_codes_string(), "\n"; }
	    
	}
	else { die "The argument --set_genetic_code need the name of the code to use to translate triplets, e.g. standard, as next argument.\n"; }
    }
    elsif ($ARGV[$i] eq '--include_missmatch_mates') {
	$include_mate_in_other_group = 'Y';
    }
    elsif ($ARGV[$i] eq '-S' || $ARGV[$i] eq '--SNPscore_cutoff') {
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
	    $score_cut_off = $ARGV[++$i];
	}
	else { die "-S/--SNPscore_cutoff need a integer value as next argument.\n"; }
    }
    elsif ($ARGV[$i] eq '--get_CDS') {
	$things_to_do{"pars_annotations"} = 'CDS';
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
	    if ($extraSQLcondition) { $extraSQLcondition .= ";"; }
	    $extraSQLcondition .= $ARGV[++$i];
	}
    }
    elsif ($ARGV[$i] eq '--get_proteins') {
	$things_to_do{"pars_annotations"} = 'protein';
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
	    if ($extraSQLcondition) { $extraSQLcondition .= ";"; }
	    $extraSQLcondition .= $ARGV[++$i];
	}
    }
    elsif ($ARGV[$i] eq '-B' || $ARGV[$i] eq '--BLASTtaxon_annotation') {
        $things_to_do{"pars_blast"}='T';
        if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
            $things_to_do{"pars_blast"}=$ARGV[++$i];
        }
	if (defined($ARGV[$i+1]) && !($ARGV[$i+1] =~ /^-/)) {
	    ++$i;
	    if ($ARGV[$i]=~ /check/i) { $overwrite = 'C'; }
	    else { $overwrite = ''; }
	    if ($ARGV[$i]=~ /overwrite/i) { $overwrite .= 'T'; }
	    else { $overwrite .= 'F'; }
	    if ( $ARGV[$i]=~ /query:gene/i ) { $protein_query = 'T'; }
	    else { $protein_query = 'F'; }
	    if ( $ARGV[$i] =~ /GB:\s*('|")(.+)('|")\s*/i ) {
		$BLAST_GBdatabase = $2;
		if (!$BLAST_GBdatabase =~ /file:/i) { $BLAST_GBdatabase = lc($BLAST_GBdatabase); }
	    }
	    if ( $ARGV[$i] =~ /annotation:acc/i ) { $BLASTtype = 'accno'; }
	    if ( $ARGV[$i] =~ /ecutoff:([0-9e.-]+)/i || $ARGV[$i] =~ /e-cutoff:([0-9e.-]+)/i || $ARGV[$i] =~ /e_cut_off:([0-9e.-]+)/i || $ARGV[$i] =~ /e-cut_off:([0-9e.-]+)/i) { $e_cut_off = $1; }
	    if ( $ARGV[$i] =~ /delimiter:(.)/i) {
		$delimiter = $1;
		if ($delimiter eq '|' || $delimiter eq '.' || $delimiter eq '(' || $delimiter eq ')' || $delimiter eq '[' || $delimiter eq ']' || $delimiter eq '{' || $delimiter eq '}' || $delimiter eq '\\'){
		    $delimiter = "\\$delimiter";
		}
	    }
	    if ( $ARGV[$i] =~ /position:([0-9]+)/ ) { $position = $1; }
	    if ( $ARGV[$i] =~ /taxa:\s*('.+'|".+"|[^;]+)/i ) {
		undef (%categories);
		my $temp = $1;
		$temp =~ s/^\s+//;
		$temp =~ s/\s+$//;
		$temp =~ s/^('|")//;
		$temp =~ s/('|")$//;
		#print $temp, "\n";
		if ($temp =~ /file:\s*([^;]+)/i) {
		    open (TAXA, '<', $1) || die "Could not open $1: $!.\n";
		    while (<TAXA>) {
			if (/^#/) { next; }
			chomp;
			my @temp = split /\t/;
			$temp[0] =~ s/^\s+//;
			$temp[0] =~ s/\s+$//;
			$temp[1] =~ s/^\s+//;
			$temp[1] =~ s/\s+$//;
			$categories{$temp[0]} = $temp[1];
		    }
		    close TAXA;
		}
		elsif ($temp =~ /,/) {
		    %categories = split /\s*,\s*/, $temp;
		}
		else { die "Unrecognized format for taxa in: $ARGV[$i].\n"; }
		if (!%categories) { die "Failed to pars taxonomic categories.\n"; }
	    }
	}
    }
    elsif ($ARGV[$i] eq '-h' || $ARGV[$i] eq '--help') {
	&help();
	exit;
    }
    elsif ($i+1 == scalar @ARGV && (($things_to_do{"read_scaffolds"} && $things_to_do{"read_scaffolds"} eq 'T') || ($things_to_do{"read_annotations"} && $things_to_do{"read_annotations"} eq 'T') || ($things_to_do{"read_bam"} && $things_to_do{"read_bam"} eq 'T') || ($things_to_do{"pars_blast"} && $things_to_do{"pars_blast"} eq 'T')) && !($ARGV[$i] =~ /^-/)) {
	if ($things_to_do{"read_scaffolds"} && $things_to_do{"read_scaffolds"} eq 'T') {
	    $things_to_do{"read_scaffolds"} = $ARGV[$i];
	}
	elsif ($things_to_do{"read_annotations"} && $things_to_do{"read_annotations"} eq 'T') {
	    $things_to_do{"read_annotations"} = $ARGV[$i];
	}
	elsif ($things_to_do{"read_bam"} && $things_to_do{"read_bam"} eq 'T') {
	    $things_to_do{"read_bam"} = $ARGV[$i];
	}
	elsif ($things_to_do{"pars_blast"} && $things_to_do{"pars_blast"} eq 'T') {
	    $things_to_do{"pars_blast"} = $ARGV[$i];
	}
    }
    else {
	die "Unrecognized argument: $ARGV[$i].\n";
    }
}
##############################################
### Open database                          ###
##############################################
my $dbh;
if ($database_name) {
    $dbh = DBI->connect("dbi:SQLite:dbname=$database_name",'','');
}
else { die "Database name required (-db//--database_name).\n"; }
#############################################
### Things that require an open database  ###
#############################################
if ($dbh) {
    if ($things_to_do{"Create database"}) {
	$last_time = time;
	if ($dbh->do("CREATE TABLE scaffolds (name TEXT PRIMARY KEY, sequence TEXT DEFAULT 'empty', taxon TEXT DEFAULT 'empty', coverage_median INTEGER, coverage_mean REAL, coverage_sd REAL, coverage_max INTEGER, coverage_min INTEGER,  GCcontent REAL, medianSNPdiversity REAL, nSNP INTEGER, npolySNP INTEGER, percdiSNP40_60 REAL, note TEXT DEFAULT 'empty')")) {
	    print "Created table scaffolds (",time-$last_time,"s).\n";
	    $last_time = time;
	}
	else { die "Could not create table scaffolds in database $database_name.\n"; }
	if ($dbh->do("CREATE TABLE annotations (seqid TEXT, source TEXT, type TEXT, start INTEGER, end INTEGER, score REAL, strand TEXT, phase INTEGER, id TEXT, name TEXT, alias TEXT, parent TEXT, target TEXT, gap TEXT, derives_from TEXT, note TEXT, dbxref TEXT, ontology_term TEXT, is_circular TEXT, extras TEXT, coverage_median INTEGER, nSNP INTEGER, GCcontent REAL, PRIMARY KEY (seqid,source,type,start,end,strand,id) FOREIGN KEY(seqid) REFERENCES scaffolds(name))")) {
	    print "Created table annotations (",time-$last_time,"s).\n";
	    $last_time = time;
	}
	else { die "Could not create table annotations in database $database_name.\n"; }
    }
###############################################################################
### Alternativs to add data to database                                     ###
###############################################################################
    if ($things_to_do{"read_scaffolds"}) {
	if ($things_to_do{"read_scaffolds"} ne 'T') {
	    if (!$taxon) { $taxon = "'empty'"; }
	    else { $taxon = "'$taxon'"; }
	    if (!$note) { $note ="'empty'"; }
	    else { $note = "'$note'"; }
	    open my $INPUT, "<$things_to_do{'read_scaffolds'}" or die "Could not open $things_to_do{'read_scaffolds'}: $!.\n";
	    my ($name, $seq);
	    my $added_seq=0;
	    print "Reading scaffolds...\n";
	    my $begin_sth = $dbh->prepare("BEGIN");
	    $begin_sth->execute();
	    my $commit_sth = $dbh->prepare("COMMIT");
	    $last_time = time;
	    while (my $row = <$INPUT>) {
		chomp $row;
		if ($row =~ s/^>//) {
		    if ($seq && $name) {
			$seq = "'" . $seq . "'";
			if (&add_to_scaffolds($dbh,\$name,\$seq,\$taxon,\$note)) {
			    ++$added_seq;
			    if (!($added_seq % 1000)) {
				$commit_sth->execute();
				$begin_sth->execute();
				print "Added $added_seq scaffolds (",time-$last_time,"s).\n";
				$last_time = time;
			    }
			}
			else { print STDERR "Could not add $name to scaffolds in $database_name.\n"; }
		    }
		    $name = "'$row'";
		    undef($seq);
		}
		else {
		    $seq .= $row;
		}
	    }
	    if ($seq && $name) {
		$seq = "'" . $seq . "'";
		if (&add_to_scaffolds($dbh,\$name,\$seq,\$taxon,\$note)) {
			++$added_seq;
		}
		else { print STDERR "Could not add $name to scaffolds in $database_name.\n"; }
	    }
	    $commit_sth->execute();
	    $begin_sth->finish();
	    $commit_sth->finish();
	    print "Added $added_seq scaffolds (",time-$last_time,"s).\n";
	    $last_time = time;
	}
	else { die "Need a name of the sequence file.\n"; }
    }
##################################################################################################################################
    if ($things_to_do{"read_annotations"}) {
     	if ($things_to_do{"read_annotations"} ne 'T') {
	    open my $INPUT, "<$things_to_do{'read_annotations'}" or die "Could not open $things_to_do{'read_annotations'}: $!.\n";
	    my $added_anno=0;
	    print "Reading annotations (expecting GFF3 file).\n";
	    my $sth = $dbh->prepare("INSERT INTO annotations VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
	    my @other_columns = ('id', 'name', 'alias', 'parent', 'target', 'gap', 'derives_from', 'note', 'dbxref', 'ontology_term', "is_circular", "extras");
	    # coverage_median INTEGER, nSNP INTEGER, GCcontent
	    my $seq_sth = $dbh->prepare("SELECT sequence FROM scaffolds WHERE name=?");
	    my $begin_sth = $dbh->prepare("BEGIN"); # for faster processing commit will be done explicitly
	    $begin_sth->execute();
	    my $commit_sth = $dbh->prepare("COMMIT");
	    $last_time = time;
	    my $prev_scaffold;
	    my $sequence;
	    my $seq_length;
	    while (my $row = <$INPUT>) {
		chomp $row;
		if (!($row=~ /^#/)) { # if not a comment
		    my @values = split /\t/, $row; # get the different values
		    if (scalar @values != 9) {      # if not nine columns
			die "Format error: $row\n"; # it is not the right format
		    }
		    if ($values[5] eq '.') { $values[5] = 'empty'; } # if no score set score to null
		    my @temp = split /=|;/, pop @values; # get attributes from last column
		    my %attributes;
		    while (my $key = shift @temp) { $attributes{lc($key)} = shift @temp; } # type value pair
		    foreach my $col (@other_columns) {
			if (defined($attributes{$col})) { $values[scalar @values] = delete $attributes{$col}; }
			else {  push @values, "empty"; }
		    }
		    foreach my $col (keys %attributes) {
			if ($values[$#values] eq 'empty') { $values[$#values] = "$col=$attributes{$col}"; }
			else { $values[$#values] .= ";$col=$attributes{$col}"; }
		    }
		    foreach (@values) { $dbh->quote($_); }
		    push @values, -1; # add zero for coverage
		    push @values, -1; # for nSNP
		    if (!$prev_scaffold || $values[0] ne $prev_scaffold) {
			$seq_sth->execute($values[0]);
			$sequence = $seq_sth->fetchrow_array();
			$prev_scaffold = $values[0];
			$seq_length = length($sequence);
		    }
		    if ($sequence && $values[3] <= $seq_length) {
			#my $anno_seq = '';
			my $start = 0;
			my $length = 0;
			if ($values[3] < $values[4]) { $start = $values[3]-1; $length = $values[4]-$values[3]+1; }
			else { $start = $values[4]-1; $length = $values[3]-$values[4]+1; }
			#my $anno_seq = substr($sequence, $start, $length);
			push @values, &sequences::GCcontent(\substr($sequence, $start, $length));
			#if (@values[$#values] == 0) { print ("$start # $length\n"); }
		    }
		    else { push @values, -1.0; }
		    if ($sth->execute(@values)) {
			++$added_anno;
			if (!($added_anno % 5000)) {
			    $commit_sth->execute();
			    $begin_sth->execute();
			    print "Added $added_anno annotations (",time-$last_time,"s).\n";
			    $last_time = time;
			}
		    }
		    else {
			print STDERR "Could not add: $row\n";
		    }
		}
		elsif ($row =~ /##FASTA/) { last; }
	    }
	    $commit_sth->execute();
	    print "Added $added_anno annotations (",time-$last_time,"s).\n";
	    $last_time = time;
	    #$dbh->do("ALTER TABLE annotations ENABLE KEYS");
	    $begin_sth->finish();
	    $commit_sth->finish();
	    $sth->finish();
	}
	else { die "Need a name of the annotation file.\n"; }
    }
#################################################################################################################################
    if ($things_to_do{'read_bam'}) {
	$score_cut_off += 33; #  to account for Sanger format Phred quality score
	if ($things_to_do{'read_bam'} ne 'T') {
	    #if (my $sam = Bio::DB::Sam->new( -bam  => $things_to_do{'read_bam'})) {
	    if (-e $things_to_do{'read_bam'} && -e "$things_to_do{'read_bam'}.bai") {
		#print "Reading bam file...\n";
		my $sth = $dbh->prepare("SELECT name FROM scaffolds");
		my $sth_insert = $dbh->prepare("UPDATE scaffolds SET coverage_median=?, coverage_mean=?, coverage_sd=?, coverage_max=?, coverage_min=?, medianSNPdiversity=?, nSNP=?, npolySNP=?, percdiSNP40_60=? WHERE name=?");
		my $sth_annotations;
		my $sth_insert_annotations;
		if ($bam_for_annotations eq 'T') {
		    $sth_annotations = $dbh->prepare("SELECT start,end,id FROM annotations WHERE seqid=?");
		    $sth_insert_annotations = $dbh->prepare("UPDATE annotations SET coverage_median=?, nSNP=?  WHERE seqid=? AND id=?");
		}
		my $begin_sth = $dbh->prepare("BEGIN");
		$begin_sth->execute();
		my $commit_sth = $dbh->prepare("COMMIT");
		my $added_stats = 0;
		my $protein_annotations = 0;
		$sth->execute();
		print "Processing sequences.\n";
		$last_time = time;
		while (my $seq_name = $sth->fetchrow_array()) {
		    #my $debug_time = time;
		    #my $depth = 0;
		    my $positions = 0;
		    my @data;
		    my @SNPdiversity;
		    my $n_pollySNP = 0;
		    my $ndiSNP40_50 = 0;
		    #for (my $i=0; $i< scalar 20; ++$i) { $diSNP[$i] = 0; }
		    my @mpileup = `samtools mpileup -r $seq_name $things_to_do{'read_bam'}`;
		    #print "Samtools time: ", time - $debug_time, "\n";
		    if (@mpileup) {
			#$debug_time = time;
			while (my $row = shift @mpileup) {
			    my @columns = split /\t/, $row;
			    ++$positions;
			    my %n_bases;
			    push @data, $columns[3];
			    my $pos = 0;
			    for (my $i=0; $i < length $columns[4]; ++$i) {
				my $base = uc(substr($columns[4],$i,1));
				if ($base eq ',' || $base eq '.') { $base = uc($columns[2]); }
				if ($base eq '^') { ++$i; }
				elsif ($base eq '+' || $base eq '-') {
				    if (substr($columns[4],$i+1) =~ /^([0-9]+)/) {
					$i += $1 + length($1);
				    }
				}
				elsif ($base eq 'A' || $base eq 'T' || $base eq 'C' || $base eq 'G') {
				    if ($pos >= length $columns[5]) { print STDERR "$i $columns[4]  --  $pos $columns[5]"; }
				    if (ord(substr($columns[5],$pos,1))  > $score_cut_off) {
					++$n_bases{$base};
				    }
				    ++$pos;
				}
				elsif ($base eq '*') { ++$pos; }
			    }
			    if ($pos != length ($columns[5])-1) { print STDERR "Error in reading char $columns[0] - $columns[1] ($pos, ", length $columns[5], ".\n"; }
			    if (scalar keys %n_bases > 1) {
				my $sum = 0;
				my $freq_add = 0;
				my $tot=0;
				foreach (keys %n_bases) { $tot += $n_bases{$_}; }
				foreach my $key (keys %n_bases) {
				    my $freq = $n_bases{$key}/$tot;
				    $sum += $freq**2;
				}
				push @SNPdiversity, 1/$sum;
				if (scalar keys %n_bases > 2) {
				    ++$n_pollySNP;
				}
				else {
				    my $high = 0;
				    foreach(keys %n_bases) { if ( $n_bases{$_} > $high ) { $high = $n_bases{$_}; } }
				    if (($high/$tot) >= 0.4 && ($high/$tot) <= 0.6) { ++$ndiSNP40_50; }
				}
			    }
			    else { push @SNPdiversity, -1; }
			}
			#print "Pars time: ", time - $debug_time, "\n";
			if ($positions) {
			    if ($bam_for_annotations eq 'T') {
				$sth_annotations->execute($seq_name);
				#$debug_time = time;
				my $n=0;
				while (my ($start,$end,$id) = $sth_annotations->fetchrow_array()) {
				    my @region_data;
				    my @region_SNPs;
				    for (my $i=$start-1; $i < $end && $i < scalar @SNPdiversity && $i < scalar @data; ++$i) {
					push @region_data, $data[$i];
					push @region_SNPs, $SNPdiversity[$i];
				    }
				    my ($average,$sd) = &get_bam_stats(\@region_data, \@region_SNPs);
				    if ($sth_insert_annotations->execute($region_data[scalar @region_data/2], scalar @region_SNPs, $seq_name, $id)) { ++$n; }
				}
				$protein_annotations += $n;
				if ($n > 500) {
				    $commit_sth->execute();
				    $begin_sth->execute();
				}
				#if ($n) { print "\tAdd $n genome annotation stats in: ", time - $debug_time, "s\n"; }
			    }
    
    			    sub get_bam_stats {
	    			my $data = shift;
    				my $SNPs = shift;
    				@{$data} = sort {$a <=> $b} @{$data};
    				@{$SNPs} = sort {$a <=> $b} @{$SNPs};
    				while (@{$SNPs} && $SNPs->[0] < 0) { shift @{$SNPs}; }
    				my $average = 0;
    				foreach (@{$data}) { $average += $_; }
    				if (scalar @{$data}) { $average /= scalar @{$data}; }
    				my $sd = 0;
    				foreach (@{$data}) { $sd += ($_-$average)*($_-$average); }
    				if (scalar @{$data}) { $sd /= scalar @{$data}; }
    				$sd = sqrt($sd);
				return ($average,$sd);
	    		    }
			    #$debug_time = time;
	    		    my ($average,$sd) = &get_bam_stats(\@data, \@SNPdiversity);
			    if (scalar @SNPdiversity - $n_pollySNP > 0) { $ndiSNP40_50 = $ndiSNP40_50/(scalar @SNPdiversity - $n_pollySNP); }
			    else { $ndiSNP40_50 = 0; }
			    my $median_SNP_div = 0;
			    if (@SNPdiversity) { $median_SNP_div = $SNPdiversity[(scalar @SNPdiversity)/2]; }
	    		    if (@data && defined($average) && defined($sd) && defined($median_SNP_div) && defined($n_pollySNP) && defined($ndiSNP40_50) && $seq_name) {
				#print "Calc scaffold stats: ", time - $debug_time, "\n";
				#$debug_time = time;
				if ($sth_insert->execute($data[scalar @data/2], $average, $sd, $data[-1], $data[1], $median_SNP_div, scalar @SNPdiversity, $n_pollySNP, $ndiSNP40_50, $seq_name)) {
				    ++$added_stats;
				    if (!($added_stats % 1000)) {
					print "Added $added_stats annotations to scaffolds ";
					if ($protein_annotations) { print " and $protein_annotations to geneome annotations "; }
					print "(",time-$last_time,"s).\n";
					$last_time = time;
					$commit_sth->execute();
					$begin_sth->execute();
				    }
				}
				else {
				    print STDERR "Could not add stats to $seq_name.\n";
				}
	    		    }
	    		    else {
	    			print STDERR "Could not add stats to $seq_name: ";
				if (!@data) { print "No data." }
				elsif (!defined($average)) { print "No average."; }
				elsif (!defined($sd)) { print "No sd."; }
				elsif (!defined($median_SNP_div)) { print "No SNPs."; }
				elsif (!defined($n_pollySNP)) { print "No pollySNP."; }
				elsif (!defined($ndiSNP40_50)) { print "No SNP 40-50."; }
				else { print "Databasing failed."; }
				print "\n";
	    		    }
			    #print "Add scaffold time: ", time - $debug_time, "\n";
	    		}
		    }
		    else { print STDERR "No pileup using: samtools mpileup -r $seq_name $things_to_do{'read_bam'}.\n"; }
		}
		$commit_sth->execute();
		$commit_sth->finish();
		$begin_sth->finish();
		if ($bam_for_annotations eq 'T') { $sth_annotations->finish(); $sth_insert_annotations->finish(); }
		$sth_insert->finish();
		$sth->finish();
		print "Added stats to $added_stats scaffolds (",time-$last_time,"s).\n";
		$last_time = time;
	    }
    	    else { die "Could not find $things_to_do{'read_bam'} and/or $things_to_do{'read_bam'}.bai.\n"; }
	}
	else { die "Need a name of the bam file.\n"; }
    }
########################################################################################
    if ($things_to_do{'readVCF'}) {
	print "Will read and summarize SNPs from $things_to_do{'readVCF'}.\n";
	if (-e $things_to_do{'readVCF'}) {
	    open VCF, '<', $things_to_do{'readVCF'} or die "Could not open $things_to_do{'readVCF'}.\n";
	    my $sth_insert = $dbh->prepare("UPDATE scaffolds SET nSNP=?, npolySNP=?, percdiSNP40_60=?, medianSNPdiversity=? WHERE name=?");
	    my $sth_annotations;
    	    my $sth_insert_annotations;
	    if ($SNP_for_annotations eq 'T') {
		$sth_annotations = $dbh->prepare("SELECT start,end,id FROM annotations WHERE seqid=?");
		$sth_insert_annotations = $dbh->prepare("UPDATE annotations SET nSNP=? WHERE seqid=? AND id=?");
	    }
	    my $begin_sth = $dbh->prepare("BEGIN");
	    $begin_sth->execute();
	    my $commit_sth = $dbh->prepare("COMMIT");
	    my $scaffold;
	    my $n_scaffolds=0;
	    my @SNPs;
	    my $added_stats = 0;
	    my $protein_annotations = 0;
	    print "Processing SNPs.\n";
	    $last_time = time;
##############
##############
		sub handleSNPs {
		    if (@SNPs && $SNP_for_annotations eq 'T') {
			$sth_annotations->execute($scaffold);
			my $n=0;
			while (my ($start,$end,$id) = $sth_annotations->fetchrow_array()) {
			    my $nSNPs = 0;
			    foreach (@SNPs) { if ($_->[0] >= $start && $_->[0] <= $end) { ++$nSNPs } }
			    if ($sth_insert_annotations->execute($nSNPs, $scaffold, $id)) { ++$n; }
			}
			$protein_annotations += $n;
			if ($n > 500) {
			    $commit_sth->execute();
			    $begin_sth->execute();
			}
		    }
		    if (@SNPs) {
			my $n_pollySNP = 0;
			my $ndiSNP40_50 = 0;
		       	my @SNP_div;
			foreach (@SNPs) {
			    if ($_->[1] && scalar @{$_->[1]} > 1) {
				++$n_pollySNP;
				my $sum=0;
				foreach (@{$_->[1]}) { $sum += $_**2; }
				push @SNP_div, 1/$sum;
			    }
			    else {
				if ($_->[1] && $_->[1]->[0] > 0.4 && $_->[1]->[0] < 0.6) { ++$ndiSNP40_50; }
			    }
			}
			my $median_SNP_div = 0;
			if (@SNP_div) {
			    @SNP_div = sort @SNP_div;
			    if (scalar @SNP_div < 2) {
				$median_SNP_div = $SNP_div[0];
			    }
			    elsif (scalar @SNP_div % 2) {
				$median_SNP_div = $SNP_div[scalar(@SNP_div)/2];
			    }
			    else { $median_SNP_div = ($SNP_div[(scalar(@SNP_div)/2)-1] + $SNP_div[scalar(@SNP_div)/2]) / 2; }
			}
			if ($sth_insert->execute( scalar @SNPs, $n_pollySNP, $ndiSNP40_50/scalar @SNPs, $median_SNP_div, $scaffold)) {
			    ++$added_stats;
			    if (!($added_stats % 500)) {
				print "Added $added_stats annotations to scaffolds ";
				if ($protein_annotations) { print " and $protein_annotations to geneome annotations "; }
				print "(",time-$last_time,"s).\n";
				$last_time = time;
				$commit_sth->execute();
				$begin_sth->execute();
			    }
			}
			else {
			    print STDERR "Could not add stats for $scaffold.\n";
			}
		    }
		    else { print "No SNPs for $scaffold.\n"; }
		}
####################################
####################################
	    while (my $row = <VCF>) {
		if ($row =~ /^#/) { next; }
		my @columns = split /\t/, $row;
		if (defined $scaffold && $columns[0] ne $scaffold) {
		    &handleSNPs();
		    undef @SNPs;
		    $scaffold = $columns[0];
		}
		elsif (!defined($scaffold)) { $scaffold = $columns[0]; }
		push @SNPs,[];
		$SNPs[$#SNPs]->[0] = $columns[1];
		$SNPs[$#SNPs]->[1] = [];
		my %info = split /=|;/, $columns[7];
		if (defined $info{AF}) {
		    @{$SNPs[$#SNPs]->[1]} = split /,/,$info{AF};
		}
		elsif (defined $info{FR}) {
                    @{$SNPs[$#SNPs]->[1]} = split /,/,$info{FR};
                }
		else { print STDERR "No alternative base frequencies for $columns[0], pos $SNPs[$#SNPs]->[0].\n"; }
##############################################
	    }
	    if (defined $scaffold) { &handleSNPs(); }
	    $commit_sth->execute();
    	    $commit_sth->finish();
	    $begin_sth->finish();
	    if ($bam_for_annotations eq 'T') { $sth_annotations->finish(); $sth_insert_annotations->finish(); }
	    $sth_insert->finish();
	    #$sth->finish();
	    print "Added stats to $added_stats scaffolds (",time-$last_time,"s).\n";
	    $last_time = time;
	}
	else { die "Could not find VCF file name.\n"; }
############
############
    }
######################################################################################
    if ($things_to_do{'pars_blast'}) {
	if ($things_to_do{'pars_blast'} ne 'T') {
	    print "Will read BLAST results from $things_to_do{'pars_blast'} to uppdate taxon annotation of scaffolds.\n";
	    if ($BLASTtype eq 'taxonomy') {
		print " Assuming first word of each match is a taxon name.";
		foreach(keys %categories) { print " $_ will be interpreted as $categories{$_}."; }
	    }
	    elsif ($BLASTtype eq 'accno') { print " Assuming accession number can be found in position $position based on $delimiter, and corresponding annotation can be found in $BLAST_GBdatabase."; }
	    if ($protein_query eq 'T') { print " Will look for scaffold name (seqid) by searching id in annotations table."; }
	    else { print " Will use query name directtly as scaffolds name."; }
	    if ($overwrite =~ /^C/) { print " Will compare new with old taxon annotation and print if they are different."; }
	    if ($overwrite =~ /T$/) { print " Will overwrite old annotations."; }
	    elsif ($overwrite =~ /F$/) { print " Will not overwrite old annotations."; }
	    print "\n";
	    my $query;
	    my @hits;
	    my $added_taxa=0;
	    my $changed_annotations=0;
	    my $taxa_without_hit=0;
	    my $taxa_with_unsignificant_hits=0;
	    open my $BLAST, "<$things_to_do{'pars_blast'}" or die "Could not open $things_to_do{'pars_blast'}: $!.\n";
	    my $sth = $dbh->prepare("UPDATE scaffolds SET taxon=? WHERE name=?");
	    my $sth_protein;
	    if ($protein_query eq 'T') { $sth_protein = $dbh->prepare("SELECT seqid FROM annotations WHERE id=?"); }
	    my $sth_check = $dbh->prepare("SELECT taxon FROM scaffolds WHERE name==?");
	    my $begin_sth = $dbh->prepare("BEGIN");
    	    $begin_sth->execute();
	    my $commit_sth = $dbh->prepare("COMMIT");
	    my %index;
	    if (($BLAST_GBdatabase =~ /^file:download(:)?/ || $BLAST_GBdatabase =~ /^file:complement:[^\:]+/) && $BLASTtype eq 'accno') {
		my $NCBIdatabase = 'protein';
		my $filename;
		if ($BLAST_GBdatabase =~ /^file:download((:)(\w+))?$/) { if ($3) { $NCBIdatabase = $3; } }
		elsif ($BLAST_GBdatabase =~ /^file:complement:([^\:]+)((:)(\w+))?$/) {
		    if ($1) { $filename = $1; }
		    else { die "Could not pars file name from '$BLAST_GBdatabase'.\n"; }
		    if ($4) { $NCBIdatabase = $4; }
		}
		else { die "Do not recognize option '$BLAST_GBdatabase'.\n"; }
		print STDERR "Will use NCBI database $NCBIdatabase.\n";
		my %accnos;
		while (<$BLAST>) {
		    if (/Sequences producing significant alignments:/) { # Reached list of hits
			while(<$BLAST>) {
			    if (/^>/ || /^ALIGNMENTS/) { last; } # hit first alignment
			    elsif (/^\s*(.+)/) {
				my $temp = $1;
				my @annotations = split /$delimiter/, $temp;
				if ($annotations[$position]) {
				    $annotations[$position] =~ s/\.[0-9]+$//;
				    ++$accnos{$annotations[$position]};
				    #print $annotations[$position], "\n";
				}
			    }
			}
		    }
		}
		close $BLAST;
		open $BLAST, "<$things_to_do{'pars_blast'}" or die "Could not re-open $things_to_do{'pars_blast'}: $!.\n";
		print STDERR "Found ", scalar keys %accnos, " accession numbers. Will try to download GenBank annotations from NCBI database $NCBIdatabase.\n";
		if ($filename) {
		    print STDERR "Checking what records are already present in $filename.\n";
		    open INPUT, '<', $filename or die "Could not open $filename to see what accnos are present: $!.\n";
		    while (my $row = <INPUT>) {
			if ($row =~ /^ACCESSION\s+(\w+)/) { delete $accnos{$1}; }
		    }
		    print STDERR "\t$filename needs to be complemented by ", scalar keys %accnos, " records.\n";
		}
		my $i=1;
		my $query;
		my $GBFILE;
		if (!$filename) {
		    $filename = "genomeDB.genbankfile.gb";
		    while (-e $filename) {
			print STDERR "$filename exist, ";
			if ($filename =~ /gb$/) { $filename = "genomeDB.genbankfile.gb_1"; }
			elsif ($filename =~ s/gb_(\d+)$/gb_/) { my $number=$1; ++$number; $filename .= "$number"; }
			else { die "ERROR in assigning file name ($filename) to GenBank annotation file.\n"; }
			print STDERR "will try $filename.\n";
		    }
		    print STDERR "Will download GenBank annotations to $filename.\n";
		    open $GBFILE,'>', $filename || die "Could not open $filename: $!.\n";
		}
		else {open $GBFILE,'>>', $filename || die "Could not open $filename: $!.\n"; }
		$BLAST_GBdatabase = "file:$filename";
		my $no_queries = 0;
		my $time = time;
		foreach my $accno (keys %accnos) {
		    if ($query) {  $query .= ",$accno"; }
		    else { $query = $accno; }
		    if (!($i % 250) || $i == scalar keys %accnos) {
			#print "$query\n";
			my $printed = 1;
			while ($printed) {
			    #print $printed, "\n";
			    if ($printed > 10) { print STDERR "Could not download GenBank annotations for $query.\n"; last; }
			    # print "Getting entry from NCBI: ", time,"\n";
			    my $get= sequences::get_entry_from_GB($query,$NCBIdatabase);
			    ++$no_queries;
			    if ($get && $get =~ /^Resource temporarily unavailable/) { print STDERR "NCBI eutils temporarily unavailable, waiting 30 seconds...\n"; ++$printed; sleep(30);}
			    if ($get && $get =~ /^Database: .+ - is not supported/) { die "The database $NCBIdatabase is not supported by NCBI.\n"; }
			    elsif ($get) {
				$printed = 0;
				print $GBFILE $get;
				print STDERR "Processed $i accession numbers.\n";
			    }
			    else { print STDERR "No results from NCBI eutils, waiting 30 seconds...\n"; ++$printed; sleep(30);}
			}
			undef $query;
		    }
		    ++$i;
		}
		close $GBFILE;
		print "Made $no_queries queries from GenBank in ", time - $time, " seconds.\n";
	    }
	    if ($BLAST_GBdatabase =~ /file:(.+)/i) {
		open (FILE, '<', $1) || die "Could not open $1: $!.\n";
		print "Indexing $1 ...\n";
		my $pos;
		while (<FILE>) {
		    if (/^LOCUS/) { $pos = tell(FILE); }
		    if ($pos && /^ACCESSION\s+(\w+)/) { $index{$1} = $pos; $pos = 0;}
		}
		print "Done. Indexed ", scalar keys %index, " sequences. Proceeding.\n";
		close FILE;
	    }
	    ### sub to process query results ###
		    sub process_query {
			my $query = shift;
			my $hits = shift;
			my $taxon;
			# Get new taxon name
		       	if ($BLASTtype eq 'taxonomy') {
			    $taxon = &process_search($hits,\%categories,$ambig_critearia);
			    if ($categories{$taxon}) { $taxon = $categories{$taxon}; }
			    elsif ( $taxon eq 'error' ) { print STDERR "Error parsing BLAST query $query.\n"; undef $taxon; }
			    elsif ( $taxon eq 'unc') {}
			    elsif (  $taxon eq 'ambi') { $taxon = 'empty'; }
			    else { undef $taxon; }
			}
			elsif ($BLASTtype eq 'accno') {
			    #print "processing $query.\n";
			    $taxon = &process_accno_search($hits,$delimiter,$position,$BLAST_GBdatabase,\%index);
			    if ( $taxon eq 'error' ) { print STDERR "Error parsing BLAST query $query.\n"; undef $taxon; }
			    elsif (  $taxon eq 'ambi') { $taxon = 'empty'; }
			}
			# Insert taxon name
			if ($taxon && $taxon ne 'unc') {
			    if ($taxon eq 'empty') { print "Taxon assignment ambigious.\n"; }
			    #print "$query - $taxon\n";
			    my $insert = 'T';
			    my $scaffold = $query;
			    if ($protein_query eq 'T') {
		    		$sth_protein->execute($query);
	    			$scaffold = $sth_protein->fetchrow_array();
    			    }
			    if ($scaffold) {
				if ($overwrite =~/F$/ || $overwrite =~ /^C/) {
				    $sth_check->execute($scaffold);
				    my $temp = $sth_check->fetchrow_array();
				    if ($overwrite =~ /^C/ && $temp && $temp ne 'empty' && !($temp =~ /^$taxon/)) {
					print "BLAST based on $query say taxon should be '$taxon' for $scaffold.\n\tIt is now: $temp.\n";
					if ($overwrite =~/T$/) { print "\tIt will be changed to: $taxon.\n"; ++$changed_annotations}
				    }
				    #if ($overwrite =~/T$/) {
					#if ($temp && $temp ne 'empty') {
					#    $insert = 'F';
					#    print "\nWill keep annotation '$temp' for $scaffold, and not set it to 'empty'.\n";
					#}
					#else { print "\tIt will be changed to: $taxon.\n"; }
				    #}
				    if (!$temp || $temp eq 'empty') { $insert = 'T'; }
				    elsif ($overwrite =~ /F$/) { $insert = 'F'; }
				}
				#print "$query - $taxon - $scaffold - $insert\n";
				if ($insert eq 'T') {
				    if ( $sth->execute($taxon,$scaffold) ) {
					#print "$query - $taxon - $scaffold - $insert\n";
					++$added_taxa;
					if (!($added_taxa % 2500)) {
					    $commit_sth->execute();
					    $begin_sth->execute();
					    print "Added $added_taxa annotations (",time-$last_time,"s).\n";
					    $last_time = time;
					}
				    }
				    else { print STDERR "Could not add taxon annotation '$taxon' to '$scaffold'.\n"; }
				}
			    }
			    else { print STDERR "Could not find scaffold for query: $query\n"; }
			}
			elsif ($hits->[0] ne 'no hits') {
			    if (!$taxon) { print STDOUT "Parsing failure for $query\n\tBest hit: $hits->[0]\n"; }
			    else { print STDERR "Unable to determine taxon ($taxon) for $query\n\tBest hit: $hits->[0]\n"; }
			}
		    }
	    while (<$BLAST>) {
		if (/Query=\s*(.+)/) { # if at new entry
		    if ($query) {
			&process_query($query,\@hits);
		    }
		    $query = $1;
		    undef @hits;
		}
	    	elsif (/Sequences producing significant alignments:/) { # Reached list of hits
		    while(<$BLAST>) {
			if (/^>/ || /^ALIGNMENTS/) { last; } # hit first alignment
			elsif (/^\s*(.+)/) {
			    my $hit = $1;
			    if (defined($e_cut_off)) {
				
				if ($hit =~ /([0-9\.-]+)\s*$/) {
				    if ($1 > $e_cut_off) { last; }
				}
				else { print STDERR "Could not find e-value for '", $hit, "'.\n"; }
			    }
			    push @hits, $hit;
			}
		    }
		    if (!@hits) { $hits[0] = 'no hits'; ++$taxa_with_unsignificant_hits; }
		}
		elsif (/\*{5} No hits found \*{5}/) { # No hits
		    $hits[0] = 'no hits';
		    ++$taxa_without_hit;
		}
	    }
	    if (@hits) { 
		&process_query($query,\@hits);
	    }
	    $commit_sth->execute();
	    $begin_sth->finish();
	    $commit_sth->finish();
	    $sth->finish();
	    $sth_check->finish();
	    if ($protein_query eq 'T') { $sth_protein->finish(); }
	    my $output_string = "Added $added_taxa annotations ";
	    if ($changed_annotations) { $output_string .= "of which $changed_annotations had previous taxonomic annotation "; }
	    $output_string .= "(" . (time-$last_time) . "s). ";
	    $output_string .= "$taxa_without_hit queries did not have any hit, and $taxa_with_unsignificant_hits did not have any hits with e-value below the cut off.\n";
	    print wrap('','',$output_string);
	}
	else { die "Need name of blast output file.\n"; }
    }
######################################################################################
### Alternative to get output from script                                          ###
######################################################################################
    if ($things_to_do{'pars_annotations'}) {
	my $count = 'N';
	if ($extraSQLcondition =~ s/^COUNT;?//i) { $count = 'Y'; }
	my $not_fragmented = 'N';
	if ($extraSQLcondition=~ /NotFragmented/) { $not_fragmented = 'Y'; }
	$extraSQLcondition = &pars_extraSQLcondition($extraSQLcondition, $dbh);
	if ($extraSQLcondition) { $extraSQLcondition = " AND $extraSQLcondition"; }
	#print "$extraSQLcondition\n";
	if ($things_to_do{'pars_annotations'} eq 'protein' || $things_to_do{'pars_annotations'} eq 'CDS') {
	    my $query;
	    if ($count eq 'Y') { $query = "SELECT COUNT(*) FROM scaffolds INNER JOIN annotations ON scaffolds.name=annotations.seqid WHERE type='gene'$extraSQLcondition"; }
	    else { $query = "SELECT scaffolds.name,scaffolds.sequence,annotations.start,annotations.end,annotations.strand,annotations.id,annotations.source FROM scaffolds INNER JOIN annotations ON scaffolds.name=annotations.seqid WHERE type='gene'$extraSQLcondition"; }
	    print STDERR "Query to select genes: $query\n";
	    my $sth = $dbh->prepare($query);
	    $sth->execute();
	    if ($count eq 'Y') {
		my $n = $sth->fetchrow_array();
		print "Number of protein coding genes: $n\n";
	    }
	    else {
		my $sth_cds = $dbh->prepare("SELECT start,end,phase FROM annotations WHERE seqid=? AND start >= ? AND end <= ? AND type='CDS' AND strand=? AND source=? ORDER BY start");
		while (my $gene_info = $sth->fetchrow_arrayref()) {
		    #print "$gene_info->{'name'}\n";
		    my $seq = '';
		    $sth_cds->execute($gene_info->[0],$gene_info->[2], $gene_info->[3], $gene_info->[4], $gene_info->[6]);
		    while (my ($start,$end,$phase) = $sth_cds->fetchrow_array()) {
			$seq .= substr($gene_info->[1],$start-1,$end-$start+1);
		    }	
		    if ($gene_info->[4] eq '-') { &sequences::revcomp(\$seq); }
		    if ($not_fragmented eq 'N' || &sequences::translateDNA(\(substr $seq, 0, 3), $genetic_code) eq 'M' && &sequences::translateDNA(\(substr $seq, -3, 3), $genetic_code) eq '*') {
			print $SEQOUT '>', $gene_info->[5],"\n";
			if ($things_to_do{'pars_annotations'} eq 'protein') {
			    print $SEQOUT &sequences::translateDNA(\$seq,$genetic_code),"\n";
			}
			else { print $SEQOUT $seq,"\n"; }
		    }
		}
	    }
	}
    }
##############################################################################################################################
    if ($things_to_do{'get_reads_of_taxon'} || $things_to_do{'get_reads_of_scaffold'}){
    	my ($argument, $file_name);
	my $condition;
	if ($things_to_do{'get_reads_of_scaffold'}) {
	    ($argument, $file_name) = split " ;###; ", $things_to_do{'get_reads_of_scaffold'};
	    if ($argument =~ s/^NOT//i) { $condition="name NOT LIKE '$argument'"; }
	    else { $condition="name LIKE '$argument'"; }
	}
	elsif ($things_to_do{'get_reads_of_taxon'}) {
	    ($argument, $file_name) = split " ;###; ", $things_to_do{'get_reads_of_taxon'};
	    if ($argument =~ s/^NOT//i) { $condition= "taxon NOT LIKE '$argument'"; }
	    else { $condition="taxon LIKE '$argument'"; }
	}
	if ( -e $file_name) { 
	    print "here\n";
	    my $sth = $dbh->prepare("SELECT name FROM scaffolds WHERE $condition");
	    $sth->execute();
	    my %reads;
	    open my $ONE, ">reads1.fastq" or die "Could not open reads1.fastq: $!\n";
	    open my $TWO, ">reads2.fastq" or die "Could not open reads2.fastq: $!\n";
	    while ( my $scaffold = $sth->fetchrow_array() ){
		print "Getting reads from $scaffold\n";
		my @input = `samtools view $file_name $scaffold`;
		foreach (@input) {
		    &pars_samtool_view_row_to_hash($_,\%reads);
		}
		&print_mates(\%reads,$ONE,$TWO);
	    }
	    if ($include_mate_in_other_group eq 'Y') {
		my %mate_scaffolds;
		foreach my $name (keys %reads) {
		    if ($reads{$name}->[0]->{'mate_scaf'}) {
			$mate_scaffolds{$reads{$name}->[0]->{'mate_scaf'}->[0]}->{$name} = 1;
		    }
		}
		foreach my $scaffold (keys %mate_scaffolds) {
		    print "Getting reads from $scaffold\n";
		    my @input = `samtools view $file_name $scaffold`;
		    foreach (@input) {
			my @columns = split "\t";
			if ($mate_scaffolds{$scaffold}->{$columns[0]}) { 
			    &pars_samtool_view_row_to_hash($_,\%reads);
			}
		    }
		    &print_mates(\%reads,$ONE,$TWO);
		}
	    }
	    close $ONE;
	    close $TWO;
	    if (%reads) {
		open my $SINGLE, ">single.reads.fastq" or die "Could not open single.reads.fastq: $!\n";
		foreach my $name (keys %reads) {
		    my $stats = shift @{$reads{$name}};
		    if ($stats->{'mate_scaf'}) { print STDERR "Mate at $stats->{'mate_scaf'}->[0].\n"; }
		    else { print STDERR "No mate.\n"; }
		    if (scalar @{$reads{$name}} > 1) { print "More than one???\n"; }
		    foreach my $seq ( @{$reads{$name}} ) {
			print $SINGLE "\@$name\n$seq\n";
		    }
		    delete $reads{$name};
		}
		close $SINGLE;
	    }
	}
	else { die "Could not open bam file $file_name.\n"; }
    }
    if ($things_to_do{'N50'}) {
	my $critearia = '';
	if ($things_to_do{'N50'} ne 'all') {
	    $critearia = " WHERE taxon ";
	    #if ($things_to_do{'N50'} =~ /null/i) { $critearia .= "is "; } 
	    if ($things_to_do{'N50'} =~ s/^not //i) { $critearia .= "NOT "; }
	    #if ($things_to_do{'N50'} =~ /null/i) { $critearia .= "NULL"; }
	    $critearia .= "LIKE '$things_to_do{'N50'}'";
	}
	my $sth = $dbh->prepare("SELECT COUNT(*),SUM(LENGTH(sequence)) FROM scaffolds$critearia");
	$sth->execute();
	my ($count,$total) = $sth->fetchrow_array();
	$sth->finish();
	if ($total) {
	    $sth = $dbh->prepare("SELECT LENGTH(sequence) FROM scaffolds$critearia ORDER BY LENGTH(sequence)");
	    my $sum = 0;
	    my $length;
	    $sth->execute;
	    while ($length = $sth->fetchrow_array()) {
		$sum += $length;
		if ($sum > $total/2) { last; }
	    }
	    $sth->finish();
	    print "Number of scaffolds: $count\n";
	    print "Total length: $total\n";
	    print "N50: $length\n";
	}
	else {
	    die "Could not determin the total sequence length for taxon $things_to_do{'N50'}.\n";
	}
    }
    # else { die "Need a name of the sequence file.\n"; }
    $dbh->disconnect();
}
else {
    die "Could not connect to database $database_name.\n";
}
print STDERR "Total time: ", time-$start_time,"s\n";

exit;


###############################
### Subroutines ###############
###############################

sub help {
print "Genome database handler.\nmartin.ryberg\@ebc.uu.se\n\n";
print "Arguments\n";
print "-B/--BLASTtaxon_annotation       Give taxonomic annotation to scaffolds based on blast file (default BLASTN 2.2.28+ output) given\n";
print "                                 as next argument. Sequences in the BLAST database must start with unambigous taxon name\n";
print "                                 matching keys in the hash \%categories (default: fungi and bact), or have easily parsable GeneBank accession\n";
print "                                 numbers. The accetion number should be in a given possition (index starting with 0; default 1) in\n";
print "                                 relation to a delimiter (| by default), e.g. position 1 in 'gb|KIL54703.1|  hypothetical...'.\n";
print "                                 Extra settings can be given as an extra argument. If the string contain: 'annotation:acc',\n";
print "                                 (or longer version of accession number) it will set the function to take an accession number, and\n";
print "                                 check the taxonomy in GenBank (require internet connection unless file in GenBank format is given);\n";
print "                                 'GB:' followed by the quoted ('|\") name of an NCBI entrez database name (eg. GB:'nucleotide'; default: protein),\n";
print "                                 will set which database to look in. If database is set to file: followed by the name of a file in full GenBank\n";
print "                                 format, the corresponding file will be used as database (eg. GB:'file:myGBannotations.txt'). By giving\n";
print "                                 'file:download' all accession numbers for matching sequences in the blast file will be downloaded to a\n";
print "                                 file before the annotation starts. If 'file:complement' is given followed by another colon (:) and the name of a\n";
print "                                 file in genbank format, accession numbers not in the given file will be downloaded and added to the file. If any\n";
print "                                 of these two latter option are used the NCBI database can be given after another colon (eg.\n";
print "                                 'file:download:nucleotide'). 'delimiter:' will set the delimiter to the next character (no whitespace in between)\n";
print "                                 as delimiter; 'position:' will set the possition of the accession number to the following integer number (no\n";
print "                                 whitespace in between); 'query:gene' will set that genes were used for the queries (or something else with an\n";
print "                                 annotation in the annotations table in the database), and not scaffolds; 'e-cut_off' followed by a value will set\n";
print "                                 the maximum e-value to consider as a significant match (ignoring other matches). 'check' will set that it is checked if\n";
print "                                 the new taxon annotations are different from the previous; 'overwrite' set that previous taxon annotation will be\n";
print "                                 overwritten. 'taxon:' followed by a quoted ('|\") string will set the \%categories hash. If the string starts with\n";
print "                                 'file:' and is followed by  a file name, the taxon translation will be read from the file. The file shold be a tab\n";
print "                                 separated text file where the keys are in the first column and taxon names in the second (one pair per row).\n";
print "                                 Otherwise the string should be comma separated with pairs of key followed by taxon name. Examples:\n";
print "                                     perl genomeDB.pl -db database.db -B  blast.file \"annotation:acc;query:gene;GB:'file:sequence.gp'\"\n";
print "                                     perl genomeDB.pl -db database.db -B  blast.file \"taxa: 'fungi, Eukaryota; Fungi, bact, Prokaryote\"\n";
print "-C/--create_database             Create the database structure.\n";
print "-db/--database_name              Give the name of the database.\n";
print "--get_CDS                        Print protein coding genes (annotation type gene, including at least one streach with annotation\n";
print "                                 type CDS) as nucleotide sequences. It is possible to give extra conditions for which sequences\n";
print "                                 to print, e.g. coverage, proteinCoverage, GC, proteinGC, taxon, GO, GOname or GOname_space for\n";
print "                                 limiting the output based on median coverage, GC content, taxon of of scaffold that the protein is\n";
print "                                 coded on, median coverage, GC content of proteins, or GO number, GO number associated with certain\n";
print "                                 name, or GO number belonging to a certain name space in the gene ontology_term. For the later two the\n";
print "                                 file go.obo (http://geneontology.org/page/download-ontology) is needed to be in the same folder as the\n";
print "                                 script is running in. The type of condition can be given as <, > (default), =, <=, >=, or LIKE follow\n";
print "                                 by a value. If like is choosen % can be used as wildcard character. Examples: 'coverage 50' or 'coverage\n";
print "                                 < 10'. Several constraints can be given separated by semicolon (;). The constraints will be additive\n"; 
print "                                 (AND). It is also possible to give SQL syntax directly. It will the be treated as added to the SQL statment\n"; 
print "                                 after an AND. If you do not want the sequences but only the number of genes, you can give COUNT as the first\n";
print "                                 condition. If NotFragmented is given, only genes that do not start or end at the first or last base of the scaffold\n";
print "                                 is counted. If printing sequences it also only print sequences for which the first codon match a start codon and\n";
print "                                 the last match a stop codon. This additional requirement when printing sequences mean that you may get a different\n";
print "                                 number of sequences printed than the value you get when counting. Example:\n";
print "                                     perl genomeDB.pl -db database.db --get_CDS \"coverage > 30; coverage<50;GC > 0.41; GC < 0.5\"\n";
print "--get_proteins                   Print protein coding genes (annotation type gene, including at least one streach with annotation\n";
print "                                 type CDS) as aminoacid sequences. It is possible to give extra conditions for which sequences\n";
print "                                 to print as for --get_CDS.\n";
print "-R/--get_reads_of_scaffold       Print the reads that match given scaffold. The scaffold name should be the next argument and a\n";
print "                                 bam file the next argument after that.\n";
print "-T/--get_reads_of_taxon          Print the reads that match scaffolds of given taxon. The taxon should be given as next argument,\n";
print "                                 % may be used as wildchard character, and a bam file should be the next argument after that.\n";
print "-h/--help                        Print this help.\n";
print "--include_missmatch_mates        Include reads that do not match taxon/scaffold if their mate match.\n";
print "--N/--note                       Set note annotation when reading scaffolds\n";
print "--N50                            Get statistics, including N50, for scaffolds belonging to given taxon. The taxon must be the\n";
print "                                 next argument, if the taxon is preceded by 'not ' it will give the statistics for scaffold\n";
print "                                 not belonging to the taxon. If all (default) is given the statistics for all scaffolds will\n";
print "                                 be given. If 'null' is given statistics will be provided for the scaffolds that do not have\n";
print "                                 any taxon annotation. % may be used as wildchard character.\n";
print "-a/--read_annotations            Read annotations from GFF3 formated file. Can only read annotations for scaffolds that are in\n";
print "                                 the database. File name given as next argument.\n";
print "-b/--read_bam                    Read coverage and SNPs for scaffolds in the database from a bam file, given as next argument. If\n";
print "                                 you want to read Median coverage and SNPs for annotations too give anno as an extra argument\n";
print "-v/--read_vcf                    Read SNPs for scaffolds in the database from a VCF file, given as next argument. If you want to read\n";
print "                                 number of SNPs for annotations as well, give anno as an extra argument\n";
print "-s/--read_scaffolds              Read scaffold data fasta file, given as next argument.\n";
print "--set_genetic_code               Set the genetic code that will be used to translate base triplets into amino acids. eg. standard (default).\n";
print "-S/--SNPscore_cutoff             Set base quality score cut-off to use for calling SNPs when reading bam file\n";
print "-t/--taxon                       Set taxon name when reading scaffolds.\n";

}

sub add_to_scaffolds {
    my $dbh = shift;
    my $name_ref = shift;
    my $seq_ref = shift;
    my $taxon_ref = shift;
    my $note_ref = shift;
    my $GCcontent = &sequences::GCcontent($seq_ref);
    return $dbh->do("INSERT INTO scaffolds (name, sequence, GCcontent, taxon, note) VALUES ($$name_ref, $$seq_ref, $GCcontent, $$taxon_ref, $$note_ref)");
}

sub avarage_and_sd {
    my $avarage=0;
    foreach (@{$_[0]}) { $avarage += $_; }
    $avarage /= scalar @{$_[0]};
    my $sd = 0;
    foreach (@{$_[0]}) { $sd += ($_-$avarage)*($_-$avarage); }
    $sd /= scalar @{$_[0]};
    $sd = sqrt($sd);
    return $avarage, $sd;
}

sub process_accno_search {
    my $hits = shift;
    my $delimiter = shift;
    my $position = shift;
    my $BLAST_GBdatabase = shift;
    my $index_ref = shift;
    my $ntop_hits = 10;
    if ($hits and @{$hits}) {
	if ($hits->[0] eq 'no hits') {
	    #print "no hits\n";
            return 'unc';
        }
	else {
	    my $organism = 'unc';
	    my @accnos;
	    for (my $i = 0; $i < scalar @{$hits} && $i < $ntop_hits; ++$i) {
		if (!$hits->[$i]) { push @accnos, "empty"; }
		else {
		    my @annotations = split /$delimiter/, $hits->[$i];
		    if ($annotations[$position]) { push @accnos, $annotations[$position]; }
		    else { push @accnos, "empty"; }
		}
	    }
	    my $accno_string = '';
	    my @search;
	    my $FILE;
	    my $file_name;
	    if ($BLAST_GBdatabase =~ /^file:(.+)/) {
		$file_name = $1;
		open ($FILE, '<', $file_name) or die "Could not open $file_name: $!.\n";
		#print "$file_name is open.\n";
		if(tell($FILE) == -1) { print "Possible file reading error, when reopening $file_name.\n"; }
	    }
	    if (@accnos) {
		foreach (@accnos) {
		    if ($_ ne "empty") {
			if ($accno_string) { $accno_string .= ",$_"; }
			else { $accno_string = $_; }
		    }
		}
		if (!(defined($FILE))) {
		    # print "Getting taxon annotations: http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$BLAST_GBdatabase&id=$accno_string&rettype=gb&retmode=text\n";
		   # print "Getting entry from NCBI: ", time,"\n";
		    my $get= sequences::get_entry_from_GB($accno_string,$BLAST_GBdatabase); #get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$BLAST_GBdatabase&id=$accno_string&rettype=gb&retmode=text");
		    if ($get && $get =~ /^Resource temporarily unavailable/) { print STDERR "NCBI eutils temporarily unavailable.\n"; }
		    elsif ($get) { @search = split /\n/, $get; }
		    else { print STDERR "No results for $accno_string in $BLAST_GBdatabase.\n"; }
		}
	    }
	    else { print STDERR "Could not pars any accnos.\n"; }
	    my @taxonomy;
	    #print "Processing annotations.\n";
	    ####### Sub to parse taxon from GB data ########
	    sub taxon_from_GB {
		my $row = shift;
		my $accno = shift;
		my $right_entry = shift;
		my $taxon = shift;
		if ($$row =~ /^ACCESSION\s+(\w+)/) {
		    #print $$row, "\n";
		    my $entry = $1;
		    if (!$accno) { print "No accno\n"; }
		    if ($accno =~ /$entry/) { $right_entry = 'T'; }#print "Found it.\n";} 
		    else { $right_entry = 'F'; }
		}
		elsif ($right_entry eq 'T' && $$row =~ /\s*ORGANISM\s+(.+)/) { $$taxon = '';}
		elsif (defined($$taxon) && $right_entry eq 'T') {
		    s/\s+//g;
		    $$taxon .= $_;
		    if (/\.$/) { $right_entry = 'L'; }
		}
		return $right_entry;
	    }
	    #################################################
	    for (my $i = 0; $i < scalar @{$hits} && $i < $ntop_hits; ++$i) {
	       	if ($accnos[$i] && $accnos[$i] ne 'empty') {
		    #print "Made it here.\n";
		    {
			my @temp = split /\s+/, $hits->[$i];
			push @taxonomy, [];
			$taxonomy[$#taxonomy]->[0] = pop @temp;
			$taxonomy[$#taxonomy]->[1] = pop @temp;
			#print "e-value: $taxonomy[$#taxonomy]->[0]; score: $taxonomy[$#taxonomy]->[1]\n";
		    }
		    my $taxon;
		    my $right_entry = 'F';
		    if ($FILE) {
			if(tell($FILE) == -1) {
			    open ($FILE, '<', $file_name) or die "Could not open $file_name: $!.\n";
			    #print "... or is it?\n";
			}
			if ($index_ref) {
			    #print "Looking in index.\n";
			    my $temp = $accnos[$i];
			    $temp =~ s/\.[0-9]+$//;
			    if ($index_ref->{$temp}) {
				seek $FILE,$index_ref->{$temp}-1,0;
				#print "Found $temp at $index_ref->{$temp}.\n";
			    }
			    else { seek $FILE,-1,2; }
			}
			else { seek $FILE, 0, 0; }
		    }
		    my $condition = 1; # Flag to stop loop
		    if (defined($FILE)) {
			#print "looking for taxonomy.\n";
			my $j=0;
			while (<$FILE>) {
			    #print ++$j, "\n";
			    $right_entry = &taxon_from_GB(\$_,$accnos[$i],$right_entry,\$taxon);
			    if ( $right_entry eq 'L' ) { last; }
			}
		    }
		    else {
			foreach (@search) {
			    $right_entry = &taxon_from_GB(\$_,$accnos[$i],$right_entry,\$taxon);
			    if ( $right_entry eq 'L' ) { last; }
			}
		    }
		    if ($taxon) { push @{$taxonomy[$#taxonomy]}, split /\s*;\s*/, $taxon }
		    #else { print STDERR "Could not find taxon annotation for $accnos[$i].\n"; }
		}
		else { print STDERR "No GB accession number found for:\n\t $hits->[$i].\n"; }
	    }
	    my $taxonstring='';
	    my $scoresum=0;
	    foreach (@taxonomy) { $scoresum+= $_->[1]; }
	    my $cut_off = 0;
	    my $part_sum = 0;
	    while ($part_sum < $scoresum/2) { $part_sum += $taxonomy[$cut_off++]->[1]; }
	    my $done = 'F';
	    if ($cut_off > 1) {
		my $pos = 2;
		while ($done eq 'F') {
		    my $check = 'F';
		    my $taxon;
		    for ( my $i=0; $i < $cut_off-1; ++$i) {
			if ($pos < scalar @{$taxonomy[$i]} || $pos < scalar @{$taxonomy[$i+1]}) {
			    $check = 'T';
			    if ($taxonomy[$i]->[$pos] && $taxonomy[$i+1]->[$pos] && $taxonomy[$i]->[$pos] ne $taxonomy[$i+1]->[$pos]) { $done = 'T'; last}
			    elsif (!$taxon && $taxonomy[$i]->[$pos]) { $taxon = $taxonomy[$i]->[$pos];}
			    elsif (!$taxon && $taxonomy[$i+1]->[$pos]) { $taxon = $taxonomy[$i+1]->[$pos];}
			}
		    }
		    if ($check eq 'F') { $done = 'T'; }
		    if ($done eq 'F' && $taxon) {
			if ($taxonstring) { $taxonstring .= "; $taxon"; }
			else { $taxonstring = $taxon; }
		    }
		    ++$pos;
		}
		if ($taxonstring) { $organism = $taxonstring; }
		else {
		    print "No taxon determined for $accno_string.\n";
		    if ($taxonomy[0] && $taxonomy[0]->[2]) {
			print "\tKingdoms (score) of hits in order";
			foreach (@taxonomy) { if ($_->[2]) { print ", ", $_->[2],' (', $_->[1], ')';} }
			print "\n";
			$organism = 'ambi';
		    }
	       	}
	    }
	    elsif (@taxonomy && $taxonomy[0]->[2]) {
		$taxonstring = $taxonomy[0]->[2];
		for (my $i=3; $i<scalar @{$taxonomy[0]}; ++$i) {
		    $taxonstring .= "; $taxonomy[0]->[$i]";
		}
		if ($taxonstring) { $organism = $taxonstring; }
	    }
	    else {
		print STDERR "Not able to determin taxon string.\n";
		if (!@taxonomy) { print STDERR "\tParsing of hits failed.\n"; }
		elsif (!$taxonomy[0]->[2]) { print STDERR "\tNo taxonomy parsed for first hit.\n"; }
		elsif ($cut_off == 0) { print STDERR "\tNot able to determin score cut off.\n"; }
	    }
	    #print $organism, "\n";
	    return $organism;
	}
    }
    else { return "error"; }
}

sub process_search {
    my $hits = shift;
    my $categories = shift;
    my $criteria = shift;
    my $ntop_hits = 10;
    if ($hits and @{$hits}) {
        if (${$hits}[0] eq 'no hits') {
            return 'unc';
        }
        else {
            my $organism = 'unc';
	    my $sum_score=0;
	    for (my $i = 0; $i < $ntop_hits && $i < scalar @{$hits}; ++$i) {
		my @temp = split /\s+/, $hits->[$i];
		my $evalue = pop @temp;
		my $score = pop @temp;
		if ($score && $score > 0) { $sum_score += $score; }
	    }
	    my $cut_off = $sum_score/2;
	    $sum_score=0;
            for (my $i = 0; $i < scalar @{$hits}; ++$i) {
                if (!$hits->[$i]) { next; }
                my $flag = 'n';
                foreach my $type (keys %{$categories}) {
                    if ($hits->[$i] =~ /^$type/) {
                        $flag = 'y';
			if ($criteria == 1) {
			    if ($organism eq 'unc') { $organism = $type; }
			    elsif ($organism ne $type) { $organism = 'ambi'; }
			}
			elsif ($criteria == 2) {
			    if ($organism eq 'unc') { $organism = $type; }
			}
			elsif ($criteria == 3) { 
			    if ($i < $ntop_hits && $sum_score < $cut_off) {
				if ($organism eq 'unc') { $organism = $type; }
				elsif ($organism ne $type) { $organism = 'ambi'; }
				my @temp = split /\s+/, $hits->[$i];
				my $evalue = pop @temp;
				$sum_score += pop @temp;
			    }
			}
                         #       $organism = &evaluate_ambiguity($organism,$hits,$i,$critearia);
			last;
		    }
		}
		if ($organism eq 'ambi') { last; }
                if ($flag eq 'n') { print STDERR "Do not recognize $hits->[$i] --- $i --- ", scalar @{$hits}, " --- ", join (', ', keys %{$categories}), "\n"; }
            }
            return $organism;
        }
    }
    else { return "error"; }
}

sub evaluate_ambiguity {
    my $type = shift; # The organism it is already classified as
    my $hits = shift; # A reference to an array of blast hits
    my $position = shift; # The index of the present hit being evaluated
    my $criterion = shift; # Which of the below criterions to use
    my $ntop_hits = shift;
    if ($criterion == 1) {
        return 'ambi';
    }
    elsif ($criterion == 2) {
        if ($position > 0) { return $type; }
        else { return 'ambi'; }
    }
    elsif ($criterion == 3) {

    }
    else { return 'ambi'; }
}

sub pars_extraSQLcondition {
    my @string = split /;/, shift;
    my $dbh = shift;
    my $return_string = '';
    foreach (@string) {
	print "###$_###\n";
	#s/('|")/\\$1/g;
	if ($return_string) { $return_string .= " AND "; }
	if ( s/^SQL://i) { $return_string .= $_; }
	elsif ($_ eq 'NotFragmented') { $return_string .= "NOT (annotations.start==1 OR annotations.end==LENGTH(scaffolds.sequence))"; }
	elsif (/(proteincoverage|coverage|proteinGC|GC|taxon|GOname_space|GOname|GO)(:| )*(>=|<=|==|!=|<|>|=|LIKE|NOT LIKE)*\s*(.+)/i) {
	    print "$_\t$1 $3 $4\n";
	    my $field = uc($1);
	    my $sign;
	    if ($3) {$sign= uc($3);}
	    my $value = $4;
	    if ($value =~ /[^0-9\.]/ && !($field eq 'GONAME' || $field eq 'GONAME_SPACE')) { $value = "'$value'"; }
	    if (!$sign) { $sign = '>'; }
	    #print "$field $sign $value\n";
	    if ($field eq 'COVERAGE') {	$return_string .= "scaffolds.coverage_median "; }
	    elsif ($field eq 'GC') { $return_string .= "scaffolds.GCcontent "; }
	    elsif ($field eq 'PROTEINCOVERAGE') { $return_string .= "annotations.coverage_median "; }
	    elsif ($field eq 'PROTEINGC') { $return_string .= "annotations.GCcontent "; }
	    elsif ($field eq 'TAXON') { $return_string .= "scaffolds.taxon "; }
	    elsif ($field eq 'GO') { $return_string .= "annotations.ontology_term LIKE \"\%GO:$value\%\""; next}
	    elsif ($field eq 'GONAME' || $field eq 'GONAME_SPACE') {
		if ($return_string) { $return_string = "($return_string"; }
		open GOTERMS, "<go.obo" or die "Need the file go.obo (http://geneontology.org/page/download-ontology) to translate GO names or name spaces into GO numbers.\n";
		my $present_number = '';
		my @GOnumbers;
		my $match = 'N';
		my %GOinDB;
		&getGOnumbersINdb(\%GOinDB,$dbh);
		while (<GOTERMS>) {
		    if ($present_number && /^\[Term\]|^\[Typedef\]/) {
			if ($match eq 'Y' && $GOinDB{$present_number}) { # only include sequences that match the value and are in the database
			    push @GOnumbers, $present_number;
			    $match = 'N';
			}
			$present_number = '';
			if (/^\[Typedef\]/) { last; }
		    }
		    elsif (/^id: (\S+)/) {
			$present_number = $1;
		    }
		    elsif ($field eq 'GONAME' && /^(name: |synonym: |is_a: .+!).*$value.*/ ) { $match = 'Y'; }
		    elsif ($field eq 'GONAME_SPACE' && /^namespace: $value/) { $match = 'Y'; }
		}
		close GOTERMS;
		if (scalar @GOnumbers < 1) { die "No matching GO numbers for '$value'.\n"; }
		else {
		    $return_string .= "(";
		    @GOnumbers = sort @GOnumbers;
		    print STDERR "Number of GO numbers added: ", scalar @GOnumbers, "\n";
		    # Since SQLite can not handle more than 1000 alternatives the following was an attempt to decrease the number
		    # of GO numbers by reducing the number of LIKE statements if all numbers for a ten or hundred are included
		    # but it seem like this is rarely the case.
		    my $first_in_series = 0;
		    for (my $i=1; $i <= scalar @GOnumbers; ++$i) {
			my ($present,$previous);
			if ($i < scalar @GOnumbers) {
			    $present = $GOnumbers[$i]; $present =~ s/^GO://;
			    $previous = $GOnumbers[$i-1]; $previous =~ s/^GO://;
			}
			if ($i < scalar @GOnumbers && $present == $previous+1) {
			    next; # if in series of number and not gone outside the array, move on
			}
		     	else {
			    my $zeroes = '0';
			    while ($first_in_series < $i) { # Do not go outside the array
				if ($GOnumbers[$first_in_series] =~ /$zeroes$/ && ($i - $first_in_series) > 10**length($zeroes)) { # If at even number and we have more numbers than to next even number
				    $zeroes .= '0'; # try with one order of magnitude larger number
				}
				else { # if not at number as even as we tried
				    if (($i - $first_in_series) < 10**length($zeroes) && length($zeroes) > 1) {
					$zeroes = substr($zeroes,1); # if we do not have numbers enough to reach next even number we must start printing
				    }
				    my $value = substr ($GOnumbers[$first_in_series],0,length($GOnumbers[$first_in_series])-length($zeroes)+1); # remove numbers at end if we have them all in the series
				    if ($first_in_series > 0) { $return_string .= " OR "; } # if not the first number we should print
				    $return_string .= "annotations.ontology_term LIKE \"\%$value\%\""; # print the number to search string
				    $first_in_series += 10**(length($zeroes)-1); # move on to as many numbers we have covered
				}
			    }
			    $first_in_series = $i;
			}
		    }
		    $return_string .= ")";
		}
		next;
	    }
	    $return_string .= "$sign ";
	    $return_string .= "$value";
	}
	else { $return_string .= $_; }
    }
    return $return_string;
}

sub getGOnumbersINdb {
    my $hash_ref = shift;
    my $dbh = shift;
    my $sth = $dbh->prepare("SELECT ontology_term FROM annotations");
    $sth->execute();
    while (my $ontology_term = $sth->fetchrow_array()) {
	my @temp = split /;/, $ontology_term;
	foreach (@temp) {
	    if (/GO:/) { ++$hash_ref->{$_}; }
	}
    }
    $sth->finish();
}
######################################
### For parsing and printing reads ###
######################################

sub pars_samtool_view_row_to_hash {
    my @columns = split "\t", shift;
    my $hash_ref = shift;
    if (!$hash_ref->{$columns[0]}) { $hash_ref->{$columns[0]}->[0] = {}; }
    if ($columns[1] & 0x0010) {
	&sequences::revcomp(\$columns[9]);
	$columns[10] = reverse $columns[10];
	--$hash_ref->{$columns[0]}->[0]->{'rev'};
    }
    else { ++$hash_ref->{$columns[0]}->[0]->{'rev'}; }
    if ($columns[6] ne '=') {
	if (!$hash_ref->{$columns[0]}->[0]->{'mate_scaf'}) { $hash_ref->{$columns[0]}->[0]->{'mate_scaf'} = []; }
	push @{$hash_ref->{$columns[0]}->[0]->{'mate_scaf'}}, $columns[6];
    }
    push @{$hash_ref->{$columns[0]}}, "$columns[9]\n+\n$columns[10]";
}

sub print_mates {
    my $hash_ref = shift;
    my $ONE = shift;
    my $TWO = shift;
    foreach my $name (keys %{$hash_ref}) {
	if (scalar @{$hash_ref->{$name}} > 2) {
	    my $stats = shift @{$hash_ref->{$name}};
	    print $ONE "\@$name\n", shift @{$hash_ref->{$name}}, "\n";
	    print $TWO "\@$name\n", shift @{$hash_ref->{$name}}, "\n";
	    if (@{$hash_ref->{$name}}) { print STDERR "More than two reads with the id $name. Only two printed to files!!!\n"; }
	    elsif ($stats->{'rev'} > 0) { print STDERR "Both reads of $name in forward direction.\n"; }
	    elsif ($stats->{'rev'} < 0) { print STDERR "Both reads of $name in reverse direction.\n"; }
	    if ($stats->{'mate_scaf'}) { print STDERR "Mate at $stats->{'mate_scaf'}->[0].\n"; }
	    delete $hash_ref->{$name};
	}
    }
}
