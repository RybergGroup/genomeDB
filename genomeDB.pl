#! /usr/bin/perl -w

use strict;
use DBI;
use Bio::DB::Sam;
use LWP::Simple;
use translation;

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
my %categories = ( 'fungi' => 'Fungi', 'bact' => 'Bacteria' ); # Give your organisms annotations as keys
my $ambig_critearia = 1; # The criteria to decide if the classification is ambiguous
                         # The function that use the criteria is evaluate_ambiguity
my $overwrite = 'CF';
my $delimiter = '\|';
my $position = 1;
my $BLAST_GBdatabase = 'protein';

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
	    if (!&translation::code_available($genetic_code)) { die "'$genetic_code' is not an available genetic code. Available options are: ", &translation::available_codes_string(), "\n"; }
	    
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
	    if ( $ARGV[$i] =~ /GB:([A-Za-z]+)/i ) { $BLAST_GBdatabase = lc($1); }
	    if ( $ARGV[$i] =~ /annotation:acc/i ) { $BLASTtype = 'accno'; }
	    if ( $ARGV[$i] =~ /delimiter:(.)/i) {
		$delimiter = $1;
		if ($delimiter eq '|' || $delimiter eq '.' || $delimiter eq '(' || $delimiter eq ')' || $delimiter eq '[' || $delimiter eq ']' || $delimiter eq '{' || $delimiter eq '}' || $delimiter eq '\\'){
		    $delimiter = "\\$delimiter";
		}
	    }
	    if ( $ARGV[$i] =~ /position:([0-9]+)/ ) { $position = $1; }
	    
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
	if ($dbh->do("CREATE TABLE scaffolds (name TEXT PRIMARY KEY, sequence TEXT DEFAULT 'empty', taxon TEXT DEFAULT 'empty', coverage_median INTEGER, coverage_mean REAL, coverage_sd REAL, coverage_max INTEGER, coverage_min INTEGER,  GCcontent REAL, medianSNPdiversity REAL, nSNP INTEGER, note TEXT DEFAULT 'empty')")) {
	    print "Created table scaffolds (",time-$last_time,"s).\n";
	    $last_time = time;
	}
	else { die "Could not create table scaffolds in database $database_name.\n"; }
	if ($dbh->do("CREATE TABLE annotations (seqid TEXT, source TEXT, type TEXT, start INTEGER, end INTEGER, score REAL, strand TEXT, phase INTEGER, id TEXT, name TEXT, alias TEXT, parent TEXT, target TEXT, gap TEXT, derives_from TEXT, note TEXT, dbxref TEXT, ontology_term TEXT, is_circular TEXT, extras TEXT, FOREIGN KEY(seqid) REFERENCES scaffolds(name))")) {
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
	    if (!$taxon) { $taxon = 'empty'; }
	    else { $taxon = "'$taxon'"; }
	    if (!$note) { $note ='empty'; }
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
	    my $sth = $dbh->prepare("INSERT INTO annotations VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
	    my @other_columns = ('id', 'name', 'alias', 'parent', 'target', 'gap', 'derives_from', 'note', 'dbxref', 'ontology_term', "is_circular", "extras");
	    my $begin_sth = $dbh->prepare("BEGIN"); # for faster processing commit will be done explicitly
	    $begin_sth->execute();
	    my $commit_sth = $dbh->prepare("COMMIT");
	    $last_time = time;
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
	if ($things_to_do{'read_bam'} ne 'T') {
	    if (my $sam = Bio::DB::Sam->new( -bam  => $things_to_do{'read_bam'})) {
		print "Reading bam file...\n";
		my $sth = $dbh->prepare("SELECT name FROM scaffolds");
		my $added_stats = 0;
		$sth->execute();
		my $sth_insert = $dbh->prepare("UPDATE scaffolds SET coverage_median=?, coverage_mean=?, coverage_sd=?, coverage_max=?, coverage_min=?, medianSNPdiversity=?, nSNP=? WHERE name=?");
		print "Processing sequences.\n";
		$last_time = time;
		while (my $seq_name = $sth->fetchrow_array()) {
		    my $depth = 0;
		    my $positions = 0;
		    my @data;
		    my @SNPdiversity;
		    my $callback = sub {
			my ($seqid,$pos,$pileup) = @_;
			$positions++;
			my $count=0;
			my %n_bases;
			foreach my $p (@$pileup) {
			    my $ali = $p->alignment;
			    if ($ali->qscore->[$p->pos - $ali->query->start + $p->indel] > $score_cut_off) {
				++$n_bases{uc(substr($ali->qseq,$p->pos - $ali->query->start + $p->indel,1))};
			    }
			    ++$depth;
			    ++$count;
			}
			push @data, $count;
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
			}
		    };
		    $sam->fast_pileup($seq_name,$callback);
		    if ($positions) {
			@data = sort {$a <=> $b} @data;
			my $average = 0;
			foreach (@data) { $average += $_; }
			$average /= scalar @data;
			my $sd = 0;
			foreach (@data) { $sd += ($_-$average)*($_-$average); }
			$sd /= scalar @data;
			$sd = sqrt($sd);
			if ($sth_insert->execute($data[scalar @data/2], $average, $sd, $data[-1], $data[1], $SNPdiversity[(scalar @SNPdiversity)/2], scalar @SNPdiversity, $seq_name)) {
			    ++$added_stats;
			    if (!($added_stats % 1000)) {
				print "Added $added_stats annotations (",time-$last_time,"s).\n";
				$last_time = time;
			    }
			}
			else {
			    print STDERR "Could not add stats to: $seq_name\n";
			}
		    }
		    else { print STDERR "No pileup for $seq_name.\n"; }
		}
		$sth_insert->finish();
		$sth->finish();
		print "Added stats to $added_stats scaffolds (",time-$last_time,"s).\n";
		$last_time = time;
	    }
    	    else { die "Could not open $things_to_do{'read_bam'}.\n"; }
	}
	else { die "Need a name of the bam file.\n"; }
    }
######################################################################################
    if ($things_to_do{'pars_blast'}) {
	if ($things_to_do{'pars_blast'} ne 'T') {
	    print "Will read BLAST results from $things_to_do{'pars_blast'}.";
	    if ($BLASTtype eq 'taxonomy') { print " Assuming first word of each match is a taxon name."; }
	    elsif ($BLASTtype eq 'accno') { print " Assuming accession number can be found in position $position based on $delimiter, and corresponding annotation can be found in $BLAST_GBdatabase."; }
	    if ($protein_query eq 'T') { print " Will look for query name in annotations table."; }
	    else { print " Will look for query name in scaffolds table."; }
	    if ($overwrite =~ /^C/) { print " Will compare new with old taxon annotation and print if they are different."; }
	    if ($overwrite =~ /T$/) { print " Will overwrite old annotations."; }
	    elsif ($overwrite =~ /F$/) { print " Will not overwrite old annotations."; }
	    print "\n";
	    my $query;
	    my @hits;
	    my $added_taxa=0;
	    open my $BLAST, "<$things_to_do{'pars_blast'}" or die "Could not open $things_to_do{'pars_blast'}: $!.\n";
	    my $sth = $dbh->prepare("UPDATE scaffolds SET taxon=? WHERE name=?");
	    my $sth_protein;
	    if ($protein_query eq 'T') { $sth_protein = $dbh->prepare("SELECT seqid FROM annotations WHERE id=?"); }
	    my $sth_check = $dbh->prepare("SELECT taxon FROM scaffolds WHERE name==?");
	    my $begin_sth = $dbh->prepare("BEGIN");
    	    $begin_sth->execute();
	    my $commit_sth = $dbh->prepare("COMMIT");
	    while (<$BLAST>) {
		if (/Query=\s*(.+)/) { # if at new entry
		    if ($query) {
			my $taxon;
			# Get new taxon name
		       	if ($BLASTtype eq 'taxonomy') {
			    $taxon = &process_search(\@hits,\%categories,$ambig_critearia);
			    if ($categories{$taxon}) { $taxon = $categories{$taxon}; }
			    elsif ( $taxon eq 'error' ) { print STDERR "Error parsing BLAST query $query.\n"; undef $taxon; }
			    else { undef $taxon; }
			}
			elsif ($BLASTtype eq 'accno') {
			    $taxon = &process_accno_search(\@hits,$delimiter,$position,$BLAST_GBdatabase);
			    if ( $taxon eq 'error' ) { print STDERR "Error parsing BLAST query $query.\n"; undef $taxon; }
			    #print "$taxon\n";
			}
			# Insert taxon name
			if ($taxon && $taxon ne 'unc') {
			    my $insert = 'T';
			    my $scaffold = $query;
			    if ($protein_query eq 'T') {
		    		$sth_protein->execute($query);
	    			$scaffold = $sth_protein->fetchrow_array();
				#print "$query - $scaffold\n"
    			    }
			    if ($scaffold) {
				#print "$query - $scaffold\n";
				if ($overwrite eq 'F' || $overwrite eq 'CT' || $overwrite eq 'CF') {
				    $sth_check->execute($scaffold);
				    my $temp = $sth_check->fetchrow_array();
				    if ($overwrite =~ /^C/ && $temp && $temp ne 'empty' && $temp ne $taxon) {
					print "BLAST based on $query say taxon should be $taxon for $scaffold.\n\tit is now $temp.\n";
					if ($overwrite =~/T$/) { print "\tIt will be changed to $taxon."; }
				    }
				    if (!$temp || $temp eq 'empty') { $insert = 'T'; }
				    elsif ($overwrite =~/F$/) { $insert = 'F'; }
				}
				#print "$query - $taxon - $scaffold\n";
				if ($insert eq 'T') {
				    if ( $sth->execute($taxon,$scaffold) ) {
					#print "$scaffold - $taxon\n";
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
			elsif ($hits[0] ne 'no hits') { print STDERR "Unable to determine taxon for $query\nBest hit: $hits[0]\n"; }
		    }
		    $query = $1;
		    undef @hits;
		}
	    	elsif (/Sequences producing significant alignments:/) { # Reached list of hits
		    while(<$BLAST>) {
			if (/^>/) { last; } # hit first alignment
			elsif (/^\s*(\S.+)/) {
			    push @hits, $1;
			}
		    }
		}
		elsif (/\*{5} No hits found \*{5}/) { # No hits
		    $hits[0] = 'no hits';
		}
	    }
	    $commit_sth->execute();
	    $begin_sth->finish();
	    $commit_sth->finish();
	    $sth->finish();
	    $sth_check->finish();
	    if ($protein_query eq 'T') { $sth_protein->finish(); }
	    print "Added $added_taxa annotations (",time-$last_time,"s).\n";
	}
	else { die "Need name of blast output file.\n"; }
    }
######################################################################################
### Alternative to get output from script                                           ###
######################################################################################
    if ($things_to_do{'pars_annotations'}) {
	my $count = 'N';
	if ($extraSQLcondition =~ s/^COUNT;?//i) { $count = 'Y'; }
	$extraSQLcondition = &pars_extraSQLcondition($extraSQLcondition, $dbh);
	if ($extraSQLcondition) { $extraSQLcondition = " AND $extraSQLcondition"; }
	#print "$extraSQLcondition\n";
	if ($things_to_do{'pars_annotations'} eq 'protein' || $things_to_do{'pars_annotations'} eq 'CDS') {
	    my $query;
	    if ($count eq 'Y') { $query = "SELECT COUNT(*) FROM scaffolds INNER JOIN annotations ON scaffolds.name=annotations.seqid WHERE type='gene'$extraSQLcondition"; }
	    else { $query = "SELECT scaffolds.name,scaffolds.sequence,annotations.start,annotations.end,annotations.strand,annotations.id FROM scaffolds INNER JOIN annotations ON scaffolds.name=annotations.seqid WHERE type='gene'$extraSQLcondition"; }
	    print STDERR "Query to select genes: $query\n";
	    my $sth = $dbh->prepare($query);
	    $sth->execute();
	    if ($count eq 'Y') {
		my $n = $sth->fetchrow_array();
		print "Number of protein coding genes: $n\n";
	    }
	    else {
		my $sth_cds = $dbh->prepare("SELECT start,end,phase FROM annotations WHERE seqid=? AND start >= ? AND end <= ? AND type='CDS' ORDER BY start");
		while (my $gene_info = $sth->fetchrow_arrayref()) {
		    #print "$gene_info->{'name'}\n";
		    my $seq = '';
		    $sth_cds->execute($gene_info->[0],$gene_info->[2], $gene_info->[3]);
		    while (my ($start,$end,$phase) = $sth_cds->fetchrow_array()) {
			$seq .= substr($gene_info->[1],$start-1,$end-$start+1);
		    }	
		    if ($gene_info->[4] eq '-') { &revcomp(\$seq); }
		    print $SEQOUT '>', $gene_info->[5],"\n";
		    if ($things_to_do{'pars_annotations'} eq 'protein') {
			print $SEQOUT &translation::translate(\$seq,$genetic_code),"\n";
		    }
		    else { print $SEQOUT $seq,"\n"; }
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
print "                                     as next argument. Sequences in the BLAST database must start with unambigous taxon name\n";
print "                                     matching keys in the hash \%categories (default: fungi and bact), or have an GeneBank accession\n";
print "                                     number. The accetion number should be in a given possition (index starting with 0; default 1) in\n";
print "                                     relation to a delimiter (| by default), e.g. position 1 in 'gb|KIL54703.1|  hypothetical...'\n";
print "                                     delimited by |. Extra settings can be given a string as an extra argument. If the string contain:\n";
print "                                     'annotation:acc', 'annotation:accno', 'annotation:accesion number' or similar will set the function\n";
print "                                     to take an accession number, and check the taxonomy in GenBank (require internet connection; 'GB:'\n";
print "                                     followed by the name of an NCBI entrez database name (default: protein), will set which database to\n";
print "                                     look in; 'delimiter:' the next character (no whitespace in between) will be used as delimiter; 'position:'\n";
print "                                     followed by an integer number (no whitespace in between) will be used to give the position; 'query:gene'\n";
print "                                     will set that it is genes (or something with an annotation in the annotations table in the database)\n";
print "                                     and not scaffolds that have been BLASTed; 'check' it will be checked if the new taxon annotation is\n";
print "                                     different from the previous; 'overwrite' previous taxon annotation will be overwritten.\n";
print "-C/--create_database             Create the database structure.\n";
print "-db/--database_name              Give the name of the database.\n";
print "--get_CDS                        Print protein coding genes (annotation type gene, including at least one streach with annotation\n";
print "                                     type CDS) as nucleotide sequences. It is possible to give extra conditions for which sequences\n";
print "                                     to print, e.g. coverage, GC, taxon, GO, GOname or GOname_space for limiting the output based on\n";
print "                                     median coverage, GC content, taxon of of scaffold that the protein is coded on, or GO number, GO\n";
print "                                     number associated with certain name, or GO number belonging to a certain name space in the gene\n";
print "                                     ontology_term. For the later two the file go.obo (http://geneontology.org/page/download-ontology)\n";
print "                                     is needed to be in the same folder as the script is running in. The type of condition can be given\n";
print "                                     as <, > (default), =, <=, >=, or LIKE followe by a value. If like is choosen % can be used as wildcard\n";
print "                                     character. Examples: 'coverage 50' or 'coverage < 10'. Several constraints can be given separated\n"; 
print "                                     by semicolon (;). The constraints will be additive (AND). It is also possible to give SQL syntax\n"; 
print "                                     directly. It will the be treated as added to the SQL statment after an AND. If you do not want the\n";
print "                                     sequences but only the number of genes, you can give COUNT as the first condition.\n";
print "--get_proteins                   Print protein coding genes (annotation type gene, including at least one streach with annotation\n";
print "                                     type CDS) as aminoacid sequences. It is possible to give extra conditions for which sequences\n";
print "                                     to print as for --get_CDS.\n";
print "-R/--get_reads_of_scaffold       Print the reads that match given scaffold. The scaffold name should be the next argument and a\n";
print "                                     bam file the next argument after that.\n";
print "-T/--get_reads_of_taxon          Print the reads that match scaffolds of given taxon. The taxon should be given as next argument,\n";
print "                                     % may be used as wildchard character, and a bam file should be the next argument after that.\n";
print "-h/--help                        Print this help.\n";
print "--include_missmatch_mates        Include reads that do not match taxon/scaffold if their mate match.\n";
print "--N/--note                       Set note annotation when reading scaffolds\n";
print "--N50                            Get statistics, including N50, for scaffolds belonging to given taxon. The taxon must be the\n";
print "                                     next argument, if the taxon is preceded by 'not ' it will give the statistics for scaffold\n";
print "                                     not belonging to the taxon. If all (default) is given the statistics for all scaffolds will\n";
print "                                     be given. If 'null' is given statistics will be provided for the scaffolds that do not have\n";
print "                                     any taxon annotation. % may be used as wildchard character.\n";
print "-a/--read_annotations            Read annotations from GFF3 formated file. Can only read annotations for scaffolds that are in\n";
print "                                     the database. File name given as next argument.\n";
print "-b/--read_bam                    Read coverage and SNPs for scaffolds in the database from a bam file given as next argument.\n";
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
    my $GCcontent = 0;
    my $length = 0;
    for (my $i=0; $i< length $$seq_ref; ++$i) {
	my $char = substr($$seq_ref, $i, 1);
	if ($char eq 'g' || $char eq 'G' || $char eq 'c' || $char eq 'C') { ++$GCcontent; ++$length; }
	elsif ($char eq 'a' || $char eq 'A' || $char eq 't' || $char eq 'T') { ++$length; }
    }
    if ($length) { $GCcontent/=$length; }
    else { $GCcontent = 0; }
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

sub revcomp {
    my $dna = shift @_; # parameters passed to the sub
    $$dna = reverse $$dna;
    $$dna =~ tr/ACGTacgt/TGCAtgca/; # translate is similar tu substitute but does the substitutions simultaneously
}

sub process_accno_search {
    my $hits = shift;
    my $delimiter = shift;
    my $position = shift;
    my $BLAST_GBdatabase = shift;
    if ($hits and @{$hits}) {
	if (${$hits}[0] eq 'no hits') {
            return 'unc';
        }
	else {
	    my $organism = 'unc';
	    for (my $i = 0; $i < scalar @{$hits}; ++$i) {
		if (!$hits->[$i]) { next; }
		my @annotations = split /$delimiter/, $hits->[$i];
		my @search;
	       	if ($BLAST_GBdatabase && $annotations[$position]) {
		    my $get= get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$BLAST_GBdatabase&id=$annotations[$position]&rettype=gb&retmode=text");
		    if ($get) { @search = split /\n/, $get; }
		    else { print STDERR "Could not get taxon annotation from http://eutils.ncbi.nlm.nih.gov for $annotations[$position] in database $BLAST_GBdatabase \n"; }
		}
		else { print "GB database and/or accession number missing.\n"; }
		my $taxon;
		foreach (@search) {
		    if (/\s*ORGANISM\s+(.+)/) { $taxon = $1; }
		    elsif ($taxon) {
			s/\s+/ /g;
			$taxon .= $_;
			if (/\.$/) {last; }
		    }
		}
		if ($taxon && $organism eq 'unc') { $organism = $taxon; last;}
	    }
	    return $organism;
	}
    }
    else { return "error"; }
}

sub process_search {
    my $hits = shift;
    my $categories = shift;
    my $critearia = shift;

    if ($hits and @{$hits}) {
        if (${$hits}[0] eq 'no hits') {
            return 'unc';
        }
        else {
            my $organism = 'unc';
            for (my $i = 0; $i < scalar @{$hits}; ++$i) {
                if (!$hits->[$i]) { next; }
                my $flag = 'n';
                foreach my $type (keys %{$categories}) {
                    if ($hits->[$i] =~ /^$type/) {
                        $flag = 'y';
                        if ($organism eq 'unc') { $organism = $type; }
                        else {
                            if ($organism ne $type) {
                                $organism = &evaluate_ambiguity($organism,$hits,$i,$critearia);
                            }
                            if ($organism eq 'ambi') { last; }
                        }
                    }
                }
                if ($flag eq 'n') { print STDERR "Do not recognize $hits->[$i] --- $i --- ", scalar @{$hits}, " --- ", keys %{$categories}, "\n"; }
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
    if ($criterion == 1) {
        return 'ambi';
    }
    elsif ($criterion == 2) {
        if ($position > 0) { return $type; }
        else { return 'ambi'; }
    }
    else { return 'ambi'; }
}

sub pars_extraSQLcondition {
    my @string = split /;/, shift;
    my $dbh = shift;
    my $return_string = '';
    foreach (@string) {
	#s/('|")/\\$1/g;
	if ($return_string) { $return_string .= " AND "; }
	if (/(coverage|GC|taxon|GOname_space|GOname|GO)(:| )*(>=|<=|==|!=|<|>|=|LIKE|NOT LIKE)*\s*(.+)/i) {
	    #print "$_\t$1 $3 $4\n";
	    my $field = uc($1);
	    my $sign;
	    if ($3) {$sign= uc($3);}
	    my $value = $4;
	    if ($value =~ /\D/ && !($field eq 'GONAME' || $field eq 'GONAME_SPACE')) { $value = "'$value'"; }
	    if (!$sign) { $sign = '>'; }
	    #print "$field $sign $value\n";
	    if ($field eq 'COVERAGE') {	$return_string .= "scaffolds.coverage_median "; }
	    elsif ($field eq 'GC') { $return_string .= "scaffolds.GCcontent "; }
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
	&revcomp(\$columns[9]);
	$columns[10] = reverse $columns[10];
	--$hash_ref->{$columns[0]}->[0]->{'rev'};
    }
    else { ++$hash_ref->{$columns[0]}->[0]->{'rev'}; }
    if ($columns[6] ne '=') {
	if (!$hash_ref->{$columns[0]}->[0]->{'mate_scaf'}) { $hash_ref->{$columns[0]}->[0]->{'mate_scaf'} = []; }
	push $hash_ref->{$columns[0]}->[0]->{'mate_scaf'}, $columns[6];
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
#}
