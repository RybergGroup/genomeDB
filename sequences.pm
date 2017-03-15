package sequences;
use strict;
use LWP::Simple;

my %genecode = ( standard => {'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG'=>'A',
		'TTA' => 'L', 'TTG' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG'=> 'L',
		'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R', 'AGG' => 'R',
		'AAA' => 'K', 'AAG' => 'K',
		'AAT' => 'N', 'AAC' => 'N',
		'ATG' => 'M',
		'GAT' => 'D', 'GAC' => 'D',
		'TTT' => 'F', 'TTC' => 'F',
		'TGT' => 'C', 'TGC' => 'C',
		'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
		'CAA' => 'Q', 'CAG' => 'Q',
		'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S',
		'GAA' => 'E', 'GAG' => 'E',
		'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
		'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G',
		'TGG' => 'W',
		'CAT' => 'H', 'CAC' => 'H',
		'TAT' => 'Y', 'TAC' => 'Y',
		'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I',
		'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',
		'TAA' => '*', 'TGA' => '*', 'TAG' => '*'},
            bacterial => {'TTT'=>'F', 'TCT'=>'S', 'TAT'=>'Y', 'TGT'=>'C',
		'TTC'=>'F', 'TCC'=>'S', 'TAC'=>'Y', 'TGC'=>'C',
		'TTA'=>'L', 'TCA'=>'S', 'TAA'=>'*', 'TGA'=>'*',
		'TTG'=>'L', 'TCG'=>'S', 'TAG'=>'*', 'TGG'=>'W',
		'CTT'=>'L', 'CCT'=>'P', 'CAT'=>'H', 'CGT'=>'R',
		'CTC'=>'L', 'CCC'=>'P', 'CAC'=>'H', 'CGC'=>'R',
		'CTA'=>'L', 'CCA'=>'P', 'CAA'=>'Q', 'CGA'=>'R',
		'CTG'=>'L', 'CCG'=>'P', 'CAG'=>'Q', 'CGG'=>'R',
		'ATT'=>'I', 'ACT'=>'T', 'AAT'=>'N', 'AGT'=>'S',
		'ATC'=>'I', 'ACC'=>'T', 'AAC'=>'N', 'AGC'=>'S',
		'ATA'=>'I', 'ACA'=>'T', 'AAA'=>'K', 'AGA'=>'R',
		'ATG'=>'M', 'ACG'=>'T', 'AAG'=>'K', 'AGG'=>'R',
		'GTT'=>'V', 'GCT'=>'A', 'GAT'=>'D', 'GGT'=>'G',
		'GTC'=>'V', 'GCC'=>'A', 'GAC'=>'D', 'GGC'=>'G',
		'GTA'=>'V', 'GCA'=>'A', 'GAA'=>'E', 'GGA'=>'G',
		'GTG'=>'V', 'GCG'=>'A', 'GAG'=>'E', 'GGG'=>'G'}
);

sub gene_codes_available {
    my $test = shift;
    if ($genecode{$test}) { return 1; }
    else { return 0; }
}

sub available_gene_codes_string {
    return join ', ', keys %genecode;
}

sub translateDNA {
    my $dna = shift @_;
    my $used_code = shift @_;
    my $start=0;
    my $seq = '';
    my $length = length $$dna;
    for (my $start=0; $start+2 < $length; $start+=3) {
        if ($genecode{$used_code}) {
            my $aa = $genecode{$used_code}->{uc(substr($$dna,$start,3))};
            if ($aa) { $seq .= $aa; }
	    else { print STDERR "No matching aa for: ", uc(substr($$dna,$start,3)), ".\n"; }
        }
        else { die "I do not know the code $used_code.\n"; }
    }
    return $seq;
}

sub revcomp {
    my $dna = shift @_; # parameters passed to the sub
    $$dna = reverse $$dna;
    $$dna =~ tr/ACGTacgt/TGCAtgca/; # translate is similar tu substitute but does the substitutions simultaneously
}

sub GCcontent {
    my $seq_ref = shift;
    my $GCcontent = 0;
    my $length = 0;
    for (my $i=0; $i< length $$seq_ref; ++$i) {
        my $char = substr($$seq_ref, $i, 1);
        if ($char eq 'g' || $char eq 'G' || $char eq 'c' || $char eq 'C') { ++$GCcontent; ++$length; }
        elsif ($char eq 'a' || $char eq 'A' || $char eq 't' || $char eq 'T') { ++$length; }
    }
    if ($length) { $GCcontent/=$length; }
    else {
	$GCcontent = 0;
	#print "$$seq_ref - $GCcontent - $length\n";
    }
    return $GCcontent;
}

sub get_entry_from_GB {
    my $database = pop; # last argument is database name
    my $accnos = shift; # first argument is comma separated string of accnos
    foreach (@_) { $accnos .= ",$_"; }
    #print "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$database&id=$accnos&rettype=gb&retmode=text\n"; # for debugging
    return get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$database&id=$accnos&rettype=gb&retmode=text");
}

1;
