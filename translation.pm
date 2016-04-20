package translation;
use strict;

my %code = ( standard => {'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG'=>'A',
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

sub code_available {
    my $test = shift;
    if ($code{$test}) { return 1; }
    else { return 0; }
}

sub available_codes_string {
    return join ', ', keys %code;
}

sub translate {
    my $dna = shift @_;
    my $used_code = shift @_;
    my $start=0;
    my $seq = '';
    my $length = length $$dna;
    for (my $start=0; $start+2 < $length; $start+=3) {
        if ($code{$used_code}) {
            my $aa = $code{$used_code}->{uc(substr($$dna,$start,3))};
            if ($aa) { $seq .= $aa; }
	    else { print STDERR "No matching aa for: ", uc(substr($$dna,$start,3)), ".\n"; }
        }
        else { die "I do not know the code $used_code.\n"; }
#if ($code eq 'bacterial') {
#            my $aa .= $code{bacterial}->{uc(substr($$dna,$start,3))};
#	    if ($aa) { $seq .= $aa; }
#	    else { print STDERR "No matching aa for: ", uc(substr($$dna,$start,3)), ".\n"; }
#        }
    }
    return $seq;
}


1;
