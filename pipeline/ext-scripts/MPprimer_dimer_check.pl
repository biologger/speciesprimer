#!/usr/bin/perl -w

# Author: Wubin Qu
# E-mail: quwubin@gmail.com
# Date: 2009-09-30

use Getopt::Std;
use List::Util qw[min max];

my %opt = ();
my $primer_file = '';
my $usage = "Usage: $0 [-h] [-f primer_file]

\tOptional:

\t-o oligo_conc, default 50 nM;
\t-d divalent_cation_conc, default 0 mM;
\t-m monovalent_cation_conc, default 50 mM;
\t-n dntp_conc, default 0 mM\n";

# Starting ionic concentration variables
my $oligo_conc = 50; #in nM
my $mg_conc=0; #in mM
my $monovalent_cation_conc=50; #in mM
my $dntp_conc=0; #in mM

if (@ARGV <= 0){
    print $usage;
    exit;
}else{
    getopts('hf:r:o:d:m:n:',\%opt); 
    die $usage if ($opt{h});
    $primer_file = $opt{f} if $opt{f};
    $oligo_conc = $opt{o} if $opt{o};
    $mg_conc = $opt{d} if $opt{d};
    $monovalent_cation_conc = $opt{m} if $opt{m};
    $dntp_conc = $opt{n} if $opt{n};
    die $usage if ((length $primer_file) == 0);
}

my (%oligo_dH, %oligo_dH_full, %oligo_dS, %oligo_dS_full, %genetic_code);
load_data();

# Primer-dimer parameters
my $pd_extensible=1;
my $pd_temperature=37;

my %oligo_dG=(
	qw(initC 0.98 	initG 0.98 
	initA 1.03 	initT 1.03), 
);
recalculate_dG();

my $primers = parse_primer_file($primer_file);

my @primer_key = keys %{$primers};

for (my $i = 0; $i < @primer_key - 1; $i++){
    for (my $j = $i+1; $j < @primer_key; $j++){
        my $primer_i = ${$primers}{$primer_key[$i]};
        my $primer_j = ${$primers}{$primer_key[$j]};
        my $min_dG = max_stable($primer_i, $primer_j);
        print "$primer_key[$i]\t$primer_key[$j]\t$min_dG\n";
    }
}

sub max_stable{
    my ($fp, $rp) = @_;
    my $rating_hash = '';
    my $align_chart = '';

    (my $dG_fr, $rating_hash) = primer_dimer($fp, $rp, 1);
    #(my $dG_ff, $rating_hash) = primer_dimer($fp, $fp, 1);
    #(my $dG_rr, $rating_hash) = primer_dimer($rp, $rp, 1);

    #return min($dG_fr, $dG_ff, $dG_rr);
    return $dG_fr;
}

sub max_stable_test{
    my ($fp, $rp) = @_;
    my $rating_hash = '';
    my $align_chart = '';
    my $cutoff = 0;

    (my $dG_fr, $rating_hash) = primer_dimer($fp, $rp, 1);
    #print "FR: $dG_fr kcal/mol\n\n";
    $align_chart = run_draw_dimer($fp, $rp, $dG_fr, $rating_hash);
    if ($dG_fr >= $cutoff ){
        print "$dG_fr: $fp\t$rp\n$align_chart";
    }

    (my $dG_ff, $rating_hash) = primer_dimer($fp, $fp, 1);
    #print "FF: $dG_ff kcal/mol\n\n";
    $align_chart = run_draw_dimer($fp, $fp, $dG_ff, $rating_hash);
    if ($dG_ff >= $cutoff){
        print "$dG_ff: $fp\t$fp\n$align_chart";
    }

    (my $dG_rr, $rating_hash) = primer_dimer($rp, $rp, 1);
    #print "RR: $dG_rr kcal/mol\n\n";
    $align_chart = run_draw_dimer($rp, $rp, $dG_rr, $rating_hash);
    if ($dG_rr >= $cutoff){
        print "$dG_rr: $rp\t$rp\n$align_chart";
    }

    return min($dG_fr, $dG_ff, $dG_rr);
}

sub run_draw_dimer{
    my ($fp, $rp, $dG, $rating_hash) = @_;
    if ($dG >= 0){
        return "No match found.\n";
    }else{
        my $align_chart = draw_dimer($fp, $rp, $rating_hash{$dG});
        return $align_chart;
    }
}

sub parse_primer_file{
    my ($primer_file) = @_;
    my %primers = ();
    open (FL, "$primer_file");
    while (my $line = <FL>){
        chomp $line;
        if ($line =~ /^>/){
            my $id = substr($line, 1);
            my $seq = <FL>;
            chomp $seq;
            $primers{$id} = $seq;
        }
    }
    return \%primers;
}

#--------------------------#
# Primer-dimer calculation #
#--------------------------#

sub primer_dimer {
	# This is my second attempt at implementing a primer-dimer system:
	# The first attempt aligned each combination together explicitly; this attempt
	# creates a binding matrix and then reads each primer combination from the
	# matrix.  It's not significantly faster, but it does have the advantage of
	# extending this algorithm to allow for unequal loops (although this has not
	# been implemented as yet) and providing the possiblity of reading all
	# primer-dimers (not just 3') by adjusting a single variable (which is used when
	# displaying primer information.
		
	my ($primer_f, $primer_r, $pd_full) = @_;
	return unless ($primer_f) && ($primer_r);
			
	my ($k, $l);
	@score=();
	%primer_hash=();
	@score_sort=();
	@bind_string=();
	%rating_hash=();
		
	# $pl = greatest length
	$pfl=length($primer_f);
	$prl=length($primer_r);
	$pl = ($pfl>$prl ? $pfl : $prl);
	
	my $rcompr = reverse(complement($primer_r));
	my $rcomprlc = lc($rcompr);
	my $fprimer_r=lc(reverse($primer_f));
	$rprimer_r=reverse($primer_r);
	
	# Scan the primers against each other:
	# The default is to only consider 3' primer-dimers, for speed concerns (5'
	# pd's will reduce the primer population, but won't cause extendible dimers) -
	# however, setting $pd_full to 1 will calculate all primer-dimers.  This is
	# currently used only when viewing individual primers, for both speed concerns
	# and because it's 3' primer-dimers that are the real problem in PCR.
	
	# create a binding array for each of the four bases
	for $l (0 .. $pfl-1) {
		my $mbase = substr($fprimer_r, $l, 1);
		$primer_hash{$mbase}[$l]=1;
		for $k (qw(a g c t)) {
			$primer_hash{$k}[$l] ||=0;
		}
	}
		
	# create the primer matrix
	my @primer_comp;
	for $k (0 .. $prl-1) {
		$primer_comp[$k]=$primer_hash{substr($rcomprlc, $k, 1)};
	}
		
	# read each combination from the matrix, calculate dG for each dimer
	my $pd_len = ($pd_full ? $pfl+$prl-1 : $pl-2);
	for $k (0 .. $pd_len) {
		$score[$k]=0;
		my $bind;
		my $score_p=0;
		
		# extensible primer short-circuit - ignore all primers that will
		# not create extensible (i.e. amplifiable) dimers
		my $start = $k>$pfl-1 ? $pfl-1 : $k;
		my $end = $k>$prl-1 ? $prl-1 : $k;
		if ($pd_extensible && !$pd_full) {
			next unless $primer_comp[0][$start] == 1;
			next unless $primer_comp[$end][$start-$k] == 1;
		}
		
		# elsif ($pd_extensible) {
			# # no point reconsidering them!
			# next if $primer_comp[0][$start] == 1 && $primer_comp[$end][$start-$k] == 1;
		# }
		
		# read the binding data
		for $l (0 .. $prl-1) {
			if (($k-$l)<$pfl) {
				$bind .= $primer_comp[$l][$k-$l] if ($k-$l)>=0;
			} else {
				# spacer
				$bind .= "2";
			}
		}
		
		# Single matched bases surrounded by mismatches are unstable,
		# so we remove them with the regexp (look ahead is needed otherwise
		# strings of consecutive match/mismatches are not caught)
		$bind =~ s/01(?=[^1])/00/gx;
		
		# Short circuit if there's nothing to bind
		next unless $bind =~ /[1]/;
		
		# Find start and end of similarity
		my ($pb_init,$pb_end);
		for $l (0 .. length($bind)-1) {
			# at first I tried finding the initiating terminal bases with
			# regexps, but that was much slower ...
			if (substr($bind, $l, 1) eq "1") {
				defined($pb_init) || ($pb_init = $l);
				$pb_end=$l;
			}
		}
				
		if (defined($pb_init)) {
			# deltaG calculation
			for $l ($pb_init .. $pb_end-1) {
				next if substr($bind, $l, 2) eq "00";
				next if substr($bind, $l, 1) eq "2";
				$score_p+=$oligo_dG{substr($primer_f, $pfl-$k+$l-1, 2).substr($rprimer_r, $l, 2)};
			}
			
			# init term corrections
			my $initterm="init" . substr($rprimer_r, $pb_init, 1);
			$score_p+= $oligo_dG{$initterm};
			
			my $endterm="init" . substr($rprimer_r, $pb_end, 1);
			$score_p+= $oligo_dG{$endterm};
			
			# add to the hash ...
			$score[$k]=sprintf("%.2f",$score_p);
			$bind_string[$k]=$bind;
			$rating_hash{$score[$k]}=$k;
		}
	}
	
	# sort the dimers to give the most stable:	
	@score_sort = sort { $a <=> $b } @score;
		
	# Returns the most stable dimer
	return ($score_sort[0], \%rating_hash);
}


sub draw_dimer {
	# This all seems a bit cumbersome!!
	#my ($primer_f, $primer_r, $pos, $FH) = @_;
	my ($primer_f, $primer_r, $pos) = @_;
	
	my $rprimer_r=reverse($primer_r);
	my $dimer_binding="";
	my $pr_space="";
	my $fspace="";
	my $rspace="";
			
	my $fspace_def = $pl-$pfl>0 ? $pl-$pfl : 0;
	$fspace=" "x($fspace_def+($pos>$pl-1?$pos-$pl+1:0));
	
	if ($pos+1>=$pfl) {
		$rspace=" "x($pl-$pos-1);
	} else {
		$rspace=$fspace;
	}
	
	$pr_space=" "x($pfl-$pos-1);
	
	for my $j (0 .. $pos) {
		next unless $j < $prl;
		if (substr($bind_string[$pos],$j,1)==1) {
			$dimer_binding=$dimer_binding."|"
		} elsif (substr($bind_string[$pos],$j,1)==0) {
			$dimer_binding=$dimer_binding."."
		} else {
			$dimer_binding=$dimer_binding." "
		}
	}
				
	#print $FH "$fspace"."5' "."$primer_f"." 3'\n".
	#print   "$fspace"."5' "."$primer_f"." 3'\n".
	my $align_chart = "$fspace"."5' "."$primer_f"." 3'\n".
		"$rspace"."   "."$pr_space"."$dimer_binding\n".
		"$rspace"."$pr_space"."3' "."$rprimer_r"." 5'\n\n";

        return $align_chart;
}	


sub complement {
	$_ = shift;
	tr/AGCTagct/TCGAtcga/;
	return $_;
}

#----------------------------------------------------------#
# Recalculate the oligo_dG hash on current salt conditions #
#----------------------------------------------------------#

sub recalculate_dG {
	# because dG = dH - TdS, and dS is dependent on the salt concentration ...
	my $temperature = shift || $pd_temperature;
	
	# Big problems if [dNTPs] > [Mg++] !!  This is a quick fix ...
	my $salt_correction;
	if ($mg_conc > $dntp_conc) {
		$salt_correction = sqrt($mg_conc - $dntp_conc);
	} else {
		$salt_correction = 0;
	}
	
	my $na_eq=($monovalent_cation_conc + 120 * $salt_correction)/1000;
	
	# the length of each NN dimer is 2, therefore the modifier is 1
	my $entropy_adjust = (0.368 * log($na_eq));
		
	foreach my $key (keys(%oligo_dH_full)) {
		next if $key =~ /init/; # the length of each monomer is 1, thus the modifier of dS is 0 and the values are precalulated
		
		my $dS = $oligo_dS_full{$key} + $entropy_adjust;
		my $dG = $oligo_dH_full{$key}-((273.15+$temperature)*($dS/1000));
		$oligo_dG{$key} = $dG;
	}
}

#--------------------------------#
# Rountine to draw primer-dimers #
#--------------------------------#
sub check_degenerate {
	$_ = shift;
	if (/[^ATGC]/i) {
		dialogue("One of your primer sequences has a degenerate or non-DNA character.  PerlPrimer cannot calculate the Tm of degenerate sequences") if shift;
		return 1;
	} 
}




sub load_data {
	#-----
	#
	# NN thermodynamics hashes (AA = 5' AA 3'/3' TT 5') derived from ...
	# 
	# Allawi HT, SantaLucia J Jr.  Thermodynamics and NMR of internal G.T mismatches in DNA.
	# 	Biochemistry. 1997 Aug 26;36(34):10581-94
	#
	# SantaLucia J Jr.  A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics.
	# 	Proc Natl Acad Sci U S A. 1998 Feb 17;95(4):1460-5. 
	# 
	# ... with mismatch dG data (AGTG = 5' AG 3'/3' TG 5') derived from ...
	# 
	# Peyret N, Seneviratne PA, Allawi HT, SantaLucia J Jr.  Nearest-neighbor thermodynamics and NMR of DNA sequences with internal A.A, C.C, G.G, and T.T mismatches.
	# 	Biochemistry. 1999 Mar 23;38(12):3468-77. 
	# 
	# Allawi HT, SantaLucia J Jr.  Nearest-neighbor thermodynamics of internal A.C mismatches in DNA: sequence dependence and pH effects.
	# 	Biochemistry. 1998 Jun 30;37(26):9435-44.
	# 
	# Allawi HT, SantaLucia J Jr.  Thermodynamics of internal C.T mismatches in DNA.
	# 	Nucleic Acids Res. 1998 Jun 1;26(11):2694-701. 
	# 
	# Allawi HT, SantaLucia J Jr.  Nearest neighbor thermodynamic parameters for internal G.A mismatches in DNA.
	# 	Biochemistry. 1998 Feb 24;37(8):2170-9.
	# 
	# Allawi HT, SantaLucia J Jr.  Thermodynamics and NMR of internal G.T mismatches in DNA.
	# 	Biochemistry. 1997 Aug 26;36(34):10581-94
	# 
	#-----
	
	#-------------------#
	# deltaH (kcal/mol) #
	#-------------------#
	
	%oligo_dH=qw(
		AA -7.9 TT -7.9 
		AT -7.2 TA -7.2 
		CA -8.5 TG -8.5 
		GT -8.4 AC -8.4 
		CT -7.8 AG -7.8 
		GA -8.2 TC -8.2 
		CG -10.6 GC -9.8 
		GG -8.0 CC -8.0 
		initC 0.1 initG 0.1 
		initA 2.3 initT 2.3
	);
	
	%oligo_dH_full=(
		qw(AATT -7.9 	TTAA -7.9 
		ATTA -7.2 	TAAT -7.2 
		CAGT -8.5 	TGAC -8.5 
		GTCA -8.4 	ACTG -8.4 
		CTGA -7.8 	AGTC -7.8 
		GACT -8.2 	TCAG -8.2 
		CGGC -10.6 	GCCG -9.8 
		GGCC -8.0 	CCGG -8.0
			
		initC 0.1 	initG 0.1 
		initA 2.3 	initT 2.3),
		
		# Like pair mismatches 
			
		qw(AATA 1.2 	ATAA 1.2
		CAGA -0.9 	AGAC -0.9
		GACA -2.9 	ACAG -2.9
		TAAA 4.7 	AAAT 4.7 
		
		ACTC 0.0 	CTCA 0.0 
		CCGC -1.5 	CGCC -1.5
		GCCC 3.6 	CCCG 3.6 
		TCAC 6.1 	CACT 6.1 
		
		AGTG -3.1 	GTGA -3.1
		CGGG -4.9 	GGGC -4.9
		GGCG -6.0 	GCGG -6.0
		TGAG 1.6 	GAGT 1.6 
		
		ATTT -2.7 	TTTA -2.7
		CTGT -5.0 	TGTC -5.0
		GTCT -2.2 	TCTG -2.2
		TTAT 0.2 	TATT 0.2  ),
		
		# G.T mismatches 
		
		qw(AGTT 1.0  	TTGA 1.0
		ATTG  -2.5 	GTTA  -2.5
		CGGT  -4.1 	TGGC  -4.1
		CTGG  -2.8 	GGTC  -2.8
		GGCT  3.3 	TCGG  3.3
		GGTT  5.8 	TTGG  5.8
		GTCG  -4.4 	GCTG  -4.4
		GTTG  4.1 	GTTG  4.1
		TGAT  -0.1 	TAGT  -0.1
		TGGT  -1.4 	TGGT  -1.4
		TTAG  -1.3 	GATT  -1.3), 
		
		# G.A mismatches 
		
		qw(AATG  -0.6 	GTAA  -0.6
		AGTA  -0.7 	ATGA  -0.7
		CAGG  -0.7 	GGAC  -0.7
		CGGA  -4.0 	AGGC  -4.0
		GACG  -0.6 	GCAG  -0.6
		GGCA  0.5 	ACGG  0.5
		TAAG  0.7 	GAAT  0.7
		TGAA  3.0 	AAGT  3.0), 
		
		# C.T mismatches 
		
		qw(ACTT  0.7 	TTCA  0.7
		ATTC  -1.2 	CTTA  -1.2
		CCGT  -0.8 	TGCC  -0.8
		CTGC  -1.5 	CGTC  -1.5
		GCCT  2.3 	TCCG  2.3 
		GTCC  5.2 	CCTG  5.2 
		TCAT  1.2 	TACT  1.2 
		TTAC  1.0 	CATT  1.0), 
		
		# A.C mismatches 
		
		qw(AATC  2.3	CTAA  2.3
		ACTA  5.3 	ATCA  5.3 
		CAGC  1.9 	CGAC  1.9 
		CCGA  0.6 	AGCC  0.6 
		GACC  5.2 	CCAG  5.2 
		GCCA  -0.7 	ACCG  -0.7
		TAAC  3.4  	CAAT  3.4 
		TCAA  7.6 	AACT  7.6),
	
	);
	
	#--------------------#
	# deltaS (cal/K.mol) #
	#--------------------#
	
	%oligo_dS=qw(
		AA -22.2 TT -22.2 
		AT -20.4 TA -21.3 
		CA -22.7 TG -22.7 
		GT -22.4 AC -22.4 
		CT -21.0 AG -21.0 
		GA -22.2 TC -22.2 
		CG -27.2 GC -24.4 
		GG -19.9 CC -19.9 
		initC -2.8 initG -2.8 
		initA 4.1 initT 4.1 
		sym -1.4
	);
	
	%oligo_dS_full=(
		qw(AATT -22.2 	TTAA -22.2 
		ATTA -20.4 	TAAT -21.3 
		CAGT -22.7 	TGAC -22.7 
		GTCA -22.4 	ACTG -22.4 
		CTGA -21.0 	AGTC -21.0 
		GACT -22.2 	TCAG -22.2 
		CGGC -27.2 	GCCG -24.4 
		GGCC -19.9 	CCGG -19.9
			
		initC -2.8 	initG -2.8 
		initA 4.1 	initT 4.1
		sym -1.4),
		
		# Like pair mismatches
			
		qw(AATA 1.7 	ATAA 1.7
		CAGA -4.2 	AGAC -4.2 
		GACA -9.8 	ACAG -9.8 
		TAAA 12.9 	AAAT 12.9 
		
		ACTC -4.4 	CTCA -4.4 
		CCGC -7.2 	CGCC -7.2 
		GCCC 8.9 	CCCG 8.9 
		TCAC 16.4 	CACT 16.4 
		
		AGTG -9.5 	GTGA -9.5 
		CGGG -15.3 	GGGC -15.3
		GGCG -15.8 	GCGG -15.8
		TGAG 3.6 	GAGT 3.6 
		
		ATTT -10.8 	TTTA -10.8
		CTGT -15.8 	TGTC -15.8
		GTCT -8.4 	TCTG -8.4 
		TTAT -1.5 	TATT -1.5),
		
		# G.T mismatches
		
		qw(AGTT 0.9 	TTGA 0.9
		ATTG  -8.3 	GTTA  -8.3
		CGGT  -11.7 	TGGC  -11.7
		CTGG  -8.0 	GGTC  -8.0
		GGCT  10.4 	TCGG  10.4
		GGTT  16.3 	TTGG  16.3
		GTCG  -12.3 	GCTG  -12.3
		GTTG  9.5 	GTTG  9.5
		TGAT  -1.7 	TAGT  -1.7
		TGGT  -6.2 	TGGT  -6.2
		TTAG  -5.3 	GATT  -5.3), 
		
		# G.A mismatches
		
		qw(AATG  -2.3 	GTAA  -2.3
		AGTA  -2.3 	ATGA  -2.3
		CAGG  -2.3 	GGAC  -2.3
		CGGA  -13.2 	AGGC  -13.2
		GACG  -1.0 	GCAG  -1.0
		GGCA  3.2 	ACGG  3.2
		TAAG  0.7 	GAAT  0.7
		TGAA  7.4 	AAGT  7.4), 
		
		# C.T mismatches
		
		qw(ACTT  0.2 	TTCA  0.2
		ATTC  -6.2 	CTTA  -6.2
		CCGT  -4.5 	TGCC  -4.5
		CTGC  -6.1 	CGTC  -6.1
		GCCT  5.4 	TCCG  5.4 
		GTCC  13.5 	CCTG  13.5
		TCAT  0.7 	TACT  0.7 
		TTAC  0.7 	CATT  0.7), 
		
		# A.C mismatches
		
		qw(AATC  4.6 	CTAA  4.6
		ACTA  14.6 	ATCA  14.6
		CAGC  3.7 	CGAC  3.7 
		CCGA  -0.6 	AGCC  -0.6
		GACC  14.2 	CCAG  14.2
		GCCA  -3.8 	ACCG  -3.8
		TAAC  8.0  	CAAT  8.0 
		TCAA  20.2 	AACT  20.2),
	
	);
	
	
	# Genetic code hash
	%genetic_code=qw(
			TTT F TTC F TTA L TTG L
			CTT L CTC L CTA L CTG L
			ATT I ATC I ATA I ATG M
			GTT V GTC V GTA V GTG V
			TCT S TCC S TCA S TCG S
			CCT P CCC P CCA P CCG P
			ACT T ACC T ACA T ACG T
			GCT A GCC A GCA A GCG A
			TAT Y TAC Y TAA * TAG *
			CAT H CAC H CAA Q CAG Q
			AAT N AAC N AAA K AAG K
			GAT D GAC D GAA E GAG E
			TGT C TGC C TGA * TGG W
			CGT R CGC R CGA R CGG R
			AGT S AGC S AGA R AGG R
			GGT G GGC G GGA G GGG G
	);
}



