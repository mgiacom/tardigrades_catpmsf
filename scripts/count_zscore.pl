use strict; 
use warnings;


my @files = <*gaps>;
my $original_dataset = $ARGV[0];
open (OUT, ">diversity.pbr_gaps");
open (OUT1, ">diversity_scores_bootstrapped_data.txt_gaps");

#### MAIN
#### calculate diversity main dataset.
my %original_data = ReadFasta($original_dataset);
my @original_taxa = keys(%original_data);
my $original_ntaxa = scalar(@original_taxa);
my $original_ncharacter= length($original_data{$original_taxa[0]});
my @original_aligemnt = MakeAligment(\@original_taxa, \%original_data);
my @original_characters = CharMatrix($original_ncharacter, $original_ntaxa, \@original_aligemnt);
my $original_data_AA_diversity = AASiteDiversity(\@original_characters, $original_ncharacter);  ## This store AA diversity original data
#### DONE

#### CALCULATE diversity bootstrapped data
my %para_boot_diversities;
for my $infile (@files)
{
    print "processing $infile \n";
    my %data = ReadFasta ($infile);

    ### taxa and characters
    my @taxa = keys (%data);
    my $ntaxa = scalar(@taxa);
    my $ncharacter= length($data{$taxa[0]});
    my @aligemnt = MakeAligment(\@taxa, \%data);
    my @characters = CharMatrix($ncharacter, $ntaxa, \@aligemnt);

    $para_boot_diversities{$infile} = AASiteDiversity(\@characters, $ncharacter);
}
#### DONE

#### calculate average diversity bootstrapped data and SD
my $total_diversity_simulated_data=0;
my @simulated_keys = keys(%para_boot_diversities);
my $number_simulations= scalar(@simulated_keys);

for my $simulation (@simulated_keys)
{
    $total_diversity_simulated_data += $para_boot_diversities{$simulation};
}
my $average_simulated_data = $total_diversity_simulated_data/$number_simulations;
#### DONE

#### Calculate SD simulated data
my $sum_of_squares=0;
for my $simulation (@simulated_keys)
{
    $sum_of_squares += ($para_boot_diversities{$simulation} - $average_simulated_data)**2
}
my $variance = $sum_of_squares / ($number_simulations-1);
my $standard_deviation_simulated_data = sqrt($variance);
#### Done

### calculate Z score
my $standard_deviate_dataset= ($original_data_AA_diversity - $average_simulated_data)/$standard_deviation_simulated_data;


#### generating ouoputs
print OUT "Test of model adequacy / pattern heterogeneity / across sites compositional heterogeneity\n";
print OUT "Diversity original data: " . $original_data_AA_diversity . "\n";
print OUT "Average diversity simulated data: " . $average_simulated_data . "\n";
print OUT "SD simulated data: " . $standard_deviation_simulated_data . "\n";
print OUT "Z-score: " . $standard_deviate_dataset . "\n";
    

my @keys_div_hash = keys(%para_boot_diversities);
print OUT1 "file \t diversity\n";
for my $file (@keys_div_hash)
{
    print OUT1 $file . "\t" . $para_boot_diversities{$file} . "\n";
}
####
    
    

###### SUBS ########
sub ReadFasta{
    my $fasta_file = shift @_;     
    open (FASTA, "<$fasta_file") || die "cannot open $fasta_file";
    
    my %seqs;
    my $current_title = "";
    my $current_seq = "";
    while (<FASTA>)
    {
	chomp;
	$_ =~ s/ //g; #remove white spases either in the name or in the sequence when sites are griuped by 5 for example
	if (m/^\s*$/) # remove empty lines
	{
	    next;
	}
	if (m/^>/) 
	{ 
	    if ($current_title ne "") 
	    { 
		$seqs{$current_title} = $current_seq;
		$current_seq = "";
	    }
	    ($current_title = $_) =~ s/>//;
	} 
	else {
	    $current_seq .= $_;
	}
    }
    $seqs{$current_title} = $current_seq;
    close FASTA;
    return %seqs;
}

sub MakeAligment{ # this just make an alignment as an array or arrays ordered following taxon array order 
 ## there are two parameters dataset (fasta) as a ref has and the list of taxa (from @taxa) // could be made better by passing only hash and calculating all here
    my @local_taxa_list = @{$_[0]};
    my %local_sequences = %{$_[1]};
    my @returned_aligemnt;
    for my $taxon (@local_taxa_list)
    {
	my @local_array_characters = split (//, $local_sequences{$taxon});
	push (@returned_aligemnt, \@local_array_characters);
    }
    return @returned_aligemnt;
}

sub CharMatrix { ## this just make an character matrix where each character is an array ordered by taxon list (again array or array - trasposition of Alignment function
    my $local_nchar= $_[0];
    my $local_ntaxa= $_[1];
    my @local_alignment = @{$_[2]};
    my @local_characters;
    
    for (my $character=0; $character < $local_nchar; $character++)
    {
	my @local_character;
	for (my $taxon =0; $taxon < $local_ntaxa; $taxon++)
	{
	    push (@local_character, $local_alignment[$taxon][$character]);
	}
	push (@local_characters, \@local_character);
    }
    return @local_characters;
}

sub AASiteDiversity {
    my @local_matrix = @{$_[0]};
    my $alignment_lenght= $_[1];
    
    my @array_of_aa_per_site;
    for my $site (@local_matrix)
    {
        my %aa_diversity_hash;
        for my $state (@{$site})
        {
            if ($state eq "-")
            {next;}
            if ($state eq "X")
            {next;}
            if ($state eq "?")
            {next;}
            $aa_diversity_hash{$state}=1;
        }
        push (@array_of_aa_per_site, scalar(keys(%aa_diversity_hash)));
    }

    my $total_aa_div=0;
    for my $site (@array_of_aa_per_site)
    {
    $total_aa_div += $site;
    }
    my $average_diversity_per_site = $total_aa_div/$alignment_lenght;
    return  $average_diversity_per_site;
}
