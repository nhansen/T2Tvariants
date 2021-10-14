#!/usr/bin/perl -w
#

use strict;

my $Usage = "Usage: pull_superpop_frequencies.pl <region> <ref>\nref is either chm13v1.0 or grch38\n";
$#ARGV==1
    or die $Usage;

my $region = $ARGV[0];
my $ref = lc($ARGV[1]);

my $samples_pops = '1000G_2504_high_coverage.sequence.index';
my $superpops = 'superpopulations.txt';

my $rh_sample_pops = read_sample_pops($samples_pops);
my $rh_superpops = read_superpops($superpops);

# print totals for each pop:
my %pop_samples = ();
my %superpop_samples = ();
foreach my $sample (keys %{$rh_sample_pops}) {
    my $pop = $rh_sample_pops->{$sample};
    my $superpop = $rh_superpops->{$pop};
    $pop_samples{$pop}++;
    $superpop_samples{$superpop}++;
}
foreach my $pop (keys %pop_samples) {
    my $no_samples = $pop_samples{$pop};
    #print "$pop\t$no_samples\n";
}
foreach my $superpop (keys %superpop_samples) {
    my $no_samples = $superpop_samples{$superpop};
    #print "$superpop\t$no_samples\n";
}

my $chrom = $region;
$chrom =~ s/:.*//;

# obtain variants/genotypes from region:

my $vcf = ($ref eq 'grch38') ? "/data/nhansen/multisamplevcfs/GRCh38_withgenos/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_$chrom.recalibrated_variants.pass.vcf.gz" :
                               "/data/nhansen/multisamplevcfs/chm13v1.0_withgenos/1kgp.$chrom.recalibrated.snp_indel.pass.vcf.gz";

my @vcf_field_headers = grab_field_headers($vcf);
check_field_headers([@vcf_field_headers], $rh_sample_pops);

my $vcf_command = "tabix $vcf $region | ";

#print "$vcf_command\n";

open GENOS, $vcf_command
    or die "Couldn\'t execute $vcf_command: $!\n";

my $first_line = 1;
while (<GENOS>) {
    #print;
    chomp;
    my @fields = split /\t/, $_;
    my %sp_ac = (); # by superpop, then allele
    my %sp_ac_hom = (); # by superpop, then allele
    my %sp_ac_het = (); # by superpop, then allele
    my %sp_an = (); # by superpop
    my $alt_field = $fields[4];
    my $no_alt_alleles = split /,/, $alt_field;
    my $info_field = $fields[7];

    my $format = $fields[8];
    if ($format !~ /^GT/) {
        die "Unexpected FORMAT field order in file $vcf: $format\n";
    }
    for (my $i=9; $i<=$#fields; $i++) {
        my $fieldname = $vcf_field_headers[$i]; 
        if ($rh_sample_pops->{$fieldname}) {
            my $pop = $rh_sample_pops->{$fieldname};
            my $superpop = $rh_superpops->{$pop};
            $sp_an{$superpop} = $sp_an{$superpop} || 0;
            my $geno = $fields[$i];
            $geno =~ s/:.*//;
            if ($geno =~ /\./) {
                #print "Skipping sample $fieldname with genotype $geno\n";
                next;
            }
            
            my @alleles = split m:[/|]:, $geno;
            my $hom = (($#alleles==1) && ($alleles[0] == $alleles[1])) ? 1 : 0;
            my $het = (($#alleles==1) && ($alleles[0] != $alleles[1])) ? 1 : 0;
            $sp_ac_hom{$superpop}->[$alleles[0]]+= 2 if ($hom);
            $sp_ac_het{$superpop}->[$alleles[0]]++ if ($het);
            $sp_ac_het{$superpop}->[$alleles[1]]++ if ($het);

            foreach my $allele (@alleles) {
                $sp_ac{$superpop}->[$allele]++;
                $sp_an{$superpop}++;
            }
        }
    }
    # write some totals:
    my @ac_info_additions = (); # extra info fields
    my @ac_hom_info_additions = (); # extra info fields
    my @ac_het_info_additions = (); # extra info fields
    my @af_info_additions = (); # extra info fields
    my @an_info_additions = (); # extra info fields
    foreach my $superpop (keys %sp_an) {
        my $ra_acs = $sp_ac{$superpop} || [];
        my $ra_het_acs = $sp_ac_het{$superpop} || [];
        my $ra_hom_acs = $sp_ac_hom{$superpop} || [];
        my $total_alleles = $sp_an{$superpop};

        my @allele_counts = ();
        my @het_allele_counts = ();
        my @hom_allele_counts = ();
        my @freqs = ();
        for (my $i=1; $i<=$no_alt_alleles; $i++) {
            my $ac = $ra_acs->[$i] || 0;
            my $het_ac = $ra_het_acs->[$i] || 0;
            my $hom_ac = $ra_hom_acs->[$i] || 0;
            push @allele_counts, $ac;
            push @het_allele_counts, $het_ac;
            push @hom_allele_counts, $hom_ac;
            my $af = ($total_alleles) ? sprintf("%9.8f", $ac/$total_alleles) : '';
            $af =~ s/0+$//;
            $af =~ s/\.$//;
            push @freqs, $af if ($af);
        }
        my $allelestring = join ',', @allele_counts;
        my $hetallelestring = join ',', @het_allele_counts;
        my $homallelestring = join ',', @hom_allele_counts;
            
        my $freqstring = join ',', @freqs;

        if ($ref eq 'grch38') { # just check that allele counts are right:
            my $grch38_total = ($info_field =~ /AN_$superpop\_unrel=(\d+)/) ? $1 : 'NA';
            my $grch38_ac = ($info_field =~ /AC_$superpop\_unrel=([\d,]+)/) ? $1 : 'NA';
            my $grch38_ac_hom = ($info_field =~ /AC_Hom_$superpop\_unrel=([\d,]+)/) ? $1 : 'NA';
            my $grch38_ac_het = ($info_field =~ /AC_Het_$superpop\_unrel=([\d,]+)/) ? $1 : 'NA';
            my $grch38_af = ($info_field =~ /AF_$superpop\_unrel=([\d.,]+)/) ? $1 : 'NA';
            if ($grch38_total != $total_alleles || $grch38_ac ne $allelestring) {
                print STDERR "Unequal estimates:\t$superpop $grch38_total/$total_alleles $grch38_ac/$allelestring $grch38_af/$freqstring\n";
            }
            if ($grch38_ac_hom ne $homallelestring || $grch38_ac_het ne $hetallelestring) {
                print STDERR "Unequal hom/het counts:\t$superpop $grch38_ac_hom/$homallelestring $grch38_ac_het/$hetallelestring\n";
            }
        }
        else {
            push @ac_info_additions, "AC_$superpop\_unrel=$allelestring";
            push @ac_hom_info_additions, "AC_Hom_$superpop\_unrel=$homallelestring";
            push @ac_het_info_additions, "AC_Het_$superpop\_unrel=$hetallelestring";
            push @af_info_additions, "AF_$superpop\_unrel=$freqstring" if ($freqstring);
            push @an_info_additions, "AN_$superpop\_unrel=$total_alleles";
        }
    }
    if ($ref ne 'grch38') {
        my $ac_info_tags = join ";", @ac_info_additions;
        my $achom_info_tags = join ";", @ac_hom_info_additions;
        my $achet_info_tags = join ";", @ac_het_info_additions;
        my $af_info_tags = (@af_info_additions) ? join ";", @af_info_additions : '';
        my $an_info_tags = join ";", @an_info_additions;
        $fields[7] .= ";$an_info_tags" if ($an_info_tags);
        $fields[7] .= ";$af_info_tags" if ($af_info_tags);
        $fields[7] .= ";$ac_info_tags;$achom_info_tags;$achet_info_tags";
    }
    my $vcf_record = join "\t", @fields[0..7];
    print "$vcf_record\n";
    $first_line = 0;
}
close GENOS;

sub read_sample_pops {
    my $file = shift;

    open SAMPLES, $file
        or die "Couldn\'t open $file: $!\n";

    my %samplepops = ();
    while (<SAMPLES>) {
        chomp;
        next if (/^#/);

        my @fields = split /\t/, $_;
        $fields[9] =~ s/\s+$//;
        $fields[10] =~ s/\s+$//;
        $samplepops{$fields[9]} = $fields[10];
    }
    close SAMPLES;

    return {%samplepops};
}

sub read_superpops {
    my $file = shift;

    open SUPERS, $file
        or die "Couldn\'t open $file: $!\n";

    my %superpops = ();
    while (<SUPERS>) {
        chomp;
        next if (/^Population/);

        my @fields = split /\t/, $_;

        $superpops{$fields[0]} = $fields[6];
    }
    close SUPERS;

    return {%superpops};
}

sub grab_field_headers {
    my $bgzipped_file = shift;

    my $grab_command = "gunzip -c $bgzipped_file | head -50000 | grep 'CHROM' | ";
    #print "$grab_command\n";

    open GRAB, $grab_command
        or die "Couldn\'t open $grab_command: $!\n";
    while (<GRAB>) {
        chomp;
        my @fields = split /\t/, $_;

        close GRAB;
        return @fields;
    }
}

sub check_field_headers {
    my $ra_field_headers = shift;
    my $rh_sample_pops = shift;

    my %missing_samples = ();
    my %confirmed_samples = ();
    for (my $sample_index = 9; $sample_index <= $#$ra_field_headers; $sample_index++) {
        my $samplename = $ra_field_headers->[$sample_index];
        if ($rh_sample_pops->{$samplename}) {
            $confirmed_samples{$samplename}++;
        }
        else {
            $missing_samples{$samplename}++;
        }
    }

    foreach my $sample2504 (keys %{$rh_sample_pops}) {
        if (!$confirmed_samples{$sample2504}) {
            #print "No sample $sample2504 in VCF file!\n";
        }
    }
 
    my $no_confirmed = keys %confirmed_samples;
    my $no_missing = keys %missing_samples;
    #print "$no_confirmed samples in 2504 file, $no_missing samples not\n";
    foreach my $missing_sample (keys %missing_samples) {
        #print "MISSING\t$missing_sample\n";
    }
}
