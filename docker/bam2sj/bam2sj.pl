#!/usr/bin/perl

use Data::Dumper;

$minmapq=255;
my $HoS = {};

my $sam_view = $ARGV[0];

open(my $fh, '<', $sam_view) or die "Cannot open input file $sam_view";

while (<$fh>) {
    if ($_ =~ /^@/) {
        next;
    }
    else {
        my @line = split /\t/,$_;
        my $chr = $line[2];
        my $start = $line[3];
        my $cigar = $line[5];
        my $read = $line[0];
        my $mqu = $line[3];

        if ($mqu >= $minmapq) {
            if ($cigar =~ /N/) { #if split read, then parse
                my @inter = split /\D/, $cigar;
                (my $t = $cigar) =~ s/\d//g;
                my @cops = split //,$t;
                #print "@inter @cops\n";
            
                for (my $i=0; $i <= $#cops; $i++) { #use nested loop to count multiple splits in a single read!
                    my $ofs = 0;
                    if ($cops[$i] eq 'N')   {
                        for (my $k=0; $k < $i; $k++) {
                            if (($cops[$k] eq 'M') || ($cops[$k] eq 'N') || ($cops[$k] eq 'D')) {
                                $ofs += $inter[$k];
                            }
                            elsif ($cops[k] eq 'I') {
                                $ofs -= $inter[$k];
                            }
                        ## ignore S and s                    
                        }
                    my $junc_start=$start + $ofs;
                    my $junc_end=$start+$ofs+$inter[$i]-1;

                    my $hk = $chr."_".$junc_start."_".$junc_end;
                    #print "$hk\n";
                    
                    if (exists $HoS{$hk}) {
                        $HoS{$hk}+=1;
                    }

                    else {
                        $HoS{$hk}=1;
                    }

                    #print "$chr\t$junc_start\t$junc_end\t$read\t$cigar\t$start\n";
                    }
                }
            }
        }
    }
}
close $fh;

while (($key, $counts) = each (%HoS)) {
    ( my $t = $key ) =~ s/_/\t/g;
    print "$t\t$counts\n";
}