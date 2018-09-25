use strict;
use warnings;

my $seq = "NNNCnnGnnnnnnnnNNnnnnn";
$seq =~ s/(^[nN]+)|([nN]+$)//g;
my $regex = enzREGEX(uc $seq);
print "$regex\t" . regex_inv($regex) . "\n";



sub enzREGEX
	{
		my $regex = "";
		my ($Seq) = @_;
		my @SeqIndv = split("",$Seq);
		my $SeqIndvL = @SeqIndv;
		for (my $i=0;$i<$SeqIndvL;$i++) {
			if ($SeqIndv[$i] =~ /[ATGC]/) {$regex = $regex.$SeqIndv[$i]}
			elsif ($SeqIndv[$i] eq "R") {$regex = $regex."[AG]"}
			elsif ($SeqIndv[$i] eq "Y") {$regex = $regex."[CT]"}
			elsif ($SeqIndv[$i] eq "S") {$regex = $regex."[CG]"}
			elsif ($SeqIndv[$i] eq "W") {$regex = $regex."[AT]"}
			elsif ($SeqIndv[$i] eq "K") {$regex = $regex."[GT]"}
			elsif ($SeqIndv[$i] eq "M") {$regex = $regex."[AC]"}
			elsif ($SeqIndv[$i] eq "B") {$regex = $regex."[CGT]"}
			elsif ($SeqIndv[$i] eq "D") {$regex = $regex."[AGT]"}
			elsif ($SeqIndv[$i] eq "H") {$regex = $regex."[ACT]"}
			elsif ($SeqIndv[$i] eq "V") {$regex = $regex."[ACG]"}
			elsif ($SeqIndv[$i] eq "N") {$regex = $regex."[ACGT]"}
		}
		return $regex;
	}

sub regex_inv
{
	my $input = $_[0];
	my @data = split("", $input);
	my @data_inv = ();
	
	my @temp = ();
	my $check = 0;
	for (my $i = scalar(@data) - 1; $i >= 0; $i-- )
	{
		if ($data[$i] =~ /\]/)
		{
			$check = 1;
			next;
		}
		elsif ($data[$i] =~ /\[/)
		{
			my @temp_sorted = sort { $a cmp $b } @temp;
			push @temp_sorted, ']';
			unshift @temp_sorted, '[';
			push @data_inv, join("", @temp_sorted);
			
			@temp = ();
			$check = 0;
		}
		elsif ($check == 1 and $data[$i] =~ /[^\[\]]/ )
		{
			unshift @temp, invert($data[$i]);
		}
		else
		{
			push @data_inv, invert($data[$i]);
		}
	}
	
	return join("", @data_inv);
}

sub invert
{
	my $input = $_[0];
	my $output = "";
	if ($input eq 'A') { $output = 'T' }
	elsif ($input eq 'T') { $output = 'A' }
	elsif ($input eq 'C') { $output = 'G' }
	elsif ($input eq 'G') { $output = 'C' }
	else { $output = $input }
	
	return $output;
}