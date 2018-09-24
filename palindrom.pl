use strict;
use warnings;

my $regex = '[ACG]AT[CGT]';

if (regex_inv($regex) eq $regex)
{
	print "OK";
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