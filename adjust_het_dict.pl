#!/usr/bin/perl 

open IN, "<$ARGV[0]"; 
open OUT, ">$ARGV[1]"; 

while ($line=<IN>) {

	if (substr($line, 0, 7) eq RESIDUE) {
		if (substr($line, 10,3) =~ m/.. /) { substr($line, 10,3) =~ s/(..) / $1/; }
		if (substr($line, 10,3) =~ m/.  /) { substr($line, 10,3) =~ s/(.)  /  $1/; }
		print OUT $line;  
	}

	elsif (substr($line, 0, 6) eq CONECT) {
	        $atom_name="";
	        $conn_1="";
	        $conn_2="";
	        $conn_3="";
	        $conn_4="";
	        $conn_5="";
       		$conn_6="";
       		$conn_7="";
        	$conn_8="";
        	$conn_9="";
        	$conn_10="";

        	$atom_name=substr($line, 11,4);
        	$conn_1=substr($line, 20, 4);
        	$conn_2=substr($line, 25, 4);
        	$conn_3=substr($line, 30, 4);
        	$conn_4=substr($line, 35, 4);
        	$conn_5=substr($line, 40, 4);
        	$conn_6=substr($line, 45, 4);
        	$conn_7=substr($line, 50, 4);
        	$conn_8=substr($line, 55, 4);
        	$conn_9=substr($line, 60, 4);
        	$conn_10=substr($line, 65, 4);
		@line=split(" ", $line); 
        	if (length($line[1]) == 4) { $atom_name=$line[1]; }
        	if (length($line[3]) == 4) { $conn_1=$line[3]; }
        	if (length($line[4]) == 4) { $conn_2=$line[4]; }
        	if (length($line[5]) == 4) { $conn_3=$line[5]; }
        	if (length($line[6]) == 4) { $conn_4=$line[6]; }
        	if (length($line[7]) == 4) { $conn_5=$line[7]; }
        	if (length($line[8]) == 4) { $conn_6=$line[8]; }
        	if (length($line[9]) == 4) { $conn_7=$line[9]; }
        	if (length($line[10]) == 4) { $conn_8=$line[10]; }
        	if (length($line[11]) == 4) { $conn_9=$line[11]; }
        	if (length($line[12]) == 4) { $conn_10=$line[12]; }
		write OUT;
	}
	else {print OUT $line;}
}

format OUT = 
@<<<<<<<<<<@<<< @>>>@<<< @<<< @<<< @<<< @<<< @<<< @<<< @<<< @<<< @<<<
$line[0], $atom_name,$line[2],$conn_1, $conn_2, $conn_3, $conn_4, $conn_5, $conn_6, $conn_7, $conn_8, $conn_9, $conn_10 
.
