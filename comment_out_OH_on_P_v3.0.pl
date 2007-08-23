#!/usr/bin/perl 

open IN, "<$ARGV[0]"; 
open OUT, ">temp"; 

while ($line=<IN>) {
	if ($line !~ m/RESIDUE/) { next; }
	$resn=substr($line, 10, 3); 
	%residue_hash=(); 
	$key_num=0; 
	while (1) { 
		if ($line =~ m/CONECT/) { 
			$key =$key_num;
			$key_num++; 
			$residue_hash{$key}=$line;
		}
		else {
			@array=split(" ", $line); 
			$key = $array[0];
			$residue_hash{$key}=$line;  
			if ($key =~ m/FORMUL/) { 
				last; 
			 }
		}
		$line =<IN>; 
	}

	@formula_split=split(" ", $residue_hash{FORMUL}); 
	$contains_phosphorus=0;
	$contains_oxygen=0;  
	for ($i=0; $i <= scalar(@formula_split); $i++) {
		$line_part=$formula_split[$i]; 
		if ($resn eq "ADP") { print $line_part."\n"; }
		if ($line_part =~ /P/ && $line_part !~ /P[A-Z]/) {
			$contains_phosphorus=1; 
		}
                if ($line_part =~ /O/ && $line_part !~ /O[A-Z]/) {
                        $contains_oxygen=1;
                }
	}
if ($resn eq "ADP") { print $contains_phosphorus.$contains_oxygen."\n"; }

	if ($contains_oxygen==1 && $contains_phosphorus==1) {
		for ($i=0; $i <= scalar(%residue_hash); $i++) {
			$line=$residue_hash{$i}; 
			$conect=substr($line, 0, 11);
	                $atom=substr($line, 11, 4);
	                $num_bonded_atoms=substr($line, 15, 5);
	                $bonded_atom_1=substr($line, 20, 4);
	                $bonded_atom_2=substr($line, 25, 4);
	                $bonded_atom_3=substr($line, 30, 4);
	                $bonded_atom_4=substr($line, 35, 4);
       		        $bonded_atom_5=substr($line, 40, 4);
               		$bonded_atom_6=substr($line, 45, 4);
                	$bonded_atom_7=substr($line, 50, 4);
                	$bonded_atom_8=substr($line, 55, 4);
                	$bonded_atom_9=substr($line, 60, 4);
                	$bonded_atom_10=substr($line, 65, 4);
			if ($atom =~ /O/ && $bonded_atom_1 && $bonded_atom_2) { 
				$start=substr($atom, 0,2); 
				$second=substr($atom,1,1); 
				if ($second eq "O" && !( $start eq "CO" || $start eq "MO" || $start eq "PO" )) {
					$start=substr($bonded_atom_1, 0,2);
                                	$second=substr($bonded_atom_1,1,1);
					if ($second eq "P" && !( $start eq "NP" )) {
						$start=substr($bonded_atom_2, 0,2);
	                                        $first=substr($bonded_atom_2,0,1);
						$second=substr($bonded_atom_2,1,1);
						if (($start eq " H") || 
						   ($first eq "H" and resn !~ m/ HG|HG2|HGB|HGC|HGI|MAC|MBO|MMC|PHG|PMB|AAS|AMS|BE7|CMH|EMC|EMT|PHF|HF3|HF5| HO|HO3/)) 						{
							print OUT $resn." ".$bonded_atom_2."\n"; 
						}
					}
				}
			}
		}
	}
}

close IN; 
close OUT; 

open IN, "<temp"; 
%to_be_commented_hash=(); 

while ($line=<IN>) {
	chomp($line); 
	$to_be_commented_hash{$line}="bad"; 
}
close IN; 

open DICT, "<$ARGV[0]"; 
open OUT, ">$ARGV[0]_mod"; 

while ($line=<DICT>) { 
	if ($line !~ m/RESIDUE/ && $line !~ m/CONECT/) {
		print OUT $line; 
		next; 
	}
	if ($line =~ m/RESIDUE/) { $resn=substr($line, 10, 3); }
	$atom=substr($line, 11, 4);
	$key=$resn." ".$atom; 
	if ($to_be_commented_hash{$key}) {
		$line="xx".$line;
	}
	print OUT $line;  
}

close DICT; 
close OUT; 			


		

