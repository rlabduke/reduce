#!/usr/bin/perl 

open IN, "<$ARGV[0]"; 
open OUT, ">temp"; 

while ($line=<IN>) {
	if ($line !~ m/RESIDUE/) { next; }
	$resn=substr($line, 10, 3); 
	%residue_hash=(); 
        %atom_name_hash=(); 
	$key_num=0;
	while (1) { 
		if ($line =~ m/CONECT/) { 
			$key =$key_num;
                        $key_num++;
			$atom=substr($line, 11, 4);
			$residue_hash{$key}=$line;
			$atom_name_hash{$atom}=$line; 
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
	$contains_carbon=0;
	$contains_oxygen=0;  
	for ($i=0; $i <= scalar(@formula_split); $i++) {
		$line_part=$formula_split[$i]; 
		if ($resn eq "ASP") { print $line_part."\n"; }
		if ($line_part =~ /^C/ && $line_part !~ /^C[A-Z]/) {
			$contains_carbon=1; 
		}
                if ($line_part =~ /^O/ && $line_part !~ /^O[A-Z]/) {
                        $contains_oxygen=1;
                }
	}
if ($resn eq "ASP") { print $contains_carbon.$contains_oxygen."\n"; }

	if ($contains_oxygen==1 && $contains_carbon==1) {
		for ($i=0; $i <= scalar(%residue_hash); $i++) {
			$line=$residue_hash{$i}; 
			chomp($line);
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
			if ($atom =~ /C/ && $bonded_atom_1 && $bonded_atom_2 && $bonded_atom_3 && !($bonded_atom_4)) {  # potential sp2 Carbon 
                                $start=substr($atom, 0,2);
                                $second=substr($atom,1,1);
				if ($second eq "C" && !( $start eq "SC" || $start eq "AC" || $start eq "TC" )) { # validate carbon
					$start=substr($bonded_atom_1, 0,2);
                                        $second=substr($bonded_atom_1,1,1);
					if ($second eq "C" && !( $start eq "SC" || $start eq "AC" || $start eq "TC" )) { # validate first bonded partner as carbon
						$start=substr($bonded_atom_2, 0,2);
                                		$second=substr($bonded_atom_2,1,1);
						if ($second eq "O" && !( $start eq "CO" || $start eq "MO" || $start eq "PO" )) { # validate second bonded partner as oxygen
							$start=substr($bonded_atom_3, 0,2);
							$second=substr($bonded_atom_3,1,1);
                                               		if ($second eq "O" && !( $start eq "CO" || $start eq "MO" || $start eq "PO" )) { # validate third bonded partner as oxygen
								$oxygen_1 = $bonded_atom_2;
								$oxygen_2 = $bonded_atom_3;
								$oxygen_1_line = $atom_name_hash{$oxygen_1};			
								$oxygen_2_line = $atom_name_hash{$oxygen_2};
						                chomp($oxygen_1_line);
								chomp($oxygen_2_line);
								$num_bonded_atoms_1=substr($oxygen_1_line, 15, 5)+0;
								$num_bonded_atoms_2=substr($oxygen_2_line, 15, 5)+0;
                                                                $num_bonded_atoms=$num_bonded_atoms_1 +  $num_bonded_atoms_2; 
       						       	        $bonded_atom_11=substr($oxygen_1_line, 20, 4);
						                $bonded_atom_12=substr($oxygen_1_line, 25, 4);
       						                $bonded_atom_21=substr($oxygen_2_line, 20, 4);
						                $bonded_atom_22=substr($oxygen_2_line, 25, 4);
								if ( ($num_bonded_atoms == 3 || $num_bonded_atoms == 4) && ($bonded_atom_11 eq $bonded_atom_21 && $bonded_atom_11 eq $atom) && ($bonded_atom_12 || $bonded_atom_22) )
								{
									if ($bonded_atom_12 && $bonded_atom_22) {
										$start=substr($bonded_atom_12, 0,2);
                                                                                $first=substr($bonded_atom_12, 0,1);
										if (($start eq " H") ||
                                                                               	   ($first eq "H" and resn !~ m/ HG|HG2|HGB|HGC|HGI|MAC|MBO|MMC|PHG|PMB|AAS|AMS|BE7|CMH|EMC|EMT|PHF|HF3|HF5| HO|HO3/))
                                                                               	{
                                                                               		$start=substr($bonded_atom_22, 0,2);
											$first=substr($bonded_atom_22, 0,1);
											if (($start eq " H") ||
                	                                                                   ($first eq "H" and resn !~ m/ HG|HG2|HGB|HGC|HGI|MAC|MBO|MMC|PHG|PMB|AAS|AMS|BE7|CMH|EMC|EMT|PHF|HF3|HF5| HO|HO3/))
        	                                                                        {
												print OUT $resn." ".$bonded_atom_12."\n";
                                                                                                print OUT $resn." ".$bonded_atom_22."\n";
											}
										}
									}
									elsif ($bonded_atom_12) {
										$start=substr($bonded_atom_12, 0,2);
										$first=substr($bonded_atom_12, 0,1);
                                                                                if (($start eq " H") ||
                                                                                   ($first eq "H" and resn !~ m/ HG|HG2|HGB|HGC|HGI|MAC|MBO|MMC|PHG|PMB|AAS|AMS|BE7|CMH|EMC|EMT|PHF|HF3|HF5| HO|HO3/))
                                                                                {
                                                                                        print OUT $resn." ".$bonded_atom_12."\n";
                                                                                }

									}
                                                                        else {
                                                                                $start=substr($bonded_atom_22, 0,2);
                                                                                $first=substr($bonded_atom_22, 0,1);
                                                                                if (($start eq " H") ||
                                                                                   ($first eq "H" and resn !~ m/ HG|HG2|HGB|HGC|HGI|MAC|MBO|MMC|PHG|PMB|AAS|AMS|BE7|CMH|EMC|EMT|PHF|HF3|HF5| HO|HO3/))
                                                                                {
                                                                                        print OUT $resn." ".$bonded_atom_22."\n";
                                                                                }
                                                                        }
								}
							} 
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


		

