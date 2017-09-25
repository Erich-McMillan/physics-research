#!/bin/csh

set sizes = {32,100,316}
set gamma = {-2.0,-2.5,-3.0,-3.5}
set rand = -239
set offset = -83

set cone = 1
set ctwo = 4
set ione = 1
set itwo = 5
set num = 1

mkdir results

while($cone < $ctwo)
	echo "hello"
	while($ione < $itwo)
		echo "   y"	
		while($num < 11)
			echo "      i"
			./GenSeq.exe $sizes[$cone] $gamma[$ione] $num $rand results/ 
			@ rand = $rand - $offset	
			@ num = $num + 1
		end
		
		@ ione = $ione + 1
		@ num = 1
	end
	@ ione = 1
	@ cone = $cone + 1
end
