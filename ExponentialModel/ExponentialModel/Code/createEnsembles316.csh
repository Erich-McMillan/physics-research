#!/bin/tcsh

## program variables
set ns = 0
set ne = 1
set es = 10000
set ga = 0
set to = .05
set se = -1
set gr = 2
set ns = 1000
set sq = "../../LagrangianSolver/results_2017_07_03/"
set st = "../results_2017_07_03"

## request user input for sq and st
echo $sq
echo $st

## script variables
set gamma = {-2.000,-2.500,-3.000,-3.500}
set gammap = {2.0,2.5,3.0,3.5}
set sizes = {316}
set nscnt = 1
set ngcnt = 1
set nqcnt = 1
set ngamm = 5
set nsize = 2
set nseqs = 1
set soffs = 2
@ nseqs = $ns + 1

## begin code
while($nscnt < $nsize)
	while($ngcnt < $ngamm)
		while($nqcnt < $nseqs)
			# get name for file
			set fname = ${sq}lagrMultiN${sizes[$nscnt]}r${gammap[$ngcnt]}Sq${nqcnt}.txt
				echo $fname
			# execute code
			pwd
			./NMLib/GENERGM.exe -N $sizes[$nscnt] -e $es -g $gamma[$ngcnt] -r $se -p $st -s $fname -n $nqcnt
			# inc seq counter
			@ nqcnt = $nqcnt + 1
			# modify random number seed
			@ se = $se - $soffs
		end
		# reset sequence counter
		@ nqcnt = 1
		# inc gamma counter
		@ ngcnt = $ngcnt + 1
	end
	# reset gamma counter
	@ ngcnt = 1
	# inc size counter
	@ nscnt = $nscnt + 1
end
