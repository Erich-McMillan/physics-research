#!/bin/tcsh

## program variables
set ns = 0
set ne = 1
set es = 1000
set ga = 0
set to = .05
set se = -1
set gr = 2
set ns = 1000
set sq = "../../../WeibinResults/seq316_dense_r-0123/"
set st = "../results_dense_2017_08_09"

## request user input for sq and st
echo $sq
echo $st

## script variables
set gamma = {0.000,1.000,2.000,3.000}
set gammap = {0.0,-1.0,-2.0,-3.0}
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
			set fname = ${sq}dbroot_N${sizes[$nscnt]}r${gammap[$ngcnt]}Sq${nqcnt}.txt
				echo $fname
			# execute code
			pwd
			./NMLib/GENERGM_dense.exe -N $sizes[$nscnt] -e $es -g $gamma[$ngcnt] -r $se -p $st -s $fname -n $nqcnt
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
