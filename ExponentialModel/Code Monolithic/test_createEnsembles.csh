#!/bin/tcsh

## program variables
set ns = 0
set ne = 1
set es = 1
set ga = 0
set to = .05
set se = -123
set gr = 2
set sq = "../../Code/"
set st = "../Results/"
set ns = 1
set fi = ""

## script variables
set gamma = {-2.000}
set gammap = {2.0}
set sizes = {4}
set nscnt = 1
set ngcnt = 1
set nqcnt = 1
set ngamm = 2
set nsize = 2
set nseqs = 1
set soffs = -82
@ nseqs = $ns + 1

# change dir to storage
cd $st

## begin code
while($nscnt < $nsize)
	# create size directory
	mkdir x$sizes[$nscnt]
	cd    x$sizes[$nscnt]

	while($ngcnt < $ngamm)
		# create gamma directory

		while($nqcnt < $nseqs)
			# create sequence directory


			# get name for file
			set fname = ${sq}degseq_N${sizes[$nscnt]}r${gammap[$ngcnt]}Sq${nqcnt}.txt
				#echo $fname
				#pwd
				#@ ns = $sizes[$nscnt]
				#@ ga = $gamma[$ngcnt]

			# execute code
			pwd
			./../../Code/NMLib/EXPMOD.exe -N $sizes[$nscnt] -e $es -g $gamma[$ngcnt] -t $to -r $se -p ./ -s $fname -n $nqcnt

			# inc seq counter
			@ nqcnt = $nqcnt + 1

			# modify random number seed
			@ se = $se - $soffs

			# change up directory
			cd ..
		end

		# reset sequence counter
		@ nqcnt = 1

		# inc gamma counter
		@ ngcnt = $ngcnt + 1

		# change up directory
		cd ..
	end

	# reset gamma counter
	@ ngcnt = 1

	# inc size counter
	@ nscnt = $nscnt + 1

	# change up directory
	cd ..
end
