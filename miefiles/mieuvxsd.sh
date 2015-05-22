#! /bin/sh
# sample Bourne shell script for associating filenames with
# logical units in fortran programs which use (call) the ioinit routine
# ******NOTE that FILES are in the directory that calls this script,
# unless explicitly indicated
# mie theory calculation, including x-ray extension, with size distribution

# input mie, rho, sizes
cp $1.uvxsd tsd4.uvxsd
# input wavelength and refractive index
cp $2.refr tsd2.refr

# empty dummy file for single particle output for diagnostic purposes
cp dummy1.sd tsd1.sd 
# empty dummy files for plotting output for up to 2 size distributions
cp dummy13.posd tsd13.posd 
cp dummy14.posd tsd14.posd 

mieuvxsd		 	# execute fortran

grep fail tsd3.prsd

# save output files
# printed output
mv tsd3.prsd $1$2.prsd

# the following is to remove .sd single particle diagnostic if not used
if test -s tsd1.sd              # does file exist and have size .gt. 0
then				# leave it
		 mv tsd1.sd $1$2.sd
else
	if test -f tsd1.sd      # does file exist and have size .eq. 0
	then
		 rm tsd1.sd	# remove it
	fi
fi


# for plotting and hazy, and qphase for kramers-kronig checks
# ibase = 12 used in main program
# maximum 2 size distributions; always use at least 1
mv tsd13.posd $1$21.posd

# the following is to remove the second .posd file if not used
if test -s tsd14.posd           # does file exist and have size .gt. 0
then				# leave it
		 mv tsd14.posd $1$22.posd
else
	if test -f tsd14.posd   # does file exist and have size .eq. 0
	then
		 rm tsd14.posd	# remove it
	fi
fi
