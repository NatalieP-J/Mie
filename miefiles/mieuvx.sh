#! /bin/sh
# sample Bourne shell script for associating filenames with
# logical units in fortran programs
# mie theory calculation, including x-ray extension

# input mie, rho, sizes
cp $1.uvx t4.uvx
# input wavelength and refractive index
cp $2.refr t2.refr

mieuvx	# execute fortran

# printed output
cp t3.pr $1$2.pr
# plotting output for plotq and  for hazy
cp t1.po $1$2.po

