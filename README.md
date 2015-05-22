# MIE
Containing files used to investigate Mie scattering for summer 2015

FILE INVENTORY
- anything with test in the name is a temporary file
- all output written to working directory, miefiles or to a subdirectory
  of /home/njones/Dropbox/Mie/NatalieResults

- addphase.pro contains the function addphase, which takes two arrays and
   two values by which to weight them
       - used in rgcombsd.pro, combsd.pro

- combsd.pro will create plots for all .phsd files in miefiles as of May 22
  2015, and will create plots that combine the phase functions of carbon 
  and silicate particles for each colour for each author.
      - rgcombsd.pro does the same but puts both colours on the same plot 

- ct.pro contains the function ct, written by Peter Martin, which finds the
  cosine of the scattering angle for a given set of galactic coordinates for
  both the direction of illumination and direction of scattered propagation
       - used in scattest.pro, angletest.pro, pointingangles.pro

- findnside.pro (remove this)

- genfiledata.pro contains the function genfiledata which takes a source
  directory, file name, file extension, header size and a list of column 
  names. For each file name, it reads in header and column data and returns
  the entirety as a dictionary
      - used in sinpar.pro, sd.pro, combsd.pro, rgcombsd.pro

- genfilelist.pro contains the function genfilelist, which takes a list of
  particle sizes, colours and types and from these creates a list of 
  filenames
	- used in sinpar.pro

- hg.pro contains the function hg, written by Peter Martin, which takes a value
  for the anisotropy parameter g and an angle theta and computes the
  corresponding Henyey-Greenstein phase function value
  	- used in plotoptions.pro, scattest.pro

- interp[draine/indices/zubko].pro are a set of historical files that were
  used to produce refractive indices at particular wavelengths of interest

- intpf.pro contains the function intpf. which integrates the phase function in 
  spherical coordinates to confirm normalization
  	- used in sd.pro, sinpar.pro, combsd.pro

- intpftog.pro contains the function intpftog, which integrates the phase
  function multiplied by cos(theta) to find anisotropy parameter g
  	- used in sd.pro, sinpar.pro, combsd.pro

- mie.pro contains the function mie, which interpolates a phase function to
  a new value

- plotoptions.pro contains the function plotoptions, which is designed to allow
  for easy plotting of phase functions
      	- used in rgcombsd.pro, combsd.pro, sd.pro, sinpar.pro

- pointingangles.pro contains the function pointingangles, which for a geomerty
  of nebula and galactic center in galactic coordinates calculates the
  scattering angle towards earth
  	 - used in rgcombsd.pro, combsd.pro, sd.pro, sinpar.pro

- prheader.pro contains the function prheader, which for a given filename and
  header size, converts the header information into a dictionary and returns it
  	 - used in genfiledata.pro, readsisd.pro, rgcombsd.pro, combsd.pro,
	  sd.pro, sinpar.pro

- readsisd.pro reads in files with .sisd file extension and plots the weighted
  size distribution against the particle radius

- rgcombsd.pro (see combsd.pro)

- sd.pro plots phase functions for both colours for a particular size
  distribution and particle type on the same plot - uses .phsd files

- sdgenfilelist.pro contains the function sdgenfilelist, which takes a list of
  particle size distributions, colours and types and from these creates a 
  list of filenames
       - used in readsisd.pro, rgcombsd.pro, combsd.pro, sd.pro

- sindpar.pro plots phase functions for both colours for a particular particle 
  size and type on the same plot - uses .pr files