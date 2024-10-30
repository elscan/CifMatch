# CifMatch
Given a list of reflections from your PXRD scan and the wavelength, find the CIF file that most closely matches your diffraction profile

There should just be one file, _FindBestMatchToReflectionList.py \n
The instructions are written in _FindBestMatchToReflectionList.py and are copied below

\#----- Instructions -----# \n
\# Above is a list, tThets
\# For your diffraction pattern, enter into tThets a list of 2-theta values, in degrees, where there are reflections
\# Above is a variable, wavelength
\# change wavelength to the wavelength your experiment was done at
\# Put this file into a folder with all the .cif files that you want to search through
\# Run the file from the anaconda prompt/command prompt/terminal (it is assumed the user knows how to do this)
\# A file, 2ThetaHits.csv should appear in the folder with all the .cif files that you want to search through
\# Have a look at 2ThetaHits.csv
