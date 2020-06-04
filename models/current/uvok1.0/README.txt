This is version 1.0 of UVOK, an interface to the Kiel version of the biogeochemical model 
embedded in the University of Victoria Earth System Climate Modle (UVic ESCM). (As you 
will have guessed, ``UVOK'' stands for the Universities of Victoria, Oxford and Kiel.) 

To run this model requires a number of components that are not distributed here: 
1) The full, updated UVic ESCM 2.9 model code available from: UVic_ESCM.2.9.updated.tar.gz (2013-8-1) 
at http://wikyonos.seos.uvic.ca/model/
2) Modifications to (1) involving extensive changes of the biogeochemical model and 
code implementing transport matrix extraction and interfacing of the biogeochemical 
model to the TMM: http://kelvin.earth.ox.ac.uk/spk/Research/TMM/Kiel_Jan_2017_updates_to_UVIC2.9.tar.gz 
or https://thredds.geomar.de/thredds/catalog/open_access/Kiel_Jan_2017_updates_to_UVIC2.9/catalog.html.

Once you have all the above components, set environment and other variables as described 
at the top of the Makefile.

Then compile with:
make tmmuvok
