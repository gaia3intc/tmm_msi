This is the TMM interface to the MITgcm 'gchem' package. It provides hooks to 
the various biogeochemical models implemented within MITgcm making it possible 
to use them without modification of the original code. Currently (as of May 27, 
2019) only pkg/dic is supported but other models (e.g., BLING, Darwin, ...) 
can be easily used by appropriately adapting external_forcing_mitgchem_dic.c.

You will need a recent copy of MITgcm for the required header and source files 
not distributed with this package.

List of required files:

namelist files:
data.gchem 
data.dic (for using pkg/dic)

The following source files are needed from the main MITgcm source tree. Of these 
only SIZE.h, PTRACERS_SIZE.h and dic_biotic_forcing.F need to be modified (as 
described below) and are supplied in src/. The rest are used without change by 
setting the path to the top level of the MITgcm source in Makefile.

pkg/dic:
alk_surfforcing.F
bio_export.F
calcite_saturation.F
car_flux.F
car_flux_omega_top.F
carbon_chem.F
dic_atmos.F
dic_biotic_diags.F
dic_biotic_forcing.F
dic_biotic_init.F
dic_diagnostics_init.F
dic_ini_atmos.F
dic_ini_forcing.F
dic_init_fixed.F
dic_init_varia.F
dic_read_co2_pickup.F
dic_readparms.F
dic_surfforcing.F
dic_surfforcing_init.F
fe_chem.F
insol.F
o2_surfforcing.F
phos_flux.F
DIC_ATMOS.h
DIC_CTRL.h
DIC_LOAD.h
DIC_OPTIONS.h
DIC_VARS.h

pkg/gchem:
gchem_check.F
gchem_forcing_sep.F
gchem_init_fixed.F
gchem_init_vari.F
gchem_readparms.F
GCHEM.h
GCHEM_FIELDS.h
GCHEM_OPTIONS.h

pkg/ptracers:
PTRACERS_FIELDS.h
PTRACERS_OPTIONS.h
PTRACERS_PARAMS.h
PTRACERS_SIZE.h

model/src:
packages_unused_msg.F

model/inc:
CPP_OPTIONS.h
DYNVARS.h
FFIELDS.h
GRID.h
PARAMS.h
SIZE.h
*_MACROS.h

eesupp/src:
get_periodic_interval.F
lef_zero.F
mds_flush.F
mds_reclen.F
mdsfindunit.F
nml_change_syntax.F
open_copy_data_file.F
print.F
utils.F
write_utils.F
different_multiple.F

eesupp/inc:
CPP_EEMACROS.h
CPP_EEOPTIONS.h
EEPARAMS.h
EESUPPORT.h

tools:
set64bitConst.csh

Other required MITgcm files:
PACKAGES_CONFIG.h (automatically generated; turn on PTRACERS, GCHEM and 
whichever biogeochemical package you want (e.g., DIC)).

MITgcm files that need to be changed and supplied here:
SIZE.h (set all s* and n* variables to 1, OL* to 0 and Nr to number of depth levels)
PTRACERS_SIZE.h (set PTRACERS_num to 6 for pkg/dic, 7 for pkg/bling etc)
dic_biotic_forcing.F (add lines to include TMM_MITGCHEM_DIAGS.h)

Also note that the following options in DIC_OPTIONS.h should *always* be undefined:
ALLOW_OLD_VIRTUALFLUX
USE_PLOAD
ALLOW_DIAGNOSTICS
USE_QSW
USE_QSW_UNDERICE
USE_ATMOSCO2

Note on building the code: To build, type 'make tmmmitgchemdic'. This generally works 
but on some machines it fails because the C preprocessor fails to automatically run on 
the *.F files. To fix this, first type: 'make smallf' followed by 'make tmmmitgchemdic'.
