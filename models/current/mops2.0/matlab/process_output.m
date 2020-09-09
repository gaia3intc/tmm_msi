% process_output.m 
% Script to convert .petsc files to .nc file for tracer and physics outputs 
% Tatsuro Tanioka
% .m files needed in the same directory are
% n7fluxes28.m for postprocessing petsc2 diagnostic fluxes to netcdf
% n7tracers28.m for postprocessing petsc2 snapshot of tracers to netcdf
% n7tracersavg28.m for postprocessing petsc2 tracers averaged over certain time period to netcdf
% n7physics.m for postprocessing petsc2 monthly salinity and temperature to netcdf
% load_pco2.m for postprocessing binary pco2 and co2airseaflux files to netcdf

n7fluxes28('fields_fluxes.nc');
n7tracers28('fields_tracers_snapshot.nc');
n7tracersavg28('fields_tracers_avg.nc');
n7physics('fields_TS.nc');
load_pco2

