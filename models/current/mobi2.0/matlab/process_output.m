% process_output.m 
% Script to convert .petsc files to .nc file for tracer and physics outputs 
% Tatsuro Tanioka
% .m files needed in the same directory are
% load_diagnostics.m for postprocessing diagnostic fluxes to netcdf
% load_output.m for postprocessing tracer snapshot concentrations 
% load_output_time_avg.m for postprocessing tracer time-averaged concentrations 
% load_output_sed.m for postprocessing sedimentary fluxes (not used in default MOBI runs) 

load_diagnostics;
load_output;
load_output_time_avg;
%load_output_sed;


