% File to convert pCO2 binary file to .nc file
% By Tatsuro Tanioka
% co2 air sea flux is in mmolC/m2/timestep (postive flux means CO2 is going into Sea from Air)

%basepath='~/TMM2/MITgcm_ECCO';
basepath='~/TMM2/MITgcm_2.8deg';
addpath(genpath('~/TMM2/tmm_matlab_code'));
oceanCarbonBasePath='~/TMM2/OceanCarbon';
%---------------------------------------------------------------------------
% Option to get pco2 based on Prognostic (useTimeVaryingPrescribed=0) or Prescibed CO2 (useTimeVaryingPrescribed=1)
useTimeVaryingPrescribedCO2=1

% Available options: 'historical', 'RCP3PD', 'RCP45', 'RCP6' and 'RCP85'
co2Scenario='historical';
atmosDataPath=fullfile(oceanCarbonBasePath,'AtmosphericCarbonData');
%-------------------------------------------------------------------------
% Making pco2 timeseries data
load(fullfile(basepath,'config_data'));
matrixPath=fullfile(basepath,matrixPath);
gridFile=fullfile(basepath,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');
load(gridFile,'nx','ny','nz','x','y','z','gridType');
load(boxFile,'izBox','nb');
load(profilesFile,'Irr');

if ~useTimeVaryingPrescribedCO2
    if exist('atm_output_time.txt', 'file') == 2;
       timeFile='atm_output_time.txt';
       [hdr,tdat]=hdrload(timeFile);
       T = tdat(:,2);
       nt=length(T);

       pCO2LogFile='pCO2atm_output.bin';
       diagData=read_binary(pCO2LogFile,[],'float64');
       pCO2atm= diagData;
   end
else
    if strcmp(co2Scenario,'historical')
        co2File='co2_HG11_extrap.dat'; % historical data from Heather Graven
    else
        co2File=[co2Scenario '_CO2_conc.txt']; % other scenarios
    end   
        [hdr,pCO2_atm]=hdrload(fullfile(atmosDataPath,co2File));
        T=pCO2_atm(:,1);
        pCO2atm=pCO2_atm(:,2);
end
if exist('pCO2atm','var') == 1; %
   ncid=netcdf.create('pco2.nc','CLOBBER');
   dimidt = netcdf.defDim(ncid,'Year',length(T));
   date_ID=netcdf.defVar(ncid,'Year','double',[dimidt]);
   pCO2_ID=netcdf.defVar(ncid,'pco2','double',[dimidt]);
   netcdf.endDef(ncid);
   netcdf.putVar(ncid,date_ID,T);
   netcdf.putVar(ncid,pCO2_ID,pCO2atm);
   netcdf.close(ncid);
end
% Making co2 airsea flux 2D data 
if exist('time_average_output_time.txt', 'file') == 2;
   avg_timeFile='time_average_output_time.txt';
   [hdravg,diagnostic_output_time]=hdrload(avg_timeFile);
   Tavg = diagnostic_output_time(:,2);
   Ib=find(izBox==1);
   nbb=length(Ib);
   diagFile='co2airseaflux_surf.bin';
   diagData=read_binary(diagFile,[],'float64');
   nt=length(diagData)/nbb;
   diagData=reshape(diagData,[nbb nt]);
   Diag=matrixToGrid(diagData,Ib,boxFile,gridFile);
   write2netcdf(['co2airseaflux.nc'],Diag,x,y,[],Tavg,'co2airseaflux','mmolC/m2/timestep')
end
