function write_fluxes28(netcdfFileName);
%function petsc2netcdft(netcdfFileName);
% Function for postprocessing petsc2 diagnostic fluxes to netcdf
% From Iris Kriest (GEOMAR), Jan-28-2020)
% Edited by Tatsuro Tanioka (Nov 2020)
%    - Added PAR [W m-2] at the top of every layer in the euphotic zone (fbgc8)
%    - Added SEDIMENT_C [mmolC/m2/timestep] (fbgc9)
%    - Added CP_PHYTO_UPTAKE [mol C:mol P] (fbgc10)
%    - Added CP_ZOO_UPTAKE [mol C:mol P] (fbgc11)

% To run this script, type in MATLAB command line
% > n7fluxes28('test.nc')

% Notes from Iris Kriest:
% - Unit of the fluxes is per ocean time step. So if using a time step length for the 
%     ocean of half a day, multiply by 2 to get the result per day
%   All fluxes are in mmol P/m3/ocean timestep, except:
%   "Nfixation" which is in mmol N/m3/ocean time step
%   "Denitrification" which is in mmol N/m3/ocean time step
%   "Sedimentation" which is in mmol P/m2/ocean timestep. Seimentation is organic POP flux into each box.
%    (i.e., through the upper box boundary). But for layer 1, this flux is burial at the sea floor.

['write_fluxes28']
		
% Options for output   
useORGCARBON = 1;  % set to 1 if option DORGCARBON was used in model comiplation, 0 otherwise
usePFT = 1;        % set to 1 if option DPFT was used, 0 otherwise
useFLEXCP = 1;     % set to 1 if option DFLEXCP was used, 0 otherwise
   
%%%%% Set toplevel path to GCMs configuration
basepath='~/TMM2/MITgcm_2.8deg';
%basepath='~/TMM2/MITgcm_ECCO';

addpath(genpath('~/TMM2/tmm_matlab_code'));

load(fullfile(basepath,'config_data'))
matrixPath=fullfile(basepath,matrixPath);
gridFile=fullfile(basepath,'grid');
load(gridFile,'nx','ny','nz','dz','dznom','x','y','z');
boxFile=fullfile(matrixPath,'Data','boxes');
load(boxFile,'Xboxnom','Yboxnom','Zboxnom','ixBox','iyBox','izBox','nb','volb');

profilesFile=fullfile(matrixPath,'Data','profile_data');
load(profilesFile,'Irr')

Ib=find(izBox(1,:)==1)';
nbb=length(Ib);

for is=1:nbb
  ibs=Ib(is);
  Ip{is}=find(Xboxnom==Xboxnom(ibs) & Yboxnom==Yboxnom(ibs));
  [zp,izp]=sort(Zboxnom(Ip{is}));
  Ip{is}=Ip{is}(izp);
end
Ir=cat(1,Ip{:});

v1=readPetscBinVec('fbgc1.petsc',-1);
v2=readPetscBinVec('fbgc2.petsc',-1);
v3=readPetscBinVec('fbgc3.petsc',-1);
v4=readPetscBinVec('fbgc4.petsc',-1);
v5=readPetscBinVec('fbgc5.petsc',-1);
v6=readPetscBinVec('fbgc6.petsc',-1);
v7=readPetscBinVec('fbgc7.petsc',-1);
v8=readPetscBinVec('fbgc8.petsc',-1);

%if useORGCARBON
%    v9=readPetscBinVec('fbgc9.petsc',-1);
%    v10=readPetscBinVec('fbgc10.petsc',-1);
%    v11=readPetscBinVec('fbgc11.petsc',-1);
%end

nb=size(v1,1);
[V1,x,y,z]=matrixToGrid(v1(Irr,:),[1:nb]',boxFile,gridFile);
[V2,x,y,z]=matrixToGrid(v2(Irr,:),[1:nb]',boxFile,gridFile);
[V3,x,y,z]=matrixToGrid(v3(Irr,:),[1:nb]',boxFile,gridFile);
[V4,x,y,z]=matrixToGrid(v4(Irr,:),[1:nb]',boxFile,gridFile);
[V5,x,y,z]=matrixToGrid(v5(Irr,:),[1:nb]',boxFile,gridFile);
[V6,x,y,z]=matrixToGrid(v6(Irr,:),[1:nb]',boxFile,gridFile);
[V7,x,y,z]=matrixToGrid(v7(Irr,:),[1:nb]',boxFile,gridFile);
[V8,x,y,z]=matrixToGrid(v8(Irr,:),[1:nb]',boxFile,gridFile);

if usePFT  
  if useORGCARBON   
    if useFLEXCP
      ['MOPS + (CARBON assumed) + ORGCARBON + PFT + FLEXCP']
      v9=readPetscBinVec('fbgc9.petsc',-1);        % photosyn for PFT 2
      v10=readPetscBinVec('fbgc10.petsc',-1);      % sediment_C
      v11=readPetscBinVec('fbgc11.petsc',-1);      % C:P for PFT 1
      v12=readPetscBinVec('fbgc12.petsc',-1);      % C:P for PFT 2
      v13=readPetscBinVec('fbgc13.petsc',-1);      % C:P for Zoo
      [V9,x,y,z]=matrixToGrid(v9(Irr,:),[1:nb]',boxFile,gridFile);
      [V10,x,y,z]=matrixToGrid(v10(Irr,:),[1:nb]',boxFile,gridFile);
      [V11,x,y,z]=matrixToGrid(v11(Irr,:),[1:nb]',boxFile,gridFile);
      [V12,x,y,z]=matrixToGrid(v12(Irr,:),[1:nb]',boxFile,gridFile);
      [V13,x,y,z]=matrixToGrid(v13(Irr,:),[1:nb]',boxFile,gridFile);
    else
      ['MOPS + (CARBON assumed) + ORGCARBON + PFT']
      v9=readPetscBinVec('fbgc9.petsc',-1);        % photosyn for PFT 2
      v10=readPetscBinVec('fbgc10.petsc',-1);      % sediment_C
      [V9,x,y,z]=matrixToGrid(v9(Irr,:),[1:nb]',boxFile,gridFile);
      [V10,x,y,z]=matrixToGrid(v10(Irr,:),[1:nb]',boxFile,gridFile);
    end    
  else
    ['MOPS + PFT']
    v9=readPetscBinVec('fbgc9.petsc',-1);          % photosyn for PFT 2
    [V9,x,y,z]=matrixToGrid(v9(Irr,:),[1:nb]',boxFile,gridFile);
  end
else
  if useORGCARBON
    if useFLEXCP
      ['MOPS + (CARBON assumed) + ORGCARBON + FLEXCP']
      v9=readPetscBinVec('fbgc9.petsc',-1);        % sediment_C
      v10=readPetscBinVec('fbgc10.petsc',-1);      % C:P for PFT
      v11=readPetscBinVec('fbgc11.petsc',-1);      % C:P for Zoo
      [V9,x,y,z]=matrixToGrid(v9(Irr,:),[1:nb]',boxFile,gridFile);
      [V10,x,y,z]=matrixToGrid(v10(Irr,:),[1:nb]',boxFile,gridFile);
      [V11,x,y,z]=matrixToGrid(v11(Irr,:),[1:nb]',boxFile,gridFile);
    else
      ['MOPS + (CARBON assumed) + ORGCARBON']
      v9=readPetscBinVec('fbgc9.petsc',-1);        % sediment_C
      [V9,x,y,z]=matrixToGrid(v9(Irr,:),[1:nb]',boxFile,gridFile);
    end
  else
    ['MOPS']
  end 
end
  
%useORGCARBON
% [V9,x,y,z]=matrixToGrid(v9(Irr,:),[1:nb]',boxFile,gridFile);
% [V10,x,y,z]=matrixToGrid(v10(Irr,:),[1:nb]',boxFile,gridFile);
% [V11,x,y,z]=matrixToGrid(v11(Irr,:),[1:nb]',boxFile,gridFile);
%

[nx,ny,nz,nt]=size(V1);

% fake time data
t=[0:nt-1]';

if isempty(Irr)
  Irr=[1:nb]';
end

% open file
[ncid, status] = mexnc('CREATE', netcdfFileName, 'noclobber' );
if status, error(mexnc('STRERROR',status)), end

% define dimensions
[lon_dimid, status] = mexnc('DEF_DIM',ncid,'Longitude',nx);
if status, error(mexnc('STRERROR',status)), end
[lat_dimid, status] = mexnc('DEF_DIM',ncid,'Latitude',ny);
if status, error(mexnc('STRERROR',status)), end
[dep_dimid, status] = mexnc('DEF_DIM',ncid,'Depth',nz);
if status, error(mexnc('STRERROR',status)), end
[time_dimid, status] = mexnc('DEF_DIM',ncid,'Time',nt);
if status, error(mexnc('STRERROR',status)), end

% Define the coordinate variables
[lon_varid, status] = mexnc('DEF_VAR', ncid, 'Longitude', nc_double, 1, lon_dimid );
if status, error(mexnc('STRERROR',status)), end
[lat_varid, status] = mexnc('DEF_VAR', ncid, 'Latitude', nc_double, 1, lat_dimid );
if status, error(mexnc('STRERROR',status)), end
[dep_varid, status] = mexnc('DEF_VAR', ncid, 'Depth', nc_double, 1, dep_dimid );
if status, error(mexnc('STRERROR',status)), end
[time_varid, status] = mexnc('DEF_VAR', ncid, 'Time', nc_double, 1, time_dimid );
if status, error(mexnc('STRERROR',status)), end

%PHOSY
[data_varid1, status] = mexnc('DEF_VAR', ncid, 'Photo', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid1, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%GRAZING
[data_varid2, status] = mexnc('DEF_VAR', ncid, 'Graze', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid2, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%SEDIMENTATION
[data_varid3, status] = mexnc('DEF_VAR', ncid, 'Sediment', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid3, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%REMINERALISATION, AEROBIC
[data_varid4, status] = mexnc('DEF_VAR', ncid, 'Remin', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid4, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%RUNOFF
[data_varid5, status] = mexnc('DEF_VAR', ncid, 'Runoff', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid5, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%NFixation
[data_varid6, status] = mexnc('DEF_VAR', ncid, 'NFix', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid6, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%Denitrification
[data_varid7, status] = mexnc('DEF_VAR', ncid, 'Denit', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid7, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%PAR
[data_varid8, status] = mexnc('DEF_VAR', ncid, 'PAR', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid8, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%if useORGCARBON
%% SEDIMENTATION of POC
%    [data_varid9, status] = mexnc('DEF_VAR', ncid, 'Sediment_C', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
%    if status, error(mexnc('STRERROR',status)), end
%    status=mexnc('put_att_double', ncid, data_varid9, 'missing_value', nc_double', 1, NaN );
%% Phytoplankton C:P uptake 
%    [data_varid10, status] = mexnc('DEF_VAR', ncid, 'CP_PHYTO_UPTAKE', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
%    if status, error(mexnc('STRERROR',status)), end
%    status=mexnc('put_att_double', ncid, data_varid10, 'missing_value', nc_double', 1, NaN );
%% Zooplankton C:P uptake 
%    [data_varid11, status] = mexnc('DEF_VAR', ncid, 'CP_ZOO_UPTAKE', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
%    if status, error(mexnc('STRERROR',status)), end
%    status=mexnc('put_att_double', ncid, data_varid11, 'missing_value', nc_double', 1, NaN );
%end

if usePFT  
  if useORGCARBON   
    if useFLEXCP
      ['MOPS + (CARBON assumed) + ORGCARBON + PFT + FLEXCP']
       % Photosynthesis for PFT2
      [data_varid9, status] = mexnc('DEF_VAR', ncid, 'Photo2', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
      if status, error(mexnc('STRERROR',status)), end
      status=mexnc('put_att_double', ncid, data_varid9, 'missing_value', nc_double', 1, NaN );
      if status, error(mexnc('STRERROR',status)), end
      % SEDIMENTATION of POC
      [data_varid10, status] = mexnc('DEF_VAR', ncid, 'Sediment_C', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
      if status, error(mexnc('STRERROR',status)), end
      status=mexnc('put_att_double', ncid, data_varid10, 'missing_value', nc_double', 1, NaN );
      % Phytoplankton 1 C:P uptake 
      [data_varid11, status] = mexnc('DEF_VAR',ncid,'CP_PFT1_UPTAKE', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
      if status, error(mexnc('STRERROR',status)), end
      status=mexnc('put_att_double', ncid, data_varid11, 'missing_value', nc_double', 1, NaN );
      % Phytoplankton 2 C:P uptake 
      [data_varid12, status] = mexnc('DEF_VAR',ncid,'CP_PFT2_UPTAKE', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
      if status, error(mexnc('STRERROR',status)), end
      status=mexnc('put_att_double', ncid, data_varid12, 'missing_value', nc_double', 1, NaN );
      % Zooplankton C:P uptake 
      [data_varid13, status] = mexnc('DEF_VAR', ncid, 'CP_ZOO_UPTAKE', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
      if status, error(mexnc('STRERROR',status)), end
      status=mexnc('put_att_double', ncid, data_varid13, 'missing_value', nc_double', 1, NaN );
    else
      ['MOPS + (CARBON assumed) + ORGCARBON + PFT']
       % Photosynthesis for PFT2
      [data_varid9, status] = mexnc('DEF_VAR', ncid, 'Photo2', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
      if status, error(mexnc('STRERROR',status)), end
      status=mexnc('put_att_double', ncid, data_varid9, 'missing_value', nc_double', 1, NaN );
      if status, error(mexnc('STRERROR',status)), end
      % SEDIMENTATION of POC
      [data_varid10, status] = mexnc('DEF_VAR', ncid, 'Sediment_C', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
      if status, error(mexnc('STRERROR',status)), end
      status=mexnc('put_att_double', ncid, data_varid10, 'missing_value', nc_double', 1, NaN );
    end    
  else
    % Photosynthesis for PFT2
    ['MOPS + PFT']
    [data_varid9, status] = mexnc('DEF_VAR', ncid, 'Photo2', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
    if status, error(mexnc('STRERROR',status)), end
    status=mexnc('put_att_double', ncid, data_varid9, 'missing_value', nc_double', 1, NaN );
    if status, error(mexnc('STRERROR',status)), end
  end
else
  if useORGCARBON
    if useFLEXCP
      ['MOPS + (CARBON assumed) + ORGCARBON + FLEXCP']
      % SEDIMENTATION of POC
      [data_varid9, status] = mexnc('DEF_VAR', ncid, 'Sediment_C', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
      if status, error(mexnc('STRERROR',status)), end
      status=mexnc('put_att_double', ncid, data_varid9, 'missing_value', nc_double', 1, NaN );
      % Phytoplankton C:P uptake 
      [data_varid10, status] = mexnc('DEF_VAR', ncid, 'CP_PHYTO_UPTAKE', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
      if status, error(mexnc('STRERROR',status)), end
      status=mexnc('put_att_double', ncid, data_varid10, 'missing_value', nc_double', 1, NaN );
      % Zooplankton C:P uptake 
      [data_varid11, status] = mexnc('DEF_VAR', ncid, 'CP_ZOO_UPTAKE', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
      if status, error(mexnc('STRERROR',status)), end
      status=mexnc('put_att_double', ncid, data_varid11, 'missing_value', nc_double', 1, NaN );
    else
      ['MOPS + (CARBON assumed) + ORGCARBON']
      % SEDIMENTATION of POC
      [data_varid9, status] = mexnc('DEF_VAR', ncid, 'Sediment_C', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
      if status, error(mexnc('STRERROR',status)), end
      status=mexnc('put_att_double', ncid, data_varid9, 'missing_value', nc_double', 1, NaN );
    end
  else
    ['MOPS']
  end 
end




status = mexnc('PUT_ATT_TEXT', ncid, lon_varid, 'units', nc_char, length('degrees_east'), 'degrees_east');
if status, error(mexnc('STRERROR',status)), end
status = mexnc('PUT_ATT_TEXT', ncid, lat_varid, 'units', nc_char, length('degrees_north'), 'degrees_north');
if status, error(mexnc('STRERROR',status)), end
status = mexnc('PUT_ATT_TEXT', ncid, dep_varid, 'units', nc_char, length('meter'), 'meter');
if status, error(mexnc('STRERROR',status)), end

status = mexnc('enddef', ncid );
if status, error(mexnc('STRERROR',status)), end

% write the coordinate variables
status = mexnc('put_var_double', ncid, lon_varid, x );
if status, error(mexnc('STRERROR',status)), end

status = mexnc('put_var_double', ncid, lat_varid, y );
if status, error(mexnc('STRERROR',status)), end

status = mexnc('put_var_double', ncid, dep_varid, z );
if status, error(mexnc('STRERROR',status)), end

status = mexnc('put_var_double', ncid, time_varid, t );
if status, error(mexnc('STRERROR',status)), end

% write the data

%first diagnostic
status = mexnc('put_var_double', ncid, data_varid1, V1 );
if status, error(mexnc('STRERROR',status)), end

%second diagnostic
status = mexnc('put_var_double', ncid, data_varid2, V2 );
if status, error(mexnc('STRERROR',status)), end

%third diagnostic
status = mexnc('put_var_double', ncid, data_varid3, V3 );
if status, error(mexnc('STRERROR',status)), end

%fourth diagnostic
status = mexnc('put_var_double', ncid, data_varid4, V4 );
if status, error(mexnc('STRERROR',status)), end

%fifth diagnostic
status = mexnc('put_var_double', ncid, data_varid5, V5 );
if status, error(mexnc('STRERROR',status)), end

%sixth diagnostic
status = mexnc('put_var_double', ncid, data_varid6, V6 );
if status, error(mexnc('STRERROR',status)), end

%seventh diagnostic
status = mexnc('put_var_double', ncid, data_varid7, V7 );
if status, error(mexnc('STRERROR',status)), end

%eighth diagnostic
status = mexnc('put_var_double', ncid, data_varid8, V8 );
if status, error(mexnc('STRERROR',status)), end

%useORGCARBON
% %nineth diagnostic
% status = mexnc('put_var_double', ncid, data_varid9, V9 );
% if status, error(mexnc('STRERROR',status)), end
%
% %tenth diagnostic
% status = mexnc('put_var_double', ncid, data_varid10, V10 );
% if status, error(mexnc('STRERROR',status)), end
%
% %eleventh diagnostic
% status = mexnc('put_var_double', ncid, data_varid11, V11 );
% if status, error(mexnc('STRERROR',status)), end
%

if usePFT  
  if useORGCARBON   
    if useFLEXCP
      ['MOPS + (CARBON assumed) + ORGCARBON + PFT + FLEXCP']
      status = mexnc('put_var_double', ncid, data_varid9, V9 );
      if status, error(mexnc('STRERROR',status)), end
      status = mexnc('put_var_double', ncid, data_varid10, V10 );
      if status, error(mexnc('STRERROR',status)), end
      status = mexnc('put_var_double', ncid, data_varid11, V11 );
      if status, error(mexnc('STRERROR',status)), end
      status = mexnc('put_var_double', ncid, data_varid12, V12 );
      if status, error(mexnc('STRERROR',status)), end
      status = mexnc('put_var_double', ncid, data_varid13, V13 );
      if status, error(mexnc('STRERROR',status)), end
    else
      ['MOPS + (CARBON assumed) + ORGCARBON + PFT']
      status = mexnc('put_var_double', ncid, data_varid9, V9 );
      if status, error(mexnc('STRERROR',status)), end
      status = mexnc('put_var_double', ncid, data_varid10, V10 );
      if status, error(mexnc('STRERROR',status)), end
    end    
  else
    ['MOPS + PFT']
      status = mexnc('put_var_double', ncid, data_varid9, V9 );
      if status, error(mexnc('STRERROR',status)), end
  end
else
  if useORGCARBON
    if useFLEXCP
      ['MOPS + (CARBON assumed) + ORGCARBON + FLEXCP']
      status = mexnc('put_var_double', ncid, data_varid9, V9 );
      if status, error(mexnc('STRERROR',status)), end
      status = mexnc('put_var_double', ncid, data_varid10, V10 );
      if status, error(mexnc('STRERROR',status)), end
      status = mexnc('put_var_double', ncid, data_varid11, V11 );
      if status, error(mexnc('STRERROR',status)), end
    else
      ['MOPS + (CARBON assumed) + ORGCARBON']
      status = mexnc('put_var_double', ncid, data_varid9, V9 );
      if status, error(mexnc('STRERROR',status)), end
    end
  else
    ['MOPS']
  end 
end

status = mexnc('close', ncid );
if status, error(mexnc('STRERROR',status)), end

