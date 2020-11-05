function petsc2netcdf(netcdfFileName);
% Function for postprocessing petsc2 tracers averaged over certain time period to netcdf
% Edited by Tatsuro Tanioka (Jan-28-2020)
% To run this script, type in MATLAB command line
% > n7tracersavg28('test.nc')
% Units: Chemical tracers in mmol m-3. PHYTO, ZOO and DETRITUS in mmol P m-3.  
% Units: PHYTOC, ZOOC, POC, and DOC in mmol C m-3.  

%%%% Options for output
useCarbon = 1; % 1 if option DCARBON used, 0 if otherwise
useOrgCarbon = 1; % 1 if option DORGCARBON used, 0 if otherwise
%%%%
% Set toplevel path to GCMs configuration
basepath='~/TMM2/MITgcm_2.8deg';
%basepath='~/TMM2/MITgcm_ECCO';
% basepath='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';
% basepath='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO_v4';
% basepath = '/work_O1/smomw069/MIT2.8';

addpath(genpath('~/TMM/tmm_matlab_code'));

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

tr1=readPetscBinVec('po4avg.petsc',-1);
tr2=readPetscBinVec('dopavg.petsc',-1);
tr3=readPetscBinVec('oxyavg.petsc',-1);
tr4=readPetscBinVec('phyavg.petsc',-1);
tr5=readPetscBinVec('zooavg.petsc',-1);
tr6=readPetscBinVec('detavg.petsc',-1);
tr7=readPetscBinVec('no3avg.petsc',-1);
if useCarbon
tr8=readPetscBinVec('dicavg.petsc',-1);
tr9=readPetscBinVec('alkavg.petsc',-1);
end
if useOrgCarbon
tr10=readPetscBinVec('docavg.petsc',-1);
tr11=readPetscBinVec('pocavg.petsc',-1);
tr12=readPetscBinVec('phycavg.petsc',-1);
tr13=readPetscBinVec('zoocavg.petsc',-1);
end
nb=size(tr1,1);
[TR1,x,y,z]=matrixToGrid(tr1(Irr,:),[1:nb]',boxFile,gridFile);
[TR2,x,y,z]=matrixToGrid(tr2(Irr,:),[1:nb]',boxFile,gridFile);
[TR3,x,y,z]=matrixToGrid(tr3(Irr,:),[1:nb]',boxFile,gridFile);
[TR4,x,y,z]=matrixToGrid(tr4(Irr,:),[1:nb]',boxFile,gridFile);
[TR5,x,y,z]=matrixToGrid(tr5(Irr,:),[1:nb]',boxFile,gridFile);
[TR6,x,y,z]=matrixToGrid(tr6(Irr,:),[1:nb]',boxFile,gridFile);
[TR7,x,y,z]=matrixToGrid(tr7(Irr,:),[1:nb]',boxFile,gridFile);
if useCarbon
[TR8,x,y,z]=matrixToGrid(tr8(Irr,:),[1:nb]',boxFile,gridFile);
[TR9,x,y,z]=matrixToGrid(tr9(Irr,:),[1:nb]',boxFile,gridFile);
end
if useOrgCarbon
[TR10,x,y,z]=matrixToGrid(tr10(Irr,:),[1:nb]',boxFile,gridFile);
[TR11,x,y,z]=matrixToGrid(tr11(Irr,:),[1:nb]',boxFile,gridFile);
[TR12,x,y,z]=matrixToGrid(tr12(Irr,:),[1:nb]',boxFile,gridFile);
[TR13,x,y,z]=matrixToGrid(tr13(Irr,:),[1:nb]',boxFile,gridFile);
end
[nx,ny,nz,nt]=size(TR1);

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

%PO4
[data_varid1, status] = mexnc('DEF_VAR', ncid, 'PO4', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid1, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%DOP
[data_varid2, status] = mexnc('DEF_VAR', ncid, 'DOP', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid2, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%OXYGEN
[data_varid3, status] = mexnc('DEF_VAR', ncid, 'OXYGEN', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid3, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%PHYTOPLANKTON
[data_varid4, status] = mexnc('DEF_VAR', ncid, 'PHYTO', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid4, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%ZOOPLANKTON
[data_varid5, status] = mexnc('DEF_VAR', ncid, 'ZOO', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid5, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%DETRITUS
[data_varid6, status] = mexnc('DEF_VAR', ncid, 'DET', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid6, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%NO3
[data_varid7, status] = mexnc('DEF_VAR', ncid, 'NO3', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid7, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

if useCarbon

%DIC
[data_varid8, status] = mexnc('DEF_VAR', ncid, 'DIC', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid8, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%ALK
[data_varid9, status] = mexnc('DEF_VAR', ncid, 'ALK', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid9, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

end

if useOrgCarbon

%DOC
[data_varid10, status] = mexnc('DEF_VAR', ncid, 'DOC', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid10, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%POC
[data_varid11, status] = mexnc('DEF_VAR', ncid, 'POC', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid11, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%PHYTOC
[data_varid12, status] = mexnc('DEF_VAR', ncid, 'PHYTOC', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid12, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end

%ZOOC
[data_varid13, status] = mexnc('DEF_VAR', ncid, 'ZOOC', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
if status, error(mexnc('STRERROR',status)), end
status=mexnc('put_att_double', ncid, data_varid13, 'missing_value', nc_double', 1, NaN );
if status, error(mexnc('STRERROR',status)), end
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

%first tracer avg
status = mexnc('put_var_double', ncid, data_varid1, TR1 );
if status, error(mexnc('STRERROR',status)), end

%second tracer avg
status = mexnc('put_var_double', ncid, data_varid2, TR2 );
if status, error(mexnc('STRERROR',status)), end

%third tracer avg
status = mexnc('put_var_double', ncid, data_varid3, TR3 );
if status, error(mexnc('STRERROR',status)), end

%fourth tracer avg
status = mexnc('put_var_double', ncid, data_varid4, TR4 );
if status, error(mexnc('STRERROR',status)), end

%fifth tracer avg
status = mexnc('put_var_double', ncid, data_varid5, TR5 );
if status, error(mexnc('STRERROR',status)), end

%sixth tracer avg
status = mexnc('put_var_double', ncid, data_varid6, TR6 );
if status, error(mexnc('STRERROR',status)), end

%seventh tracer avg
status = mexnc('put_var_double', ncid, data_varid7, TR7 );
if status, error(mexnc('STRERROR',status)), end


if useCarbon
%eighth tracer avg
status = mexnc('put_var_double', ncid, data_varid8, TR8 );
if status, error(mexnc('STRERROR',status)), end
%nineth tracer avg
status = mexnc('put_var_double', ncid, data_varid9, TR9 );
if status, error(mexnc('STRERROR',status)), end
end
if useOrgCarbon
%tenth tracer avg
status = mexnc('put_var_double', ncid, data_varid10, TR10 );
if status, error(mexnc('STRERROR',status)), end
%eleventh tracer avg
status = mexnc('put_var_double', ncid, data_varid11, TR11 );
if status, error(mexnc('STRERROR',status)), end
%twelveth tracer avg
status = mexnc('put_var_double', ncid, data_varid12, TR12 );
if status, error(mexnc('STRERROR',status)), end
%thirteenth tracer avg
status = mexnc('put_var_double', ncid, data_varid13, TR13 );
if status, error(mexnc('STRERROR',status)), end

end

% closing
status = mexnc('close', ncid );
if status, error(mexnc('STRERROR',status)), end

