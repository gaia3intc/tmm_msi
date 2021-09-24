function write_ts28(netcdfFileName);
%function petsc2netcdf(netcdfFileName);
%%%% Functions for making monthly mean Temperature and Salinity into net cdf
% By Tatsuro Tanioka (2020-01-28)
% To run this script, type in MATLAB command line
% > n7physics('test.nc')

['write_ts28']

base_path='~/TMM2/MITgcm_2.8deg';
%base_path='~/TMM2/MITgcm_ECCO';
addpath(genpath('~/TMM2/tmm_matlab_code'));

load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','x','y','z','gridType')

load(profilesFile,'Irr')

% Convert temperature and salinity files to netCDF
tempNames = {'Ts_00','Ts_01','Ts_02','Ts_03','Ts_04','Ts_05','Ts_06','Ts_07','Ts_08','Ts_09','Ts_10','Ts_11'};
saltNames = {'Ss_00','Ss_01','Ss_02','Ss_03','Ss_04','Ss_05','Ss_06','Ss_07','Ss_08','Ss_09','Ss_10','Ss_11'};
numTS = length(tempNames);
for its=1:numTS;
  fn_temp=[tempNames{its}];
  fn_salt=[saltNames{its}];
  tr_temp=readPetscBinVec(fn_temp,1,-1);
  tr_salt=readPetscBinVec(fn_salt,1,-1);
  TR_temp=matrixToGrid(tr_temp(Irr,end),[],boxFile,gridFile);
  TR_salt=matrixToGrid(tr_salt(Irr,end),[],boxFile,gridFile);
  TEMP(its,:,:,:) = TR_temp;
  SALT(its,:,:,:) = TR_salt;
end
TEMP_out = permute(TEMP,[2 3 4 1]);
SALT_out = permute(SALT,[2 3 4 1]);
t= [1:12]';

[nx,ny,nz,nt]=size(TEMP_out);

% Calculate mean annual TEMP and SALINITY for each i,j,k
TEMP_MEAN = nanmean(TEMP_out,4);
SALT_MEAN = nanmean(SALT_out,4);

% open file
[ncid, status] = mexnc('CREATE', netcdfFileName, 'noclobber' );
% define dimensions
[lon_dimid, status] = mexnc('DEF_DIM',ncid,'Longitude',nx);
[lat_dimid, status] = mexnc('DEF_DIM',ncid,'Latitude',ny);
[dep_dimid, status] = mexnc('DEF_DIM',ncid,'Depth',nz);
[time_dimid, status] = mexnc('DEF_DIM',ncid,'Time',nt);
% Define the coordinate variables
[lon_varid, status] = mexnc('DEF_VAR', ncid, 'Longitude', nc_double, 1, lon_dimid );
[lat_varid, status] = mexnc('DEF_VAR', ncid, 'Latitude', nc_double, 1, lat_dimid );
[dep_varid, status] = mexnc('DEF_VAR', ncid, 'Depth', nc_double, 1, dep_dimid );
[time_varid, status] = mexnc('DEF_VAR', ncid, 'Time', nc_double, 1, time_dimid );

% Temperature
[data_varid1, status] = mexnc('DEF_VAR', ncid, 'TEMP', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
status=mexnc('put_att_double', ncid, data_varid1, 'missing_value', nc_double', 1, NaN );

% Salinity
[data_varid2, status] = mexnc('DEF_VAR', ncid, 'SALT', nc_double, 4, [time_dimid dep_dimid lat_dimid lon_dimid] );
status=mexnc('put_att_double', ncid, data_varid2, 'missing_value', nc_double', 1, NaN );

% Mean annual Temperature
[data_varid3, status] = mexnc('DEF_VAR', ncid, 'TEMP_MEAN', nc_double, 3, [dep_dimid lat_dimid lon_dimid] );
status=mexnc('put_att_double', ncid, data_varid3, 'missing_value', nc_double', 1, NaN );

% Mean annual Salinity
[data_varid4, status] = mexnc('DEF_VAR', ncid, 'SALT_MEAN', nc_double, 3, [dep_dimid lat_dimid lon_dimid] );
status=mexnc('put_att_double', ncid, data_varid4, 'missing_value', nc_double', 1, NaN );

status = mexnc('PUT_ATT_TEXT', ncid, lon_varid, 'units', nc_char, length('degrees_east'), 'degrees_east');
status = mexnc('PUT_ATT_TEXT', ncid, lat_varid, 'units', nc_char, length('degrees_north'), 'degrees_north');
status = mexnc('PUT_ATT_TEXT', ncid, dep_varid, 'units', nc_char, length('meter'), 'meter');

status = mexnc('enddef', ncid );

% write the coordinate variables
status = mexnc('put_var_double', ncid, lon_varid, x );
status = mexnc('put_var_double', ncid, lat_varid, y );
status = mexnc('put_var_double', ncid, dep_varid, z );
status = mexnc('put_var_double', ncid, time_varid, t );

% write the data

%first tracer
status = mexnc('put_var_double', ncid, data_varid1, TEMP_out );

%second tracer
status = mexnc('put_var_double', ncid, data_varid2, SALT_out );

%third tracer
status = mexnc('put_var_double', ncid, data_varid3, TEMP_MEAN );

%fourth tracer
status = mexnc('put_var_double', ncid, data_varid4, SALT_MEAN );

% closing
status = mexnc('close', ncid );












