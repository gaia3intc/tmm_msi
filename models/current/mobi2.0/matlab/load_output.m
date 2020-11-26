% base_path='/data2/spk/UVic_OSU_Matrix/LGM_WindPerturbation_Experiments/no_embm_awind2/picdefault';
base_path='~/TMM2/UVicOSUpicdefault';
addpath(genpath('~/TMM2/tmm_matlab_code'));

load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','x','y','z','gridType')

load(profilesFile,'Irr')

timeFile='output_time.txt';
[hdr,tdat]=hdrload(timeFile);
Tavg = tdat(:,2)';

trNames=readtable('MOBI_tracer_names.txt','ReadVariableNames',0);
trNames=table2cell(trNames);

numTr=length(trNames);

if strcmp(gridType,'llc_v4')
  load(fullfile(base_path,'llc_v4_grid'))
  gcmfaces_global
end

for itr=1:numTr
  varName=upper(trNames{itr})
  fn=[trNames{itr} '.petsc'];
%  tr=readPetscBinVec(fn,1,-1);
  tr=readPetscBinVec(fn,-1);
%  TR=matrixToGrid(tr(Irr,end),[],boxFile,gridFile);
  TR=matrixToGrid(tr(Irr,:),[],boxFile,gridFile);

  if strcmp(gridType,'llc_v4')
    error('LLC_v4 not supported!')
%    varName=[varName '_plot'];
%	tmp=gcmfaces(TR);
%	[x,y,TRplot]=convert2pcol(mygrid.XC,mygrid.YC,tmp);
%	[n1,n2]=size(TRplot);
%	eval([trNames{itr} '=zeros([n1 n2 nz]);']);
%	for iz=1:nz
%	  eval(['[x,y,' varName '(:,:,iz)]=convert2pcol(mygrid.XC,mygrid.YC,tmp(:,:,iz));']);
%	end
  else
	eval([varName '=TR;']);
  end
%  save(varName,varName,'x','y','z')
  save(varName,varName,'x','y','z','Tavg')
%  write2netcdf([varName '.nc'],TR,x,y,z,[],upper(varName))
  write2netcdf([varName '.nc'],TR,x,y,z,Tavg,upper(varName))
end
