% MSAC_SPLITWAVE - Create a shear-wave split using the FORTRAN
%                  code sacsplitwave
function [trN,trE,trZ]=msac_splitwave(fast,tlag,spol,noise)

wdir = '/tmp/msac_splitwave_wdir' ;

exec = '/Users/ja17375/Ext_Programs/seismic/bin/sacsplitwave' ;

% Build up the command line 
cmd_str = [ ...
   'rm -rf ' wdir '; ' ...
   'mkdir ' wdir '; ' ...
   'cd ' wdir '; ' ...
   exec ...
   sprintf(' -op %7.2f %6.3f',fast, tlag) ...
   sprintf(' -spol %7.2f ',spol) ...
   sprintf(' -noise %9.7f',noise)] ;

% run executable
[status,result] = system(cmd_str) ;
if status
   error(['MSAC_SPLITWAVE: System call failed:' result] )
end

% ingest data
[trE] = msac_read([wdir '/SWAV.BHE']) ;
[trN] = msac_read([wdir '/SWAV.BHN']) ;
[trZ] = msac_read([wdir '/SWAV.BHZ']) ;

end