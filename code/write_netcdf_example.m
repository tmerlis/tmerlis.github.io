% created April 1, 2010 by Tim Merlis; uses native matlab netcdf support
% see http://www.mathworks.com/access/helpdesk/help/techdoc/ref/netcdf.html

% -------------------
% get vars from input netcdf
% -------------------
infilename = 'input_file.nc';

% open input netcdf
ncid  = netcdf.open(infilename,'NC_NOWRITE');

% get dimension variables
varid = netcdf.inqVarID(ncid, 'lat');
lat = netcdf.getVar(ncid, varid);
varid = netcdf.inqVarID(ncid, 'lon');
lon = netcdf.getVar(ncid, varid);
varid = netcdf.inqVarID(ncid, 'time');
time = netcdf.getVar(ncid, varid);

varid = netcdf.inqVarID(ncid, 't_surf');
t_surf = netcdf.getVar(ncid, varid);

% close input netcdf
netcdf.close(ncid);

% -------------------
% write output netcdf
% -------------------
outfilename = 'test_output.nc';

% create netcdf
ncid = netcdf.create(outfilename,'NC_NOCLOBBER');

% define dimensions
timedimid = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
latdimid = netcdf.defDim(ncid,'lat',length(lat));
londimid = netcdf.defDim(ncid,'lon',length(lon));

% define variables
timevarid = netcdf.defVar(ncid,'time','double',timedimid);
latvarid = netcdf.defVar(ncid,'lat','double',latdimid);
lonvarid = netcdf.defVar(ncid,'lon','double',londimid);
tsurfvarid = netcdf.defVar(ncid,'t_surf','double',[londimid latdimid timedimid]);

% attribute
netcdf.putAtt(ncid,timevarid,'dimension','days since tidal locking');
netcdf.putAtt(ncid,latvarid,'dimension','deg');
netcdf.putAtt(ncid,lonvarid,'dimension','deg');
netcdf.putAtt(ncid,tsurfvarid,'dimension','deg K');

% leave define mode and enter data mode
netcdf.endDef(ncid);

% write data
netcdf.putVar(ncid,timevarid,0,length(time),time);
netcdf.putVar(ncid,latvarid,lat);
netcdf.putVar(ncid,lonvarid,lon);
netcdf.putVar(ncid,tsurfvarid,t_surf);

% close output netcdf
netcdf.close(ncid);
