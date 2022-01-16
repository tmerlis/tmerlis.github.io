% created August 16, 2012 by Tim Merlis; uses native matlab netcdf support
% see http://www.mathworks.com/access/helpdesk/help/techdoc/ref/netcdf.html

infilename = ['in.nc'];
outfilename = ['out.nc'];

% create output netcdf
ncid_out = netcdf.create(outfilename,'NC_NOCLOBBER');

% open input netcdf
ncid_in  = netcdf.open(infilename,'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_in);

% first loop over dims 
for d=0:ndims-1
  % get dim from input file
  [dimname, dimlen] = netcdf.inqDim(ncid_in,d);

  % make dim in output file
  if strcmp(dimname,'time')
    dimid = netcdf.defDim(ncid_out,dimname,netcdf.getConstant('NC_UNLIMITED'));
  else
    dimid = netcdf.defDim(ncid_out,dimname,dimlen);
  end
end

% next loop over vars
for v=0:nvars-1
  [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_in,v);

  out_varid = netcdf.defVar(ncid_out,varname,xtype,dimids);
  for attnum = 0:natts-1
    attname = netcdf.inqAttName(ncid_in,v,attnum);
    netcdf.copyAtt(ncid_in,v,attname,ncid_out,out_varid);
  end
end

% leave define mode and enter data mode
netcdf.endDef(ncid_out);

for v=0:nvars-1
  [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_in,v);
  var = netcdf.getVar(ncid_in,v);

  out_varid = netcdf.inqVarID(ncid_out,varname);

  if strcmp(varname,'time')
    netcdf.putVar(ncid_out,out_varid,0,length(var),var);
  else
    netcdf.putVar(ncid_out,out_varid,var);
  end
end

% close input and output netcdf
netcdf.close(ncid_in);
netcdf.close(ncid_out);

