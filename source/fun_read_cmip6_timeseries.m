function dat = fun_read_cmip6_timeseries(str, varname, base)
% fun_read_cmip6_timeseries
%
% Read CMIP6 variable-per-file and return 3D array [lat,lon,time].
% Works with CDO-regridded rectilinear files. Handles common dim orders.
%
% Inputs:
%   str     - muffingen GCM struct (uses str(1).path and str(1).exp)
%   varname - CMIP6 variable id (e.g., 'uas','vas','sfcWind')
%   base    - filename base used in EXAMPLE (e.g., '1deg' or Amon_* suffix)
%
% Output:
%   dat     - double [lat,lon,time]

ncfile = [str(1).path '/' str(1).exp '/' base '_' varname '.nc'];
ncid = netcdf.open(ncfile,'nowrite');

% Read variable
varid = netcdf.inqVarID(ncid,varname);
raw = double(netcdf.getVar(ncid,varid));

% Coordinates
try, latid = netcdf.inqVarID(ncid,'lat'); catch, latid = netcdf.inqVarID(ncid,'latitude'); end
try, lonid = netcdf.inqVarID(ncid,'lon'); catch, lonid = netcdf.inqVarID(ncid,'longitude'); end
latv = double(netcdf.getVar(ncid,latid)); latv = latv(:);
lonv = double(netcdf.getVar(ncid,lonid)); lonv = lonv(:);

% Time length (if present)
nt = 1;
% Identify dims attached to var
[~,~,dimids,~] = netcdf.inqVar(ncid,varid);
for k = 1:numel(dimids)
    [dname, dlen] = netcdf.inqDim(ncid, dimids(k)); %#ok<ASGLU>
    if strcmp(dname,'time') || strcmp(dname,'Time')
        nt = dlen; break;
    end
end

netcdf.close(ncid);

% Shape to [lat,lon,time]
s = size(raw);
nlat = numel(latv); nlon = numel(lonv);
raw = squeeze(raw);
s = size(raw);

if numel(s) == 2
    % two-dimensional; add singleton time
    if isequal(s,[nlat nlon])
        dat = reshape(raw,[nlat nlon 1]);
    elseif isequal(s,[nlon nlat])
        dat = reshape(raw',[nlat nlon 1]);
    else
        error('Unexpected 2D shape for %s in %s', varname, ncfile);
    end
elseif numel(s) == 3
    % attempt to permute to [lat lon time]
    if isequal(s,[nlat nlon nt])
        dat = raw;
    elseif isequal(s,[nlon nlat nt])
        dat = permute(raw,[2 1 3]);
    elseif isequal(s,[nt nlat nlon])
        dat = permute(raw,[2 3 1]);
    elseif isequal(s,[nt nlon nlat])
        dat = permute(raw,[3 2 1]);
    else
        error('Unexpected 3D shape for %s in %s', varname, ncfile);
    end
else
    % higher dims not expected after CDO regridding
    error('Unsupported dimensionality for %s in %s', varname, ncfile);
end

