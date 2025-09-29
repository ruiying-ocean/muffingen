function fld = fun_read_cmip6_annual_mean(str, varname, base)
% fun_read_cmip6_annual_mean
%
% Read CMIP6 variable-per-file and return a [lat,lon] 2D annual mean.
% Uses CDO-regridded rectilinear files; averages over 'time' if present.
%
% Inputs:
%   str     - muffingen GCM struct (uses str(1).path and str(1).exp)
%   varname - CMIP6 variable id (e.g., 'tauu','tauv','uas','vas','sfcWind')
%   base    - filename base used in EXAMPLE (e.g., '1deg' or Amon_* suffix)
%
% Output:
%   fld     - 2D double [lat,lon], annual mean

ncfile = [str(1).path '/' str(1).exp '/' base '_' varname '.nc'];
ncid = netcdf.open(ncfile,'nowrite');

% Read variable
varid = netcdf.inqVarID(ncid,varname);
dat = double(netcdf.getVar(ncid,varid));

% Coordinates
try, latid = netcdf.inqVarID(ncid,'lat'); catch, latid = netcdf.inqVarID(ncid,'latitude'); end
try, lonid = netcdf.inqVarID(ncid,'lon'); catch, lonid = netcdf.inqVarID(ncid,'longitude'); end
latv = double(netcdf.getVar(ncid,latid)); latv = latv(:);
lonv = double(netcdf.getVar(ncid,lonid)); lonv = lonv(:);

% Detect time dimension by name and average over it
% Identify dims attached to var
[~,~,dimids,~] = netcdf.inqVar(ncid,varid);
tdim = [];
for k = 1:numel(dimids)
    [dname, ~] = netcdf.inqDim(ncid, dimids(k));
    if strcmp(dname,'time') || strcmp(dname,'Time')
        tdim = k; break;
    end
end

if ~isempty(tdim) && ndims(dat) >= tdim
    ord = 1:ndims(dat);
    ord = [tdim setdiff(ord,tdim)];
    datp = permute(dat, ord);
    datp = squeeze(nanmean(datp,1));
    dat = datp;
else
    dat = squeeze(dat);
end

% Ensure 2D and orient to [lat,lon]
s = size(dat);
nlat = numel(latv); nlon = numel(lonv);
if numel(s) ~= 2
    dat = squeeze(dat);
    s = size(dat);
end
if isequal(s, [nlon nlat])
    dat = dat';
elseif isequal(s, [nlat nlon])
    % ok
else
    if isequal(fliplr(s), [nlat nlon])
        dat = dat';
    end
end

netcdf.close(ncid);
fld = dat;

