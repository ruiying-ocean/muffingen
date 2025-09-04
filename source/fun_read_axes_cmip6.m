function [axisloncm,axislonce,axislatcm,axislatce,axislonam,axislonae,axislatam,axislatae] = fun_read_axes_cmip6(str)
%
%%

% *********************************************************************** %
% *** READ AND RETURN AXES INFO (CMIP6) ********************************* %
% *********************************************************************** %
%
% str input KEY:
% str(2).nc == par_nc_axes_name (e.g., deptho_1deg)
% str(1).path == par_pathin;  str(1).exp == par_expid
%
% Approach:
% - Read lat/lon from CDO-regridded rectilinear file
% - Construct edge coordinates by extending midpoints by half a grid step
% - Assumes all CMIP6 data preprocessed with CDO to 1-degree rectilinear grid
%
% *********************************************************************** %

% Read ocean grid axes from the bathymetry file
nc_path = [str(1).path '/' str(1).exp '/' str(2).nc '.nc'];

try
    ncid = netcdf.open(nc_path,'nowrite');
    
    % Read latitude
    try
        vlat = netcdf.inqVarID(ncid,'lat');
    catch
        vlat = netcdf.inqVarID(ncid,'latitude');
    end
    latv = double(netcdf.getVar(ncid,vlat));
    
    % Read longitude
    try
        vlon = netcdf.inqVarID(ncid,'lon');
    catch
        vlon = netcdf.inqVarID(ncid,'longitude');
    end
    lonv = double(netcdf.getVar(ncid,vlon));
    
    netcdf.close(ncid);
    
    % Verify rectilinear (should be 1D after CDO regridding)
    if ~(isvector(latv) && isvector(lonv))
        error('CMIP6 data should be CDO-regridded to rectilinear grid first');
    end
    
    % Ensure column vectors
    latv = latv(:);
    lonv = lonv(:);
    
    % Keep latitude in ascending order (standard NetCDF convention)
    % The data will be flipped in fun_read_topomask_cmip6 to match muffingen's expectation
    if numel(latv) > 1 && latv(2) < latv(1)
        latv = flipud(latv);
    end
    
    % Construct edge coordinates from midpoints
    % For longitude
    if numel(lonv) > 1
        dlon = diff(lonv);
        % Handle wraparound for global grids
        if abs(lonv(end) - lonv(1) + dlon(1)) < 2  % roughly 360 degrees
            lon_edges = [lonv(1)-dlon(1)/2; (lonv(1:end-1)+lonv(2:end))/2; lonv(end)+dlon(end)/2];
        else
            lon_edges = [lonv(1)-dlon(1)/2; (lonv(1:end-1)+lonv(2:end))/2; lonv(end)+dlon(end)/2];
        end
    else
        lon_edges = [lonv(1)-0.5; lonv(1)+0.5];
    end
    
    % For latitude
    if numel(latv) > 1
        dlat = diff(latv);
        lat_edges = [latv(1)-dlat(1)/2; (latv(1:end-1)+latv(2:end))/2; latv(end)+dlat(end)/2];
    else
        lat_edges = [latv(1)-0.5; latv(1)+0.5];
    end
    
    % Set ocean grid coordinates
    axisloncm = lonv;
    axislatcm = latv;
    axislonce = lon_edges;
    axislatce = lat_edges;
    
catch
    error('Could not read coordinates from CDO-regridded file: %s', nc_path);
end

% Read atmospheric grid coordinates if available
% Try wind files first
axislonam = [];
axislatam = [];

% Check for wind files to get atmospheric grid
wind_files = {'tauu', 'tauv', 'uas', 'vas'};
for i = 1:length(wind_files)
    if numel(str) >= 5 && ~isempty(str(5).nc)
        wind_path = [str(1).path '/' str(1).exp '/' wind_files{i} '_' str(5).nc '.nc'];
    elseif numel(str) >= 3 && ~isempty(str(3).nc)
        wind_path = [str(1).path '/' str(1).exp '/' wind_files{i} '_' str(3).nc '.nc'];
    else
        continue;
    end
    
    if exist(wind_path,'file') == 2
        try
            ncid = netcdf.open(wind_path,'nowrite');
            
            % Read atmospheric lat/lon
            try
                varid = netcdf.inqVarID(ncid,'lat');
            catch
                varid = netcdf.inqVarID(ncid,'latitude');
            end
            latv_a = double(netcdf.getVar(ncid,varid));
            
            try
                varid = netcdf.inqVarID(ncid,'lon');
            catch
                varid = netcdf.inqVarID(ncid,'longitude');
            end
            lonv_a = double(netcdf.getVar(ncid,varid));
            
            netcdf.close(ncid);
            
            % Should also be 1D after CDO regridding
            if isvector(latv_a) && isvector(lonv_a)
                axislatam = latv_a(:);
                axislonam = lonv_a(:);
                break; % Found valid atmospheric coordinates
            end
            
        catch
            try, netcdf.close(ncid); end %#ok<TRYNC>
            continue;
        end
    end
end

% If no atmospheric grid found, use ocean grid
if isempty(axislatam) || isempty(axislonam)
    axislatam = axislatcm;
    axislonam = axisloncm;
end

% Construct atmospheric edge coordinates
if numel(axislonam) > 1
    dlon_a = diff(axislonam);
    axislonae = [axislonam(1)-dlon_a(1)/2; (axislonam(1:end-1)+axislonam(2:end))/2; axislonam(end)+dlon_a(end)/2];
else
    axislonae = [axislonam(1)-0.5; axislonam(1)+0.5];
end

if numel(axislatam) > 1
    dlat_a = diff(axislatam);
    axislatae = [axislatam(1)-dlat_a(1)/2; (axislatam(1:end-1)+axislatam(2:end))/2; axislatam(end)+dlat_a(end)/2];
else
    axislatae = [axislatam(1)-0.5; axislatam(1)+0.5];
end

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %