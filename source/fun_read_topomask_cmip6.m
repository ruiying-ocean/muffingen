function [topo,mask]  = fun_read_topomask_cmip6(str)
%
%%

% *********************************************************************** %
% *** READ AND RETURN CMIP6 TOPOGRAPHY AND MASK ************************ %
% *********************************************************************** %
%
% str input KEY:
% str(2).nc == par_nc_topo_name  (e.g., deptho_1deg)
% str(4).nc == par_nc_mask_name  (e.g., sftof_1deg)
% str(1).path == par_pathin;  str(1).exp == par_expid
%
% Behavior:
% - Reads deptho (m, positive depth) from CDO-regridded rectilinear file
% - If sftof/sftlf exists, build mask from percent (>50% ocean)
% - Else derive mask from deptho > 0
% - Assumes data is already on rectilinear 1-degree grid from CDO preprocessing
%
% *********************************************************************** %

% --- Read deptho bathymetry ---
nc_topo = [str(1).path '/' str(1).exp '/' str(2).nc '.nc'];
ncid = netcdf.open(nc_topo,'nowrite');
varid = netcdf.inqVarID(ncid,'deptho'); 
topo = double(netcdf.getVar(ncid,varid));

% Handle fill/missing values
try
    fv = netcdf.getAtt(ncid,varid,'_FillValue');
    topo(topo>=fv*0.9) = NaN;
catch
    topo(topo>1.0e10) = NaN;
end

% Read coordinates (should be 1D vectors from CDO regridding)
try
    vlat = netcdf.inqVarID(ncid,'lat'); 
    lat = double(netcdf.getVar(ncid,vlat));
catch
    try
        vlat = netcdf.inqVarID(ncid,'latitude'); 
        lat = double(netcdf.getVar(ncid,vlat));
    catch
        error('Could not find latitude coordinate in CDO-regridded file');
    end
end

try
    vlon = netcdf.inqVarID(ncid,'lon'); 
    lon = double(netcdf.getVar(ncid,vlon));
catch
    try
        vlon = netcdf.inqVarID(ncid,'longitude'); 
        lon = double(netcdf.getVar(ncid,vlon));
    catch
        error('Could not find longitude coordinate in CDO-regridded file');
    end
end
netcdf.close(ncid);

% Verify rectilinear grid (should be 1D after CDO processing)
if ~(isvector(lat) && isvector(lon))
    error('CMIP6 data should be CDO-regridded to rectilinear grid first');
end

% Ensure column vectors
lat = lat(:); lon = lon(:);
nlat = length(lat); nlon = length(lon);

% Check data orientation and transpose if needed
sz = size(topo);
if sz(1) == nlon && sz(2) == nlat
    topo = topo'; % transpose [lon,lat] -> [lat,lon]
elseif sz(1) == nlat && sz(2) == nlon
    % already [lat,lon], keep as is
else
    error('Topography dimensions [%d,%d] do not match coordinates lat=%d, lon=%d', sz(1), sz(2), nlat, nlon);
end

% Always flip data to match muffingen's convention
% Follow same pattern as other GCM functions (foam, cesm, rockee)
topo = flipud(topo);

% Set land points to NaN (positive depths are ocean)
topo(topo <= 0.0) = NaN;

% Initialize mask
mask = [];

% Try to read ocean/land fraction mask if available
if numel(str) >= 4 && ~isempty(str(4).nc)
    nc_mask = [str(1).path '/' str(1).exp '/' str(4).nc '.nc'];
    if exist(nc_mask,'file') == 2
        try
            ncid = netcdf.open(nc_mask,'nowrite');
            
            % Try sftof (ocean fraction) first, then sftlf (land fraction)
            varname = 'sftof';
            try
                vid = netcdf.inqVarID(ncid,'sftof');
            catch
                vid = netcdf.inqVarID(ncid,'sftlf');
                varname = 'sftlf';
            end
            
            sft = double(netcdf.getVar(ncid,vid));
            
            % Handle fill/missing values
            try
                mv = netcdf.getAtt(ncid,vid,'_FillValue');
                sft(abs(sft-mv)<abs(mv)*0.1 | sft>=mv*0.9) = NaN;
            catch
                sft(sft>1.0e10) = NaN;
            end
            netcdf.close(ncid);
            
            % If a time dimension exists, collapse to annual/overall mean
            szm = size(sft);
            if numel(szm) > 2
                % Try common layouts and reduce over the non-(lat,lon) dim
                if isequal(szm, [nlat, nlon, szm(3)])
                    % [lat,lon,time] -> mean over time
                    sft = mean(sft, 3, 'omitnan');
                elseif isequal(szm, [nlon, nlat, szm(3)])
                    % [lon,lat,time] -> transpose to [lat,lon] then mean
                    sft = permute(sft, [2 1 3]);
                    sft = mean(sft, 3, 'omitnan');
                elseif numel(szm) == 3 && szm(1) ~= nlat && szm(1) ~= nlon && szm(2) == nlat && szm(3) == nlon
                    % [time,lat,lon]
                    sft = squeeze(mean(sft, 1, 'omitnan'));
                elseif numel(szm) == 3 && szm(1) ~= nlat && szm(1) ~= nlon && szm(2) == nlon && szm(3) == nlat
                    % [time,lon,lat]
                    sft = squeeze(mean(sft, 1, 'omitnan'));
                    sft = sft';
                else
                    error('Mask dimensions %s not compatible with lat=%d, lon=%d', mat2str(szm), nlat, nlon);
                end
                szm = size(sft);
            end

            % Check for corrupted data (all values identical)
            finite_vals = sft(isfinite(sft));
            if ~isempty(finite_vals) && all(abs(finite_vals - finite_vals(1)) < 1e-6)
                fprintf('Warning: %s data corrupted (all values = %.1f%%), using bathymetry-derived mask\n', varname, finite_vals(1));
                sft = [];
            end
            
            if ~isempty(sft)
                % Orient mask field to [lat,lon] 
                szm = size(sft);
                if szm(1) == nlon && szm(2) == nlat
                    sft = sft'; % transpose [lon,lat] -> [lat,lon]
                elseif szm(1) == nlat && szm(2) == nlon
                    % already [lat,lon]
                else
                    error('Mask dimensions [%d,%d] do not match coordinates lat=%d, lon=%d', szm(1), szm(2), nlat, nlon);
                end
                
                % Apply same flip as topography for consistent orientation
                sft = flipud(sft);
                
                % Handle missing values
                if strcmp(varname,'sftof')
                    sft(~isfinite(sft)) = 0.0; % missing ocean fraction = land
                else
                    sft(~isfinite(sft)) = 100.0; % missing land fraction = land
                end
                
                % Convert percentage to fraction
                if max(sft(:),[],'omitnan') > 1.0
                    sft = sft/100.0;
                end
                
                % Create binary mask (1=ocean, 0=land)
                if strcmp(varname,'sftof')
                    mask = double(sft >= 0.5); % ocean fraction >= 50%
                else
                    mask = double(sft < 0.5);  % land fraction < 50% = ocean
                end
                
                % Set land bathymetry to NaN
                topo(mask==0) = NaN;
            end
            
        catch ME
            try, netcdf.close(ncid); end %#ok<TRYNC>
            fprintf('Warning: Could not read mask file %s (%s), using bathymetry-derived mask\n', nc_mask, ME.message);
        end
    end
end

% If no mask was created, derive from bathymetry
if isempty(mask)
    mask = double(~isnan(topo));
end

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
