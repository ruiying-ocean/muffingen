function [grid_albd]  = fun_read_albd_cmip6(str)
%
%%

% *********************************************************************** %
% *** READ AND RETURN PLANETARY ALBEDO (CMIP6) ************************* %
% *********************************************************************** %
%
% str input KEY:
% str(3).nc == par_nc_atmos_name base (e.g., 1deg)
% Variables used: rsut (TOA outgoing SW), rsdt (TOA incident SW)
%
% *********************************************************************** %

% Read incident and outgoing shortwave at TOA
solin  = read_cmip6_2d_ann('rsdt', str(3).nc, 1e-12, str);
fsutoa = read_cmip6_2d_ann('rsut', str(3).nc, [],     str);

% Calculate albedo; reflected/incident with guarding
grid_albd = fsutoa ./ solin;
% Physically bound albedo
grid_albd(~isfinite(grid_albd)) = NaN;
grid_albd = max(0.0, min(1.0, grid_albd));

end

% -----------------------------------------------------------------------
% Local function: read a CMIP6 variable and return [LAT,LON] annual mean
% -----------------------------------------------------------------------
function fld = read_cmip6_2d_ann(varname, base, clamp_eps, str)
    ncfile = [str(1).path '/' str(1).exp '/' varname '_' base '.nc'];
    ncid = netcdf.open(ncfile,'nowrite');
    varid = netcdf.inqVarID(ncid,varname);
    [~, ~, dimids, ~] = netcdf.inqVar(ncid, varid);
    % Identify dimension positions by name
    timeDim = []; latDim = []; lonDim = [];
    for d = 1:numel(dimids)
        [dname, ~] = netcdf.inqDim(ncid, dimids(d));
        dl = lower(dname);
        if strcmp(dl,'time')
            timeDim = d;
        elseif strcmp(dl,'lat') || strcmp(dl,'latitude')
            latDim = d;
        elseif strcmp(dl,'lon') || strcmp(dl,'longitude')
            lonDim = d;
        end
    end
    dat = netcdf.getVar(ncid,varid);
    netcdf.close(ncid);
    dat = double(dat);
    % Bring to [time, lat, lon] or [lat, lon] for robust averaging
    nd = ndims(dat);
    % Build permutation to order as [time lat lon (others...)]
    front = [];
    if ~isempty(timeDim), front(end+1) = timeDim; end %#ok<AGROW>
    if ~isempty(latDim),  front(end+1) = latDim;  end %#ok<AGROW>
    if ~isempty(lonDim),  front(end+1) = lonDim;  end %#ok<AGROW>
    rest = setdiff(1:nd, front, 'stable');
    perm = [front rest];
    dat = permute(dat, perm);
    % Squeeze trailing singleton dims
    dat = squeeze(dat);
    % Average time if present
    if ndims(dat) == 3
        dat = squeeze(nanmean(dat, 1)); % now [lat, lon] or [lon, lat]
    else
        while ndims(dat) > 2
            dat = squeeze(nanmean(dat, 1));
        end
    end
    % Ensure [lat, lon]
    if size(dat,1) < size(dat,2) && ~isempty(latDim) && ~isempty(lonDim)
        % [lat,lon] common case for 1deg grids
    else
        dat = dat';
    end
    % Handle invalids and clamp if denominator
    dat(~isfinite(dat)) = NaN;
    if ~isempty(clamp_eps)
        dat(dat <= 0) = clamp_eps;
    end
    fld = dat;
end
