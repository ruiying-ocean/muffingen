function [] = make_grid_winds_cmip6(gilone,gilate,gimask,golonm,golone,golatm,golate,gomask,str,optplots)
% make_grid_winds_cmip6
%
%   *********************************************************
%   *** RE-GRID WIND PRODUCTS (CMIP6 variable-per-file)   ***
%   *********************************************************
%
%   str input KEY:
%   str(1).nc == par_nc_axes_name
%   str(2).nc == par_nc_topo_name
%   str(3).nc == par_nc_atmos_name   (base; e.g., Amon_SOURCE_EXP_MEMBER_GRID_TIMERANGE)
%   str(4).nc == par_nc_ocean_name
%   str(5).nc == par_nc_coupl_name   (use same base for tauu/tauv)
%
%   Input variables (CMIP6): tauu, tauv (Pa); uas, vas (m/s); sfcWind (m/s)
%

% determine output grid size (remember: [rows columns])
[jmax, imax] = size(gomask);

% (moved to separate helper function fun_read_cmip6_annual_mean.m)

% --- Wind stress (tauu/tauv): expected sign downward positive; follow CESM sign correction ---
try
    giwtauu = fun_read_cmip6_annual_mean(str,'tauu', str(5).nc);
    giwtauv = fun_read_cmip6_annual_mean(str,'tauv', str(5).nc);
catch
    % Fallback: use tauuo/tauvo (ocean) if available
    giwtauu = fun_read_cmip6_annual_mean(str,'tauuo', str(5).nc);
    giwtauv = fun_read_cmip6_annual_mean(str,'tauvo', str(5).nc);
end
% re-orientate to match other readers' convention before regrid
% note: make_regrid_2d expects input as [lon,lat] transposed; we supply post-flip later
% Sign: CMIP6 tauu/tauv follow CF 'downward, eastward/northward' stress.
% Do NOT negate here to match HadCM3 pattern and OUTPUT.EXAMPLES.

% --- Wind velocity (uas/vas) and optional sfcWind ---
giwvelu = fun_read_cmip6_annual_mean(str,'uas', str(3).nc);
giwvelv = fun_read_cmip6_annual_mean(str,'vas', str(3).nc);

% Wind speed: compute both metrics to mirror CESM
% 1) uvaa: magnitude from annual-mean components
giwspd_uvaa = (giwvelu.^2 + giwvelv.^2).^0.5;

% 2) uvma: monthly magnitude average if time series available; else fallback
try
    ts_u = fun_read_cmip6_timeseries(str,'uas', str(3).nc);   % [lat lon time]
    ts_v = fun_read_cmip6_timeseries(str,'vas', str(3).nc);
    nt = size(ts_u,3);
    if nt > 1
        acc = zeros(size(ts_u(:,:,1)));
        for t = 1:nt
            acc = acc + sqrt(ts_u(:,:,t).^2 + ts_v(:,:,t).^2)/nt;
        end
        giwspd_uvma = acc;
    else
        giwspd_uvma = giwspd_uvaa;
    end
catch
    giwspd_uvma = giwspd_uvaa;
end

% 3) wsma: REQUIRE direct sfcWind monthly average; do not fallback
try
    ts_ws = fun_read_cmip6_timeseries(str,'sfcWind', str(3).nc);
    giwspd_wsma = mean(ts_ws,3,'omitnan');
catch
    error('CMIP6 requires sfcWind file for wind speed (wsma). Provide file: %s', ['sfcWind_' str(3).nc '.nc']);
end

% Build a p-grid-like ocean mask (1 ocean, NaN land) on input grid
gimaskp = gimask;
gimaskp(gimaskp~=1.0) = NaN;

% Create GENIE u and v edge axes
golatue = golate;
golatve = [golatm 90.0];
% longitude edges for u/v positions
golonue = [golonm golonm(end)+360.0/imax];
golonve = golone;
% GENIE NaN mask for plotting
gm = gomask; gm(gm==0) = NaN;

% --- Regrid wind stress to GOLDSTEIN u/v grids ---
if (optplots), plot_2dgridded(flipud(giwtauu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_u.IN'],['wind stress in -- u']); end
if (optplots), plot_2dgridded(flipud(giwtauv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_v.IN'],['wind stress in -- v']); end

% u: taux(y) at u grid (do not pre-mask by input land-sea mask)
[gowtauuu,~] = make_regrid_2d(gilone,gilate,(giwtauu)',golonue,golatue,false);
gowtauuu(isnan(gowtauuu)) = 0.0; gowtauuu = gowtauuu'; gowtauuu = flipud(gowtauuu); gowtauuu = gomask.*gowtauuu;
[gowtauvu,~] = make_regrid_2d(gilone,gilate,(giwtauv)',golonue,golatue,false);
gowtauvu(isnan(gowtauvu)) = 0.0; gowtauvu = gowtauvu'; gowtauvu = flipud(gowtauvu); gowtauvu = gomask.*gowtauvu;
if (optplots), plot_2dgridded(flipud(gm.*gowtauuu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_xATu.out'],['wind stress out -- x @ u']); end
if (optplots), plot_2dgridded(flipud(gm.*gowtauvu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_yATu.out'],['wind stress out -- y @ u']); end

% v: taux(y) at v grid
[gowtauuv,~] = make_regrid_2d(gilone,gilate,(giwtauu)',golonve,golatve,false);
gowtauuv(isnan(gowtauuv)) = 0.0; gowtauuv = gowtauuv'; gowtauuv = flipud(gowtauuv); gowtauuv = gomask.*gowtauuv;
[gowtauvv,~] = make_regrid_2d(gilone,gilate,(giwtauv)',golonve,golatve,false);
gowtauvv(isnan(gowtauvv)) = 0.0; gowtauvv = gowtauvv'; gowtauvv = flipud(gowtauvv); gowtauvv = gomask.*gowtauvv;
if (optplots), plot_2dgridded(flipud(gm.*gowtauuv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_xATv.out'],['wind stress out -- x @ v']); end
if (optplots), plot_2dgridded(flipud(gm.*gowtauvv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_yATv.out'],['wind stress out -- y @ v']); end

% --- Regrid wind velocity (uas/vas) to GOLDSTEIN c-grid ---
if (optplots), plot_2dgridded(flipud(giwvelu),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_x.IN'],['wind velocity in -- x']); end
[gowvelu,~] = make_regrid_2d(gilone,gilate,giwvelu',golone,golate,false);
gowvelu(isnan(gowvelu)) = 0.0; gowvelu = gowvelu'; gowvelu = flipud(gowvelu);
if (optplots), plot_2dgridded(flipud(gowvelu),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_x.OUT'],['wind velocity out -- x']); end

if (optplots), plot_2dgridded(flipud(giwvelv),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_y.IN'],['wind velocity in -- y']); end
[gowvelv,~] = make_regrid_2d(gilone,gilate,giwvelv',golone,golate,false);
gowvelv(isnan(gowvelv)) = 0.0; gowvelv = gowvelv'; gowvelv = flipud(gowvelv);
if (optplots), plot_2dgridded(flipud(gowvelv),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_y.OUT'],['wind velocity out -- y']); end

% --- Regrid wind speed (scalar) to GOLDSTEIN c-grid; save masked/unmasked variants used elsewhere ---
if (optplots), plot_2dgridded(flipud(giwspd_uvaa),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd_uvaa.IN'],['wind speed in']); end
[gowspdall,~] = make_regrid_2d(gilone,gilate,giwspd_uvaa',golone,golate,false);
gowspdall = gowspdall'; gowspdall = flipud(gowspdall);
wspeed_uvaa = gomask.*gowspdall; wspeed_uvaa(isnan(wspeed_uvaa)) = 0.0;

if (optplots), plot_2dgridded(flipud(giwspd_uvma),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd_uvma.IN'],['wind speed in']); end
[gowspdall,~] = make_regrid_2d(gilone,gilate,giwspd_uvma',golone,golate,false);
gowspdall = gowspdall'; gowspdall = flipud(gowspdall);
wspeed_uvma = gomask.*gowspdall; wspeed_uvma(isnan(wspeed_uvma)) = 0.0;

if (optplots), plot_2dgridded(flipud(giwspd_wsma),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd_wsma.IN'],['wind speed in']); end
[gowspdall,~] = make_regrid_2d(gilone,gilate,giwspd_wsma',golone,golate,false);
gowspdall = gowspdall'; gowspdall = flipud(gowspdall);
wspeed_wsma = gomask.*gowspdall; wspeed_wsma(isnan(wspeed_wsma)) = 0.0;

% Land (ENTS) regridding like HadCM3: land mask and land-only fields
glmask = abs(gomask - 1);
gimaskpl = ones(size(gimask));
gimaskpl(gimask==1) = NaN; % land=1, ocean=NaN on input grid

% uvaa land (do not pre-mask input; mask after regrid)
[gowspdl,~] = make_regrid_2d(gilone,gilate,(giwspd_uvaa)',golone,golate,false);
gowspdl = gowspdl'; gowspdl = flipud(gowspdl);
wspeed_uvaal = glmask.*gowspdl; wspeed_uvaal(isnan(wspeed_uvaal)) = 0.0;

% uvma land (do not pre-mask input; mask after regrid)
[gowspdl_uvma,~] = make_regrid_2d(gilone,gilate,(giwspd_uvma)',golone,golate,false);
gowspdl_uvma = gowspdl_uvma'; gowspdl_uvma = flipud(gowspdl_uvma);
wspeed_uvmal = glmask.*gowspdl_uvma; wspeed_uvmal(isnan(wspeed_uvmal)) = 0.0;

% wsma land (from sfcWind) (do not pre-mask input; mask after regrid)
[gowspdl_wsma,~] = make_regrid_2d(gilone,gilate,(giwspd_wsma)',golone,golate,false);
gowspdl_wsma = gowspdl_wsma'; gowspdl_wsma = flipud(gowspdl_wsma);
wspeed_wsmal = glmask.*gowspdl_wsma; wspeed_wsmal(isnan(wspeed_wsmal)) = 0.0;

% --- Copy to output arrays and save (match CESM file conventions) ---
wstress(:,:,1) = flipud(gowtauuu); % g_taux_u
wstress(:,:,2) = flipud(gowtauuv); % g_taux_v
wstress(:,:,3) = flipud(gowtauvu); % g_tauy_u
wstress(:,:,4) = flipud(gowtauvv); % g_tauy_v
wvelocity(:,:,1) = flipud(gowvelu);
wvelocity(:,:,2) = flipud(gowvelv);

% Save to ASCII like CESM reader with messages
outname = [str(2).dir '/' str(2).exp '.taux_u.dat'];
c = wstress(:,:,1); b = permute(c,[2 1]); a = reshape(b,[imax*jmax 1]);
save(outname,'a','-ascii'); fprintf('       - Written tau u (u point) data to %s\n',outname);

outname = [str(2).dir '/' str(2).exp '.taux_v.dat'];
c = wstress(:,:,2); b = permute(c,[2 1]); a = reshape(b,[imax*jmax 1]);
save(outname,'a','-ascii'); fprintf('       - Written tau u (v point) data to %s\n',outname);

outname = [str(2).dir '/' str(2).exp '.tauy_u.dat'];
c = wstress(:,:,3); b = permute(c,[2 1]); a = reshape(b,[imax*jmax 1]);
save(outname,'a','-ascii'); fprintf('       - Written tau v (u point) data to %s\n',outname);

outname = [str(2).dir '/' str(2).exp '.tauy_v.dat'];
c = wstress(:,:,4); b = permute(c,[2 1]); a = reshape(b,[imax*jmax 1]);
save(outname,'a','-ascii'); fprintf('       - Written tau v (v point) data to %s\n',outname);

outname = [str(2).dir '/' str(2).exp '.wvelx.dat'];
c = wvelocity(:,:,1); b = permute(c,[2 1]); a = reshape(b,[imax*jmax 1]);
save(outname,'a','-ascii'); fprintf('       - Written u wind speed data to %s\n',outname);

outname = [str(2).dir '/' str(2).exp '.wvely.dat'];
c = wvelocity(:,:,2); b = permute(c,[2 1]); a = reshape(b,[imax*jmax 1]);
save(outname,'a','-ascii'); fprintf('       - Written v wind speed data to %s\n',outname);

outname = [str(2).dir '/' str(2).exp '.windspeed_uvaa.dat'];
a = wspeed_uvaa; save(outname,'a','-ascii');

outname = [str(2).dir '/' str(2).exp '.windspeed_uvma.dat'];
a = wspeed_uvma; save(outname,'a','-ascii');

outname = [str(2).dir '/' str(2).exp '.windspeed_wsma.dat'];
a = wspeed_wsma; save(outname,'a','-ascii'); fprintf('       - Written BIOGEM windspeed data to %s\n',outname);

% FULL grid (unmasked) windspeed variants
[gowspdall_uvaa,~] = make_regrid_2d(gilone,gilate,giwspd_uvaa',golone,golate,false); gowspdall_uvaa = gowspdall_uvaa'; gowspdall_uvaa = flipud(gowspdall_uvaa);
outname = [str(2).dir '/' str(2).exp '.windspeed_uvaa_all.dat']; a = gowspdall_uvaa; save(outname,'a','-ascii');
[gowspdall_uvma,~] = make_regrid_2d(gilone,gilate,giwspd_uvma',golone,golate,false); gowspdall_uvma = gowspdall_uvma'; gowspdall_uvma = flipud(gowspdall_uvma);
outname = [str(2).dir '/' str(2).exp '.windspeed_uvma_all.dat']; a = gowspdall_uvma; save(outname,'a','-ascii');
[gowspdall_wsma,~] = make_regrid_2d(gilone,gilate,giwspd_wsma',golone,golate,false); gowspdall_wsma = gowspdall_wsma'; gowspdall_wsma = flipud(gowspdall_wsma);
outname = [str(2).dir '/' str(2).exp '.windspeed_wsma_all.dat']; a = gowspdall_wsma; save(outname,'a','-ascii');

% ENTS combined ocean+land fields
outname = [str(2).dir '/' str(2).exp '.windspeed_uvaa_ents.dat']; a = wspeed_uvaa + wspeed_uvaal; save(outname,'a','-ascii');
outname = [str(2).dir '/' str(2).exp '.windspeed_uvma_ents.dat']; a = wspeed_uvma + wspeed_uvmal; save(outname,'a','-ascii');
outname = [str(2).dir '/' str(2).exp '.windspeed_wsma_ents.dat']; a = wspeed_wsma + wspeed_wsmal; save(outname,'a','-ascii');
