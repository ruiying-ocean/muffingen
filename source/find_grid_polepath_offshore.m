function [n_paths,v_paths,grid_paths,info] = find_grid_polepath_offshore(n_paths,v_paths,grid_mask,grid_k1,grid_islands,i_poles,v_de,par_Dmin)
% FIND_GRID_POLEPATH_OFFSHORE Route the S polar island path through deep water.
%
% find_grid_paths follows the coast of the S polar island, which can put
% most of the path over shallow shelf. GOLDSTEIN's island integral weights
% the wind stress by 1/depth and adds bottom-pressure terms between cells
% of different depth, so a shelf path can make the polar island
% streamfunction (the ACC) ill-conditioned. This replaces the S polar path
% with an eastward zonal path along the southernmost row that is entirely
% ocean at least PAR_DMIN (m) deep, provided all land south of that row
% belongs to the S polar island. Otherwise the coastal path is kept.
% A zonal path in open water is a valid island integral contour: the ocean
% psi points between it and the coast are solved as ordinary interior points.
%
% N_PATHS, V_PATHS are as returned by find_grid_paths (V_PATHS rows are
% [direction i j] with j counted from the S, as in GENIE). GRID_MASK,
% GRID_K1 and GRID_ISLANDS use muffingen orientation (row 1 = N).
% V_DE holds the depth of the bottom of each level, as used by fun_conv_k1.

[jmax,imax] = size(grid_mask);
kmax = length(v_de)-1;
depth = abs(v_de(1:kmax));
info.moved = false;
info.j = [];
info.n_before = [];
info.n_shelf_before = [];
grid_paths = [];
% identify the S polar island: land on the S row renumbered as the last
% island (find_grid_islands_update), or a land-free S pole (i_poles == -2);
% either way its path is the last one created by find_grid_paths
s_label = unique(grid_islands(jmax,grid_mask(jmax,:) == 0));
if isscalar(s_label) && s_label == max(grid_islands(:))
elseif isempty(s_label) && any(i_poles == -2)
    s_label = NaN;
else
    disp('       ! S polar island not identified uniquely; keeping coastal path.');
    return;
end
p_end = sum(n_paths);
p_start = p_end - n_paths(end) + 1;
old = v_paths(p_start:p_end,:);
if min(old(:,3)) > 3
    disp('       ! Last path does not reach the S polar rows; keeping coastal path.');
    return;
end
info.n_before = n_paths(end);
old_k1 = grid_k1(sub2ind([jmax imax],jmax-old(:,3)+1,old(:,2)));
info.n_shelf_before = sum(depth(old_k1) < par_Dmin);
% search northward from the S pole for a deep, all-ocean row
for j = 1:floor(jmax/2)
    row = jmax - j + 1;
    south = (row+1):jmax;
    south_land = (grid_mask(south,:) == 0);
    if any(grid_islands(south,:) ~= s_label & south_land, 'all')
        disp('       ! Land south of the candidate row is not part of the S polar island; keeping coastal path.');
        return;
    end
    if all(grid_mask(row,:) == 1) && all(depth(grid_k1(row,:)) >= par_Dmin)
        new = [ones(imax,1) (1:imax)' repmat(j,imax,1)];
        v_paths = [v_paths(1:p_start-1,:); new];
        n_paths(end) = imax;
        info.moved = true;
        info.j = j;
        break;
    end
end
if ~info.moved
    disp(['       ! No all-ocean row at least ' num2str(par_Dmin) ...
        ' m deep found near the S pole; keeping coastal path.']);
    return;
end
% rebuild the 2D paths map (as in find_grid_paths)
grid_paths = zeros(jmax,imax);
for p = (n_paths(1)+1):sum(n_paths)
    grid_paths(v_paths(p,3),v_paths(p,2)) = v_paths(p,1);
end
grid_paths = flipud(grid_paths);
end
