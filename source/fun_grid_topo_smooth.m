function [grid_topo] = fun_grid_topo_smooth(grid_topo,grid_mask,n_pass)
% FUN_GRID_TOPO_SMOOTH Smooth ocean bathymetry before it is cut into levels.
%
% Each pass replaces every ocean cell by the mean of itself and its N, S, E
% and W ocean neighbours (periodic in longitude), as the default bathymetry
% smoothing of CLIMBER-X (Willeit et al., 2022, GMD 15, 5905), which removes
% one-cell cliffs that destabilise the barotropic solve. Unlike CLIMBER-X,
% land cells are left out of the average, so coasts and narrow gateways are
% not made shallower by the adjacent land. The polar rows are not smoothed.
%
% GRID_TOPO uses muffingen orientation (row 1 = N); GRID_MASK is 1 for ocean.

[jmax,imax] = size(grid_topo);
ocean = (grid_mask == 1);
for pass = 1:n_pass
    old = grid_topo;
    for j = 2:jmax-1
        for i = 1:imax
            if ~ocean(j,i), continue; end
            im = mod(i-2,imax)+1;
            ip = mod(i,imax)+1;
            stencil = [j i; j im; j ip; j-1 i; j+1 i];
            total = 0.0;
            n = 0;
            for s = 1:5
                if ocean(stencil(s,1),stencil(s,2))
                    total = total + old(stencil(s,1),stencil(s,2));
                    n = n + 1;
                end
            end
            grid_topo(j,i) = total/n;
        end
    end
end
end
