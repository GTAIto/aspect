function [F_rg, xg, yg] = interpolate_soln_to_regular_grid(F, x, y, conn, dx, dy, display_cell_message)
% INTERPOLATE_SOLN_TO_REGULAR_GRID  Interpolate an ASPECT nodal field onto a
% uniform rectilinear grid using linear scattered interpolation.
%
% USAGE
%   [F_rg, xg, yg] = interpolate_soln_to_regular_grid(F, x, y, conn)
%   [F_rg, xg, yg] = interpolate_soln_to_regular_grid(F, x, y, conn, dx, dy)
%
% INPUTS
%   F    - nodal field: [Nnodes x 1] scalar  OR  [Nnodes x ncomp] vector/tensor
%   x, y - node coordinates [Nnodes x 1], metres
%   conn - element connectivity [Ncells x 4]
%   dx   - (optional) regular grid x-spacing [m].  Default: 2nd-smallest cell width.
%   dy   - (optional) regular grid y-spacing [m].  Default: 2nd-smallest cell height.
%   display_cell_message = true to display original cell dimensions and new
%   grid size
%
% OUTPUTS
%   F_rg - interpolated field on regular grid:
%            [ny x nx]        for a scalar field (ncomp == 1)
%            [ny x nx x ncomp] for a vector/tensor field
%          Points outside the convex hull of the mesh are NaN.
%   xg   - regular grid x-coordinates [1 x nx], metres
%   yg   - regular grid y-coordinates [ny x 1], metres
%
% The full range of original cell sizes is printed to the command window so
% the user can verify or override dx/dy.
%
% EXAMPLE
%   [phi_rg, xg, yg] = interpolate_soln_to_regular_grid(phi, x, y, conn);
%   [uf_rg,  xg, yg] = interpolate_soln_to_regular_grid(vel, x, y, conn, 500, 500);

if nargin < 7,  display_cell_message = false;  end

%% --- Compute cell sizes from connectivity ----------------------------------
x1=x(conn(:,1)); x2=x(conn(:,2)); x3=x(conn(:,3)); x4=x(conn(:,4));
y1=y(conn(:,1)); y2=y(conn(:,2)); y3=y(conn(:,3)); y4=y(conn(:,4));
cell_w = max([x1,x2,x3,x4],[],2) - min([x1,x2,x3,x4],[],2);
cell_h = max([y1,y2,y3,y4],[],2) - min([y1,y2,y3,y4],[],2);

w_sorted = sort(unique(round(cell_w)));
h_sorted = sort(unique(round(cell_h)));

%% --- Choose default dx, dy if not supplied ---------------------------------
if nargin < 5 || isempty(dx)
    dx = w_sorted(min(2, numel(w_sorted)));   % 2nd-smallest width
end
if nargin < 6 || isempty(dy)
    dy = h_sorted(min(2, numel(h_sorted)));   % 2nd-smallest height
end

%% --- Build regular grid ----------------------------------------------------
x_min = min(x);  x_max = max(x);
y_min = min(y);  y_max = max(y);

xg = x_min : dx : x_max;          % [1 x nx]
yg = (y_min : dy : y_max)';       % [ny x 1]
nx = numel(xg);
ny = numel(yg);

if display_cell_message
    fprintf('Interpolation to regular grid: Original cell dimensions: %.1f–%.1f m in x/ %.1f–%.1f m in y\n', ...
        w_sorted(1), w_sorted(end), h_sorted(1), h_sorted(end));
    fprintf('  Regular grid spacing: dx = %.1f m, dy = %.1f m  →  %d × %d (nx × ny)\n', ...
        dx, dy, nx, ny);
end
[Xq, Yq] = meshgrid(xg, yg);

%% --- Deduplicate nodes (shared nodes appear once per element in some formats)
[~, ia, ic] = unique([x(:), y(:)], 'rows', 'stable');
xu = x(ia);  yu = y(ia);

%% --- Interpolate each component --------------------------------------------
ncomp = size(F, 2);
F_rg  = nan(ny, nx, ncomp);

for k = 1:ncomp
    Fk = accumarray(ic, double(F(:,k)), [numel(ia), 1], @mean);
    I  = scatteredInterpolant(xu, yu, Fk, 'linear', 'none');
    F_rg(:,:,k) = I(Xq, Yq);
end

if ncomp == 1
    F_rg = squeeze(F_rg);   % return [ny x nx] for scalar fields
end
end
