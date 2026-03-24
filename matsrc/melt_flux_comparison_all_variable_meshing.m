% MELT_FLUX_COMPARISON  Animate melt balance at a mid-ocean ridge through time.
%
% Ridge axis at x = 0 (left boundary). Melt rises upward and exits through
% the top within a "channel zone" bounded on the right by the streamline
% with the greatest starting x that still reaches a minimum arc length.
%
% channel ZONE DETECTION (replaces uy_thresh approach)
%   For each top-boundary cell (scanned left to right), a backward
%   streamline is traced from (cx, y_max).  The rightmost starting x
%   whose streamline arc length >= sl_min_length defines x_channel.
%   This avoids the problem that uy_thresh often fails to anchor
%   streamlines that do not reach the base of the melting zone.
%
% FIGURE LAYOUT
%   subplot(211) — growing time series (updated each step):
%       blue: melt flux leaving top of channel zone  [m^2/yr]
%       red : integrated d(phi)/dt below channel zone [m^2/yr]
%
%   subplot(212) — cross-section for current step:
%       patch: d(phi)/dt  (blue = freezing, red = melting)
%       black line: backward streamline bounding the right edge of channel zone
%
% Press any key in the command window to advance to the next timestep.
%
% FIELD NAME DEFAULTS  (edit if your output uses different names):
%   velocity     — melt velocity field   [Npts x 3]
%   porosity     — melt fraction         [Npts x 1]
%   melting_rate — d(phi)/dt             [Npts x 1]
%
% REQUIRES: read_aspect_paraview_output.m on the MATLAB path.
clear;

%% ---- User settings --------------------------------------------------------
% NOTE: ASPECT outputs velocity in m/yr and coordinates in m.
%       Flux units are m^2/yr = km^2/Myr (exact equivalence, no rescaling needed).
%       dporosity_dt is in yr^-1 (divided by timestep in years in postprocess).
vel_field        = 'u_f';            % melt velocity field name (m/yr)
phi_field        = 'porosity';       % porosity (melt fraction) field name
dmdt_field       = 'dporosity_dt';   % d(phi)/dt field name (yr^-1)
bwd_speed_thresh = 0.1;             % m/yr: min speed to continue backward streamline
sl_ds_frac       = 1/500;            % streamline step size as fraction of domain height
sl_min_length_km = 50;              % minimum streamline arc length [km]
                                     % for a starting x to count as inside the channel zone
dmdt_thresh      = 1e-9;            % yr^-1: only cells with dmdt_c > this are included
                                     % in area_left and gen_below (defines bottom of melting zone)
seconds_per_year = 60*60*24*365.25;

%% ---- Load data ------------------------------------------------------------
pvd_file = 'solution.pvd'
%pvd_file = strtrim(input('Enter path to solution.pvd [solution.pvd]: ', 's'));
%if isempty(pvd_file), pvd_file = 'solution.pvd'; end

data   = read_aspect_paraview_output(pvd_file);
nsteps = numel(data.times);
vm     = iscell(data.x);       % true when AMR changed mesh between timesteps
t_yr   = data.times;           % ASPECT writes time in years

% List all available fields to help identify correct field names
fprintf('Available fields in data:\n');
fns_all = fieldnames(data);
fprintf('  %s\n', fns_all{:});

for fn = {vel_field, phi_field, dmdt_field}
    if ~isfield(data, fn{1})
        error('Required field "%s" not found. See available fields listed above.', fn{1});
    end
end

%% ---- Output arrays --------------------------------------------------------
flux_top  = nan(nsteps, 1);   % (1) m^2/yr  melt flux through top
gen_below = nan(nsteps, 1);   % (2) m^2/yr  integrated d(phi)/dt in source region


%% ---- Initialise figure ----------------------------------------------------
hfig = figure(3);  clf(hfig);

%% ---- Main loop ------------------------------------------------------------
for s = 1:nsteps
    fprintf('Step %d/%d  (t = %.4g yr)\n', s, nsteps, t_yr(s));

    % --- Extract fields for this timestep ---
    if vm
        x    = data.x{s};          y    = data.y{s};
        conn = data.connectivity{s};
        vel  = data.(vel_field){s};
        phi  = data.(phi_field){s};
        dmdt = data.(dmdt_field){s};
    else
        x    = data.x;              y    = data.y;
        conn = data.connectivity;
        vel  = slice_field(data.(vel_field),  s);
        phi  = slice_field(data.(phi_field),  s);
        dmdt = slice_field(data.(dmdt_field), s);
    end
    dmdt = dmdt .* seconds_per_year;
    

    u_x = vel(:,1);   u_y = vel(:,2);

    % --- Cell geometry ---
    cx    = mean(x(conn), 2);
    cy    = mean(y(conn), 2);
    carea = quad_cell_areas(x, y, conn);

    phi_c  = mean(phi(conn),  2);
    uy_c   = mean(u_y(conn),  2);
    ux_c   = mean(u_x(conn),  2);
    dmdt_c = mean(dmdt(conn), 2);

    % ---- Cell geometry for top boundary ------------------------------------
    y_max  = max(y);   y_min = min(y);
    % y_sorted = unique(y);
    % tol    = (y_max - y_sorted(end-1)) * 0.1;   % 10% of top cell height — adapts to local AMR resolution
    % is_top = max(y(conn), [], 2) >= y_max - tol;
    is_top = max(y(conn), [], 2) == y_max;
    % figure(99); cla; hold on;
    % for k = find(is_top)'
    %     ni = conn(k, [1 2 3 4 1]);
    %     plot(x(ni)/1e3, y(ni)/1e3, 'r-');
    % end
    % title('Top boundary cells'); xlabel('x (km)'); ylabel('y (km)');

    tol = 0;   % is_top uses exact equality so tol=0 is consistent
    fw  = top_face_widths(x, y, conn, is_top, y_max, tol);

    domain_height  = y_max - y_min;
    ds             = domain_height * sl_ds_frac;
    sl_min_length  = sl_min_length_km * 1e3;   % convert km to m

    % ---- (1) Find x_channel by scanning top cells left to right -------------
    % For each top cell, trace a backward streamline and measure its arc
    % length.  x_channel is the largest cx where arc length >= sl_min_length.
    top_idx  = find(is_top);
    top_cx   = cx(top_idx);
    [top_cx_sorted, sort_ord] = sort(top_cx);   % left to right
    top_idx_sorted = top_idx(sort_ord);

    x_channel = NaN;
    sl       = zeros(0, 2);   % streamline for final plot

    for ki = 1:numel(top_idx_sorted)
        xi = top_cx_sorted(ki);
        sl_i = backward_streamline([xi, y_max], x, y, u_x, u_y, ...
                                    bwd_speed_thresh, ds, y_min);
        arc_i = streamline_length(sl_i);
        if arc_i >= sl_min_length
            x_channel = xi;
            sl       = sl_i;   % keep this streamline for plotting
        end
    end

    fprintf('  x_channel = %.1f km  |  sl_min_length = %.1f km\n', ...
            x_channel/1e3, sl_min_length/1e3);

    % ---- (2) Melt flux through top of channel zone --------------------------
    % channel zone: top cells with centroid x <= x_channel (and positive uy).
    gen_below(s) = 0;
    if isfinite(x_channel)
        channel      = is_top & cx <= x_channel & uy_c > 0;
        flux_top(s) = sum(phi_c(channel) .* uy_c(channel) .* fw(channel));

        % ---- (3) Integrated d(phi)/dt left of bounding streamline ----------
        [gen_below(s), area_sl] = integrate_melt_generation(x, y, conn, carea, dmdt_c, sl, dmdt_thresh);
    else
        flux_top(s)  = 0;
        area_sl      = 0;
        fprintf('  WARNING: no streamline reached sl_min_length — channel zone undefined.\n');
    end

    mean_dmdt = gen_below(s) / area_sl;
    fprintf('  Fluxes: Melt Rise = %.4g  |  Melt Gen = %.4g  km^2/Myr\n', flux_top(s), gen_below(s));
    fprintf('  Area of mantle supplying melt = %.4g km^2  |  Mean dporosity/dt = %.4g yr^-1\n', ...
            area_sl/1e6, mean_dmdt);

    % ======================================================================
    % SUBPLOT (211): growing time series
    % ======================================================================
    ax1 = subplot(2,1,1);
    cla(ax1);  hold(ax1, 'on');

    idx = 1:s;
    plot(ax1, t_yr(idx), flux_top(idx),  'b-o', ...
         'MarkerSize', 3, 'LineWidth', 0.5, 'DisplayName', 'Top flux');
    plot(ax1, t_yr(idx), gen_below(idx), 'r-o', ...
         'MarkerSize', 3, 'LineWidth', 0.5, 'DisplayName', 'Gen. below streamline');

    % Highlight current step
    plot(ax1, t_yr(s), flux_top(s),  'bs', 'MarkerSize', 7, ...
         'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
    plot(ax1, t_yr(s), gen_below(s), 'rs', 'MarkerSize', 7, ...
         'MarkerFaceColor', 'r', 'HandleVisibility', 'off');

    xlim(ax1, [t_yr(1), t_yr(end)]);
    xlabel(ax1, 'Time (yr)');
    ylabel(ax1, 'Flux (km^2 Myr^{-1})');
    legend(ax1, 'Location', 'best');
    title(ax1, 'Melt balance: channel zone flux (blue) vs source generation (red)');
    grid(ax1, 'on');

    % ======================================================================
    % SUBPLOT (212): cross-section of d(phi)/dt + streamline
    % ======================================================================
    nc = 128;
    cmap = [linspace(0,1,nc/2)', linspace(0,1,nc/2)', ones(nc/2,1); ...
             ones(nc/2,1), linspace(1,0,nc/2)', linspace(1,0,nc/2)'];
    cmap=turbo;
    climits = [0 1e-14] * seconds_per_year;

    ax2 = subplot(212);
    cla(ax2);  hold(ax2, 'on');

    x_km = x / 1e3;   y_km = y / 1e3;

    valid = all(isfinite(dmdt(conn)), 2);
    patch(ax2, 'Faces', conn(valid,:), 'Vertices', [x_km, y_km], ...
          'FaceVertexCData', dmdt, 'FaceColor', 'interp', 'EdgeColor', 'none');

    colormap(ax2, cmap);
    caxis(ax2, climits);

    cb = colorbar(ax2, 'Location', 'eastoutside');
    cb.Label.String = 'd\phi/dt  (yr^{-1})';

    if size(sl, 1) >= 2
        plot(ax2, sl(:,1)/1e3, sl(:,2)/1e3, 'w-',  'LineWidth', 0.5);
        plot(ax2, sl(1,1)/1e3, sl(1,2)/1e3, 'wv',  ...
             'MarkerSize', 6, 'MarkerFaceColor', 'w');   % start marker at top
    end

    axis(ax2, 'equal', 'tight');
    xlabel(ax2, 'Distance (km)');
    ylabel(ax2, 'Depth (km)');
    if isfinite(x_channel)
        title(ax2, sprintf('d\\phi/dt  at  t = %.4g yr  |  x_{channel} = %.1f km  (sl_{min} = %.0f km)', ...
              t_yr(s), x_channel/1e3, sl_min_length_km));
    else
        title(ax2, sprintf('d\\phi/dt  at  t = %.4g yr  |  channel zone undefined', t_yr(s)));
    end

    drawnow;

    if s < nsteps
        fprintf('  Press any key for next step (%d/%d)...\n', s+1, nsteps);
        pause;
    end
end

fprintf('\nDone. All %d steps processed.\n', nsteps);
if (0)
    % --- Solid velocity x-component at y=yprof ---
    yprof=50e3; 
    if vm
        vel_s = data.velocity{s};
    else
        vel_s = slice_field(data.velocity, s);
    end
    vsx      = vel_s(:,1);          % x-component at all nodes
    top_node = y >= y_max;          % logical mask for top nodes
    xn_top   = x(top_node);
    vsx_top  = vsx(top_node);
    [xn_top, ord] = sort(xn_top);
    vsx_top  = vsx_top(ord);
    
    figure(4); clf;
    plot(xn_top/1e3, vsx_top, 'b.-');
    xlabel('x (km)');  ylabel('solid v_x  (m/yr)');
    title(sprintf('Solid horizontal velocity at top boundary  (t = %.4g yr)', t_yr(s)));
    grid on;
end;
%% ==========================================================================
function v = slice_field(arr, s)
% Extract timestep s from a stacked array (no-op for single-step data).
    nd = ndims(arr);
    if nd == 2 && size(arr,2) > 1
        v = arr(:, s);        % scalar [Npts x Nsteps]
    elseif nd == 3
        v = arr(:, :, s);     % vector [Npts x Ncomp x Nsteps]
    else
        v = arr;              % single timestep — no trailing dimension
    end
end

%% ==========================================================================
function A = quad_cell_areas(x, y, conn)
% Shoelace area of each quad cell (vectorised).
    x1=x(conn(:,1)); x2=x(conn(:,2)); x3=x(conn(:,3)); x4=x(conn(:,4));
    y1=y(conn(:,1)); y2=y(conn(:,2)); y3=y(conn(:,3)); y4=y(conn(:,4));
    A = 0.5 * abs((x1.*y2-x2.*y1) + (x2.*y3-x3.*y2) + ...
                  (x3.*y4-x4.*y3) + (x4.*y1-x1.*y4));
end

%% ==========================================================================
function fw = top_face_widths(x, y, conn, is_top, y_max, tol)
% Width of each top-boundary cell's upper face.
    top_node = y >= y_max - tol;
    n  = size(conn, 1);
    fw = zeros(n, 1);
    for k = find(is_top)'
        ni = conn(k,:);
        tn = top_node(ni);
        xi = x(ni);
        if sum(tn) >= 2
            fw(k) = max(xi(tn)) - min(xi(tn));
        else
            fw(k) = max(xi) - min(xi);
        end
    end
end

%% ==========================================================================
function L = streamline_length(sl)
% Total arc length of a streamline (sum of segment lengths).
    if size(sl, 1) < 2
        L = 0;
        return;
    end
    diffs = diff(sl, 1, 1);
    L = sum(hypot(diffs(:,1), diffs(:,2)));
end

%% ==========================================================================
function sl = backward_streamline(pt0, x, y, u_x, u_y, thresh, ds, y_min)
% Trace melt streamline BACKWARD (against velocity) from pt0 downward
% using a 4th-order Runge-Kutta (RK4) arc-length integrator.
%
% The ODE is:  d(pt)/ds = -vel(pt) / |vel(pt)|   (unit backward direction)
% RK4 evaluates vel at four intermediate points per step for high accuracy.
% Velocity is interpolated by nearest node.
    max_pts = 10000;
    sl = nan(max_pts, 2);
    sl(1,:) = pt0;

    for k = 1:max_pts-1
        pt = sl(k,:);
        if pt(2) < y_min + ds,  sl = sl(1:k,:);  break;  end

        % RK4 stages: each returns the unit backward-direction vector
        d1 = bwd_dir(pt,            x, y, u_x, u_y, thresh);
        if isempty(d1),  sl = sl(1:k,:);  break;  end

        d2 = bwd_dir(pt + ds/2*d1, x, y, u_x, u_y, thresh);
        if isempty(d2),  d2 = d1;  end

        d3 = bwd_dir(pt + ds/2*d2, x, y, u_x, u_y, thresh);
        if isempty(d3),  d3 = d2;  end

        d4 = bwd_dir(pt + ds*d3,   x, y, u_x, u_y, thresh);
        if isempty(d4),  d4 = d3;  end

        sl(k+1,:) = pt + ds/6 * (d1 + 2*d2 + 2*d3 + d4);

        if k == max_pts-1,  sl = sl(1:k,:);  end
    end
    sl = sl(all(isfinite(sl),2), :);
end
% ==========================================================================
function d = bwd_dir(pt, x, y, u_x, u_y, thresh)
% Unit vector opposite to melt velocity at pt (nearest-node interpolation).
% Returns [] if speed is below threshold (flow effectively zero).
    [~, nn] = min((x - pt(1)).^2 + (y - pt(2)).^2);
    ux = u_x(nn);   uy = u_y(nn);
    spd = hypot(ux, uy);
    if spd < thresh
        d = [];
    else
        d = -[ux, uy] / spd;
    end
end

%% ==========================================================================
function [gen, area_left] = integrate_melt_generation(x, y, conn, carea, dmdt_c, sl, dmdt_thresh)
% Integrate dmdt_c over all area to the left of the streamline where
% dmdt_c > dmdt_thresh. Returns gen (integral of dmdt_c weighted by area) and area_left (total
% area of contributing cells to the left of the streamline, in m^2).
%
% Classification per cell:
%   Cells with dmdt_c <= dmdt_thresh are excluded entirely.
%   All 4 nodes left  → full cell area contributes.
%   All 4 nodes right → zero contribution.
%   Mixed (1-3 left)  → partial area computed by Sutherland-Hodgman
%                        polygon clipping of the quad against the
%                        boundary x = x_sl(y), approximated linearly
%                        between each pair of adjacent vertices.
%
% This correctly handles variably-sized cells and the fact that the
% streamline cuts through FEM cells rather than following cell edges.
    gen = 0;  area_left = 0;
    if size(sl,1) < 2,  return;  end

    % Build sorted, deduplicated streamline interpolant
    [sl_y, ord] = sort(sl(:,2));  sl_x = sl(ord,1);
    [sl_y, ui]  = unique(sl_y);   sl_x = sl_x(ui);

    if numel(sl_y) < 2
        % Degenerate: streamline collapsed to a point — use single x threshold
        left      = mean(x(conn), 2) < sl_x(1) & dmdt_c > dmdt_thresh;
        gen       = sum(dmdt_c(left) .* carea(left));
        area_left = sum(carea(left));
        return;
    end

    % Streamline x at every mesh node (y clamped to streamline range)
    y_cl       = min(max(y, sl_y(1)), sl_y(end)); %nodes below sl_y(1) are set to sl_y(1);  
    x_sl_node  = interp1(sl_y, sl_x, y_cl, 'linear');  % [Nnodes x 1]

    % Per-node and per-cell left-of-streamline status
    node_left  = x < x_sl_node;           % [Nnodes x 1] logical
    n_left     = sum(node_left(conn), 2);  % [Ncells x 1] integer 0..4

    % Fully-left cells: whole cell area (only where dmdt exceeds threshold)
    full_left  = n_left == 4 & dmdt_c > dmdt_thresh;
    gen        = sum(dmdt_c(full_left) .* carea(full_left));
    area_left  = sum(carea(full_left));

    % Cut cells: partial area via polygon clipping (only where dmdt exceeds threshold)
    x_sl_conn  = x_sl_node(conn);         % [Ncells x 4] x_sl at each cell's nodes
    cut_idx    = find(n_left > 0 & n_left < 4 & dmdt_c > dmdt_thresh);
    for ki = 1:numel(cut_idx)
        c          = cut_idx(ki);
        ni         = conn(c,:);
        a          = clip_quad_left(x(ni), y(ni), x_sl_conn(c,:)');
        gen        = gen       + dmdt_c(c) * a;
        area_left  = area_left + a;
    end
end

%% ==========================================================================
function area = clip_quad_left(px, py, px_sl)
% Area of the portion of a quadrilateral lying to the left of the
% streamline, i.e. where x < x_sl(y).
%
% Uses the Sutherland-Hodgman algorithm with signed distance
%   d(i) = px(i) - px_sl(i)
% where px_sl(i) is the streamline x interpolated at vertex i's y.
% "Inside" means d < 0  (vertex is left of the streamline).
% Edge-boundary intersections use the standard S-H parameter
%   t = d_i / (d_i - d_j),   point = p_i + t*(p_j - p_i).
%
% The streamline within each edge is approximated as linear between the
% two endpoint values of px_sl — accurate when cells are small relative
% to the streamline curvature (typical for ASPECT FEM meshes).

    d  = px(:) - px_sl(:);   % signed distance, [4 x 1]
    vx = px(:);  vy = py(:);

    out_x = zeros(1,8);  out_y = zeros(1,8);  nout = 0;

    for i = 1:4
        j    = mod(i, 4) + 1;
        in_i = d(i) < 0;
        in_j = d(j) < 0;

        if in_i
            nout = nout + 1;
            out_x(nout) = vx(i);
            out_y(nout) = vy(i);
        end

        if in_i ~= in_j
            t    = d(i) / (d(i) - d(j));
            nout = nout + 1;
            out_x(nout) = vx(i) + t*(vx(j) - vx(i));
            out_y(nout) = vy(i) + t*(vy(j) - vy(i));
        end
    end

    if nout < 3,  area = 0;  return;  end

    ox   = out_x(1:nout);
    oy   = out_y(1:nout);
    area = 0.5 * abs(sum(ox .* oy([2:end,1]) - ox([2:end,1]) .* oy));
end
