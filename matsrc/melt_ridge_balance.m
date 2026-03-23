% MELT_RIDGE_BALANCE  Animate melt balance at a mid-ocean ridge through time.
%
% Ridge axis at x = 0 (left boundary). Melt rises upward and exits through
% the top within a "feeder zone" (cells where u_y > 0 at the surface).
%
% FIGURE LAYOUT
%   subplot(211) — growing time series (updated each step):
%       blue: melt flux leaving top of feeder zone  [m^2/s]
%       red : integrated d(phi)/dt below feeder zone [m^2/s]
%
%   subplot(212) — cross-section for current step:
%       patch: d(phi)/dt  (blue = freezing, red = melting)
%       black line: backward streamline bounding the right edge of feeder zone
%
% Press any key in the command window to advance to the next timestep.
%
% FIELD NAME DEFAULTS  (edit if your output uses different names):
%   velocity     — melt velocity field   [Npts x 3]
%   porosity     — melt fraction         [Npts x 1]
%   melting_rate — d(phi)/dt             [Npts x 1]
%
% REQUIRES: read_aspect_output.m on the MATLAB path.
clear;

%% ---- User settings --------------------------------------------------------
% NOTE: ASPECT outputs velocity in m/yr and coordinates in m.
%       Flux units in this script are therefore m^2/yr.
%       dporosity_dt is in yr^-1 (divided by timestep in years in postprocess).
vel_field   = 'fluid_velocity'; % melt velocity field name (m/yr)
phi_field   = 'porosity';       % porosity (melt fraction) field name
dmdt_field  = 'dporosity_dt';   % d(phi)/dt field name (yr^-1)
uy_thresh   = 1e-4;             % m/yr: threshold for "upward melt" at top boundary
sl_ds_frac  = 1/500;            % streamline step size as fraction of domain height

%% ---- Load data ------------------------------------------------------------
pvd_file = strtrim(input('Enter path to solution.pvd [solution.pvd]: ', 's'));
if isempty(pvd_file), pvd_file = 'solution.pvd'; end

data   = read_aspect_output(pvd_file);
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
flux_top  = nan(nsteps, 1);   % (1) m^2/s  melt flux through top
gen_below = nan(nsteps, 1);   % (3) m^2/s  integrated d(phi)/dt in source region

%% ---- Colormap: blue-white-red for d(phi)/dt -------------------------------
nc = 128;
cmap_bwr = [linspace(0,1,nc/2)', linspace(0,1,nc/2)', ones(nc/2,1); ...
             ones(nc/2,1), linspace(1,0,nc/2)', linspace(1,0,nc/2)'];

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

    u_x = vel(:,1);   u_y = vel(:,2);

    % --- Cell geometry ---
    cx    = mean(x(conn), 2);
    cy    = mean(y(conn), 2);
    carea = quad_cell_areas(x, y, conn);

    phi_c  = mean(phi(conn),  2);
    uy_c   = mean(u_y(conn),  2);
    ux_c   = mean(u_x(conn),  2);
    dmdt_c = mean(dmdt(conn), 2);

    % ---- (1) Melt flux through top boundary --------------------------------
    y_max  = max(y);   y_min = min(y);
    tol    = (y_max - y_min) * 1e-4;
    is_top = max(y(conn), [], 2) >= y_max - tol;
    fw     = top_face_widths(x, y, conn, is_top, y_max, tol);

    % Diagnostic: show top-boundary uy range to help tune uy_thresh
    uy_top = uy_c(is_top);
    fprintf('  Top-boundary u_y: min=%.3g  max=%.3g  mean=%.3g  (uy_thresh=%.3g) m/yr\n', ...
            min(uy_top), max(uy_top), mean(uy_top), uy_thresh);

    upward      = is_top & uy_c > uy_thresh;
    flux_top(s) = sum(phi_c(upward) .* uy_c(upward) .* fw(upward));

    % ---- (2) Backward streamline from rightmost feeder cell ----------------
    sl         = zeros(0, 2);
    x_feeder   = NaN;
    gen_below(s) = 0;

    if any(upward)
        x_feeder = max(cx(upward));
        ds       = (y_max - y_min) * sl_ds_frac;
        sl       = backward_streamline([x_feeder, y_max], x, y, u_x, u_y, ...
                                        uy_thresh, ds, y_min);

        % ---- (3) Integrated d(phi)/dt left of streamline ------------------
        left         = cells_left_of_path(cx, cy, sl);
        mask         = left & dmdt_c > 0;
        gen_below(s) = sum(dmdt_c(mask) .* carea(mask));
    end

    fprintf('  Flux = %.4g  |  Gen = %.4g  m^2/yr\n', flux_top(s), gen_below(s));

    % ======================================================================
    % SUBPLOT (211): growing time series
    % ======================================================================
    ax1 = subplot(2,1,1);
    cla(ax1);  hold(ax1, 'on');

    % Plot all computed steps up to and including s
    idx = 1:s;
    hl1 = plot(ax1, t_yr(idx), flux_top(idx),  'b-o', ...
               'MarkerSize', 3, 'LineWidth', 1.5, 'DisplayName', 'Top flux');
    hl2 = plot(ax1, t_yr(idx), gen_below(idx), 'r-o', ...
               'MarkerSize', 3, 'LineWidth', 1.5, 'DisplayName', 'Gen. below streamline');

    % Highlight the current step with filled markers
    plot(ax1, t_yr(s), flux_top(s),  'bs', 'MarkerSize', 9, ...
         'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
    plot(ax1, t_yr(s), gen_below(s), 'rs', 'MarkerSize', 9, ...
         'MarkerFaceColor', 'r', 'HandleVisibility', 'off');

    % Fix x-axis to full time range so the plot doesn't rescale
    xlim(ax1, [t_yr(1), t_yr(end)]);

    xlabel(ax1, 'Time (yr)');
    ylabel(ax1, 'Flux (m^2 yr^{-1})');
    legend(ax1, 'Location', 'best');
    title(ax1, 'Melt balance: feeder zone flux (blue) vs source generation (red)');
    grid(ax1, 'on');

    % ======================================================================
    % SUBPLOT (212): cross-section of d(phi)/dt + streamline
    % ======================================================================
    ax2 = subplot(2,1,2);
    cla(ax2);  hold(ax2, 'on');

    x_km = x / 1e3;   y_km = y / 1e3;

    % Plot d(phi)/dt field (per-node, interpolated across cells)
    valid = all(isfinite(dmdt(conn)), 2);
    patch(ax2, 'Faces', conn(valid,:), 'Vertices', [x_km, y_km], ...
          'FaceVertexCData', dmdt, 'FaceColor', 'interp', 'EdgeColor', 'none');

    colormap(ax2, cmap_bwr);

    % Symmetric colour limits centred on zero
    clim_val = max(abs(dmdt));
    if clim_val == 0 || ~isfinite(clim_val),  clim_val = 1;  end
    caxis(ax2, [-clim_val, clim_val]);

    cb = colorbar(ax2, 'Location', 'eastoutside');
    cb.Label.String = 'd\phi/dt  (yr^{-1})';

    % Overlay the feeder-zone boundary streamline
    if size(sl, 1) >= 2
        plot(ax2, sl(:,1)/1e3, sl(:,2)/1e3, 'k-',  'LineWidth', 2);
        plot(ax2, sl(1,1)/1e3, sl(1,2)/1e3, 'kv',  ...
             'MarkerSize', 8, 'MarkerFaceColor', 'k');   % start marker at top
    end

    axis(ax2, 'equal', 'tight');
    xlabel(ax2, 'Distance (km)');
    ylabel(ax2, 'Depth (km)');
    if isfinite(x_feeder)
        title(ax2, sprintf('d\\phi/dt  at  t = %.4g yr  |  x_{feeder} = %.1f km', ...
              t_yr(s), x_feeder/1e3));
    else
        title(ax2, sprintf('d\\phi/dt  at  t = %.4g yr  |  no upward melt', t_yr(s)));
    end

    drawnow;

    if s < nsteps
        fprintf('  Press any key for next step (%d/%d)...\n', s+1, nsteps);
        pause;
    end
end

fprintf('\nDone. All %d steps processed.\n', nsteps);

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
function left = cells_left_of_path(cx, cy, sl)
% True for cells whose centroid lies LEFT (smaller x) of the streamline.
% Streamline x is linearly interpolated at each centroid's y.
    left = false(numel(cx), 1);
    if size(sl,1) < 2,  return;  end
    [sl_y, ord] = sort(sl(:,2));
    sl_x = sl(ord, 1);
    [sl_y, ui] = unique(sl_y);
    sl_x = sl_x(ui);
    % After deduplication the streamline may collapse to a single point
    % (e.g. velocity was zero immediately below the start).
    % Fall back to a simple x-threshold comparison in that case.
    if numel(sl_y) < 2
        left = cx < sl_x(1);
        return;
    end
    x_sl = interp1(sl_y, sl_x, min(max(cy, sl_y(1)), sl_y(end)), 'linear');
    left = cx < x_sl;
end
