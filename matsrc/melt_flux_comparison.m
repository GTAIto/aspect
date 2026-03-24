% MELT_FLUX_COMPARISON  Animate melt balance at a mid-ocean ridge through time.
%
% Ridge axis at x = 0 (left boundary). Melt rises upward and exits through
% the top within a "channel zone" bounded on the right by the streamline
% with the greatest starting x that still reaches a minimum arc length.
%
% All fields are first interpolated onto a uniform rectilinear grid
% (via interpolate_soln_to_regular_grid) so that streamline tracing,
% flux integration, and plotting all use consistent, evenly-spaced data.
%
% FIGURE LAYOUT
%   subplot(211) — growing time series (updated each step):
%       blue: melt flux leaving top of channel zone  [m^2/yr]
%       red : integrated d(phi)/dt below channel zone [m^2/yr]
%
%   subplot(212) — cross-section for current step:
%       imagesc: d(phi)/dt on regular grid
%       white line: backward streamline bounding the right edge of channel zone
%
% REQUIRES: read_aspect_paraview_output.m, interpolate_soln_to_regular_grid.m
clear;

%% ---- User settings --------------------------------------------------------
% NOTE: ASPECT outputs velocity in m/yr and coordinates in m.
%       Flux units are m^2/yr = km^2/Myr (exact equivalence).
%       dporosity_dt is in yr^-1 (divided by timestep in years in postprocess).
bwd_speed_thresh = 0.1;              % m/yr: min speed to continue backward streamline
sl_ds_frac       = 1/500;            % streamline step size as fraction of domain height
sl_min_length_km = 50;               % minimum streamline arc length [km]
dmdt_thresh      = 1e-9;             % yr^-1: only grid cells with dmdt > this contribute
                                     % to area_left and gen_below
dx              = 250;               % m: regular grid x-spacing ([] = auto)
dy              = 250;               % m: regular grid y-spacing ([] = auto)
pause_to_see    = false;
seconds_per_year = 60*60*24*365.25;


%% ---- Load data ------------------------------------------------------------
pvd_file = 'solution.pvd'
%pvd_file = strtrim(input('Enter path to solution.pvd [solution.pvd]: ', 's'));
%if isempty(pvd_file), pvd_file = 'solution.pvd'; end

t_yr_all = aspect_pvd_times(pvd_file);
n_all    = numel(t_yr_all);

fprintf('Found %d timestep(s):\n', n_all);
fprintf('  %6s  %s\n', 'Index', 'Simulation time (yr)');
fprintf('  %6s  %s\n', '-----', '--------------------');
for i = 1:n_all
    fprintf('  %6d  %g\n', i, t_yr_all(i));
end

sel_str = strtrim(input( ...
    sprintf('\nEnter timestep indices to process (e.g. "all", "1", "2:5", "1:2:end") [all]: '), 's'));
if isempty(sel_str) || strcmpi(sel_str, 'all')
    sel = 1:n_all;
else
    sel_str = strrep(sel_str, 'end', num2str(n_all));
    sel = round(eval(['[' sel_str ']']));
end
t_yr   = t_yr_all(sel);
nsteps = numel(sel);
t_kyr  = t_yr / 1e3;

%% ---- Output arrays --------------------------------------------------------
flux_top  = nan(nsteps, 1);   % m^2/yr  melt flux through top of channel zone
gen_below = nan(nsteps, 1);   % m^2/yr  integrated d(phi)/dt in source region

%% ---- Initialise figure ----------------------------------------------------
hfig = figure(3);  clf(hfig);

%% ---- Main loop ------------------------------------------------------------
for s = 1:nsteps
    fprintf('Step %d/%d  (index %d,  t = %.4g kyr)\n', s, nsteps, sel(s), t_kyr(s));

    step = read_aspect_paraview_output(pvd_file, sel(s));
    x_vg = step.x;
    y_vg = step.y;
    conn = step.connectivity;

    % --- Interpolate all fields onto regular grid ---------------------------
    [phi_rg, x, y] = interpolate_soln_to_regular_grid(step.porosity, x_vg, y_vg, conn, dx, dy);
    [u_f,  ~, ~] = interpolate_soln_to_regular_grid(step.u_f, x_vg, y_vg, conn, dx, dy);
    [dmdt, ~, ~] = interpolate_soln_to_regular_grid(step.dporosity_dt, x_vg, y_vg, conn, dx, dy);
    dmdt = dmdt .* seconds_per_year;
    nx   = numel(x);
    ny = numel(y);

    u_fx = u_f(:,:,1);   % [ny x nx]  horizontal melt velocity
    u_fy = u_f(:,:,2);   % [ny x nx]  vertical   melt velocity

    y_max = y(end);   y_min = y(1);
    domain_height = y_max - y_min;
    ds            = domain_height * sl_ds_frac;
    sl_min_length = sl_min_length_km * 1e3;

    % ---- (1) Find x_channel by scanning top row left to right -------------
    % For each grid column, trace a backward streamline from (x(ki), y_max).
    % x_channel = largest x(ki) whose streamline arc length >= sl_min_length.
    x_channel = NaN;
    sl        = zeros(0, 2);

    for ki = 1:nx
        sl_i  = backward_streamline_rg([x(ki), y_max], x, y, u_fx, u_fy, ...
                                        bwd_speed_thresh, ds, y_min);
        arc_i = streamline_length(sl_i);
        if arc_i >= sl_min_length
            x_channel = x(ki);
            sl        = sl_i;
        end
    end

    fprintf('  x_channel = %.1f km  |  sl_min_length = %.1f km\n', ...
            x_channel/1e3, sl_min_length/1e3);

    % ---- (2) Melt flux through top of channel zone ------------------------
    % Top row of regular grid: iy = ny.
    % Channel zone: columns with x <= x_channel and upward melt velocity.
    gen_below(s) = 0;
    area_sl      = 0;

    if isfinite(x_channel)
        channel_mask = (x <= x_channel) & (u_fy(ny,:) > 0);
        flux_top(s)  = sum(phi_rg(ny, channel_mask) .* u_fy(ny, channel_mask)) * dx;

        % ---- (3) Integrated d(phi)/dt left of streamline ------------------
        [gen_below(s), area_sl] = integrate_melt_generation_rg( ...
                                      x, y, dmdt, sl, dmdt_thresh, dx, dy);
    else
        flux_top(s) = 0;
        fprintf('  WARNING: no streamline reached sl_min_length — channel zone undefined.\n');
    end

    mean_dmdt = gen_below(s) / max(area_sl, eps);
    fprintf('  Fluxes: Melt Rise = %.4g  |  Melt Gen = %.4g  km^2/Myr\n', flux_top(s), gen_below(s));
    fprintf('  Area supplying melt = %.4g km^2  |  Mean dporosity/dt = %.4g yr^-1\n', ...
            area_sl/1e6, mean_dmdt);

    % ======================================================================
    % Flux time series
    % ======================================================================
    ax1 = subplot(211);
    cla(ax1);  hold(ax1, 'on');

    idx = 1:s;
    plot(ax1, t_kyr(idx), flux_top(idx),  'b-o', ...
         'MarkerSize', 3, 'LineWidth', 0.5, 'DisplayName', 'Top flux');
    plot(ax1, t_kyr(idx), gen_below(idx), 'r-o', ...
         'MarkerSize', 3, 'LineWidth', 0.5, 'DisplayName', 'Gen. below streamline');
    % plot(ax1, t_kyr(s), flux_top(s),  'bs', 'MarkerSize', 7, ...
    %      'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
    % plot(ax1, t_kyr(s), gen_below(s), 'rs', 'MarkerSize', 7, ...
    %      'MarkerFaceColor', 'r', 'HandleVisibility', 'off');

    xlim(ax1, [t_kyr(1), t_kyr(end)]);
    xlabel(ax1, 'Time (kyr)');
    ylabel(ax1, 'Flux (km^2 Myr^{-1})');
    legend(ax1, 'Location', 'best');
    grid(ax1, 'on');

    % ======================================================================
    % Crustal thickness time series
    % ======================================================================
    [velocity, ~, ~] = interpolate_soln_to_regular_grid(step.velocity, x_vg, y_vg, conn, dx, dy);
    U_half=velocity(1,end,1)*1e3;  %half spreading rate in km/Myr
    ax1 = subplot(212);
    cla(ax1);  hold(ax1, 'on');

    idx = 1:s;
    plot(ax1, t_kyr(idx), flux_top(idx)/U_half,  'b-o', ...
         'MarkerSize', 3, 'LineWidth', 0.5, 'DisplayName', 'Top flux');
    plot(ax1, t_kyr(idx), gen_below(idx)/U_half, 'r-o', ...
         'MarkerSize', 3, 'LineWidth', 0.5, 'DisplayName', 'Gen. below streamline');
    % plot(ax1, t_kyr(s), flux_top(s)/U_half,  'bs', 'MarkerSize', 7, ...
    %      'MarkerFaceColor', 'b', 'HandleVisibility', 'off');
    % plot(ax1, t_kyr(s), gen_below(s)/U_half, 'rs', 'MarkerSize', 7, ...
    %      'MarkerFaceColor', 'r', 'HandleVisibility', 'off');

    xlim(ax1, [t_kyr(1), t_kyr(end)]);
    ylim([0 8])
    xlabel(ax1, 'Time (kyr)');
    ylabel(ax1, 'Crustal thickness (km)');
    legend(ax1, 'Location', 'best');
    grid(ax1, 'on');
    if (pause_to_see)
    % ======================================================================
    % SUBPLOT (212): d(phi)/dt cross-section + streamline
    % ======================================================================
    climits = [0 1e-14] * seconds_per_year;

    ax2 = subplot(212);
    cla(ax2);  hold(ax2, 'on');

    imagesc(ax2, x/1e3, y/1e3, dmdt);
    set(ax2, 'YDir', 'normal');
    colormap(ax2, turbo);
    caxis(ax2, climits);

    cb = colorbar(ax2, 'Location', 'eastoutside');
    cb.Label.String = 'd\phi/dt  (yr^{-1})';

    if size(sl, 1) >= 2
        plot(ax2, sl(:,1)/1e3, sl(:,2)/1e3, 'w-',  'LineWidth', 1);
        plot(ax2, sl(1,1)/1e3, sl(1,2)/1e3, 'wv',  ...
             'MarkerSize', 6, 'MarkerFaceColor', 'w');
    end

    axis(ax2, 'equal', 'tight');
    xlabel(ax2, 'Distance (km)');
    ylabel(ax2, 'Depth (km)');
    if isfinite(x_channel)
        title(ax2, sprintf('d\\phi/dt  at  t = %.4g kyr  |  x_{channel} = %.1f km  (sl_{min} = %.0f km)', ...
              t_kyr(s), x_channel/1e3, sl_min_length_km));
    else
        title(ax2, sprintf('d\\phi/dt  at  t = %.4g yr  |  channel zone undefined', t_yr(s)));
    end

    drawnow;

    if s < nsteps
        fprintf('  Press any key for next step (%d/%d)...\n', s+1, nsteps);
        pause;
    end
    end %end no pause
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
    vsx      = vel_s(:,1);
    top_node = y_vg >= y_max;
    xn_top   = x_vg(top_node);
    vsx_top  = vsx(top_node);
    [xn_top, ord] = sort(xn_top);
    vsx_top  = vsx_top(ord);

    figure(4); clf;
    plot(xn_top/1e3, vsx_top, 'b.-');
    xlabel('x (km)');  ylabel('solid v_x  (m/yr)');
    title(sprintf('Solid horizontal velocity at top boundary  (t = %.4g yr)', t_yr(s)));
    grid on;
end

%% ==========================================================================
function v = slice_field(arr, s)
% Extract timestep s from a stacked array (no-op for single-step data).
    nd = ndims(arr);
    if nd == 2 && size(arr,2) > 1
        v = arr(:, s);
    elseif nd == 3
        v = arr(:, :, s);
    else
        v = arr;
    end
end

%% ==========================================================================
function L = streamline_length(sl)
% Total arc length of a streamline (sum of segment lengths).
    if size(sl, 1) < 2,  L = 0;  return;  end
    diffs = diff(sl, 1, 1);
    L = sum(hypot(diffs(:,1), diffs(:,2)));
end

%% ==========================================================================
function sl = backward_streamline_rg(pt0, x, y, ux_rg, uy_rg, thresh, ds, y_min)
% Trace melt streamline BACKWARD from pt0 using RK4, interpolating velocity
% bilinearly on the regular grid.
    max_pts = 10000;
    sl = nan(max_pts, 2);
    sl(1,:) = pt0;

    for k = 1:max_pts-1
        pt = sl(k,:);
        if pt(2) < y_min + ds,  sl = sl(1:k,:);  break;  end

        d1 = bwd_dir_rg(pt,            x, y, ux_rg, uy_rg, thresh);
        if isempty(d1),  sl = sl(1:k,:);  break;  end

        d2 = bwd_dir_rg(pt + ds/2*d1, x, y, ux_rg, uy_rg, thresh);
        if isempty(d2),  d2 = d1;  end

        d3 = bwd_dir_rg(pt + ds/2*d2, x, y, ux_rg, uy_rg, thresh);
        if isempty(d3),  d3 = d2;  end

        d4 = bwd_dir_rg(pt + ds*d3,   x, y, ux_rg, uy_rg, thresh);
        if isempty(d4),  d4 = d3;  end

        sl(k+1,:) = pt + ds/6 * (d1 + 2*d2 + 2*d3 + d4);

        if k == max_pts-1,  sl = sl(1:k,:);  end
    end
    sl = sl(all(isfinite(sl),2), :);
end

%% ==========================================================================
function d = bwd_dir_rg(pt, x, y, ux_rg, uy_rg, thresh)
% Unit backward-direction vector at pt via bilinear interpolation on the
% regular grid.  Returns [] if speed < thresh or pt is outside domain.
    ux = interp2(x(:)', y(:), ux_rg, pt(1), pt(2), 'linear', NaN);
    uy = interp2(x(:)', y(:), uy_rg, pt(1), pt(2), 'linear', NaN);
    if isnan(ux) || isnan(uy),  d = [];  return;  end
    spd = hypot(ux, uy);
    if spd < thresh
        d = [];
    else
        d = -[ux, uy] / spd;
    end
end

%% ==========================================================================
function [gen, area_left] = integrate_melt_generation_rg(x, y, dmdt_rg, sl, dmdt_thresh, dx, dy)
% Integrate dmdt_rg * dx * dy over all regular-grid cells that are:
%   (a) to the left of the bounding streamline, and
%   (b) dmdt_rg > dmdt_thresh
% The streamline x is interpolated at each grid row's y-value.
    gen = 0;  area_left = 0;
    cell_area = dx * dy;
    if size(sl,1) < 2,  return;  end

    [sl_y, ord] = sort(sl(:,2));  sl_x = sl(ord,1);
    [sl_y, ui]  = unique(sl_y);   sl_x = sl_x(ui);

    if numel(sl_y) < 2
        left_mask = x(:)' < sl_x(1);                       % [1 x nx]
        left_mask = repmat(left_mask, numel(y), 1);         % [ny x nx]
        mask      = left_mask & dmdt_rg > dmdt_thresh;
        gen       = sum(dmdt_rg(mask)) * cell_area;
        area_left = sum(mask(:))       * cell_area;
        return;
    end

    % Streamline x at each grid row (y clamped to streamline range)
    y_cl     = min(max(y(:), sl_y(1)), sl_y(end));
    x_sl_row = interp1(sl_y, sl_x, y_cl);   % [ny x 1]

    % Left-of-streamline: broadcast [ny x 1] against [1 x nx] → [ny x nx]
    left_mask = x(:)' <= x_sl_row(:);
    mask      = left_mask & dmdt_rg > dmdt_thresh;

    gen       = sum(dmdt_rg(mask)) * cell_area;
    area_left = sum(mask(:))       * cell_area;
end
