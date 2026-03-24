% PLOT_ASPECT_FIELD  Plot an ASPECT 2-D scalar field for one or more timesteps.
%
% Run this script from the MATLAB command window or editor (F5).
% Prompts for the .pvd file and field name.
% Timesteps are shown in sequence (press any key to advance).
% Coordinates are converted from metres to km.
%
% Fields are interpolated onto a uniform rectilinear grid before plotting
% (via interpolate_soln_to_regular_grid).  Set dx/dy below to override the
% automatic grid spacing (default: 2nd-smallest cell dimension).

% --- Prompt for PVD file ---
% pvd_file = strtrim(input('Enter path to solution.pvd [solution.pvd]: ', 's'));
% if isempty(pvd_file)
pvd_file = 'solution.pvd';
% end

% --- Get timestep list and let user select which to plot ---
t_yr_all = aspect_pvd_times(pvd_file);
n_all    = numel(t_yr_all);

fprintf('Found %d timestep(s):\n', n_all);
fprintf('  %6s  %s\n', 'Index', 'Simulation time (yr)');
fprintf('  %6s  %s\n', '-----', '--------------------');
for i = 1:n_all
    fprintf('  %6d  %g\n', i, t_yr_all(i));
end

sel_str = strtrim(input( ...
    sprintf('\nEnter timestep indices to plot (e.g. "all", "1", "2:5", "1:2:end") [all]: '), 's'));
if isempty(sel_str) || strcmpi(sel_str, 'all')
    sel = 1:n_all;
else
    sel_str = strrep(sel_str, 'end', num2str(n_all));
    sel = round(eval(['[' sel_str ']']));
end
t_yr   = t_yr_all(sel);
nsteps = numel(sel);

% --- Load step 1 to discover available fields ---
step1 = read_aspect_paraview_output(pvd_file, 1);
fns = fieldnames(step1);
fprintf('Available fields in data:\n');
fprintf('  %s\n', fns{:});

fieldname = strtrim(input('Enter field name to plot [T]: ', 's'));
if isempty(fieldname)
    fieldname = 'T';
end

if ~isfield(step1, fieldname)
    error('Field "%s" not found. Use one of the field names shown above.', fieldname);
end

% --- User settings ----------------------------------------------------------
dx             = 250;         % m: regular grid x-spacing ([] = auto)
dy             = 250;         % m: regular grid y-spacing ([] = auto)
n_streamlines  = 20;          % number of streamlines
sl_start_y     = 15000;       % seed y-level in metres
sl_vel_field   = 'u_f';       % velocity field ('velocity' = solid, 'u_f' = fluid)
% ---------------------------------------------------------------------------

do_streamlines = isfield(step1, sl_vel_field);
if ~do_streamlines
    warning('plot_aspect_field:noVelocity', ...
        'Velocity field "%s" not found — streamlines will be skipped.', sl_vel_field);
end
clear step1;

for s = 1:nsteps
    t_val = t_yr(s);
    fprintf('Step %d/%d  (index %d,  t = %.4g yr)\n', s, nsteps, sel(s), t_val);

    step     = read_aspect_paraview_output(pvd_file, sel(s));
    x_m      = step.x;
    y_m      = step.y;
    conn     = step.connectivity;
    field_vg = step.(fieldname);
    if do_streamlines; u_f_vg = step.(sl_vel_field); else; u_f_vg = []; end

    % If field has multiple components, plot its magnitude
    if size(field_vg, 2) > 1
        plot_label = sprintf('|%s|  (m/yr)', fieldname);
        field_vg   = sqrt(sum(field_vg.^2, 2));
    else
        plot_label = fieldname;
        field_vg   = field_vg(:);
    end

    if all(isnan(field_vg))
        warning('plot_aspect_field:fieldAllNaN', ...
            'Field "%s" is all NaN at step %d (t = %.4g yr).', ...
            plot_label, s, t_val);
    end

    % --- Interpolate field onto regular grid --------------------------------
    [field, x, y] = interpolate_soln_to_regular_grid( ...
                        field_vg, x_m, y_m, conn, dx, dy);

    % --- Plot ---------------------------------------------------------------
    figure(1); clf;
    subplot(211)

    imagesc(x/1e3, y/1e3, field);
    set(gca, 'YDir', 'normal');
    colormap(jet);

    % --- Streamlines --------------------------------------------------------
    if do_streamlines && ~isempty(u_f_vg)
        [vel_rg, ~, ~] = interpolate_soln_to_regular_grid( ...
                             u_f_vg, x_m, y_m, conn, dx, dy);
        Ug = vel_rg(:,:,1);
        Vg = vel_rg(:,:,2);

        sx = linspace(min(x_m), max(x_m), n_streamlines);
        sy = repmat(sl_start_y, 1, n_streamlines);

        hold on
        h_sl = streamline(x/1e3, y/1e3, Ug, Vg, sx/1e3, sy/1e3);
        set(h_sl, 'Color', 'w', 'LineWidth', 0.5);
        hold off
    end

    axis equal tight
    xlabel('Distance (km)')
    ylabel('Depth (km)')
    title(sprintf('%s   t = %.4g yr', plot_label, t_val))

    cb = colorbar('Location', 'eastoutside');
    cb.Label.String = plot_label;

    drawnow

    if nsteps > 1 && s < nsteps
        fprintf('Step %d/%d  (t = %.4g yr) — press any key for next\n', ...
                s, nsteps, t_val);
        pause
    end
end
