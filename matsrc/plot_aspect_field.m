% PLOT_ASPECT_FIELD  Plot an ASPECT 2-D scalar field for one or more timesteps.
%
% Run this script from the MATLAB command window or editor (F5).
% Prompts for the .pvd file and field name.
% If multiple timesteps are loaded and the mesh is consistent they are shown
% in sequence (press any key to advance).
% If the mesh changes across timesteps (AMR), each step uses its own mesh.
% Coordinates are converted from metres to km.

% --- Prompt for PVD file ---
% pvd_file = strtrim(input('Enter path to solution.pvd [solution.pvd]: ', 's'));
% if isempty(pvd_file)
pvd_file = 'solution.pvd';
% end

% --- Load data (also prompts for timestep selection) ---
data = read_aspect_paraview_output(pvd_file);

fns = fieldnames(data);
fprintf('Available fields in data:\n');
fprintf('  %s\n', fns{:});

fieldname = strtrim(input('Enter field name to plot [T]: ', 's'));
if isempty(fieldname)
    fieldname = 'T';
end

if ~isfield(data, fieldname)
    error('Field "%s" not found. Use one of the field names shown above.', fieldname);
end

% --- Streamline settings (edit these to adjust behaviour) ---
n_streamlines  = 20;          % number of streamlines
sl_start_y     = 15000;       % seed y-level in metres
sl_vel_field   = 'u_f';       % velocity field to use ('velocity' = solid, 'u_f' = fluid)
n_grid_x       = 300;         % interpolation grid resolution in x
n_grid_y       = 200;         % interpolation grid resolution in y

do_streamlines = isfield(data, sl_vel_field);
if ~do_streamlines
    warning('plot_aspect_field:noVelocity', ...
        'Velocity field "%s" not found — streamlines will be skipped.', sl_vel_field);
end

varying_mesh = iscell(data.x);   % true when AMR changed mesh across steps
field_all    = data.(fieldname);
nsteps       = numel(data.times);

for s = 1:nsteps
    t_val = data.times(min(s, numel(data.times)));

    if varying_mesh
        x_km  = data.x{s} / 1e3;
        y_km  = data.y{s} / 1e3;
        conn  = data.connectivity{s};
        field = field_all{s};                          % [Npts x Ncomp] or [Npts x 1]
        if do_streamlines; vel = data.(sl_vel_field){s}; else; vel = []; end
    else
        x_km = data.x / 1e3;
        y_km = data.y / 1e3;
        conn = data.connectivity;
        if ndims(field_all) == 3
            field = field_all(:,:,s);                  % vector, multi-step
        elseif nsteps > 1
            field = field_all(:, s);                   % scalar, multi-step
        else
            field = field_all;                         % single step (scalar or vector)
        end
        if do_streamlines
            v_all = data.(sl_vel_field);
            if ndims(v_all) == 3
                vel = v_all(:,:,s);
            else
                vel = v_all;
            end
        else
            vel = [];
        end
    end

    % If field has multiple components, plot its magnitude
    if size(field, 2) > 1
        plot_label = sprintf('|%s|  (m/yr)', fieldname);
        field      = sqrt(sum(field.^2, 2));
    else
        plot_label = fieldname;
        field      = field(:);
    end

    % Warn if field is entirely NaN for this step (absent in source VTU)
    if all(isnan(field))
        warning('plot_aspect_field:fieldAllNaN', ...
            'Field "%s" is all NaN at step %d (t = %.4g yr) — it may not have been written at this timestep.', ...
            plot_label, s, t_val);
    end

    % Filter cells that have any non-finite vertex value
    valid_cells = all(isfinite(field(conn)), 2);
    fprintf('  Step %d/%d: plotting %d / %d cells (%d excluded: non-finite)\n', ...
            s, nsteps, sum(valid_cells), numel(valid_cells), sum(~valid_cells));

    figure(1); clf;
    subplot(211)
    patch('Faces',          conn(valid_cells, :), ...
          'Vertices',        [x_km, y_km],         ...
          'FaceVertexCData', field,                 ...
          'FaceColor',       'interp',              ...
          'EdgeColor',       'none',                ...
          'LineWidth',       0.5)
    colormap(jet)

    % --- Streamlines ---
    if do_streamlines && ~isempty(vel)
        x_m = x_km * 1e3;   y_m = y_km * 1e3;
        u_x = vel(:, 1);     u_y = vel(:, 2);

        % Interpolate scattered velocity onto a regular grid
        xg = linspace(min(x_m), max(x_m), n_grid_x);
        yg = linspace(min(y_m), max(y_m), n_grid_y);
        [Xg, Yg] = meshgrid(xg, yg);

        Fu = scatteredInterpolant(x_m, y_m, u_x, 'linear', 'none');
        Fv = scatteredInterpolant(x_m, y_m, u_y, 'linear', 'none');
        Ug = Fu(Xg, Yg);
        Vg = Fv(Xg, Yg);

        % Seed points at sl_start_y, evenly spaced across x
        sx = linspace(min(x_m), max(x_m), n_streamlines);
        sy = repmat(sl_start_y, 1, n_streamlines);

        hold on
        h_sl = streamline(xg/1e3, yg/1e3, Ug, Vg, sx/1e3, sy/1e3);
        set(h_sl, 'Color', 'w', 'LineWidth', 0.5);
        hold off
    end

    axis equal tight
    xlabel('Distance (km)')
    ylabel('Depth (km)')
    title(sprintf('%s   t = %.4g yr', plot_label, t_val))

    cb = colorbar('Location', 'eastoutside');
    cb.Label.String = plot_label;
    %clim([-1e-14 1e-14])

    drawnow

    if nsteps > 1 && s < nsteps
        fprintf('Step %d/%d  (t = %.4g yr) — press any key for next\n', ...
                s, nsteps, t_val);
        pause
    end
end
