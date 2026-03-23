% PLOT_ASPECT_FIELD  Plot an ASPECT 2-D scalar field for one or more timesteps.
%
% Run this script from the MATLAB command window or editor (F5).
% Prompts for the .pvd file and field name.
% If multiple timesteps are loaded and the mesh is consistent they are shown
% in sequence (press any key to advance).
% If the mesh changes across timesteps (AMR), each step uses its own mesh.
% Coordinates are converted from metres to km.

% --- Prompt for PVD file ---
pvd_file = strtrim(input('Enter path to solution.pvd [solution.pvd]: ', 's'));
if isempty(pvd_file)
    pvd_file = 'solution.pvd';
end

% --- Load data (also prompts for timestep selection) ---
data = read_aspect_output(pvd_file);

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

varying_mesh = iscell(data.x);   % true when AMR changed mesh across steps
field_all    = data.(fieldname);

if varying_mesh
    nsteps = numel(data.times);
else
    if isvector(field_all)
        field_all = field_all(:);
    end
    nsteps = size(field_all, 2);
end

for s = 1:nsteps
    t_val = data.times(min(s, numel(data.times)));

    if varying_mesh
        x_km   = data.x{s}   / 1e3;
        y_km   = data.y{s}   / 1e3;
        conn   = data.connectivity{s};
        field  = field_all{s};
    else
        x_km   = data.x / 1e3;
        y_km   = data.y / 1e3;
        conn   = data.connectivity;
        field  = field_all(:, s);
    end

    % Warn if field is entirely NaN for this step (absent in source VTU)
    if all(isnan(field))
        warning('plot_aspect_field:fieldAllNaN', ...
            'Field "%s" is all NaN at step %d (t = %.4g yr) — it may not have been written at this timestep.', ...
            fieldname, s, t_val);
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

    axis equal tight
    xlabel('Distance (km)')
    ylabel('Depth (km)')
    title(sprintf('%s   t = %.4g yr', fieldname, t_val))

    cb = colorbar('Location', 'eastoutside');
    cb.Label.String = fieldname;

    drawnow

    if nsteps > 1 && s < nsteps
        fprintf('Step %d/%d  (t = %.4g yr) — press any key for next\n', ...
                s, nsteps, t_val);
        pause
    end
end
