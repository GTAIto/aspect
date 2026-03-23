% PLOT_ASPECT_T  Plot ASPECT 2-D temperature field in figure(1) subplot(2,1,1).
%
% Run this script from the MATLAB command window or editor (F5).
% It will prompt for the .pvd file path and the timestep(s) to display.
% Converts spatial coordinates from metres to km.
% Plots temperature on the unstructured quad mesh using patch().
% If multiple timesteps are loaded, each is shown in sequence
% (press any key to advance).
% Time is in years, distances in meters

% --- Prompt for PVD file ---
pvd_file = strtrim(input('Enter path to solution.pvd [solution.pvd]: ', 's'));
if isempty(pvd_file)
    pvd_file = 'solution.pvd';
end

% --- Load data (also prompts for timestep selection) ---
data = read_aspect_output(pvd_file);

x_km = data.x / 1e3;   % m -> km
y_km = data.y / 1e3;

fns = fieldnames(data);
fprintf('Available fields in data:\n');
fprintf('  %s\n', fns{:});

if ~isfield(data, 'T')
    error('Field "T" not found. Use one of the field names shown above.');
end

% T is [Npts x Nsteps] for multiple steps, [Npts x 1] for one
T_all = data.T;
if isvector(T_all)
    T_all = T_all(:);
end
nsteps = size(T_all, 2);

for s = 1:nsteps
    T     = T_all(:, s);
    t_val = data.times(min(s, numel(data.times)));

    % Filter out cells that have any NaN-temperature vertex (outside domain)
    valid_cells = all(isfinite(T(data.connectivity)), 2);
    fprintf('  Plotting %d / %d cells (%d excluded: NaN T)\n', ...
            sum(valid_cells), numel(valid_cells), sum(~valid_cells));

    figure(1); clf;
    subplot(211)
    patch('Faces',          data.connectivity(valid_cells, :), ...
          'Vertices',        [x_km, y_km],                     ...
          'FaceVertexCData', T,                                 ...
          'FaceColor',       'interp',                          ...
          'EdgeColor',       'none',                               ...
          'LineWidth',       0.5)
    colormap(jet)

    axis equal tight
    xlabel('Distance (km)')
    ylabel('Depth (km)')
    title(sprintf('Temperature (K)   t = %.4g yr', t_val))

    cb = colorbar('Location', 'eastoutside');
    cb.Label.String = 'T  (K)';

    drawnow

    if nsteps > 1 && s < nsteps
        fprintf('Step %d/%d  (t = %.4g yr) — press any key for next\n', ...
                s, nsteps, t_val);
        pause
    end
end
