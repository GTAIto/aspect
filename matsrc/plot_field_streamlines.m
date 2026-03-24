    
    field_all    = data.(fieldname);
    if varying_mesh
        nsteps = numel(data.times);
    else if isvector(field_all)
        field_all = field_all(:);
    end
    nsteps = size(field_all, 2);
    end

        if varying_mesh
        x_km   = data.x{s}   / 1e3;
        y_km   = data.y{s}   / 1e3;
        conn   = data.connectivity{s};
        field  = field_all{s};
        if do_streamlines; vel = data.(sl_vel_field){s}; else; vel = []; end
    else
        x_km   = data.x / 1e3;
        y_km   = data.y / 1e3;
        conn   = data.connectivity;
        field  = field_all(:, s);
        if do_streamlines
            v_all = data.(sl_vel_field);
            if ndims(v_all) == 3
                vel = v_all(:, :, s);   % [Npts x Ncomp] for this step
            else
                vel = v_all;            % single timestep
            end
        else
            vel = [];
        end
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

    % Warn if field is entirely NaN for this step (absent in source VTU)
    if all(isnan(field))
        warning('plot_aspect_field:fieldAllNaN', ...
            'Field "%s" is all NaN at step %d (t = %.4g yr) — it may not have been written at this timestep.', ...
            fieldname, s, t_val);
    end

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
        set(h_sl, 'Color', 'w', 'LineWidth', 1);
        hold off
    end

    axis equal tight
    xlabel('Distance (km)')
    ylabel('Depth (km)')
    title(sprintf('%s   t = %.4g yr', fieldname, t_val))

    cb = colorbar('Location', 'eastoutside');
    cb.Label.String = fieldname;
    clim([-1e-14 1e-14])

    drawnow