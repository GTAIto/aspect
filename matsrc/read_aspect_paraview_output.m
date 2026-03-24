function varargout = read_aspect_paraview_output(pvd_file, output_mat_file)
%READ_ASPECT_PARAVIEW_OUTPUT  Read ASPECT Paraview VTU output and save to a .mat file.
%
% USAGE:
%   read_aspect_paraview_output('solution.pvd')
%   read_aspect_paraview_output('path/to/solution.pvd', 'output.mat')
%
% INPUT:
%   pvd_file        - path to the .pvd index file (e.g. 'output/solution.pvd')
%   output_mat_file - (optional) output filename; defaults to <pvd_name>.mat
%                     in the same directory as pvd_file
%
% After parsing the PVD file the script prints the available timesteps and
% prompts you to select which ones to save.  Enter one of:
%   all            - read every timestep (default if you press Enter)
%   5              - read only timestep index 5
%   1 3 5          - read timestep indices 1, 3, and 5
%   2:10           - read indices 2 through 10
%   end            - read only the last timestep
%   1:2:end        - every other timestep (MATLAB colon notation supported)
%
% OUTPUT .mat file contains:
%   times      [N x 1] double  - simulation time at each step
%   timesteps  [1 x N] struct  - one struct per timestep with fields:
%       .coordinates  [Npts  x 3]   node (x,y,z) coordinates
%       .connectivity [Ncells x K]  1-based node indices (K nodes per cell)
%       .cell_types   [Ncells x 1]  VTK cell type codes
%       .<name>       [Npts  x C]   C-component data field named <name>
%                                   e.g. velocity [Npts x 3],
%                                        temperature [Npts x 1],
%                                        melting_rate [Npts x 1]
%   data.times        [Nsteps x 1]  simulation times (s)
%   data.x            [Npts x 1]    node x-coordinates (m)
%   data.y            [Npts x 1]    node y-coordinates (m)
%   data.z            [Npts x 1]    node z-coordinates (m, 0 for 2-D)
%   data.connectivity [Ncells x K]  1-based node indices
%   data.cell_types   [Ncells x 1]  VTK cell type codes
%   data.temperature  [Npts x Nsteps]        (scalar field example)
%   data.velocity     [Npts x Ndim x Nsteps] (vector field example)
%   ...one field per Paraview output variable, using the names from the VTU files...
%
% For a single selected timestep the trailing Nsteps dimension is absent.
%
% REQUIREMENTS:
%   MATLAB R2015b+ (uses java.util.Base64 for base64 decoding)
%   Java available (built-in to MATLAB) for zlib decompression
%
% NOTES:
%   - Handles both serial (.vtu) and parallel (.pvtu + .vtu pieces) output
%   - Handles ASCII, binary (base64), and zlib-compressed binary VTU formats
%   - ASPECT default output is zlib-compressed binary, which this handles
%   - Connectivity is converted from 0-based (VTK) to 1-based (MATLAB)
%   - Uses -v7.3 mat format to support files larger than 2 GB
%   - Field names are taken directly from the VTU file (with non-alphanumeric
%     characters replaced by underscores to form valid MATLAB field names)

    narginchk(1, 2);
    returning_data = (nargout > 0);
    if ~returning_data && nargin < 2
        [d, n] = fileparts(pvd_file);
        output_mat_file = fullfile(d, [n '.mat']);
    end

    fprintf('Reading %s\n', pvd_file);
    [times, rel_files] = parse_pvd(pvd_file);
    pvd_dir = fileparts(pvd_file);
    if isempty(pvd_dir), pvd_dir = '.'; end

    n_steps = numel(times);
    fprintf('Found %d timestep(s):\n', n_steps);
    fprintf('  %6s  %s\n', 'Index', 'Simulation time');
    fprintf('  %6s  %s\n', '-----', '---------------');
    for i = 1:n_steps
        fprintf('  %6d  %g\n', i, times(i));
    end

    % --- Prompt user for which timesteps to read ---
    sel_str = strtrim(input( ...
        sprintf('\nEnter timestep indices to read (e.g. "all", "1", "2:5", "1:2:end") [all]: '), ...
        's'));
    if isempty(sel_str) || strcmpi(sel_str, 'all')
        sel = 1:n_steps;
    else
        % Evaluate as MATLAB expression with 'end' replaced by n_steps
        sel_str = regexprep(sel_str, '\bend\b', num2str(n_steps));
        try
            sel = eval(['[' sel_str ']']);
        catch
            error('read_aspect_paraview_output:badSelection', ...
                  'Could not parse selection "%s". Use MATLAB index notation.', sel_str);
        end
        sel = round(sel);
        if any(sel < 1) || any(sel > n_steps)
            error('read_aspect_paraview_output:outOfRange', ...
                  'Index out of range. Valid range is 1 to %d.', n_steps);
        end
    end
    fprintf('Reading %d timestep(s).\n', numel(sel));

    timesteps = struct([]);
    times     = times(sel);   % trim times vector to match selection
    for ii = 1:numel(sel)
        i = sel(ii);
        fprintf('  Step %d/%d  (index %d)  t = %g\n', ii, numel(sel), i, times(ii));
        fpath = fullfile(pvd_dir, strrep(rel_files{i}, '/', filesep));
        [~, ~, ext] = fileparts(fpath);
        if strcmpi(ext, '.pvtu')
            ts = read_pvtu(fpath);
        else
            ts = read_vtu(fpath);
        end
        if isempty(timesteps)
            timesteps = ts;
        else
            % MATLAB struct arrays require identical fields on every element.
            % Early timesteps (e.g. t=0) may lack fields that appear later
            % (e.g. melting_rate), so pad missing fields with [] before appending.
            for f = fieldnames(ts)'
                if ~isfield(timesteps, f{1})
                    [timesteps.(f{1})] = deal([]);
                end
            end
            for f = fieldnames(timesteps)'
                if ~isfield(ts, f{1})
                    ts.(f{1}) = [];
                end
            end
            timesteps(end+1) = ts; %#ok<AGROW>
        end
    end

    % --- Flatten into top-level variables for easy workspace access ---
    out.times = times;                                     % [Nsteps x 1]

    % Check whether the mesh is the same for every selected timestep
    npts_all       = arrayfun(@(s) size(s.coordinates, 1), timesteps);
    mesh_consistent = numel(unique(npts_all)) == 1;

    if mesh_consistent
        % Common case: mesh fixed across steps
        out.x            = timesteps(1).coordinates(:, 1);    % [Npts x 1]
        out.y            = timesteps(1).coordinates(:, 2);    % [Npts x 1]
        out.z            = timesteps(1).coordinates(:, 3);    % [Npts x 1]
        out.connectivity = timesteps(1).connectivity;         % [Ncells x K]
        out.cell_types   = timesteps(1).cell_types;           % [Ncells x 1]
    else
        % AMR changed the mesh: store geometry per step as cell arrays
        fprintf('  Note: mesh size differs across timesteps (%s nodes) — ', ...
                num2str(npts_all));
        fprintf('x/y/z/connectivity stored as cell arrays.\n');
        out.x            = arrayfun(@(s) s.coordinates(:,1), timesteps, 'UniformOutput', false);
        out.y            = arrayfun(@(s) s.coordinates(:,2), timesteps, 'UniformOutput', false);
        out.z            = arrayfun(@(s) s.coordinates(:,3), timesteps, 'UniformOutput', false);
        out.connectivity = arrayfun(@(s) s.connectivity,     timesteps, 'UniformOutput', false);
        out.cell_types   = arrayfun(@(s) s.cell_types,       timesteps, 'UniformOutput', false);
    end

    % Data fields
    % Consistent mesh : scalar -> [Npts x Nsteps], vector -> [Npts x Ncomp x Nsteps]
    % Varying mesh    : any field -> {[Npts_i x Ncomp], ...} cell array (one per step)
    skip = {'coordinates', 'connectivity', 'cell_types'};
    fn_list = fieldnames(timesteps(1));
    nsteps  = numel(timesteps);
    for k = 1:numel(fn_list)
        fn = fn_list{k};
        if ismember(fn, skip), continue; end
        slices = arrayfun(@(s) s.(fn), timesteps, 'UniformOutput', false);

        % Warn and pad if this field is absent in some timesteps (empty []).
        % This happens when a field (e.g. melting_rate) is only written after
        % the simulation starts generating melt.
        empty_mask = cellfun(@isempty, slices);
        if any(empty_mask)
            warning('read_aspect_paraview_output:missingField', ...
                'Field "%s" absent in %d of %d timestep(s) — those steps will be NaN.', ...
                fn, sum(empty_mask), nsteps);
            ref = find(~empty_mask, 1);
            if ~isempty(ref)
                % Consistent mesh: all non-empty slices have the same size
                if mesh_consistent
                    fill = nan(size(slices{ref}));
                    slices(empty_mask) = {fill};
                else
                    % Varying mesh: each empty step needs its own Npts
                    for ei = find(empty_mask)'
                        npts_ei = size(timesteps(ei).coordinates, 1);
                        ncomp   = size(slices{ref}, 2);
                        slices{ei} = nan(npts_ei, ncomp);
                    end
                end
            end
        end

        if nsteps == 1
            out.(fn) = slices{1};                      % [Npts x Ncomp]
        elseif mesh_consistent
            stacked  = cat(3, slices{:});              % [Npts x Ncomp x Nsteps]
            out.(fn) = squeeze(stacked);               % drops size-1 Ncomp for scalars
        else
            out.(fn) = slices;                         % {[Npts_i x Ncomp], ...}
        end
    end

    if returning_data
        % Return the struct directly — no file written
        varargout{1} = out;
    else
        % Save each field as a top-level workspace variable
        fprintf('Saving to %s\n', output_mat_file);
        save(output_mat_file, '-struct', 'out', '-v7.3');
        fprintf('Done.  Variables saved:\n');
        fns = fieldnames(out);
        for k = 1:numel(fns)
            sz = size(out.(fns{k}));
            fprintf('  %-25s  [%s]\n', fns{k}, num2str(sz, '%d x '));
        end
        fprintf(['\nTo plot on the unstructured quad mesh use patch(), e.g.:\n' ...
                 '  load(''%s'')\n' ...
                 '  patch(''Faces'',connectivity,''Vertices'',[x y]/1e3,...\n' ...
                 '        ''FaceVertexCData'',temperature,''FaceColor'',''interp'',' ...
                 '''EdgeColor'',''none'')\n' ...
                 '  colorbar, axis equal tight\n'], output_mat_file);
    end

end

%% =========================================================================
function [times, files] = parse_pvd(pvd_file)
% Parse a .pvd XML file and return simulation times and VTU file paths.
    dom   = xmlread(pvd_file);
    nodes = dom.getElementsByTagName('DataSet');
    n     = nodes.getLength();
    times = zeros(n, 1);
    files = cell(n, 1);
    for i = 0:n-1
        node       = nodes.item(i);
        times(i+1) = str2double(char(node.getAttribute('timestep')));
        files{i+1} = char(node.getAttribute('file'));
    end
end

%% =========================================================================
function ts = read_pvtu(pvtu_file)
% Read a parallel .pvtu file by reading all its .vtu piece files and merging.
    pvtu_dir = fileparts(pvtu_file);
    dom      = xmlread(pvtu_file);
    pieces   = dom.getElementsByTagName('Piece');
    n_pieces = pieces.getLength();
    parts    = cell(n_pieces, 1);
    for i = 0:n_pieces-1
        src        = char(pieces.item(i).getAttribute('Source'));
        parts{i+1} = read_vtu(fullfile(pvtu_dir, src));
    end
    ts = merge_pieces(parts);
end

%% =========================================================================
function ts = merge_pieces(parts)
% Concatenate data from multiple VTU pieces into one struct.
    ts   = parts{1};
    skip = {'coordinates', 'connectivity', 'cell_types'};
    for i = 2:numel(parts)
        d  = parts{i};
        n0 = size(ts.coordinates, 1);   % point count so far (before appending)

        % Connectivity indices must be offset by the current point count
        ts.connectivity = [ts.connectivity; d.connectivity + n0];
        ts.cell_types   = [ts.cell_types;   d.cell_types];
        ts.coordinates  = [ts.coordinates;  d.coordinates];

        % Append all data fields
        fn_list = fieldnames(d);
        for k = 1:numel(fn_list)
            fn = fn_list{k};
            if ismember(fn, skip), continue; end
            if isfield(ts, fn)
                ts.(fn) = [ts.(fn); d.(fn)];
            else
                % Field missing in earlier pieces — pad with NaN
                nc = size(d.(fn), 2);
                ts.(fn) = [NaN(n0, nc); d.(fn)];
            end
        end
    end
end

%% =========================================================================
function ts = read_vtu(vtu_file)
% Read a single .vtu (UnstructuredGrid) file and return a struct.
% Handles files with multiple <Piece> blocks (produced when ASPECT groups
% output from several MPI processes into one file via "Number of grouped
% files" < number of MPI processes).
    dom = xmlread(vtu_file);

    % --- Global file attributes ---
    vtkNode    = dom.getElementsByTagName('VTKFile').item(0);
    byte_order = char(vtkNode.getAttribute('byte_order'));
    hdr_type   = char(vtkNode.getAttribute('header_type'));
    compressor = char(vtkNode.getAttribute('compressor'));
    if isempty(hdr_type), hdr_type = 'UInt32'; end
    is_le    = strcmpi(byte_order, 'LittleEndian');
    use_zlib = ~isempty(regexpi(compressor, 'zlib'));

    % --- Iterate over all <Piece> blocks in this VTU file ---
    piece_nodes = dom.getElementsByTagName('Piece');
    n_vtu_pieces = piece_nodes.getLength();
    parts = cell(n_vtu_pieces, 1);

    for pi = 0:n_vtu_pieces-1
        piece = piece_nodes.item(pi);
        ps = struct();

        % Node coordinates
        pts_node = piece.getElementsByTagName('Points').item(0);
        da       = pts_node.getElementsByTagName('DataArray').item(0);
        coords   = read_da(da, is_le, hdr_type, use_zlib);
        ps.coordinates = reshape(coords, 3, [])';   % [Npts x 3]

        % Cell topology
        raw_conn = []; offsets = []; cell_types_vec = [];
        cells_node = piece.getElementsByTagName('Cells').item(0);
        das = cells_node.getElementsByTagName('DataArray');
        for i = 0:das.getLength()-1
            da   = das.item(i);
            name = char(da.getAttribute('Name'));
            vals = read_da(da, is_le, hdr_type, use_zlib);
            switch name
                case 'connectivity', raw_conn       = vals + 1;
                case 'offsets',      offsets         = vals;
                case 'types',        cell_types_vec  = vals;
            end
        end
        n_cells        = numel(offsets);
        nodes_per_cell = offsets(1);
        ps.connectivity = reshape(raw_conn, nodes_per_cell, n_cells)';
        ps.cell_types   = cell_types_vec;

        % PointData
        pd = piece.getElementsByTagName('PointData').item(0);
        if ~isempty(pd)
            ps = read_data_section(ps, pd, is_le, hdr_type, use_zlib, '');
        end

        % CellData
        cd = piece.getElementsByTagName('CellData').item(0);
        if ~isempty(cd)
            ps = read_data_section(ps, cd, is_le, hdr_type, use_zlib, 'cell_');
        end

        parts{pi+1} = ps;
    end

    % Merge all pieces within this VTU file
    if n_vtu_pieces == 1
        ts = parts{1};
    else
        ts = merge_pieces(parts);
    end

    % --- Diagnostics ---
    fn_list = fieldnames(ts);
    for k = 1:numel(fn_list)
        fn = fn_list{k};
        if isnumeric(ts.(fn))
            nnan = sum(isnan(ts.(fn)(:)));
            if nnan > 0
                fprintf('      WARNING: %d NaN in field "%s"\n', nnan, fn);
            end
        end
    end
end

%% =========================================================================
function ts = read_data_section(ts, section_node, is_le, hdr_type, use_zlib, prefix)
% Read all DataArrays in a PointData or CellData section into struct fields.
    das = section_node.getElementsByTagName('DataArray');
    for i = 0:das.getLength()-1
        da     = das.item(i);
        name   = char(da.getAttribute('Name'));
        n_comp = str2double(char(da.getAttribute('NumberOfComponents')));
        if isnan(n_comp) || n_comp < 1, n_comp = 1; end
        vals = read_da(da, is_le, hdr_type, use_zlib);
        fld  = valid_fieldname([prefix name]);
        ts.(fld) = reshape(vals, n_comp, [])';   % [Npts x Ncomp]
    end
end

%% =========================================================================
function vals = read_da(da_node, is_le, hdr_type, use_zlib)
% Read one DataArray XML node and return its values as a double column vector.
    fmt      = char(da_node.getAttribute('format'));
    type_str = char(da_node.getAttribute('type'));
    text     = strtrim(char(da_node.getTextContent()));

    if strcmpi(fmt, 'ascii')
        vals = double(str2num(text)); %#ok<ST2NM>
        return;
    end

    if strcmpi(fmt, 'appended')
        error('read_aspect_paraview_output:appended', ...
              'Appended format VTU not supported. Set ASPECT output format to binary or ascii.');
    end

    % Binary format: inline base64-encoded data, optionally zlib-compressed
    raw = b64decode(text);

    if use_zlib
        data = vtk_zlib_decomp(raw, hdr_type);
    else
        hb   = hdr_nbytes(hdr_type);
        data = raw(hb+1:end);          % strip the leading byte-count header
    end

    vals = bytes_to_double(data, type_str, is_le);
end

%% =========================================================================
function data = vtk_zlib_decomp(raw, hdr_type)
% Parse VTK zlib block-compression header and decompress all blocks.
%
% VTK zlib header layout (each field is hdr_nbytes wide):
%   [nblocks] [uncompressed_block_size] [last_partial_block_size]
%   [compressed_size_0] [compressed_size_1] ... [compressed_size_{N-1}]
    hb = hdr_nbytes(hdr_type);
    if strcmpi(hdr_type, 'UInt64'), ht = 'uint64'; else, ht = 'uint32'; end

    nblocks = double(typecast(raw(1:hb), ht));
    % fields 2 and 3 (uncompressed sizes) are not needed for decompression

    comp_sizes = zeros(1, nblocks);
    for b = 1:nblocks
        o = (3 + b - 1) * hb;
        comp_sizes(b) = double(typecast(raw(o+1 : o+hb), ht));
    end

    header_end = (3 + nblocks) * hb;
    pos   = header_end + 1;
    parts = cell(nblocks, 1);
    for b = 1:nblocks
        block    = raw(pos : pos + comp_sizes(b) - 1);
        parts{b} = zlib_inflate(block);
        pos      = pos + comp_sizes(b);
    end
    data = vertcat(parts{:});
end

%% =========================================================================
function out = zlib_inflate(compressed)
% Decompress one zlib-format block using Java InflaterInputStream.
% Uses InterruptibleStreamCopier (built into MATLAB) to avoid needing
% javaArray, which doesn't accept Java primitive type names.
    import java.io.ByteArrayInputStream
    import java.io.ByteArrayOutputStream
    import java.util.zip.InflaterInputStream
    import com.mathworks.mlwidgets.io.InterruptibleStreamCopier

    % Java byte[] is signed int8; reinterpret MATLAB uint8 bits as int8
    signed = typecast(compressed(:), 'int8');

    bais = ByteArrayInputStream(signed);
    zis  = InflaterInputStream(bais);
    baos = ByteArrayOutputStream();

    isc = InterruptibleStreamCopier.getInterruptibleStreamCopier();
    isc.copyStream(zis, baos);
    zis.close();

    % Convert Java signed byte[] back to MATLAB uint8
    out = typecast(int8(baos.toByteArray()), 'uint8');
end

%% =========================================================================
function bytes = b64decode(str)
% Pure-MATLAB base64 decoder that handles deal.II / VTK output correctly.
%
% deal.II writes the compression header and the compressed data as two
% SEPARATE base64 strings, each with its own '=' padding, concatenated
% in the XML text.  E.g.:  "...header_b64...==...data_b64...="
% Java decoders fail when they see '=' in the middle of the stream.
% This decoder splits on '=' boundaries and decodes each chunk separately.

    str = char(str);
    str = str(~isspace(str));   % strip all whitespace

    % Find each chunk: maximal run of base64 chars followed by '=' padding.
    % This naturally splits the two deal.II blobs at the '=' boundary.
    tokens = regexp(str, '([A-Za-z0-9+/]+)(=*)', 'tokens');

    parts = cell(numel(tokens), 1);
    for i = 1:numel(tokens)
        parts{i} = b64_chunk(tokens{i}{1});
    end
    bytes = vertcat(parts{:});
end

%% -------------------------------------------------------------------------
function bytes = b64_chunk(str)
% Decode one padding-free base64 chunk to a uint8 column vector.
    n = numel(str);
    if n == 0, bytes = zeros(0, 1, 'uint8'); return; end

    % Build 256-element lookup table (index = char code + 1)
    lut = zeros(1, 256, 'uint32');
    b64 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
    lut(double(b64) + 1) = 0:63;
    vals = lut(double(str) + 1);   % [1 x n] uint32, each 0-63

    % Pad to a multiple of 4 with zeros (the extra bits will be discarded)
    r = mod(n, 4);
    if r > 0
        vals(end + 1 : end + (4 - r)) = 0;
    end

    % Vectorised decode: every 4 six-bit values -> 3 bytes
    vals    = reshape(vals, 4, []);                          % [4 x ngroups]
    combined = bitshift(vals(1,:), 18) + bitshift(vals(2,:), 12) + ...
               bitshift(vals(3,:),  6) + vals(4,:);          % [1 x ngroups] uint32

    raw = [uint8(bitshift(combined, -16)); ...
           uint8(bitand(bitshift(combined, -8), 255)); ...
           uint8(bitand(combined, 255))];                    % [3 x ngroups]

    % Flatten column-major (b0,b1,b2 of group1, b0,b1,b2 of group2, ...)
    % and keep only the true byte count (padding zeros add spurious bytes)
    bytes = raw(:);
    bytes = bytes(1 : floor(n * 3 / 4));
end

%% =========================================================================
function n = hdr_nbytes(hdr_type)
% Number of bytes in a VTK data-array header field.
    if strcmpi(hdr_type, 'UInt64'), n = 8; else, n = 4; end
end

%% =========================================================================
function vals = bytes_to_double(bytes, type_str, is_le)
% Reinterpret a uint8 byte array as the given VTK type, then return double.
    bytes = bytes(:);
    switch lower(type_str)
        case 'float32', v = typecast(bytes, 'single');
        case 'float64', v = typecast(bytes, 'double');
        case 'int8',    v = typecast(bytes, 'int8');
        case 'int16',   v = typecast(bytes, 'int16');
        case 'int32',   v = typecast(bytes, 'int32');
        case 'int64',   v = typecast(bytes, 'int64');
        case 'uint8',   v = bytes;
        case 'uint16',  v = typecast(bytes, 'uint16');
        case 'uint32',  v = typecast(bytes, 'uint32');
        case 'uint64',  v = typecast(bytes, 'uint64');
        otherwise
            warning('read_aspect_paraview_output:unknownType', ...
                    'Unknown VTK type "%s"; treating as float64.', type_str);
            v = typecast(bytes, 'double');
    end
    if ~is_le
        v = swapbytes(v);
    end
    vals = double(v(:));
end

%% =========================================================================
function name = valid_fieldname(str)
% Convert an arbitrary string to a valid MATLAB struct field name.
    name = regexprep(str, '[^a-zA-Z0-9_]', '_');
    if ~isempty(name) && ~isempty(regexp(name(1), '\d', 'once'))
        name = ['f_' name];
    end
    if isempty(name), name = 'unnamed_field'; end
end
