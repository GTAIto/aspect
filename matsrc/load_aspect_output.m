function data = load_aspect_output(pvd_file)
%LOAD_ASPECT_OUTPUT  Read ASPECT Paraview output and return as a struct.
%
% USAGE:
%   data = load_aspect_output('solution.pvd')
%
% Prompts for timestep selection, reads the selected steps, and returns
% a struct with fields (coordinates in metres, as written by ASPECT):
%
%   data.times        [Nsteps x 1]  simulation times (s)
%   data.x            [Npts x 1]    node x-coordinates (m)
%   data.y            [Npts x 1]    node y-coordinates (m)
%   data.z            [Npts x 1]    node z-coordinates (m, 0 for 2-D)
%   data.connectivity [Ncells x K]  1-based node indices
%   data.cell_types   [Ncells x 1]  VTK cell type codes
%   data.temperature  [Npts x Nsteps]        (scalar field example)
%   data.velocity     [Npts x Ndim x Nsteps] (vector field example)
%   ...one field per Paraview output variable...
%
% For a single selected timestep the trailing Nsteps dimension is absent.
%
% See also: read_aspect_output, plot_aspect_T

    data = read_aspect_output(pvd_file);
end
