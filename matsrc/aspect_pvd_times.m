function times = aspect_pvd_times(pvd_file)
%ASPECT_PVD_TIMES  Return simulation times from an ASPECT .pvd index file.
%
%   times = aspect_pvd_times(pvd_file)
%
%   Parses the PVD XML and returns a column vector of simulation times
%   without loading any VTU data.
    dom   = xmlread(pvd_file);
    nodes = dom.getElementsByTagName('DataSet');
    n     = nodes.getLength();
    times = zeros(n, 1);
    for i = 0:n-1
        times(i+1) = str2double(char(nodes.item(i).getAttribute('timestep')));
    end
end
