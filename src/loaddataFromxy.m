function data_on_mesh = loaddataFromxy(md, filename)
	% load the data
	fId = fopen(filename);
	xydata = textscan(fId, '%f %f %f', 'Delimiter', '\t');
	fclose(fId);
	xdata = xydata{1};
	ydata = xydata{2};
	data = xydata{3};
	% construct grid data
	[x, ~, xid] = unique(xdata, 'sorted');
	[y, ~, yid] = unique(ydata, 'sorted');

	% create Cartesian grid
	[X, Y] = meshgrid(x, y);
	data_on_grid = NaN(size(X));
	for i = 1:length(xid)
		data_on_grid(yid(i), xid(i)) = data(i); 
	end

	% project to mesh
	data_on_mesh = InterpFromGridToMesh(x, y, data_on_grid, md.mesh.x, md.mesh.y, NaN);

	% fix the boundary nodes
	flags = isnan(data_on_mesh);
	pos1  = find(flags); pos2 = find(~flags);
	data_on_mesh(pos1) = griddata(md.mesh.x(pos2), md.mesh.y(pos2), data_on_mesh(pos2), md.mesh.x(pos1), md.mesh.y(pos1), 'nearest');
end

