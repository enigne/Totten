function hout = interpICESat2ATL1415(X,Y)

filename_h = '/totten_1/ModelData/Greenland/ICESat2_ATL1415/ATL14_GL_0314_100m_002_01.nc';
filename_dh = '/totten_1/ModelData/Greenland/ICESat2_ATL1415/ATL15_GL_0314_01km_002_01.nc';

xh = ncread(filename_h,'x');
yh = ncread(filename_h,'y');

xdh = ncread(filename_dh,'delta_h/x');
ydh = ncread(filename_dh,'delta_h/y');
tdh = ncread(filename_dh,'delta_h/time');

offset=2;

xmin=min(X(:)); xmax=max(X(:));

posxh=find(xh<=xmax);
id1xh=max(1,find(xh>=xmin,1)-offset);
id2xh=min(numel(xh),posxh(end)+offset);

posxdh=find(xdh<=xmax);
id1xdh=max(1,find(xdh>=xmin,1)-offset);
id2xdh=min(numel(xdh),posxdh(end)+offset);

ymin=min(Y(:)); ymax=max(Y(:));

posyh=find(yh>=ymin);
id1yh=max(1,find(yh<=ymax,1)-offset);
id2yh=min(numel(yh),posyh(end)+offset);

posydh=find(ydh>=ymin);
id1ydh=max(1,find(ydh<=ymax,1)-offset);
id2ydh=min(numel(ydh),posydh(end)+offset);

xh = xh(id1xh:id2xh);
yh = yh(id1yh:id2yh);

xdh = xdh(id1xdh:id2xdh);
ydh = ydh(id1ydh:id2ydh);

h = double(ncread(filename_h,'h',[id1xh id1yh],[id2xh-id1xh+1 id2yh-id1yh+1],[1 1]));
dh = double(ncread(filename_dh,'delta_h/delta_h',[id1xdh id1ydh 1],[id2xdh-id1xdh+1 id2ydh-id1ydh+1 length(tdh)],[1 1 1]));
geoid = interpBedmachineGreenland(X,Y,'geoid');

href = InterpFromGrid(xh,yh,h',X,Y) - geoid;
%hout(end+1) = 2020.0; %Data is the 2020 DEM

tref = datetime('2018-01-01-00-00-00','format','yyyy-MM-dd-hh-mm-ss');
ths = date2decyear(datenum(tref + days(tdh)));

hout = ones([length(Y)+1,length(tdh)]);
for ii = 1:length(tdh)
	hout(1:end-1,ii) = InterpFromGrid(xdh,ydh,dh(:,:,ii)',X,Y) + href;
	hout(end,ii) = ths(ii);
end
