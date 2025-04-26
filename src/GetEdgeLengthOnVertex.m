function edgeLength=GetEdgeLengthOnVertex(elements, x, y)
%GetEdgeLength - compute the length of edges for each vertices
% current version is an approximation by averaging over all the connected
% edges of the vertex

%get number of nodes
nx = length(x);

nodec = NodeConnectivity(elements, nx);

%initialization
x1=x(elements(:,1)); x2=x(elements(:,2));x3=x(elements(:,3));
y1=y(elements(:,1)); y2=y(elements(:,2));y3=y(elements(:,3));

elementH = 1/3*(sqrt((x1-x2).^2+(y1-y2).^2) + ...
    sqrt((x1-x3).^2+(y1-y3).^2) + ...
    sqrt((x2-x3).^2+(y2-y3).^2));

edgeLength = zeros(nx, 1);
for i = 1: nx
    edgeLength(i) = mean(elementH(nodec(i,1:nodec(i,end))));
end
