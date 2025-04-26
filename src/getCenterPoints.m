function [rx, ry, r] = getCenterPoints(elements, x, y)
% compute the center of mass, and radius is the max distance to the vertices
x1 = x(elements(:,1));
x2 = x(elements(:,2));
x3 = x(elements(:,3));
y1 = y(elements(:,1));
y2 = y(elements(:,2));
y3 = y(elements(:,3));

% Circumcenter coordinates
rx = 1/3*(x1+x2+x3);
ry = 1/3*(y1+y2+y3);

% Calculate radius
r = max([sqrt((x1 - rx).^2 + (y1 - ry).^2), sqrt((x2 - rx).^2 + (y2 - ry).^2),sqrt((x3 - rx).^2 + (y3 - ry).^2)], [], 2);
