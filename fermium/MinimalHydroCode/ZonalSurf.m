function h=ZonalSurf(tri,x,y,z,zonalval,varargin)
% h=ZonalSurf(tri,x,y,z,zonalval,varargin)

ax = newplot;

h = patch('faces',tri,'vertices',[x(:) y(:) z(:)],...
    'FaceVertexCData',zonalval','FaceColor','flat','CDataMapping','scaled', ...
      'edgecolor',get(ax,'defaultsurfaceedgecolor'));
if nargin>5
    caxis(varargin{:});
end
% if ~ishold, view(3), grid on, end