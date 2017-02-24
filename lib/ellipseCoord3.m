function [V]=ellipseCoord3(e,t)

V=[e.radii(1).*cos(t(:)) e.radii(2).*sin(t(:)) zeros(numel(t),1)];
V=(e.axes*V')';
V=V+e.centre(ones(numel(t),1),:);