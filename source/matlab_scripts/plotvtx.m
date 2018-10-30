function plotvtx(gridfile,int)
% PLOTVTX
%
% plotvtx(gridfile,int)
%
% This function plots a .vtx file (velocity grid file in the format
% compatible with the FMST codes by Nick Rawlinson).
%
% Interpolation onto a new grid can be done through parameter 'int' (i.e.
% if int=32, a new grid will be defined with nodes every 1/32 of a degree
% in latitude and longitude).
%
% Erica Galetti, April 2014
% erica.galetti@ed.ac.uk

% Read vtx file
[lats,longs,grid] = readvtx(gridfile);

north=lats(1);
south=lats(end);
west=longs(1);
east=longs(end);


% Calculate interpolated grid
if nargin == 2
    Y=repmat(lats.',1,length(longs));
    X=repmat(longs,length(lats),1);
    newlats=north:-1/int:south;
    newlongs=west:1/int:east;
    Yi=repmat(newlats.',1,length(newlongs));
    Xi=repmat(newlongs,length(newlats),1);
end

% Plot
figure; imagesc(longs,lats,grid)
set(gca,'ydir','normal')
axis equal
xlabel('Longitude')
ylabel('Latitude')
colorbar

% Plot interpolated grid
if nargin == 2
    grid_int=interp2(X,Y,grid,Xi,Yi,'linear');
    figure; imagesc(newlongs,newlats,grid_int)
    set(gca,'ydir','normal')
    axis equal
    xlabel('Longitude')
    ylabel('Latitude')
    colorbar
end

end

