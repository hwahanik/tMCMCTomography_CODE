function [lats,longs,grid] = readvtx(gridfile)
% READVTX
%
% [lats,longs,grid] = readvtx(gridfile)
%
% This function reads a .vtx file (velocity grid file in the format
% compatible with the FMST codes by Nick Rawlinson).
% 
% Erica Galetti, April 2014
% erica.galetti@ed.ac.uk

fid=fopen(gridfile);

% Read number of nodes in latitude and longitude
info = textscan(fid,'%f%f',1);
Yd=info{1,1}(1);
Xd=info{1,2}(1);

% Read north and west corner of area
info = textscan(fid,'%f%f',1);
north=info{1,1}(1);
west=info{1,2}(1);

% Read node spacing in latitude and longitude
info = textscan(fid,'%f%f',1);
Ygs=info{1,1}(1);
Xgs=info{1,2}(1);

% Calculate south and east corner of area
south=north-(Yd-1).*Ygs;
east=west+(Xd-1).*Xgs;

% Calculate latitudes and longitudes
lats = north:-Ygs:south;       % latitudes
longs = west:Xgs:east;         % longitudes

% Get grid
grid=zeros(Yd,Xd);
for xx = 0:Xd+1
    for yy = 0:Yd+1
        velvals = textscan(fid,'%f',1);
        if xx~=0 && yy~=0 && xx~=Xd+1 && yy~=Yd+1
            grid(yy,xx) = velvals{1,1}(1);
        end
    end
end
fclose(fid);

end

