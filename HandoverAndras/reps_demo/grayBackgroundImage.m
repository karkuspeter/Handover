function grayBackgroundImage(xx, yy, zz, resolution, maxColor)

if nargin < 4
    maxColor = 1;
    resolution = 10;
elseif nargin < 5
    maxColor = 1;
end

map = 1-bsxfun(@times, ones(resolution, 3), linspace(0, maxColor, resolution)');


minz = min(min(zz));
maxz = max(max(zz));

zzDiscrete = round((zz - minz)/(maxz-minz)*resolution);
image(xx, yy,zzDiscrete)

colormap(map)
