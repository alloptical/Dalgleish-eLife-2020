function [xyz_um_transformed,xyz_um] = etlMagTransform(xyz)

% transforms pixel co-ordinates into µm. ETL magnification changes with
% each z-depth so also need to supply plane information. This is only
% relevant for our Bruker all-optical system.
%
% Input = numROIs * 3 matrix where each row is xyz co-ordinate in pixels
% (xy) and plane number (z; 1 - 4).

fovSizePx           = 512;
fovSizeUm           = 712;
planeSpacing        = 33.3;
etlMagChangePerUm   = 0.000445;
umPerPx             = fovSizeUm/fovSizePx;

% xyz coordinates
xyz_um_transformed = xyz .* [umPerPx umPerPx planeSpacing] - [0 0 planeSpacing];
xyz_um_transformed(:,3) = round(xyz_um_transformed(:,3),1);
xyz_um = xyz_um_transformed;

% adjust for etl magnification: centre, scale, then add back offset
fovSizeUm = fovSizePx * umPerPx;
toCentre = [.5 .5 0] * fovSizeUm;
xyz_um_transformed = xyz_um_transformed - toCentre;
scaleFactors = 1 - (etlMagChangePerUm * xyz_um_transformed(:,3));
xyz_um_transformed(:,1:2) = xyz_um_transformed(:,1:2) .* scaleFactors;
xyz_um_transformed = xyz_um_transformed + toCentre;

end