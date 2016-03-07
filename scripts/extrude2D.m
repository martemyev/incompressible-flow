% This script extrude a given 2D media properties described in a binary file
% to 3D media copying the properties along a Y-axis, i.e. repeating the
% properties in XZ plane

fid = fopen('../build/spe10_bottom.bin','r');
if fid == 0
    disp 'Error opening file'
    return
end
K = fread(fid,[100 100],'single');
fclose(fid);

fid = fopen('../build/spe10_3D_50.bin','w');
if fid == 0
    disp 'Error opening file'
    return
end
for i = 1:size(K, 1)
    for j = 1:50
        fwrite(fid, K(i,:), 'single');
    end
end
fclose(fid);
