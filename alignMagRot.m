function M = alignMagRot(magaxis)
% aligns magaxis along z, rotaxis along x

rot_calc = vrrotvec(magaxis,[0,0,1]);
M = rotaxi2mat(rot_calc(1:3),rot_calc(4));