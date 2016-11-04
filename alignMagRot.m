function M = alignMagRot(magaxis,rotaxis)
% aligns magaxis along z, rotaxis along x

rot_calc = vrrotvec(rotaxis,[1,0,0]);
rot_mat = rotaxi2mat(rot_calc(1:3),rot_calc(4));

wrong_rot_axis = rot_mat*[1 0 0]';
angle = atan2(norm(cross(magaxis,wrong_rot_axis)),dot(magaxis,wrong_rot_axis))
rot_mat2 = rotaxi2mat([1,0,0],angle);

M = rot_mat2*rot_mat;