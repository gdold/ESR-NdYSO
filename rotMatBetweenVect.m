function init_rot_mat = rotMatBetweenVect(from,to)

init_rot_calc = vrrotvec(from,to);
init_rot_mat = rotaxi2mat(init_rot_calc(1:3),init_rot_calc(4));