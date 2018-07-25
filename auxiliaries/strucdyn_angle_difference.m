function adif=strucdyn_angle_difference(before,after)

a=after; b=before; 
adif=atan2(sin(a-b), cos(a-b));