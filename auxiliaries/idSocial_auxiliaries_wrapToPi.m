function anglePiPi = idSocial_auxiliaries_wrapToPi(angle_rad)

anglePiPi = angle_rad - 2*pi*floor( (angle_rad+pi)/(2*pi));