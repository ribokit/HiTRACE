function G = get_gaussian(x, center, width )

G = exp( - ( x - center ).^2 / (2*width^2) );