function I=confidence_interval(x, stddev, alpha)

dx1 = stddev * normal_invcdf(alpha/2);
dx2 = stddev * normal_invcdf(1-alpha/2);
I = [x+dx1, x+dx2];
