N=10;
x = rand(2,N);

[y,dy,H]=rosenbrock(x);

dx = 1e-8;
[y1, dy1] = rosenbrock(x + repmat([dx;0], 1, N));
[y2, dy2] = rosenbrock(x + repmat([0;dx], 1, N));

clc
[(y1-y)/dx; dy(1,:)]
[(y2-y)/dx; dy(2,:)]



[reshape(H(1,:,:),2,[]); (dy1 - dy)/dx]
[reshape(H(2,:,:),2,[]); (dy2 - dy)/dx]
