n=10000;

hold off;

x=linspace(0,1,n);
dx=x(2)-x(1);
[r_i_k, sigma_k]=kl_solve_1d_exp(x, 1, 0.3, 5);
plot(x,r_i_k(:,1:3 ))
abs(r_i_k'*binfun(@times,[dx/2; repmat(dx,n-2,1); dx/2], r_i_k))
hold on;

DX=10;
x=linspace(DX*1.3,DX*2.3,n);
dx=x(2)-x(1);
[r_i_k, sigma_k]=kl_solve_1d_exp(x, 1, 0.3*DX, 7);
plot(x,r_i_k(:,1:3 ))

abs(r_i_k'*binfun(@times,[dx/2; repmat(dx,n-2,1); dx/2], r_i_k))


DX=10;
x=linspace(DX*0,DX*1,n);
dx=x(2)-x(1);
[r_i_k, sigma_k]=kl_solve_1d_exp(x, 1, 0.3*DX, 7);
plot(x,r_i_k(:,1:3 ))

abs(r_i_k'*binfun(@times,[dx/2; repmat(dx,n-2,1); dx/2], r_i_k))

