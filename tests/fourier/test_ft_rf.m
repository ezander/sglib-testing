clear
clf

a=2; b=3.5;
func = @(x)(exp(-abs(x)));
TB = kd_fourier(func, [b-a], 'max_funcs', 50);


xi = linspace(a,b,100);
[x1,x2]=meshgrid(xi);
C_ex=func(x1-x2);

u_k = trig_basis_eval(TB, xi);
C = u_k.'*u_k;

surf(x1,x2,C_ex); hold all;
surf(x1,x2,C); hold off;
view(135,0)


%%
clf
func = @(x)(exp(-(x/0.2).^2));
[A_k0, wp_k] = fourier_series_expand(func, -1, 1, 35, 'symmetry', 'even');

xi = linspace(0,1,100);
[x1,x2]=meshgrid(xi);
C_ex=func(x1-x2);
C = reshape(trig_eval(A_k0, wp_k, x1(:)'-x2(:)')',size(x1));
surf(x1,x2,C_ex); hold all;
surf(x1,x2,C); hold off;
view(135,0)


%%
hold on;
M = 35;
[~, TB] = fourier_series_expand(func, 0, 2, M);
r_k = trig_basis_eval(TB, xi);

a_k=reshape(sqrt(A_k0),[],1);
s_k=[a_k(1); a_k(2:end); a_k(2:end)];

norm(r_k'*diag(s_k.^2) * r_k-C)/norm(C)
norm(r_k'*diag(s_k.^2) * r_k-C_ex)/norm(C_ex)
clf
u_k = diag(s_k) * r_k;
surf(x1,x2,u_k.'*u_k); hold all;
surf(x1,x2,C_ex); hold off;
view(135, 0)



%%
clf
warning('off', 'MATLAB:interp1:ppGriddedInterpolant')

rand_seed(101);
M = 10;
xx = linspace(0,1,M+1);
yy = randn(M,1); yy(end+1)=yy(1);
pp = interp1(xx,yy,'pchip', 'pp');
func=@(x)(ppval(pp,x));

clf
x = linspace(0,1);
plot(x,func(x))

hold all;

[A_k, wp_k, x_i] = fourier_series_expand(func, 0, 1, 10);
plot(x,trig_eval(A_k, wp_k, x))
plot(xx,yy,'xg')
