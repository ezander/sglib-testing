function test_ft_rf

% a=2; b=3.5;
a=0; b=1.5;
func = @(x)(exp(-abs(x)));
[sigma_k, wp_k] = kd_fourier(func, [b-a], 'max_funcs', 15);


xi = linspace(a,b,100);
[x1,x2]=meshgrid(xi);
C_ex=func(x1-x2);

u_k = trig_basis_eval(sigma_k, wp_k, xi)
C = u_k.'*u_k;

surf(x1,x2,C_ex); hold all;
surf(x1,x2,C); hold off;
view(135,0)


%%
clf
func = @(x)(exp(-abs(x)));
func = @(x)(exp(-(x/0.2).^2));
[A_k0, wp_k, x_i] = fourier_series_expand(func, -1, 1, 35, 'symmetry', 'even');
%assert_equals(trig_eval(A_k, wp_k, x_i), func(x_i))

xi = linspace(0,1,100);
[x1,x2]=meshgrid(xi);
C_ex=func(x1-x2);
C = reshape(trig_eval(A_k0, wp_k, x1(:)'-x2(:)')',size(x1));
%C_ex=C;
surf(x1,x2,C_ex); hold all;
surf(x1,x2,C); hold off;
view(135,0)


%%
hold on;
M = 35;
[~, wp_k] = fourier_series_expand(func, 0, 2, M);
r_k = trig_basis_eval([], wp_k, xi)
%r_k = binfun(@minus, r_k, mean(r_k,2));
M2=floor(M/2);
k=(1:M2)';
%s_k = (1:size(r_k,1)).^(-0.5); % rand(M,1)+0.2;

% beta = 1.4;
% s_k = 0.6*[1.5; beta.^-k; beta.^-k];
% 
% alpha=0.9; delta=-0.03;
% s_k = 0.5*[1.45; (k+delta).^-alpha; (k+delta).^-alpha];

a_k=reshape(sqrt(abs(A_k0)),[],1);
a_k=reshape(sqrt(A_k0),[],1);
s_k=[a_k(1); a_k(2:end); a_k(2:end)];

%eee=1.3;
%s_k=[a_k(1); eee*a_k(2:end); eee*a_k(2:end)];
%s_k=s_k/sqrt(s_k(1).^2 + 0.5*sum(abs(s_k(2:end)).^2));


%s_k=sqrt(diag(pinv(r_k)'*C*pinv(r_k)));
norm(r_k'*diag(s_k.^2) * r_k-C)/norm(C)
norm(r_k'*diag(s_k.^2) * r_k-C_ex)/norm(C_ex)
clf
u_k = diag(s_k) * r_k;
surf(x1,x2,u_k.'*u_k); hold all;
surf(x1,x2,C_ex); hold off;
view(135, 0)
%zlim([0,1])
%plot(xi, r_k(1:5,:))



%%
clf
warning('off', 'MATLAB:interp1:ppGriddedInterpolant')

rand_seed(101);
M = 10;
xx = linspace(0,1,M+1); %xx=xx(1:end-1);
yy = randn(M,1); yy(end+1)=yy(1);
pp = interp1(xx,yy,'pchip', 'pp');
func=@(x)(ppval(pp,x));

clf
x = linspace(0,1);
plot(x,func(x))

hold all;

[A_k, wp_k, x_i] = fourier_series_expand(func, 0, 1, 10);
%plot(x_i, func(x_i), x_i, trig_eval(A_k, wp_k, x_i))
plot(x,trig_eval(A_k, wp_k, x))
plot(xx,yy,'xg')


%s_j_k = 


