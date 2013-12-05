function test_kl_evs

%test2
%plot_lmk1
%plot_lmk2
plot_lmk3


function plot_lmk1
% fig 2.1 page 22 (lemaitre & knio spectral methods)

sig = 1;
N = 10;
b=1;
[pos, els]=create_mesh_1d(0, 1, 100);
[r_i_k, sigma_k] = kl_evp_analytic(pos, els, N, sig, b);

clf
subplot(2,2,1);
plot(pos, binfun(@times, r_i_k, sigma_k));

subplot(2,2,2);
semilogy(1:N, sigma_k.^2,'x');
grid on;

[r_i_k, sigma_k] = kl_evp_numerically(pos, els, N, sig, b);
r_i_k = binfun(@times, r_i_k, sign(r_i_k(52,:)));

subplot(2,2,3);
plot(pos, binfun(@times, r_i_k, sigma_k));

subplot(2,2,4);
semilogy(1:N, sigma_k.^2,'x');
grid on;




function plot_lmk2
% fig 2.2 page 23

clf
sig = 1;
N = 20;
[pos, els]=create_mesh_1d(0, 1, 2);
hold off;
subplot(1,2,1);
for b=linspace(0.1, 1, 10)
    [r_i_k, sigma_k] = kl_evp_analytic(pos, els, N, sig, b);
    loglog(1:N, sigma_k.^2, '-'); hold all;
    xlim([1,N])
end
subplot(1,2,2);
for b=linspace(1, 10, 10)
    [r_i_k, sigma_k] = kl_evp_analytic(pos, els, N, sig, b);
    loglog(1:N, sigma_k.^2, '-'); hold all;
    xlim([1,N])
end
hold off



function plot_lmk3
% fig 2.3 page 25
clf
sig = 1;
N = 6;
[pos, els]=create_mesh_1d(0, 1, 40);
b=1;
[r_i_k, sigma_k] = kl_evp_analytic(pos, els, N, sig, b);
K = @(x,y)(sig^2 * exp(-abs(x-y)/b));

[x,y]=meshgrid(pos);
surf(x,y,K(x,y));
surf(x,y,r_i_k*diag(sigma_k.^2)*r_i_k');
surf(x,y,K(x,y)-r_i_k*diag(sigma_k.^2)*r_i_k');
    


function test2

N=20;
sig = 1;
b = 1;

[pos, els]=create_mesh_1d(0, 1, 100);

[ra_i_k, sigmaa_k] = kl_evp_analytic(pos, els, N, sig, b)
[rn_i_k, sigman_k] = kl_evp_numerically(pos, els, N, sig, b)

format short g
[sigmaa_k; sigman_k; sigmaa_k-sigman_k]'
1;

function [r_i_k, sigma_k] = kl_evp_analytic(pos, els, N, sig, b)
% LeMaitre p. 22, 2.20
K = @(x,y)(sig^2 * exp(-abs(x-y)/b));

M=ceil(N/2)*2;
w1=fzeros_mult( @(w)(1-b*w.*tan(w/2)), 0, M*pi );
w2=fzeros_mult( @(w)(b*w+tan(w/2)), 0, M*pi, 'N', 1000000 );

n=0:(M/2)-1;
assert(all(n*2<=w1/pi & w1/pi<=n*2+1));
assert(all(n*2+1<=w2/pi & w2/pi<=n*2+2));

w=sort([w1, w2]);
w=w(1:N);

sigma_k = sqrt(sig^2 * (2*b) ./ (1 + (w*b).^2));

x = pos';

r_i_k = zeros(size(pos,2), N);

ind1 = 2:2:N;
ind2 = 1:2:N;

d = 1./sqrt(0.5 + 0.5 * sin(w(ind2))./w(ind2));
r_i_k(:,ind2) = binfun(@times, cos((x-0.5)*w(ind2)), d);

d = 1./sqrt(0.5 - 0.5 * sin(w(ind1))./w(ind1));
r_i_k(:,ind1) = binfun(@times, sin((x-0.5)*w(ind1)), d);



function [r_i_k, sigma_k] = kl_evp_numerically(pos, els, N, sig, b)
G_N=mass_matrix(pos,els);
C=covariance_matrix(pos,...
    funcreate(@exponential_covariance, @funarg, @funarg, b, sig))
[r_i_k,sigma_k]=kl_solve_evp( C, G_N, N)

