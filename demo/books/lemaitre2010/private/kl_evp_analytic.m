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



