function fig2_7_pce_num_terms
% Fig 2.7 page 35 (lemaitre & knio spectral methods)

[N,p] = meshgrid(0:15, 0:5);
P1 = nan(size(N));

for i=1:numel(N)
    P1(i) = multiindex_size(N(i), p(i));
end

clf
mesh(p, N, P1, 'EdgeColor', 'k', 'FaceColor', 'none');
logaxis([], 'z');
view([-25, 30])
xlabel('p');
ylabel('N');
zlabel('P+1');
zlim([10^-2, 10^6]);
zTicks = 10.^(0:5);
set(gca,'ZTick',zTicks);
set(gca,'ZTickLabel',num2str(zTicks(:)));
