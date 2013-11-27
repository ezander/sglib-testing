n=10000;
x=rand(n,2);

Z=sum(sum(x.^2,2)<=1);

for i=1:10
    if i>4
        fprintf('Favorable outcomes: %d of %d\n', Z, n);
        fprintf('Estimate for 4*Z/n: %g\n', 4*Z/n);
    end
end

% loop over all cells
  % for each cell do whatever