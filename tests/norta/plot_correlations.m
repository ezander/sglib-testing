function plot_correlations(xi)
% PLOT_CORRELATIONS Plots densities and correlations between RVs given by samples

num_rv=size(xi,1);
for i=1:num_rv
    for j=1:num_rv
        multiplot;
        if i==j
            plot_density(xi(i,:));
        else
            scatter(xi(j,:), xi(i,:), 1);
        end
        axis square
    end
end
