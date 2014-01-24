function demo_density_estimation
% DEMO_DENSITY_ESTIMATION Demo showing density estimation functions.
%  Show kernel density plots for different distributions and different
%  numbers of samples.


%   Elmar Zander
%   Copyright 2013, Inst. of Scientific Computing, TU Braunschweig
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version. 
%   See the GNU General Public License for more details. You should have
%   received a copy of the GNU General Public License along with this
%   program.  If not, see <http://www.gnu.org/licenses/>.

clf
dists = {
    gendist_create('beta', {0.5, 0.7});
    gendist_create('lognormal', {1.5, 0.5});
    gendist_create('uniform', {0.5, 2.5});
    gendist_create('normal', {1.5, 0.5})};
%number_of_samples = [100, 1000, 10000, 100000];
number_of_samples = [100, 300, 1000, 3000, 10000];

n=1;
for dist_ = dists'
    dist = dist_{1};
    for N = number_of_samples
        subplot(length(dists), length(number_of_samples), n);
        hold all;
        % Generate the samples for the density plot
        x_samples = gendist_sample([1, N], dist);
        kernel_density(x_samples);
        % Generate values for plotting the PDF (somewhat enlarged range of the
        % range given by the samples)
        x_vals = linspace(min(x_samples), max(x_samples));
        x_vals = (x_vals - mean(x_vals)) * 1.2 + mean(x_vals);
        plot(x_vals, gendist_pdf(x_vals, dist));
        grid on;
        hold off;
        drawnow;
        n=n+1;
    end
end