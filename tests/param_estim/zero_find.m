function [ind, dir]=zero_find(y)
ind1 = find(sign(y(1:end-1)) ~= sign(y(2:end)));
ind = ind1 + y(ind1) ./ (y(ind1) - y(ind1+1));
dir = sign(y(ind1+1)-y(ind1));
