function [ind, minmax]=peak_find(y)
dy = diff(y);
[ind, dir] = zero_find(dy);
ind = ind+1;
minmax = dir;
