function y=do_measurement(funcs, u)
m = length(funcs);
n = size(u,2);
y = zeros(m,n);
for i=1:m
    y(i,:) = funcall(funcs{i}, u);
end


