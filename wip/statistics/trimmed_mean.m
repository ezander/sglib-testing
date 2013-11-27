function m=trimmed_mean(x,p)

N=length(x);
if ~issorted(x)
    x=sort(x);
end
p = min(1,max(0, p));
i1= min(N,max(1,round(0.5*p*N)));
i2= min(N,max(i1,N-round(0.5*p*N)));
m = mean(x(i1:i2));
