function test_romberg

clc

f=@sin;
for n=0:3
    for m=0:4
        %[n,m,log(abs(romberg2(f, n, m)-(cos(0)-cos(1))))]
        [n,m,log(abs(romberg(f, n, m)-romberg2(f, n, m)))]
    end
end
fprintf('\n\n\n');
for m=0:4
    [x,w]=romberg_rule(0, m);
    [x,inda,indb]=unique(x);
    w=accumarray(indb', w);
    fprintf('\n\n\n%d %d \n', m, length(w));
    rats(x)
    rats(w')
end
1

function I=romberg2(f,n,m)
[x,w]=romberg_rule(n,m);
I=f(x)*w;

function [x,w]=romberg_rule(n,m)
x = 0;
w = 1;
assert(m>=0 && n>=0)
if m==0
    if n==0
        x=[1,0];
        w=[0.5; 0.5];
    else
        hn = 0.5^n;
        k = 1:2^(n-1);

        [x,w] = romberg_rule(n-1,0);
        x = [x,  hn*(2*k-1)];
        w = [0.5*w;  hn*ones(size(k))'];
    end
else
    [x1,w1] = romberg_rule(n+1, m-1);
    [x2,w2] = romberg_rule(n, m-1);
    x = [x1, x2];
    w = [4^m * w1; -w2]/(4^m-1);
end
[x,inda,indb]=unique(x);
w=accumarray(indb', w);


function I=romberg(f,n,m)
assert(m>=0 && n>=0)
if m==0
    if n==0
        I = 0.5 * (f(1)+f(0));
    else
        hn = 0.5^n;
        k = 1:2^(n-1);
        I = 0.5 * romberg(f,n-1,0) + sum(hn * f((2*k-1)*hn));
    end
else
    I1 = romberg(f, n+1, m-1);
    I2 = romberg(f, n, m-1);
%    I = I1 + (I1 - I2)/(4^m-1);
    I = (4^m*I1 - I2)/(4^m-1);
end
