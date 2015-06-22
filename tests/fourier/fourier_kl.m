function fourier_kl
l=1.3;
%func = @(s)exp(-abs(s/l))
func = @(s)exp(-s/l)


s=2;
N=21;
x=linspace(-s,s,N+1);
x=x(1:end-1);
y = func(x);
plot(x, y)

a=fft(func(x))
w_N = exp(-2*pi*i/N)
k = (0:N-1)';
W = w_N.^( -k*k') / N;
b=a * W;
norm(b-y)

b=real(b);
hold all; plot(x,b,'r-.'); hold off

%%
m=1;
W3 = [W(1,:); W(2:2+m,:)+W(N:-1:N-m,:)];
norm(imag(W3))
b2=real(a(1:m+2))*W3;


%%
b2 = a * exp(2*pi*i*k*k'/N)/N
plot(x,y,x,b,x,b2)
b2 = a * (cos(2*pi*k*k'/N)/N + i*sin(2*pi*k*k'/N)/N)
b2 = real(a) * (cos(2*pi*k*k'/N)/N) - imag(a)*sin(2*pi*k*k'/N)/N
ar = real(a); ai = imag(a);
b2 = ar * (cos(2*pi*k*(x-s)/(2*s))/N) - ai*sin(2*pi*k*(x-s)/(2*s))/N
C = cos(2*pi*k*(x-s)/(2*s))/N;
S = sin(2*pi*k*(x-s)/(2*s))/N;
b2 = ar * C - ai * S
plot(x,y,x,b,x,b2)

%%
m = ceil(N/2);
ai2 = 2*ai(2:m)
S2 = S(2:m,:);
norm(ai*S-ai2*S2)

ar2 = [ar(1) 2*ar(2:m)]
C2 = C(1:m,:)
norm(ar*C-ar2*C2)

b2 = ar2 * C2 - ai2 * S2



plot(x,y,x,b,x,b2)
%%
A_k=[];
w_k=[];
p_k=[];


b1 = ar2 * C2 -ai2 * S2
A_k=[A_k; [ar(1) 2*ar(2:m) ar(m+1)]/N];
w_k=[w_k; 2*pi*(0:m)'/(2*s)];
p_k=[p_k; 2*pi*(0:m)'*(-s)/(2*s)+pi/2];

A_k=[A_k, -2*ai(2:m)/N];
w_k=[w_k; 2*pi*(1:m-1)'/(2*s)];
p_k=[p_k; 2*pi*(1:m-1)'*(-s)/(2*s)];

b2 = sin_eval(A_k, [w_k, p_k], x)
b2-b1

[A_k, wp_k, x_i]=myfft(func, -s, s, N);
b3 = sin_eval(A_k, [w_k, p_k], x)
b3-b1


%%
hold off
k=[0;1;2;N-1;N]
plot(x,cos(2*pi*k*(x-s)/(2*s)))

%%
clf
plot(x,real(a(2)*W(2,:)+a(N)*W(N,:)), x,real(2*a(2)*W(2,:)) ); hold all;
plot(x,imag(a(2)*W(2,:)+a(N)*W(N,:)), x,imag(2*a(2)*W(2,:)) ); hold off;

a(2)
a(N)

%%
m=10
a(2:2+m)-conj(a(N:-1:N-m))

%%
%a(1)
%a(2:10)
%%a(21:-1:11)

m=10
aa=a(2:2+m)+a(N:-1:N-m)
k2 = (0:m)';
W2 = cos(2*pi*k2*(x-s)/(2*s))/N
y./(aa*W2)
plot(x,W2(2,:))

%%
plot(x,W(2,:)+W(N+2-2,:),x,2*cos(2*pi*(x-s)/(2*s))/N)
plot(x,W(3,:)+W(N+2-3,:),x,2*cos(2*pi*2*(x-s)/(2*s))/N)
plot(x,W(4,:)+W(N+2-4,:),x,2*cos(2*pi*3*(x-s)/(2*s))/N)
%plot(x,W(:,3)+W(:,N+2-3))


om_n = 2*pi*(0:N-1)/N
a


%a-conj(fliplr(a))



