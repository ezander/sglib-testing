clc
clear

n = 30
n=5
sys = 'h'
m = 5
s1= 0.005.*rand(n,1);
sigma1 = sort(s1);
s2 = 0.002.*rand(n,1);
sigma2 = sort(s2);
s3= 0.01.*rand(n,1);
sigma3= sort(s3);
%sigma4= 0.1.*rand(n,1);
%sigma5= 1000000.*rand(n,1);
V_p = {sys, multiindex(m, 1)};


h1_mean = 0.2;
%h1_std  = 0.011;
b1_mean = 0.03;
%b1_std  = 0.005;
b2_mean = 0.1;
%b2_std  = 0.01;

winkel_mean = 0;
winkel_std = 0.1;
E_mean = 200*(10^9);
E_std = 1000000;

for k=1:n
    h1_std= sigma1(k);
    for i=1:n
        b1_std= sigma2(i);
        for j=1:n
            b2_std= sigma3(j);
            
            clear theta_i p_i p_i_alpha
            
            p_i_alpha = [[h1_mean; b1_mean; b2_mean; winkel_mean; E_mean],  ...
                diag([h1_std, b1_std, b2_std, winkel_std, E_std])];
            
            
            %% Monte Carlo
            
            N = 10000;
            theta_i = gpcgerm_sample(V_p, N);
            p_i = gpc_evaluate(p_i_alpha, V_p, theta_i);
            %x_i= beam(p_i);
            x_i(k,i,j,:)= beam(p_i);
        end
    end
end
%figure(1)
%plot(p_i(1,:), p_i(2,:), '.');
plot3(p_i(1,:), p_i(2,:), p_i(3,:), '.');
axis equal
grid on

%figure(110)
%subplot(2,1,1)
%empirical_density(x_i(1,:));

%mean(x_i,2)


nio = x_i< -2.75*10^(-3);
x_nio = x_i(nio);

%subplot(2,1,2)
%empirical_density(x_nio(1,:));
%kde(x_i(1,:));

%%

%x_mean = gpc_integrate(@beam, V_p, 4, 'grid', 'smolyak')

%% GPC

[p, w] = gpc_integrate([], V_p, 4, 'gpc_coeffs', p_i_alpha, 'grid', 'full_tensor');
%figure(111)
%plot(p(1,:), p(2,:), '.');

N=length(w);

wmax = -2.5*10^(-4);
nio(k)=0;
for i=1:N
    delta = beam(p(:,i));
    nio(k) = nio(k) + (delta<wmax) * w(i);
end
%nio(k);

if ~false
    figure(k)
    disp(k)
    plot([1:N],beam(p),'.',[1:N],wmax,'-')
    xlabel('Durchbiegung im letzten Knoten ueber die resultierenden 1024 Stuetzstellen');
    ylabel('Durchbiegeung in Meter');
    str= sprintf('GPC Simulation mit Std(b1)=%f n.i.O=%f',h1_std,mean(nio(:)));
    title(str);
    legend('Durchbiegung im letzten Knoten','NiO Grenzwert');
    
    
    %plot([1:n],sigma,'r*-',[1:n],nio,'b*-')
    %legend('Standardabweichung','Prozentanteil n.i.O.');
    %}
end
