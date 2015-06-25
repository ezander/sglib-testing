function show_smolyak_rules

N = 20;
d = 3;
animate_3d = false;

[x, w] = smolyak_grid(d, N, @pseudo_int_rule)

clf
ipos = (w>0);
ineg = (w<0);
if d==2
    plot(x(1,ipos), x(2,ipos), 'bo', ...
        x(1,ineg), x(2,ineg), 'ro')
else
    plot3(x(1,ipos), x(2,ipos), x(3,ipos), 'bo', ...
        x(1,ineg), x(2,ineg), x(3,ineg), 'ro')
    view(230, 5)
end
axis square
grid on

if animate_3d && d==3
    for angle = 200:1:(200+360*10)
        view(angle, 5);
        drawnow;
    end
end
    

function [x,w] = pseudo_int_rule(n)
x = repmat(n,1,n);
w = ones(n,1);

