clear
m=4; k=1;
clf
hold all
l = {};
T0 = 20;
for d=[1, 2, 3, 4]
    D=d*d-k*m;
    state = undamped_spring_init('solver', 'numerical', 'd', d, 'T', T0);
    [u, solve_info, state] = undamped_spring_solve(state, [m; k]-1.5);
    plot(solve_info.t,solve_info.ut(1,:))
    l{end+1} = strvarexpand('d=$d$, D=$D$');
    state = undamped_spring_init('T', T0, 'solver', 'direct', 'd', d);
    u2 = undamped_spring_solve(state, [m; k]-1.5);
    [m, k, d, k/m - (d/m)^2]
end
legend(l{:})
grid on
for T=linspace(0.1,T0,10)
    for d=[] %[0, 1, 2, 3]
        state = undamped_spring_init('T', T, 'solver', 'numerical', 'd', d);
        u1 = undamped_spring_solve(state, [m; k]);
        line([T], [u1(1)], 'Marker', 'x')
        state = undamped_spring_init('T', T, 'solver', 'direct', 'd', d);
        u2 = undamped_spring_solve(state, [m; k]);
        line([T], [u2(1)], 'Marker', 'o')
    end
    hold off
    drawnow
end
