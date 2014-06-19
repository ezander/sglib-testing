clear
clf
state = undamped_spring_init()

m=1; k=1; 
[u, solve_info, state] = undamped_spring_solve(state, [0.3*m; k])
plot(solve_info.t,solve_info.u)
