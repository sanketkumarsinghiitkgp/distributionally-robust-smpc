rng('default')  % For reproducibility
H = [0 1; 0 -1];
h = [1.2; 1.2];
G = [1; -1];
g = [6; 6];
A = [ 1 1; 0 1 ];
B = [0.5;1];
P = diag([0.1 1]);
Q = P;
R = 0.1;
K = -dlqr(A,B,P,R,0);
eps_x = 0.4; %decreased from 0.4 to see the effect 4/11 violations
eps_u = 0.1;
T = 10;
w_hist = [];
z_cur = [2;1];
x_cur = [2;1];
e_cur = x_cur-z_cur;
x = [x_cur'];
w = [];
u = [];
v = [];
z = [];
theta = 1e-3;
for i=1:10
    w = [w; w_generator()];
end
for k=1:10
    w_cur = w_generator()';
    w = [w; w_cur'];
    v_cur = find_v_drsmpc(P,Q,R,w,A,B,K,z_cur,e_cur,eps_x,eps_u,H,h,G,g,T,theta);
    u = [u; (K*e_cur+v_cur)'];
    v = [v; v_cur'];
    [new_z, new_e, new_u] = compute_next_state(A,B,z_cur,v_cur,e_cur,w_cur,K); %new_u must be removed, it is useless.
    z_cur = new_z;
    e_cur = new_e;
    x_cur = z_cur+e_cur;
    x = [x; x_cur'];
    z = [z; z_cur'];
end
plot_trajectory(x,'Trajectory');
prob_constraint_satisfied_x = find_prob_constraint_satisfied(H,h,x);
prob_constraint_satisfied_u = find_prob_constraint_satisfied(G,g,u);
fprintf("P(Hx<=h) is %f and 1-eps_x is %f\n",prob_constraint_satisfied_x,1-eps_x);
fprintf("P(Gu<=g) is %f and 1-eps_u is %f\n",prob_constraint_satisfied_u,1-eps_u);
v