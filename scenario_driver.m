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
eps_x = 0.4;
eps_u = 0.1;
T = 10;
sigma = 1;
w_support = Polyhedron('A',[0 0],'b',[0]);
w_hist = [];
z_cur = [2;0.5];
x_cur = [2;0.5];
e_cur = x_cur-z_cur;
x = [x_cur'];
w = [];
u = [];
v = [];
z = [];
theta = 1e-3;
%verify if 100 samples are sufficient
beta = 0.95;
n_x = size(x_cur);
n_x = n_x(1);
n_u = size(G);
n_u = n_u(2);
M_x = (2/eps_x)*(log(1/beta)+n_x);
M_u = (2/eps_u)*(log(1/beta)+n_u);
M = max(M_x, M_u);
for i=1:max(M,T+1)
    w = [w; w_generator()];
end
%w size determines the number of e samples generated at each instance.

%ugly hack please replace
[eta_x_history, ~] = find_eta_list_scenario(A,B,K,eps_x,eps_u,H,G,e_cur,T,w(1:M_x,:),w_support,sigma,theta);
[~, eta_u_history] = find_eta_list_scenario(A,B,K,eps_x,eps_u,H,G,e_cur,T,w(1:M_u,:),w_support,sigma,theta);

X_2_ul = [];
X_2_ll = [];
U_ul = [];
U_ll = [];
for i=1:T+1
    temp = h-eta_x_history{i}';
    X_2_ul = [X_2_ul temp(1)];
    X_2_ll = [X_2_ll -temp(2)];
end
for i=1:T
    temp = g-eta_u_history{i}';
    U_ul = [U_ul temp(1)];
    U_ll = [U_ll -temp(2)];
end
figure;
hold on
plot(X_2_ul,'DisplayName','z_2 upper limits');
plot(X_2_ll,'DisplayName','z_2 lower limits');
plot(U_ul,'DisplayName','v upper limits');
plot(U_ll,'DisplayName','v lower limits');
legend
hold off
for k=1:10
    w_cur = w_generator()';
    w = [w; w_cur'];
    v_cur = find_v_drsmpc_precompute_eta(P,Q,R,A,B,H,h,G,g,x_cur,T,eta_x_history,eta_u_history); %we set z_k = e_k for the purpose of this computation.
    %%
    %debug
    %v_cur = 0;
    %%
    u = [u; (K*e_cur+v_cur)'];
    v = [v; v_cur'];
    [new_z, new_e, new_u] = compute_next_state(A,B,z_cur,v_cur,e_cur,w_cur,K); %new_u must be removed, it is useless.
    z_cur = new_z;
    e_cur = new_e;
    x_cur = z_cur+e_cur;
    x = [x; x_cur'];
    z = [z; z_cur'];
end
figure;
plot_trajectory(x,num2str(theta));
prob_constraint_satisfied_x = find_prob_constraint_satisfied(H,h,x);
prob_constraint_satisfied_u = find_prob_constraint_satisfied(G,g,u);
fprintf("P(Hx<=h) is %f and 1-eps_x is %f\n",prob_constraint_satisfied_x,1-eps_x);
fprintf("P(Gu<=g) is %f and 1-eps_u is %f\n",prob_constraint_satisfied_u,1-eps_u);
v