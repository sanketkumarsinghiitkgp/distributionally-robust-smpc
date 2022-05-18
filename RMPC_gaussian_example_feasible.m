rng('default')  % For reproducibility
A = [1 1;0 1];
B = [0;1];
Q = [20 0;0 20];
R = 1;
P = Q;
H = [1 0; -1 0; 0 1; 0 -1];
h = [4;4;1;1];
G = [1;-1];
g = [1;1];
w_support = Polyhedron('A',[1 0;-1 0;0 1;0 -1],'b',[0.2;0.2;0.2;0.2]);
W = w_support;
X = Polyhedron('A',H,'b',h);
U = Polyhedron('A',G,'b',g);
S0 = Polyhedron('Ae',[1 0; 0 1], 'be',[0;0]);
x0 = [2.4;-0.8];
z0 = x0;
e0 = x0-z0;
n_x = size(A,2);
n_u = size(B,2);
%S_cur = W;
K = [-0.4 -1.2];
num_iter = 15;
beta_scenario = 0.95
sigma_multiplier = 0.01;
eps_x = 0.1;
eps_u = 0.1;
sigma = 1; %(doubtful)
M = 100;
T = 6;
theta = 1e-5;
M_orig = max(M,T+1);
M_x = 0;
M_u = 0;
p = size(H,1);
q = size(G,1);
M_x = ceil((2/(eps_x/p))*(log(1/beta_scenario)+n_x)); %eps_x/p or eps_x
M_u = ceil((2/(eps_u/q))*(log(1/beta_scenario)+n_u));
M = max(M, M_x);
M = max(M, M_u);

%% Generate w samples
w = [];
w_samples = {};
num_samples = max(M,T+1);
%if (1 == 1)
mu = [0;0];
Sigma = sigma_multiplier*[0.01 0; 0 1];
AA = w_support.A;
bb = w_support.b;
for i=1:num_samples
    w_cur = [];
    for j=1:T+1
        w_cur = [w_cur; w_generator("rmvnrnd",mu,Sigma,1,AA,bb)];
    end
    w_samples{i} = w_cur;
end
for i = 1:M+num_iter
    w = [w; w_generator("rmvnrnd",mu,Sigma,1,AA,bb)];
end
eta_x_history_precompute_eta_2 = find_eta_list_2(A, B, K, H, T+1, M_orig, e0, w_samples,sigma,theta,eps_x,w_support);
eta_u_history_precompute_eta_2 = find_eta_list_2(A, B, K, G, T, M_orig,  e0, w_samples,sigma,theta,eps_u,w_support);
eta_x_history_scenario_2 = find_eta_list_scenario_2(A, B, K, H, T+1, M_x, w_samples);
eta_u_history_scenario_2 = find_eta_list_scenario_2(A, B, K, G, T, M_u, w_samples);

Acl = A+B*K;
S_cur = S0;
S = [S_cur];
for i=1:num_iter+T
    S_cur = Acl*S_cur+W;
    S = [S, S_cur];
end
X_1_ul = [];
X_1_ll = [];
X_2_ul = [];
X_2_ll = [];
U_ul = [];
U_ll = [];
N = T;
for i=1:N+1
    Z = X-S(i);
    V = U-K*S(i);
    X_1_ul = [X_1_ul Z.b(5)]; %fishy method, relies on MPT preserving the order of inequalities.
    X_1_ll = [X_1_ll -Z.b(6)];
    X_2_ul = [X_2_ul Z.b(7)];
    X_2_ll = [X_2_ll -Z.b(8)];
    U_ul = [U_ul V.b(3)];
    U_ll = [U_ll -V.b(4)];
end
h_eta_x_robust = [];
for i=2:N+1
    temp = [X_1_ul(i), -X_1_ll(i), X_2_ul(i), -X_2_ll(i)];
    h_eta_x_robust= [h_eta_x_robust;temp'];
end
g_eta_u_robust = [];
for i=1:N
    temp = [U_ul(i), -U_ll(i)];
    g_eta_u_robust= [g_eta_u_robust;temp'];
end


HH_diag = {};

for i=1:N-1
    HH_diag = {HH_diag{:},H};
end
HH_diag = {HH_diag{:},H}; % since H_bar = H
HH = blkdiag(HH_diag{:});

GG_diag = {};
for i=1:N
    GG_diag = {GG_diag{:},G};
end
GG = blkdiag(GG_diag{:});
TT = [];
A_power = 1;
for i=1:N
    A_power = A_power*A;
    
    TT = [TT;A_power];
end
SS = [];
zero_block = zeros(size(B));
for i = 1:N
    cur = [];
    for j=1:N
        temp = zero_block;
        if (j<i)
            temp = TT(i-j)*B;
        elseif (j==i)
            temp = B;
        end
        cur = [cur temp];
    end
    SS = [SS;cur];
end
h_eta_x_precompute_eta_2 = [];
for i=2:N+1
    h_eta_x_precompute_eta_2 = [h_eta_x_precompute_eta_2;h-eta_x_history_precompute_eta_2{i}'];
end
g_eta_u_precompute_eta_2 = [];
for i=1:N
    g_eta_u_precompute_eta_2 = [g_eta_u_precompute_eta_2;g-eta_u_history_precompute_eta_2{i}'];
end

h_eta_x_scenario_2 = [];
for i=2:N+1
    h_eta_x_scenario_2 = [h_eta_x_scenario_2;h-eta_x_history_scenario_2{i}'];
end
g_eta_u_scenario_2 = [];
for i=1:N
    g_eta_u_scenario_2 = [g_eta_u_scenario_2;g-eta_u_history_scenario_2{i}'];
end

Polytope_A = [HH*SS HH*TT;GG zeros(size(GG,1),size(H,2));zeros(size(H,1),size(GG,2)) H];
Polytope_b_precompute_eta_2 = [h_eta_x_precompute_eta_2;g_eta_u_precompute_eta_2;h-eta_x_history_precompute_eta_2{1}'];
Polytope_b_scenario_2 = [h_eta_x_scenario_2;g_eta_u_scenario_2;h-eta_x_history_scenario_2{1}'];
temp0 = [X_1_ul(1), -X_1_ll(1), X_2_ul(1), -X_2_ll(1)];
Polytope_b_robust = [h_eta_x_robust;g_eta_u_robust;temp0'];

P_precompute_eta_2 = polytope(Polytope_A,Polytope_b_precompute_eta_2);
dims = size(Polytope_A,2);
Q_precompute_eta_2 = projection(P_precompute_eta_2,dims-1:dims);

P_scenario_2 = polytope(Polytope_A,Polytope_b_scenario_2);
dims = size(Polytope_A,2);
Q_scenario_2 = projection(P_scenario_2,dims-1:dims);

P_robust = polytope(Polytope_A,Polytope_b_robust);
dims = size(Polytope_A,2);
Q_robust = projection(P_robust,dims-1:dims);

hold on
plot(Q_precompute_eta_2,'r');
plot(Q_scenario_2,'b');
plot(Q_robust,'y');
hold off
xlabel('x0(1)');
ylabel('x0(2)');
title('Feasible values of x0');
legend({'DRSMPC','Scenario Approach','RMPC'})
alpha(0.5);
saveas(gcf,'RMPC_feasible.png')

figure;
Polytope_A = [HH*SS HH*TT;GG zeros(size(GG,1),size(H,2))];
Polytope_b_precompute_eta_2 = [h_eta_x_precompute_eta_2;g_eta_u_precompute_eta_2];
Polytope_b_scenario_2 = [h_eta_x_scenario_2;g_eta_u_scenario_2];
temp0 = [X_1_ul(1), -X_1_ll(1), X_2_ul(1), -X_2_ll(1)];
Polytope_b_robust = [h_eta_x_robust;g_eta_u_robust];

P_precompute_eta_2 = polytope(Polytope_A,Polytope_b_precompute_eta_2);
dims = size(Polytope_A,2);
Q_precompute_eta_2 = projection(P_precompute_eta_2,dims-1:dims);

P_scenario_2 = polytope(Polytope_A,Polytope_b_scenario_2);
dims = size(Polytope_A,2);
Q_scenario_2 = projection(P_scenario_2,dims-1:dims);

P_robust = polytope(Polytope_A,Polytope_b_robust);
dims = size(Polytope_A,2);
Q_robust = projection(P_robust,dims-1:dims);
hold on
plot(Q_precompute_eta_2,'r');
plot(Q_scenario_2,'b');
plot(Q_robust,'y');
hold off
xlabel('x0(1)');
ylabel('x0(2)');
title('Feasible values of x0');
legend({'DRSMPC','Scenario Approach','RMPC'})
alpha(0.5);
saveas(gcf,'RMPC_feasible.png')