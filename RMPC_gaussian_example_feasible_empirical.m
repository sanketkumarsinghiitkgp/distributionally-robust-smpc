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
feasible_points_list = containers.Map('KeyType','char','ValueType','any');
X = Polyhedron('A',H,'b',h);
W = Polyhedron('A',[1 0;-1 0;0 1;0 -1],'b',[0.2;0.2;0.2;0.2]);
w_support = W;
U = Polyhedron('A',G,'b',g);
num_iter = 1;
%num_iter = 1;
n_x = size(A,2);
n_u = size(B,2);
S0 = Polyhedron('Ae',[1 0; 0 1], 'be',[0;0]);
%S_cur = W;
K = [-0.4 -1.2];

eps_x = 0.1;
eps_u = 0.1;
sigma = 1; %(doubtful)
M = 100;
T = 6;
method_list = ["precompute_eta_2_drsmpc" "robust_mpc" "scenario_precompute_eta_2_smpc"];
%method_list = ["robust_mpc"];
beta_scenario = 0.95;
theta = 1e-5;
w_distribution = "RMPC_example";
sigma_multiplier = 0.01;
M_orig = max(M,T+1);
M_x = 0;
M_u = 0;
p = size(H,1);
q = size(G,1);
M_x = ceil((2/(eps_x/p))*(log(1/beta_scenario)+n_x)); %eps_x/p or eps_x
M_u = ceil((2/(eps_u/q))*(log(1/beta_scenario)+n_u));
M = max(M, M_x);
M = max(M, M_u);
method_name= containers.Map('KeyType','char','ValueType','any');
method_name("scenario_precompute_eta_2_smpc") = 'Scenario approach'
method_name("precompute_eta_2_drsmpc") = 'DRSMPC';
method_name("robust_mpc") = 'RMPC';
%% Generate w samples
w = [];
w_samples = {};
num_samples = max(M,T+1);
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
for method=method_list
    feasible_points = [];
    figure;
    hold on;
    for x0_1 = -6:1:6
        for x0_2 = -3:0.25:3
            x0 = [x0_1;x0_2];
            z0 = x0;
            e0 = x0-z0;
            if (isProblemFeasible(x0,z0,e0,A,B,K,H,h,G,g,P,Q,R,X,W,U,M_x,M_u,M_orig,T,sigma,theta,eps_x,eps_u,w_support,S0,w_samples,method))
                feasible_points = [feasible_points;x0'];
                scatter(x0(1),x0(2),'r');
            end
        end
    end
    feasible_points_list(method) = feasible_points;
    hold off
    xlabel('x0(1)');
    ylabel('x0(2)');
    title("Feasible values of x0 for "+method_name(method));
    saveas(gcf,method+"feasibility_empirical.png");
end