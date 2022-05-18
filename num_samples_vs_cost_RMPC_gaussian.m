format longE
rng('default');
cost_list = containers.Map('KeyType','char','ValueType','any');
prob_constraint_satisfied_x_list_list = containers.Map('KeyType','char','ValueType','any');
prob_constraint_satisfied_u_list_list = containers.Map('KeyType','char','ValueType','any');
A = [1 1;0 1];
B = [0;1];
Q = [20 0;0 20];
R = 1;
P = Q;
H = [1 0; -1 0; 0 1; 0 -1];
h = [4;4;1;1];
G = [1;-1];
g = [1;1];
X = Polyhedron('A',H,'b',h);
W = Polyhedron('A',[1 0;-1 0;0 1;0 -1],'b',[0.2;0.2;0.2;0.2]);
U = Polyhedron('A',G,'b',g);
c_theta = 1e4;

x0 = [2.4;-0.8];
z0 = x0;
e0 = x0-z0;
eps_x = 0.1;
eps_u = 0.1;
sigma = 1;
S0 = Polyhedron('Ae',eye(2), 'be',[0;0]);
w_support = W;
K = -dlqr(A,B,P,R,0);
num_iter = 20;
M = 100;
T = 6;
beta_scenario = 0.95;
theta = 1e-5;


%%%CHANGED
w_distribution = "RMPC_gaussian_example";
%w_distribution = "Two_mass_spring_triangle";
%%%%%
constriction = 1;

sigma_multiplier = 0.01;
c_theta=1e4;
method_list = ["online_drsmpc","precompute_eta_2_drsmpc", "scenario_precompute_eta_2_smpc"];

for method=method_list
    prob_constraint_satisfied_u_list_list(method)=[];
    prob_constraint_satisfied_x_list_list(method)=[];
    cost_list(method) = [];
end
cntt =0;
wbar = waitbar(0,"Starting computation");
%theta_list = 2.691681018e-3;%2.691681015e-3 is feasible,2.69168102e-3 isn't, for M = 100;
%theta_list = 1.6e-4:1e-6:1.69e-4; %1.6e-4 feasible, 1.7e-4 isn't
T_list =6:2:12;
for T= T_list
    waitbar(cntt/size(T_list,2),wbar,"Computing for T = "+num2str(T));
    [x_list_cur,u_list_cur,z_list_cur,e_list_cur,v_list_cur,X_1_ul_list_cur,X_1_ll_list_cur,X_2_ul_list_cur,X_2_ll_list_cur,U_ul_list_cur,U_ll_list_cur,prob_constraint_satisfied_x_list_cur,prob_constraint_satisfied_u_list_cur] = compare_approaches(A,B,K,G,g/constriction,H,h/constriction,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier,c_theta);
    for method=method_list
        prob_constraint_satisfied_x_list_list(method) = [ prob_constraint_satisfied_x_list_list(method) prob_constraint_satisfied_x_list_cur(method) ];
        prob_constraint_satisfied_u_list_list(method) = [ prob_constraint_satisfied_u_list_list(method) prob_constraint_satisfied_u_list_cur(method) ];
        cost_list(method) = [ cost_list(method) find_cost(z_list_cur(method),v_list_cur(method),P,Q,R) ];
    end
    cntt=cntt+1;
end
delete(wbar)
method_linetype= containers.Map('KeyType','char','ValueType','any');
method_linetype("scenario_precompute_eta_2_smpc") = '-.';
method_linetype("precompute_eta_2_drsmpc") = '-';
method_linetype("robust_mpc") = ':';
method_linetype("online_drsmpc") = '-';
method_linetype("online_drsmpc_2") = ':';
method_name= containers.Map('KeyType','char','ValueType','any');
method_name("scenario_precompute_eta_2_smpc") = 'Scenario approach';
method_name("precompute_eta_2_drsmpc") = 'DRSMPC';
method_name("robust_mpc") = 'RMPC';
method_name("online_drsmpc") = 'Online DRSMPC';
method_name("online_drsmpc_2") = 'Online DRSMPC 2';

figure;
hold on
for method=method_list
    plot(T_list,cost_list(method),method_linetype(method),'DisplayName',method_name(method));
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'T'; % xlabel
plt.YLabel = 'Cost'; %ylabel
plt.Title = 'Cost vs Prediction Horizon Length (T)'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'cost_vs_T_100_RMPC.png')
