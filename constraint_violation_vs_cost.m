format longE
rng('default');
cost_list = containers.Map('KeyType','char','ValueType','any');
prob_constraint_satisfied_x_list_list = containers.Map('KeyType','char','ValueType','any');
prob_constraint_satisfied_u_list_list = containers.Map('KeyType','char','ValueType','any');
   

spring_k = 1;
m_1 = 0.5;
m_2 = 2;
A = [1, 0, 0.1, 0;0, 1, 0, 0.1;-spring_k/m_1, 0.1*spring_k/m_1, 1, 0;spring_k/m_2, -0.1*spring_k/m_2, 0, 1];
B = [0;0;0.1/m_1;0];
B_w = [1;0.5;0.3;0.4];
x0 =[0.2;1;-0.1;0.1];
z0 = x0;
e0 = x0-z0;
Q = 5*eye(4);
R = 1;
P = 5*eye(4);
eps_x = 0.2;
eps_u = 0.2;
sigma = 1;
H = [0,0,1,0;0,0,-1,0;0,0,0,1;0,0,0,-1];
h = [0.38;0.38;0.38;0.38];
G = [1;-1];
g = [1.6;1.6];
X = Polyhedron('A',H,'b',h);
W = Polyhedron('A',[0 0 0 0],'b',[1]);
S0 = Polyhedron('Ae',eye(4), 'be',[0;0;0;0]);
w_support = W;
U = Polyhedron('A',G,'b',g);
K = -dlqr(A,B,P,R,0);
num_iter = 20;
M = 100;
T = 6;
beta_scenario = 0.95;

%%%CHANGED
w_distribution = "Two_mass_spring";
%w_distribution = "Two_mass_spring_triangle";
%%%%%
constriction = 1.0085;

sigma_multiplier = 0.0001;
c_theta=1e4;
method_list = ["online_drsmpc"];

for method=method_list
    prob_constraint_satisfied_u_list_list(method)=[];
    prob_constraint_satisfied_x_list_list(method)=[];
    cost_list(method) = [];
end
cntt =0;
wbar = waitbar(0,"Starting computation");
%theta_list = 2.691681018e-3;%2.691681015e-3 is feasible,2.69168102e-3 isn't, for M = 100;
%theta_list = 1.6e-4:1e-6:1.69e-4; %1.6e-4 feasible, 1.7e-4 isn't
theta = 6.66e-4; %6.66e-4 is feasible.
constriction_list = 1.033333:0.0000001:1.0333339;
for constriction=constriction_list
    waitbar(cntt/size(constriction_list,2),wbar,"Computing for constriction = "+num2str(constriction));
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
    [cost_list(method),sortIdx] = sort(cost_list(method),'ascend');
    yy = prob_constraint_satisfied_u_list_list(method);
    yy = yy(sortIdx);
    prob_constraint_satisfied_u_list_list(method)=yy;
    plot(cost_list(method),prob_constraint_satisfied_u_list_list(method),'DisplayName',method_name(method));
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraint Tightening for Velocity of Mass 1'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_prob_constraint_satisfied_u.png')



figure;
hold on
for method=method_list
    [cost_list(method),sortIdx] = sort(cost_list(method),'ascend');
    yy = prob_constraint_satisfied_x_list_list(method);
    yy = yy(sortIdx);
    prob_constraint_satisfied_x_list_list(method)=yy;
    plot(cost_list(method),prob_constraint_satisfied_x_list_list(method),'DisplayName',method_name(method));
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Cost'; % xlabel
plt.YLabel = 'Frequency of constraint violation (fraction)'; %ylabel
plt.Title = 'Violation frequency vs cost'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_prob_constraint_satisfied_x.png')