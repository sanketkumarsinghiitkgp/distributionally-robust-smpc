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
X = Polyhedron('A',H,'b',h);
W = Polyhedron('A',[1 0;-1 0;0 1;0 -1],'b',[0.2;0.2;0.2;0.2]);
U = Polyhedron('A',G,'b',g);
c_theta = 1e4;
num_iter = 31;
%num_iter = 1;
x0 = [2.4;-0.8];
%x0 = [-3.74504;0.0532443];
%x0 = [10;10];
z0 = x0;
e0 = x0-z0;
S0 = Polyhedron('Ae',[1 0; 0 1], 'be',[0;0]);
%S_cur = W;
K = [-0.4 -1.2];

eps_x = 0.1;
eps_u = 0.1;
sigma = 1; %(doubtful)
M = 50;
T = 6;
method_list = ["precompute_eta_2_drsmpc", "scenario_precompute_eta_2_smpc", "robust_mpc"];
%method_list = ["online_drsmpc_2"];
beta_scenario = 0.95;
theta = 1e-5;
w_distribution = "RMPC_gaussian_example";
sigma_multiplier = 0.01;
[x_list,u_list,z_list,e_list,v_list,X_1_ul_list,X_1_ll_list,X_2_ul_list,X_2_ll_list,U_ul_list,U_ll_list] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier,c_theta);
num_realizations = 10;
method_cost = containers.Map('KeyType','char','ValueType','any');
for method=method_list
    method_cost(method)=0;
end
%for i=1:num_realizations
%    [x_list,u_list,z_list,e_list,v_list,X_1_ul_list,X_1_ll_list,X_2_ul_list,X_2_ll_list,U_ul_list,U_ll_list] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier);
%    for method=method_list
%        method_cost(method) = method_cost(method) + find_cost(z_list(method),v_list(method),P,Q,R);
%    end
%end
%for method=method_list
%    method_cost(method) = method_cost(method)/num_realizations;
%end

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
%% Plot Trajectory
figure;
hold on
for method=method_list
    plot_trajectory(x_list(method),method);
end
hold off
legend('Interpreter','none');
xlabel('x_1');
ylabel('x_2');
title('Trajectories');

%% Plot Z_1 limits
figure;
hold on
for method=method_list
    temp = X_1_ul_list(method);
    plot(1:size(X_1_ul_list(method),2)-1,temp(2:end),method_linetype(method),'DisplayName',method_name(method)+" upper limit");
    temp = X_1_ll_list(method);
    plot(1:size(X_1_ll_list(method),2)-1,temp(2:end),method_linetype(method),'DisplayName',method_name(method)+" lower limit");
end
plt = Plot();
plt.BoxDim = [5, 3];
hold off
legend('Interpreter','none');
xlabel('k');
ylabel('limits');
title('z(1) limits');
saveas(gcf,'RMPC_z1.png')
%% Plot Z_2 limits
figure;
hold on
for method=method_list
    temp = X_2_ul_list(method);
    plot(1:size(X_2_ul_list(method),2)-1,temp(2:end),method_linetype(method),'DisplayName',method_name(method)+" upper limit");
    temp = X_2_ll_list(method);
    plot(1:size(X_2_ll_list(method),2)-1,temp(2:end),method_linetype(method),'DisplayName',method_name(method)+" lower limit");
end
hold off
legend('Interpreter','none');
xlabel('k');
ylabel('limits');
title('z(2) limits');
plt = Plot();
plt.BoxDim = [5, 3];

saveas(gcf,'RMPC_z2.png')
%% Plot U limits
figure;
hold on
for method=method_list
    plot(0:size(U_ul_list(method),2)-1,U_ul_list(method),method_linetype(method),'DisplayName',method_name(method)+" upper limit");
    plot(0:size(U_ll_list(method),2)-1,U_ll_list(method),method_linetype(method),'DisplayName',method_name(method)+" lower limit");
end
hold off
legend('Interpreter','none');
xlabel('k');
ylabel('limits');
title('v limits');
plt = Plot();
plt.BoxDim = [5, 3];

saveas(gcf,'RMPC_v.png')


%don't do this. Just say for the most part drsmpc is less conservative.
%num_realizations = 10;
%method_cost = containers.Map('KeyType','char','ValueType','any');
%for method=method_list
%    method_cost(method)=0;
%end
%for i=1:num_realizations
%    [x_list,u_list,z_list,e_list,v_list,X_1_ul_list,X_1_ll_list,X_2_ul_list,X_2_ll_list,U_ul_list,U_ll_list] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier);
%    for method=method_list
%        method_cost(method) = method_cost(method) + find_cost(z_list(method),v_list(method),P,Q,R);
%    end
%end
%for method=method_list
%    method_cost(method) = method_cost(method)/num_realizations;
%end
