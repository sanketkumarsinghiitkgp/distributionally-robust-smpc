rng('default')  % For reproducibility
A = [1 1;0 1];
B = [0.5;1];
Q = [0.1 0;0 1];
R = 0.1;
P = Q;
H = [1 0; -1 0; 0 1; 0 -1];
h = [1e20;1e20;1.2;1.2];
G = [1;-1];
g = [6;6];
X = Polyhedron('A',H,'b',h);
W = Polyhedron('A',[0 0],'b',1);
U = Polyhedron('A',G,'b',g);
num_iter = 31;
x0 = [2;0.5];
z0 = x0;
e0 = x0-z0;
S0 = Polyhedron('Ae',[1 0; 0 1], 'be',[0;0]);
%S_cur = W;
K = -dlqr(A,B,P,R,0);
eps_x = 0.4;
eps_u = 0.1;
sigma = 1; %(doubtful)
M = 100;
T = 6;
method_list = ["precompute_eta_2_drsmpc" "scenario_precompute_eta_2_smpc"];
beta_scenario = 0.95;
theta = 1e-5;
sigma_multiplier = 0.01
w_distribution = "SMPC_reachability_truncated_example";
[x_list,u_list,z_list,e_list,v_list,X_1_ul_list,X_1_ll_list,X_2_ul_list,X_2_ll_list,U_ul_list,U_ll_list] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier);
%% Plot Trajectory

method_linetype= containers.Map('KeyType','char','ValueType','any');
method_linetype("scenario_precompute_eta_2_smpc") = '-.'
method_linetype("precompute_eta_2_drsmpc") = '-';
method_linetype("robust_mpc") = ':';
method_linetype("online_drsmpc") = '-';
method_name= containers.Map('KeyType','char','ValueType','any');
method_name("scenario_precompute_eta_2_smpc") = 'Scenario approach'
method_name("precompute_eta_2_drsmpc") = 'DRSMPC';
method_name("robust_mpc") = 'RMPC';
method_name("online_drsmpc") = 'Online DRSMPC';
figure;
hold on
for method=method_list
    plot_trajectory(x_list(method),method);
end
%plot_trajectory(x_list("robust_mpc"),'robust mpc');
%plot_trajectory(x_list("precompute_eta_drsmpc"),'precompute eta drsmpc');
%plot_trajectory(x_list("scenario_precompute_eta_smpc"),'scenario precompute eta smpc');
%plot_trajectory(x_list("scenario_precompute_eta_2_smpc"),'scenario precompute eta 2 smpc');
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
%plot(X_1_ul_list('robust_mpc'),'DisplayName','robust upper limits');
%plot(X_1_ll_list('robust_mpc'),'DisplayName','robust lower limits');
%plot(X_1_ul_list('precompute_eta_drsmpc'),'DisplayName','precompute eta drsmpc upper limits');
%plot(X_1_ll_list('precompute_eta_drsmpc'),'DisplayName','precompute eta drsmpc lower limits');
%plot(X_1_ul_list('scenario_precompute_eta_smpc'),'DisplayName','scenario precompute eta smpc upper limits');
%plot(X_1_ll_list('scenario_precompute_eta_smpc'),'DisplayName','scenario precompute eta smpc lower limits');
%plot(X_1_ul_list('scenario_precompute_eta_2_smpc'),'DisplayName','scenario precompute eta 2 smpc upper limits');
%plot(X_1_ll_list('scenario_precompute_eta_2_smpc'),'DisplayName','scenario precompute eta 2 smpc lower limits');
hold off
legend('Interpreter','none');
xlabel('i');
ylabel('limits');
title('z(1) limits');
saveas(gcf,'SMPC_reachability_z1.png')
%% Plot Z_2 limits
figure;
hold on
for method=method_list
    temp = X_2_ul_list(method);
    plot(1:size(X_2_ul_list(method),2)-1,temp(2:end),method_linetype(method),'DisplayName',method_name(method)+" upper limit");
    temp = X_2_ll_list(method);
    plot(1:size(X_2_ll_list(method),2)-1,temp(2:end),method_linetype(method),'DisplayName',method_name(method)+" lower limit");
end

%plot(X_2_ul_list('robust_mpc'),'DisplayName','robust upper limits');
%plot(X_2_ll_list('robust_mpc'),'DisplayName','robust lower limits');
%plot(X_2_ul_list('precompute_eta_drsmpc'),'DisplayName','precompute eta drsmpc upper limits');
%plot(X_2_ll_list('precompute_eta_drsmpc'),'DisplayName','precompute eta drsmpc lower limits');
%plot(X_2_ul_list('scenario_precompute_eta_smpc'),'DisplayName','scenario precompute eta smpc upper limits');
%plot(X_2_ll_list('scenario_precompute_eta_smpc'),'DisplayName','scenario precompute eta smpc lower limits');
%plot(X_2_ul_list('scenario_precompute_eta_2_smpc'),'DisplayName','scenario precompute eta 2 smpc upper limits');
%plot(X_2_ll_list('scenario_precompute_eta_2_smpc'),'DisplayName','scenario precompute eta 2 smpc lower limits');
hold off
legend('Interpreter','none');
xlabel('i');
ylabel('limits');
title('z(2) limits');
plt = Plot();
plt.BoxDim = [5, 3];

saveas(gcf,'SMPC_reachability_z2.png')
%% Plot U limits
figure;
hold on
for method=method_list
    plot(0:size(U_ul_list(method),2)-1,U_ul_list(method),method_linetype(method),'DisplayName',method_name(method)+" upper limit");
    plot(0:size(U_ll_list(method),2)-1,U_ll_list(method),method_linetype(method),'DisplayName',method_name(method)+" lower limit");
end
%plot(U_ul_list('robust_mpc'),'DisplayName','robust upper limits');
%plot(U_ll_list('robust_mpc'),'DisplayName','robust lower limits');
%plot(U_ul_list('precompute_eta_drsmpc'),'DisplayName','precompute eta drsmpc upper limits');
%plot(U_ll_list('precompute_eta_drsmpc'),'DisplayName','precompute eta drsmpc lower limits');
%plot(U_ul_list('scenario_precompute_eta_smpc'),'DisplayName','scenario precompute eta smpc upper limits');
%plot(U_ll_list('scenario_precompute_eta_smpc'),'DisplayName','scenario precompute eta smpc lower limits');
%plot(U_ul_list('scenario_precompute_eta_2_smpc'),'DisplayName','scenario precompute eta 2 smpc upper limits');
%plot(U_ll_list('scenario_precompute_eta_2_smpc'),'DisplayName','scenario precompute eta 2 smpc lower limits');
hold off
legend('Interpreter','none');
xlabel('i');
ylabel('limits');
title('v limits');
plt = Plot();
plt.BoxDim = [5, 3];

saveas(gcf,'SMPC_reachability_v.png')


num_realizations = 10;
method_cost = containers.Map('KeyType','char','ValueType','any');
for method=method_list
    method_cost(method)=0;
end
for i=1:num_realizations
    [x_list,u_list,z_list,e_list,v_list,X_1_ul_list,X_1_ll_list,X_2_ul_list,X_2_ll_list,U_ul_list,U_ll_list] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier);
    for method=method_list
        method_cost(method) = method_cost(method) + find_cost(z_list(method),v_list(method),P,Q,R);
    end
end
for method=method_list
    method_cost(method) = method_cost(method)/num_realizations;
end
%cost is the same for both