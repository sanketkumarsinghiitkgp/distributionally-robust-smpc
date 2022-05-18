x_list = containers.Map('KeyType','char','ValueType','any');
z_list = containers.Map('KeyType','char','ValueType','any');
u_list = containers.Map('KeyType','char','ValueType','any');
e_list = containers.Map('KeyType','char','ValueType','any');
v_list = containers.Map('KeyType','char','ValueType','any');
X_1_ul_list = containers.Map('KeyType','char','ValueType','any');
X_1_ll_list = containers.Map('KeyType','char','ValueType','any');
X_2_ul_list = containers.Map('KeyType','char','ValueType','any');
X_2_ll_list = containers.Map('KeyType','char','ValueType','any');
U_ul_list = containers.Map('KeyType','char','ValueType','any');
U_ll_list = containers.Map('KeyType','char','ValueType','any');
prob_constraint_satisfied_x_list = containers.Map('KeyType','char','ValueType','any');
prob_constraint_satisfied_u_list = containers.Map('KeyType','char','ValueType','any');

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
num_iter = 100;
num_samples = 20;
M = 100;
T = 6;
beta_scenario = 0.95;
theta = 1e-5;


%%%CHANGED
%w_distribution = "Two_mass_spring";
w_distribution = "Two_mass_spring_triangle";
%%%


sigma_multiplier = 0.0001; %5 TIMES THE PREVIOUS SIGMA_MULTIPLIER
c_theta=1e4;
method_list = ["precompute_eta_2_drsmpc", "online_drsmpc"];
for method=method_list
    X_1_ul_list(method) = 0;
    X_1_ll_list(method) = 0;
    X_2_ul_list(method) = 0;
    X_2_ll_list(method) = 0;
    U_ul_list(method) = 0;
    U_ll_list(method) = 0;
    x_list(method) = 0;
    u_list(method) = 0;
    z_list(method) = 0;
    e_list(method) = 0;
    v_list(method) = 0;
end
for i=1:num_samples
    [x_list_cur,u_list_cur,z_list_cur,e_list_cur,v_list_cur,X_1_ul_list_cur,X_1_ll_list_cur,X_2_ul_list_cur,X_2_ll_list_cur,U_ul_list_cur,U_ll_list_cur,~,~] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier,c_theta);
    for method=method_list
        X_1_ul_list(method) = X_1_ul_list(method)+X_1_ul_list_cur(method);
        X_1_ll_list(method) = X_1_ll_list(method)+X_1_ll_list_cur(method);
        X_2_ul_list(method) = X_2_ul_list(method)+X_2_ul_list_cur(method);
        X_2_ll_list(method) = X_2_ll_list(method)+X_2_ll_list_cur(method);
        U_ul_list(method) = U_ul_list(method)+U_ul_list_cur(method);
        U_ll_list(method) = U_ll_list(method)+U_ll_list_cur(method);
        x_list(method) = x_list(method)+x_list_cur(method);
        u_list(method) = u_list(method)+u_list_cur(method);
        z_list(method) = z_list(method)+z_list_cur(method);
        e_list(method) = e_list(method)+e_list_cur(method);
        v_list(method) = v_list(method)+v_list_cur(method);
    end
end
for method=method_list
    X_1_ul_list(method) = X_1_ul_list(method)/num_samples;
    X_1_ll_list(method) = X_1_ll_list(method)/num_samples;
    X_2_ul_list(method) = X_2_ul_list(method)/num_samples;
    X_2_ll_list(method) = X_2_ll_list(method)/num_samples;
    U_ul_list(method) = U_ul_list(method)/num_samples;
    U_ll_list(method) = U_ll_list(method)/num_samples;
    x_list(method) = x_list(method)/num_samples;
    u_list(method) = u_list(method)/num_samples;
    z_list(method) = z_list(method)/num_samples;
    e_list(method) = e_list(method)/num_samples;
    v_list(method) = v_list(method)/num_samples;
end

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
    if method ~= "online_drsmpc"
        plot(0:(size(X_1_ul_list(method),2)-1),X_1_ul_list(method),method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        plot(0:(size(X_1_ll_list(method),2)-1),X_1_ll_list(method),method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    end
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraint Tightening for Velocity of Mass 1'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_vm1_ensemble_triangle.png')

figure;
hold on
for method=method_list
    if method ~= "online_drsmpc"
        plot(0:(size(X_2_ul_list(method),2)-1),X_2_ul_list(method),method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        plot(0:(size(X_2_ll_list(method),2)-1),X_2_ll_list(method),method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    end
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraint Tightening for Velocity of Mass 2'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_vm2_ensemble_triangle.png')

figure;
hold on
for method=method_list
    if method ~= "online_drsmpc"
        plot(0:(size(U_ul_list(method),2)-1),U_ul_list(method),method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        plot(0:(size(U_ll_list(method),2)-1),U_ll_list(method),method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    end
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraint Tightening for Control Effort'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_ce_ensemble_triangle.png')

figure;
hold on
for method=method_list
    temp = x_list(method);
    plot(0:(size(temp,1)-1),temp(:,1),method_linetype(method),'DisplayName',strrep(method_name(method),'_',' '));
end

hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time'; % xlabel
plt.YLabel = 'Position'; %ylabel
plt.Title = 'Position of Mass 1 vs Time'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_pm1_vs_k_ensemble_triangle.png')

figure;
hold on
for method=method_list
    temp = x_list(method);
    plot(0:(size(temp,1)-1),temp(:,2),method_linetype(method),'DisplayName',strrep(method_name(method),'_',' '));
end

hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time'; % xlabel
plt.YLabel = 'Position'; %ylabel
plt.Title = 'Position of Mass 2 vs Time'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_pm2_vs_k_ensemble_triangle.png')

figure;
hold on
for method=method_list
    temp = x_list(method);
    plot(0:(size(temp,1)-1),temp(:,3),method_linetype(method),'DisplayName',strrep(method_name(method),'_',' '));
end

hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time'; % xlabel
plt.YLabel = 'Velocity'; %ylabel
plt.Title = 'Velocity of Mass 1 vs Time'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_vm1_vs_k_ensemble_triangle.png')




figure;
hold on
for method=method_list
    temp = x_list(method);
    plot(0:(size(temp,1)-1),temp(:,4),method_linetype(method),'DisplayName',strrep(method_name(method),'_',' '));
end

hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time'; % xlabel
plt.YLabel = 'Velocity'; %ylabel
plt.Title = 'Velocity of Mass 2 vs Time'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_vm2_vs_k_ensemble_triangle.png')


figure;
hold on
for method=method_list
    temp = u_list(method);
    plot(0:(size(temp,1)-1),temp(:,1),method_linetype(method),'DisplayName',strrep(method_name(method),'_',' '));
end

hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'Time'; % xlabel
plt.YLabel = 'Control Effort'; %ylabel
plt.Title = 'Control Effort vs Time'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_ce_vs_k_ensemble_triangle.png')

%%------

figure;
hold on
for method=method_list
    if method == 'online_drsmpc' || method == 'online_drsmpc_2'
        temp = U_ul_list(method);
        temp = temp(2,:);
        plot(0:(size(temp,2)-1),temp,method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        temp = U_ll_list(method);
        temp = temp(2,:);
        plot(0:(size(temp,2)-1),temp,method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
        
    elseif method =='precompute_eta_2_drsmpc'
        temp = U_ul_list(method);
        plot(0:(num_iter-1),temp(6)*ones(1,num_iter),method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        temp = U_ll_list(method);
        plot(0:(num_iter-1),temp(6)*ones(1,num_iter),method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    end
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'k'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraints on v_5 vs Time (k)'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_ce_6_ensemble_triangle.png')

figure;
hold on
for method=method_list
    if method == 'online_drsmpc' || method == 'online_drsmpc_2'
        temp = X_1_ul_list(method);
        temp = temp(2,:);
        plot(0:(size(temp,2)-1),temp,method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        temp = X_1_ll_list(method);
        temp = temp(2,:);
        plot(0:(size(temp,2)-1),temp,method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    elseif method =='precompute_eta_2_drsmpc'
        temp = X_1_ul_list(method);
        plot(0:(num_iter-1),temp(6)*ones(1,num_iter),method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        temp = X_1_ll_list(method);
        plot(0:(num_iter-1),temp(6)*ones(1,num_iter),method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    end
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'k'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraints on [z_5]_3 vs Time (k)'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_vm1_6_ensemble_triangle.png')

figure;
hold on
for method=method_list
    if method == 'online_drsmpc' || method == 'online_drsmpc_2'
        temp = X_2_ul_list(method);
        temp = temp(2,:);
        plot(0:(size(temp,2)-1),temp,method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        temp = X_2_ll_list(method);
        temp = temp(2,:);
        plot(0:(size(temp,2)-1),temp,method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    elseif method =='precompute_eta_2_drsmpc'
        temp = X_2_ul_list(method);
        plot(0:(num_iter-1),temp(6)*ones(1,num_iter),method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        temp = X_2_ll_list(method);
        plot(0:(num_iter-1),temp(6)*ones(1,num_iter),method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    end
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'k'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraints on [z_5]_4 vs Time (k)'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_vm2_6_ensemble_triangle.png')






figure;
hold on
for method=method_list
    if method == 'online_drsmpc' || method == 'online_drsmpc_2'
        temp = U_ul_list(method);
        temp = temp(1,:);
        plot(0:(size(temp,2)-1),temp,method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        temp = U_ll_list(method);
        temp = temp(1,:);
        plot(0:(size(temp,2)-1),temp,method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
        
    elseif method =='precompute_eta_2_drsmpc'
        temp = U_ul_list(method);
        plot(0:(num_iter-1),temp(6)*ones(1,num_iter),method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        temp = U_ll_list(method);
        plot(0:(num_iter-1),temp(6)*ones(1,num_iter),method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    end
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'k'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraints on v_0 vs Time (k)'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_ce_1_ensemble_triangle.png')

figure;
hold on
for method=method_list
    if method == 'online_drsmpc' || method == 'online_drsmpc_2'
        temp = X_1_ul_list(method);
        temp = temp(1,:);
        plot(0:(size(temp,2)-1),temp,method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        temp = X_1_ll_list(method);
        temp = temp(1,:);
        plot(0:(size(temp,2)-1),temp,method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    elseif method =='precompute_eta_2_drsmpc'
        temp = X_1_ul_list(method);
        plot(0:(num_iter-1),temp(6)*ones(1,num_iter),method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        temp = X_1_ll_list(method);
        plot(0:(num_iter-1),temp(6)*ones(1,num_iter),method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    end
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'k'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraints on [z_1]_3 vs Time (k)'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_vm1_2_ensemble_triangle.png')

figure;
hold on
for method=method_list
    if method == 'online_drsmpc' || method == 'online_drsmpc_2'
        temp = X_2_ul_list(method);
        temp = temp(1,:);
        plot(0:(size(temp,2)-1),temp,method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        temp = X_2_ll_list(method);
        temp = temp(1,:);
        plot(0:(size(temp,2)-1),temp,method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    elseif method =='precompute_eta_2_drsmpc'
        temp = X_2_ul_list(method);
        plot(0:(num_iter-1),temp(6)*ones(1,num_iter),method_linetype(method),'DisplayName',strrep(method_name(method)+" upper",'_',' '));
        temp = X_2_ll_list(method);
        plot(0:(num_iter-1),temp(6)*ones(1,num_iter),method_linetype(method),'DisplayName',strrep(method_name(method)+" lower",'_',' '));
    end
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'k'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraints on [z_1]_4 vs Time (k)'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'Two_mass_spring_system_vm2_2_ensemble_triangle.png')
