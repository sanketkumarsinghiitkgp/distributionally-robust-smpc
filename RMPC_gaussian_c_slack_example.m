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
num_iter = 31;
x0 = [2.4;-0.8];
z0 = x0;
e0 = x0-z0;
S0 = Polyhedron('Ae',[1 0; 0 1], 'be',[0;0]);
%S_cur = W;
K = [-0.4 -1.2];

eps_x = 0.1;
eps_u = 0.1;
sigma = 1; %(doubtful)
M = 100;
T = 15;
method_list_c_theta = ["online_drsmpc_2"];
beta_scenario = 0.95;
sigma_multiplier = 0.01;
w_distribution = "RMPC_gaussian_example";
theta_list = [1e-3,1e-4,1e-5];
x_list_theta = containers.Map('KeyType','double','ValueType','any');
z_list_theta = containers.Map('KeyType','double','ValueType','any');
u_list_theta = containers.Map('KeyType','double','ValueType','any');
e_list_theta = containers.Map('KeyType','double','ValueType','any');
v_list_theta = containers.Map('KeyType','double','ValueType','any');
X_1_ul_list_theta = containers.Map('KeyType','double','ValueType','any');
X_1_ll_list_theta = containers.Map('KeyType','double','ValueType','any');
X_2_ul_list_theta = containers.Map('KeyType','double','ValueType','any');
X_2_ll_list_theta = containers.Map('KeyType','double','ValueType','any');
U_ul_list_theta = containers.Map('KeyType','double','ValueType','any');
U_ll_list_theta = containers.Map('KeyType','double','ValueType','any');
prob_constraint_satisfied_x_list_theta = containers.Map('KeyType','double','ValueType','any');
prob_constraint_satisfied_u_list_theta = containers.Map('KeyType','double','ValueType','any');
theta=1e-5;
c_theta_list = [1e-4,1,1e4];
for c_theta = c_theta_list
    [x_list,u_list,z_list,e_list,v_list,X_1_ul_list,X_1_ll_list,X_2_ul_list,X_2_ll_list,U_ul_list,U_ll_list] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list_c_theta, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier,c_theta);
    x_list_theta(c_theta) = x_list;
    u_list_theta(c_theta) = u_list;
    z_list_theta(c_theta) = z_list;
    e_list_theta(c_theta) = e_list;
    v_list_theta(c_theta) = v_list;
    X_1_ul_list_theta(c_theta) = X_1_ul_list;
    X_1_ll_list_theta(c_theta) = X_1_ll_list;
    X_2_ul_list_theta(c_theta) = X_2_ul_list;
    X_2_ll_list_theta(c_theta) = X_2_ll_list;
    U_ul_list_theta(c_theta) = U_ul_list;
    U_ll_list_theta(c_theta) = U_ll_list;
end
method_list = [ "online_drsmpc"];
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
c_theta_linetype = containers.Map('KeyType','double','ValueType','any');
c_theta_linetype(c_theta_list(1)) = ':';
c_theta_linetype(c_theta_list(2)) = '-.';
c_theta_linetype(c_theta_list(3)) = '-';
c_theta = 1e4;
[x_list,u_list,z_list,e_list,v_list,X_1_ul_list,X_1_ll_list,X_2_ul_list,X_2_ll_list,U_ul_list,U_ll_list] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier,c_theta);
%% Plot Z_1 limits
figure;
%frame_h = get(handle(gcf),'JavaFrame');
%set(frame_h,'Maximized',1);
hold on
%for method=method_list
%    plot(X_1_ul_list(method),'DisplayName',strrep(method+" upper",'_',' '));
%    plot(X_1_ll_list(method),'DisplayName',strrep(method+" lower",'_',' '));
%end

for c_theta = c_theta_list
    temp = X_1_ul_list_theta(c_theta);
    plot(0:(size(temp(method_list_c_theta(1)),2)-1),temp(method_list_c_theta(1)),c_theta_linetype(c_theta),'DisplayName',strrep("upper limit for c_{\Theta} = "+c_theta,'_',' '));
    temp = X_1_ll_list_theta(c_theta);
    plot(0:(size(temp(method_list_c_theta(1)),2)-1),temp(method_list_c_theta(1)),c_theta_linetype(c_theta),'DisplayName',strrep("lower limit for c_{\Theta} = "+c_theta,'_',' '));
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'z(1) limits'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'RMPC_gaussian_c_theta_example_Z_1.png')
%% Plot Z_2 limits
figure;
%frame_h = get(handle(gcf),'JavaFrame');
%set(frame_h,'Maximized',1);
hold on
%for method=method_list
%    plot(X_2_ul_list(method),'DisplayName',strrep(method+" upper",'_',' '));
%    plot(X_2_ll_list(method),'DisplayName',strrep(method+" lower",'_',' '));
%end
for c_theta = c_theta_list
    temp = X_2_ul_list_theta(c_theta);
    plot(0:(size(temp(method_list_c_theta(1)),2)-1),temp(method_list_c_theta(1)),c_theta_linetype(c_theta),'DisplayName',strrep("upper limit for c_{\Theta} = "+c_theta,'_',' '));
    temp = X_2_ll_list_theta(c_theta);
    plot(0:(size(temp(method_list_c_theta(1)),2)-1),temp(method_list_c_theta(1)),c_theta_linetype(c_theta),'DisplayName',strrep("lower limit for c_{\Theta} = "+c_theta,'_',' '));
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'z(2) limits'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'RMPC_gaussian_c_theta_example_Z_2.png')
%% Plot U limits

figure;
%frame_h = get(handle(gcf),'JavaFrame');
%set(frame_h,'Maximized',1);
hold on
%for method=method_list
%    plot(U_ul_list(method),'DisplayName',strrep(method+" upper",'_',' '));
%    plot(U_ll_list(method),'DisplayName',strrep(method+" lower",'_',' '));
%end
for c_theta = c_theta_list
    temp = U_ul_list_theta(c_theta);
    plot(0:(size(temp(method_list_c_theta(1)),2)-1),temp(method_list_c_theta(1)),c_theta_linetype(c_theta),'DisplayName',strrep("upper limit for c_{\Theta} = "+c_theta,'_',' '));
    temp = U_ll_list_theta(c_theta);
    plot(0:(size(temp(method_list_c_theta(1)),2)-1),temp(method_list_c_theta(1)),c_theta_linetype(c_theta),'DisplayName',strrep("lower limit for c_{\Theta} = "+c_theta,'_',' '));
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'v limits'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'RMPC_gaussian_c_theta_example_U.png')

figure
hold on
for c_theta=c_theta_list
temp = x_list_theta(c_theta);
temp(method_list_c_theta(1))
plot_trajectory(temp(method_list_c_theta(1)),"c_theta="+c_theta);
x = input("abc")
end
hold off
legend