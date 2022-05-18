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
T = 6;
method_list_theta = ["precompute_eta_2_drsmpc"];
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
for theta = theta_list
    [x_list,u_list,z_list,e_list,v_list,X_1_ul_list,X_1_ll_list,X_2_ul_list,X_2_ll_list,U_ul_list,U_ll_list] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list_theta, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier);
    x_list_theta(theta) = x_list;
    u_list_theta(theta) = u_list;
    z_list_theta(theta) = z_list;
    e_list_theta(theta) = e_list;
    v_list_theta(theta) = v_list;
    X_1_ul_list_theta(theta) = X_1_ul_list;
    X_1_ll_list_theta(theta) = X_1_ll_list;
    X_2_ul_list_theta(theta) = X_2_ul_list;
    X_2_ll_list_theta(theta) = X_2_ll_list;
    U_ul_list_theta(theta) = U_ul_list;
    U_ll_list_theta(theta) = U_ll_list;
end
method_list = [ ];
method_linetype= containers.Map('KeyType','char','ValueType','any');
method_linetype("scenario_precompute_eta_2_smpc") = '-.';
method_linetype("precompute_eta_2_drsmpc") = '-';
method_linetype("robust_mpc") = ':';
method_name= containers.Map('KeyType','char','ValueType','any');
method_name("scenario_precompute_eta_2_smpc") = 'Scenario approach'
method_name("precompute_eta_2_drsmpc") = 'DRSMPC';
method_name("robust_mpc") = 'RMPC';
theta_linetype = containers.Map('KeyType','double','ValueType','any');
theta_linetype(theta_list(1)) = ':';
theta_linetype(theta_list(2)) = '-.';
theta_linetype(theta_list(3)) = '-';
%[x_list,u_list,z_list,e_list,v_list,X_1_ul_list,X_1_ll_list,X_2_ul_list,X_2_ll_list,U_ul_list,U_ll_list] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list, beta_scenario,theta,X,W,U,S0,w_distribution);
%% Plot Z_1 limits
figure;
%frame_h = get(handle(gcf),'JavaFrame');
%set(frame_h,'Maximized',1);
hold on
%for method=method_list
%    plot(X_1_ul_list(method),'DisplayName',strrep(method+" upper",'_',' '));
%    plot(X_1_ll_list(method),'DisplayName',strrep(method+" lower",'_',' '));
%end

for theta = theta_list
    temp = X_1_ul_list_theta(theta);
    temp = temp(method_list_theta(1));
    plot(1:(size(temp,2)-1),temp(2:end),theta_linetype(theta),'DisplayName',strrep("upper limit for \theta = "+theta,'_',' '));
    temp = X_1_ll_list_theta(theta);
    temp = temp(method_list_theta(1));
    plot(1:(size(temp,2)-1),temp(2:end),theta_linetype(theta),'DisplayName',strrep("lower limit for \theta = "+theta,'_',' '));
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'z(1) limits'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'RMPC_gaussian_theta_example_Z_1.png')
%% Plot Z_2 limits
figure;
%frame_h = get(handle(gcf),'JavaFrame');
%set(frame_h,'Maximized',1);
hold on
%for method=method_list
%    plot(X_2_ul_list(method),'DisplayName',strrep(method+" upper",'_',' '));
%    plot(X_2_ll_list(method),'DisplayName',strrep(method+" lower",'_',' '));
%end
for theta = theta_list
    temp = X_2_ul_list_theta(theta);
    temp = temp(method_list_theta(1));
    plot(1:(size(temp,2)-1),temp(2:end),theta_linetype(theta),'DisplayName',strrep("upper limit for \theta = "+theta,'_',' '));
    temp = X_2_ll_list_theta(theta);
    temp = temp(method_list_theta(1));
    plot(1:(size(temp,2)-1),temp(2:end),theta_linetype(theta),'DisplayName',strrep("lower limit for \theta = "+theta,'_',' '));
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'z(2) limits'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'RMPC_gaussian_theta_example_Z_2.png')
%% Plot U limits

figure;
%frame_h = get(handle(gcf),'JavaFrame');
%set(frame_h,'Maximized',1);
hold on
%for method=method_list
%    plot(U_ul_list(method),'DisplayName',strrep(method+" upper",'_',' '));
%    plot(U_ll_list(method),'DisplayName',strrep(method+" lower",'_',' '));
%end
for theta = theta_list
    temp = U_ul_list_theta(theta);
    theta
    plot(0:(size(temp(method_list_theta(1)),2)-1),temp(method_list_theta(1)),theta_linetype(theta),'DisplayName',strrep("upper limit for \theta = "+theta,'_',' '));
    temp = U_ll_list_theta(theta);
    plot(0:(size(temp(method_list_theta(1)),2)-1),temp(method_list_theta(1)),theta_linetype(theta),'DisplayName',strrep("lower limit for \theta = "+theta,'_',' '));
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'v limits'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'RMPC_gaussian_theta_example_U.png')