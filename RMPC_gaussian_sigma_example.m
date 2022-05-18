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
method_list_sigma_multiplier = ["precompute_eta_2_drsmpc"];
beta_scenario = 0.95;
theta = 1e-5;
w_distribution = "RMPC_gaussian_example";
sigma_multiplier_list = [1e-3,1e-2,1e-1];
sigma_linetype = containers.Map('KeyType','double','ValueType','any');
sigma_linetype(sigma_multiplier_list(1)) = ':';
sigma_linetype(sigma_multiplier_list(2)) = '-.';
sigma_linetype(sigma_multiplier_list(3)) = '--';
x_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
z_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
u_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
e_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
v_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
X_1_ul_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
X_1_ll_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
X_2_ul_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
X_2_ll_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
U_ul_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
U_ll_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
prob_constraint_satisfied_x_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
prob_constraint_satisfied_u_list_sigma_multiplier = containers.Map('KeyType','double','ValueType','any');
for sigma_multiplier = sigma_multiplier_list
    [x_list,u_list,z_list,e_list,v_list,X_1_ul_list,X_1_ll_list,X_2_ul_list,X_2_ll_list,U_ul_list,U_ll_list] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list_sigma_multiplier, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier);
    x_list_sigma_multiplier(sigma_multiplier) = x_list;
    u_list_sigma_multiplier(sigma_multiplier) = u_list;
    z_list_sigma_multiplier(sigma_multiplier) = z_list;
    e_list_sigma_multiplier(sigma_multiplier) = e_list;
    v_list_sigma_multiplier(sigma_multiplier) = v_list;
    X_1_ul_list_sigma_multiplier(sigma_multiplier) = X_1_ul_list;
    X_1_ll_list_sigma_multiplier(sigma_multiplier) = X_1_ll_list;
    X_2_ul_list_sigma_multiplier(sigma_multiplier) = X_2_ul_list;
    X_2_ll_list_sigma_multiplier(sigma_multiplier) = X_2_ll_list;
    U_ul_list_sigma_multiplier(sigma_multiplier) = U_ul_list;
    U_ll_list_sigma_multiplier(sigma_multiplier) = U_ll_list;
end
method_list = [ "robust_mpc"];
sigma_multiplier = 0.01;
[x_list,u_list,z_list,e_list,v_list,X_1_ul_list,X_1_ll_list,X_2_ul_list,X_2_ll_list,U_ul_list,U_ll_list] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,W,num_iter,M,T,method_list, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier);
%% Plot Z_1 limits
figure;
%frame_h = get(handle(gcf),'JavaFrame');
%set(frame_h,'Maximized',1);
hold on
for method=method_list
    temp = X_1_ul_list(method);
    plot(1:(size(X_1_ul_list(method),2)-1),temp(2:end),'DisplayName','RMPC upper limit');
    temp = X_1_ll_list(method);
    plot(1:(size(X_1_ll_list(method),2)-1),temp(2:end),'DisplayName','RMPC lower limit');
end
for sigma_multiplier = sigma_multiplier_list
    temp = X_1_ul_list_sigma_multiplier(sigma_multiplier);
    temp = temp(method_list_sigma_multiplier(1));
    plot(1:(size(temp,2)-1),temp(2:end),sigma_linetype(sigma_multiplier),'DisplayName',"DRSMPC upper limit for"+" \sigma_m = "+sigma_multiplier);
    temp = X_1_ll_list_sigma_multiplier(sigma_multiplier);
    temp = temp(method_list_sigma_multiplier(1));
    plot(1:(size(temp,2)-1),temp(2:end),sigma_linetype(sigma_multiplier),'DisplayName',"DRSMPC lower limit for"+" \sigma_m = "+sigma_multiplier);
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'z(1) limits'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'RMPC_gaussian_sigma_example_Z_1.png')
%% Plot Z_2 limits
figure;
%frame_h = get(handle(gcf),'JavaFrame');
%set(frame_h,'Maximized',1);
hold on
for method=method_list
    temp = X_2_ul_list(method);
    plot(1:(size(X_2_ul_list(method),2)-1),temp(2:end),'DisplayName','RMPC upper limit');
    temp = X_2_ll_list(method);
    plot(1:(size(X_2_ll_list(method),2)-1),temp(2:end),'DisplayName','RMPC lower limit');
end
for sigma_multiplier = sigma_multiplier_list
    temp = X_2_ul_list_sigma_multiplier(sigma_multiplier);
    temp = temp(method_list_sigma_multiplier(1));
    plot(1:(size(temp,2)-1),temp(2:end),sigma_linetype(sigma_multiplier),'DisplayName',"DRSMPC upper limit for"+" \sigma_m = "+sigma_multiplier);
    temp = X_2_ll_list_sigma_multiplier(sigma_multiplier);
    temp = temp(method_list_sigma_multiplier(1));
    plot(1:(size(temp,2)-1),temp(2:end),sigma_linetype(sigma_multiplier),'DisplayName',"DRSMPC lower limit for"+" \sigma_m = "+sigma_multiplier);
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'z(2) limits'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'RMPC_gaussian_sigma_example_Z_2.png')
%% Plot U limits

figure;
%frame_h = get(handle(gcf),'JavaFrame');
%set(frame_h,'Maximized',1);
hold on
for method=method_list
    plot(0:(size(U_ul_list(method),2)-1),U_ul_list(method),'DisplayName','RMPC upper limit');
    plot(0:(size(U_ll_list(method),2)-1),U_ll_list(method),'DisplayName','RMPC lower limit');
end
for sigma_multiplier = sigma_multiplier_list
    temp = U_ul_list_sigma_multiplier(sigma_multiplier);
    plot(0:(size(temp(method_list_sigma_multiplier(1)),2)-1),temp(method_list_sigma_multiplier(1)),sigma_linetype(sigma_multiplier),'DisplayName',"DRSMPC upper limit for"+" \sigma_m = "+sigma_multiplier);
    temp = U_ll_list_sigma_multiplier(sigma_multiplier);
    plot(0:(size(temp(method_list_sigma_multiplier(1)),2)-1),temp(method_list_sigma_multiplier(1)),sigma_linetype(sigma_multiplier),'DisplayName',"DRSMPC lower limit for"+" \sigma_m = "+sigma_multiplier);
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'v limits'; % plot title
plt.BoxDim = [5, 3];
saveas(gcf,'RMPC_gaussian_sigma_example_U.png')