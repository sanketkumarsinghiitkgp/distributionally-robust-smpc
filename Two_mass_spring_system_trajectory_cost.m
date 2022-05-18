spring_k = 1;
m_1 = 0.5;
m_2 = 2;
A = [1, 0, 0.1, 0;0, 1, 0, 0.1;-spring_k/m_1, 0.1*spring_k/m_1, 1, 0;spring_k/m_2, -0.1*spring_k/m_2, 0, 1];
B = [0;0;0.1/m_1;0];
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
num_iter = 31;
M = 100;
T = 6;
beta_scenario = 0.95;
theta = 1e-5;
w_distribution = "Two_mass_spring";
sigma_multiplier = 0.0001;

method_list = ["precompute_eta_2_drsmpc", "scenario_precompute_eta_2_smpc"];
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

figure;
hold on
for method=method_list
    plot(X_1_ul_list(method),'DisplayName',strrep(method+" upper",'_',' '));
    plot(X_1_ll_list(method),'DisplayName',strrep(method+" lower",'_',' '));
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraint Tightening for Velocity of Mass 1'; % plot title
plt.BoxDim = [10, 6];


figure;
hold on
for method=method_list
    plot(X_2_ul_list(method),'DisplayName',strrep(method+" upper",'_',' '));
    plot(X_2_ll_list(method),'DisplayName',strrep(method+" lower",'_',' '));
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraint Tightening for Velocity of Mass 2'; % plot title
plt.BoxDim = [10, 6];


figure;
hold on
for method=method_list
    plot(U_ul_list(method),'DisplayName',strrep(method+" upper",'_',' '));
    plot(U_ll_list(method),'DisplayName',strrep(method+" lower",'_',' '));
end
hold off
legend('Interpreter','none');
plt = Plot(); % create a Plot object and grab the current figure
plt.XLabel = 'i'; % xlabel
plt.YLabel = 'limits'; %ylabel
plt.Title = 'Constraint Tightening for Control Effort'; % plot title
plt.BoxDim = [10, 6];

for method=method_list
    "cost of method "+method+" is "+num2str(method_cost(method))
end