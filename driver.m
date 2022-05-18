
rng('default')  % For reproducibility
H = [0 1; 0 -1];
h = [1.2; 1.2];
G = [1; -1];
g = [6; 6];
A = [ 1 1; 0 1 ];
B = [0.5;1];
K = [-0.4 -1.2];
z = [10; 0.5];
x0 = [10;0.5];
eps_x = 0.4;
eps_u = 0.1;
N = 10;
sigma = 1;
w_support = Polyhedron('A',[0 0],'b',[0]);
[z_history,e_history, u_history, radii_x_history, radii_u_history, eta_x_history, eta_u_history] = find_eta_list(A,B,K,eps_x, eps_u, H,G,z,N,w_support,sigma);
eta_x_history
v = sdpvar(N,1,'full');
P = diag([0.1 1]);
R = 0.1;
z0 = sdpvar(2,1,'full');
z_var = [ z0 ];
for i=1:N
    z_last = z_var(:,end); 
    z_var = [z_var A*z_last+B*v(i,:)'];
end
objective = z_var(:,N+1)'*P*z_var(:,N+1);
for i=1:N
    zz = z_var(:,i);
    objective = objective + (zz'*P*zz+v(i,:)*R*v(i,:)');
end
constraints = [z_var(:,1) == x0]
for i=1:N+1
    zz = z_var(:,i);
    constraints = [constraints, H*zz<=h-eta_x_history{i}'] % eta_x_history length should be increased
end
for i=1:N
    constraints = [constraints G*v(i,:)' <= g-eta_u_history{i}'];
end
sol = optimize(constraints,objective,sdpsettings('solver', 'mosek'));
solution = 0;
if sol.problem == 0
    solution = value(v);
else
    display('Something went wrong');
    sol.info;
    yalmiperror(sol.problem)
end
solution