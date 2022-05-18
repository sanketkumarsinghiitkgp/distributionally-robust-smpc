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
S0 = Polyhedron('Ae',[1 0; 0 1], 'be',e0);
%S_cur = W;
K = [-0.4 -1.2];

eps_x = 0.1;
eps_u = 0.1;
sigma = 1; %(doubtful)
M = 100;
T = 6;
eta_scenario = 0.95;
theta = 1e-5;
sigma_multiplier = 0.01;


x = [x0'];
u = [];
v = [];
z = [z0'];
e = [e0'];
z_cur = z0;
x_cur = x0;
e_cur = e0;
for k = 1:num_iter
    w_cur = w(M+k,:)';
    
    v_cur = find_v_drsmpc_precompute_eta(P,Q,R,A,B,H,h,G,g,x_cur,T,eta_x_history,eta_u_history);
    if(size(v_cur,2)==2)
        "Infeasibility"
        abc=input("abc");
        v_cur = 0;
    end
    u = [u; (K*e_cur+v_cur)'];
    v = [v;(v_cur)'];
    [new_z, new_e, ~] = compute_next_state(A,B,z_cur,v_cur,e_cur,w_cur,K); %new_u must be removed, it is useless.
    z_cur = new_z;
    e_cur = new_e;
    x_cur = z_cur+e_cur;
    x = [x; x_cur'];
    z = [z; z_cur'];
    e = [e; e_cur'];
end