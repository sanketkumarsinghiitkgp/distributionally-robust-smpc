rng('default')  % For reproducibility
A = [1 1;0 1];
B = [0;1];
Q = [20 0;0 20];
R = 1;
X = Polyhedron('A',[1 0;-1 0;0 1;0 -1],'b',[4;4;1;1]);
W = Polyhedron('A',[1 0;-1 0;0 1;0 -1],'b',[0.2;0.2;0.2;0.2]);
U = Polyhedron('A',[1;-1],'b',[0.7;0.7]);
sim_time = 31;
x_cur = [2.4;-0.8];
z_cur = x_cur;
e_cur = x_cur-z_cur;
S_cur = Polyhedron('Ae',[1 0; 0 1], 'be',[0;0]);
%S_cur = W;
K = [-0.4 -1.2];
Acl = A+B*K;
x =[x_cur'];
u=[];
T = 10
S = [S_cur];
for i=1:sim_time+T
    S_cur = Acl*S_cur+W;
    S = [S, S_cur];
end
X_1_ul = [];
X_1_ll = [];
X_2_ul = [];
X_2_ll = [];
U_ul = [];
U_ll = [];
for i=1:sim_time+T+1
    Z = X-S(i);
    V = U-K*S(i);
    X_1_ul = [X_1_ul Z.b(5)]; %fishy method, relies on MPT preserving the order of inequalities.
    X_1_ll = [X_1_ll -Z.b(6)];
    X_2_ul = [X_2_ul Z.b(7)];
    X_2_ll = [X_2_ll -Z.b(8)];
    U_ul = [U_ul V.b(3)];
    U_ll = [U_ll -V.b(4)];
end
figure;
hold on
plot(X_1_ul,'DisplayName','z_1 upper limits');
plot(X_1_ll,'DisplayName','z_1 lower limits');
plot(X_2_ul,'DisplayName','z_2 upper limits');
plot(X_2_ll,'DisplayName','z_2 lower limits');
plot(U_ul,'DisplayName','v upper limits');
plot(U_ll,'DisplayName','v lower limits');
legend
hold off

for k = 1:sim_time
    v_cur = find_v_rmpc(X,S(k:end),A,B,K,U,x_cur,T,Q,R);
    u=[u v_cur+K*e_cur];
    w_cur = unifrnd(-0.2,0.2,[2 1]);
    z_cur = A*z_cur+B*v_cur;
    e_cur = (Acl)*e_cur+w_cur;
    x_cur = z_cur+e_cur;
    x = [x;x_cur'];
    S_cur = Acl*S_cur+W;
end
figure;
plot_trajectory(x,'State Trajectory');
legend
figure;
hold on
plot(x(:,1),'DisplayName','x_1');
plot(x(:,2),'DisplayName','x_2');
legend
hold off
figure;
plot(u,'DisplayName','u');
legend
