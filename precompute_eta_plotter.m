rng('default')  % For reproducibility
theta_list = [1e-3];
eta_x_history_list = containers.Map('KeyType','double','ValueType','any');
eta_u_history_list = containers.Map('KeyType','double','ValueType','any');
w_list = containers.Map('KeyType','double','ValueType','any');
x_list = containers.Map('KeyType','double','ValueType','any');
z_list = containers.Map('KeyType','double','ValueType','any');
u_list = containers.Map('KeyType','double','ValueType','any');
e_list = containers.Map('KeyType','double','ValueType','any');
v_list = containers.Map('KeyType','double','ValueType','any');
for theta=theta_list
    H = [0 1; 0 -1];
    h = [1.2; 1.2];
    G = [1; -1];
    g = [6; 6];
    %G=[0;0];
    %g = [100000000000; 1000000000000];
    A = [ 1 1; 0 1 ];
    B = [0.5;1];
    P = diag([0.1 1]);
    Q = P;
    R = 0.1;
    K = -dlqr(A,B,P,R,0);
    eps_x = 0.4; %decreased from 0.4 to see the effect 4/11 violations
    eps_u = 0.1;
    T = 10;
    sigma = 1;
    w_support = Polyhedron('A',[0 0],'b',[0]);
    w_hist = [];
    z_cur = [2;1];
    x_cur = [2;1];
    e_cur = x_cur-z_cur;
    e = [e_cur'];
    x = [x_cur'];
    w = [];
    u = [];
    v = [];
    z = [z_cur'];
    for i=1:T+1
        w = [w; w_generator()];
    end
    N = T;
    [radii_x_history, radii_u_history, eta_x_history, eta_u_history] = find_eta_list(A,B,K,eps_x,eps_u,H,G,e_cur,T,w,w_support,sigma,theta);
    for k=1:10
        w_cur = w_generator()';
        w = [w; w_cur'];
        v_cur = find_v_drsmpc_precompute_eta(P,Q  ,R,A,B,H,h,G,g,x_cur,N,eta_x_history,eta_u_history); %we set z_k = e_k for the purpose of this computation.
        u = [u; (K*e_cur+v_cur)'];
        v = [v; v_cur'];
        [new_z, new_e, new_u] = compute_next_state(A,B,z_cur,v_cur,e_cur,w_cur,K); %new_u must be removed, it is useless.
        z_cur = new_z;
        e_cur = new_e;
        x_cur = z_cur+e_cur;
        x = [x; x_cur'];
        z = [z; z_cur'];
        e = [e;e_cur'];
    end
    prob_constraint_satisfied_x = find_prob_constraint_satisfied(H,h,x);
    prob_constraint_satisfied_u = find_prob_constraint_satisfied(G,g,u);
    fprintf("P(Hx<=h) is %f and 1-eps_x is %f\n",prob_constraint_satisfied_x,1-eps_x);
    fprintf("P(Gu<=g) is %f and 1-eps_u is %f\n",prob_constraint_satisfied_u,1-eps_u);
    v
    eta_x_history_list(theta)=eta_x_history;
    eta_u_history_list(theta)=eta_u_history;
    w_list(theta)=w;
    x_list(theta)=x;
    z_list(theta)=z;
    u_list(theta)=u;
    e_list(theta)=e;
    v_list(theta)=v;
end
hold on
for theta=theta_list
    plot_trajectory(x_list(theta),num2str(theta));
end
legend
hold off
figure;
hold on
for theta = theta_list
plot(u_list(theta),'DisplayName',num2str(theta));
end
hold off
legend
