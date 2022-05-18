
for outer=1:20
    A = [ 1 1; 0 1 ];
    B = [0;1];
    K = [-0.4 -1.2];
    a = [1; 1];
    z = [0; 0];
    v = 0;
    e = [2; 2];
    mu = [0; 0];
    sigma = 1;
    SIG = [sigma 0; 0 sigma];
    AA = [0 1; 1 0 ; 0 -1 ; -1 0];
    BB = [1; 1; 1; 1];
    N = 100;
    eta_history = zeros(N,1);
    w_history = []; %O(N2) method, optimize if this is causes ineffeciency
    e_history = [];
    radius_history = zeros(N,1);
    for i = 1:N
        w = rmvnrnd(mu,SIG,1,AA,BB);
        w = w';
        w_history = [w_history w(1)+w(2)];
        e_history = [e_history e(1)+e(2)];
        e_list = [];
        sz = size(w_history);
        [z, e, u] = compute_next_state(A,B,z,v,e,w,K);
    end
    plot(e_history);
    hold on
end
hold off
