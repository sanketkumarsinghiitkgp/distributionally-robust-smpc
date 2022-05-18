 function [fleq_eta_list,eta_history] = fleq_eta_list_generator(theta,px)
    A = [ 1 1; 0 1 ];
    B = [0;1];
    K = [-0.4 -1.2];
    a = [1; 1];
    z = [0; 0];
    v = 0;
    e = [0; 0];
    mu = [0; 0];
    sigma = 1;
    support_a = 1;
    N = 100;
    SIG = [sigma 0; 0 sigma];
    AA = [0 1; 1 0 ; 0 -1 ; -1 0];
    BB = support_a*[1; 1; 1; 1];
    w_history = []; %O(N2) method, optimize if this is causes ineffeciency
    z_history = [];
    e_history = [];
    u_history = [];
    radius_history = [];
    eta_history = [];
    fleq_eta_list = [];
    num_samples = 100;
    for i = 1:N
        next_e_list = [];
        for j=1:num_samples
            temp_w = rmvnrnd(mu,SIG,1,AA,BB)';
            next_e_list = [next_e_list; ((A+B*K)*e+temp_w)' ];
        end
        w = rmvnrnd(mu,SIG,1,AA,BB);
        w = w';
        w_history = [w_history; w'];
        e_list = [];
        sz = size(w_history);
        for j = 1:sz(1)
            e_list = [e_list; ((A+B*K)*e+w_history(j,:))'];
        end
        retval = find_eta(a , e_list , 1-px, sigma,theta);
        radius = retval.radius;
        z_history = [z_history [z(1); z(2)]];
        e_history = [e_history [e(1); e(2)]];
        radius_history = [radius_history radius];
        eta_history = [eta_history retval.solution];
        [z, e, u] = compute_next_state(A,B,z,v,e,w,K);
        u_history = [u_history u];
        fleq_eta = 0;
        size(next_e_list)
        for j=1:num_samples
            if (next_e_list(j,1)+next_e_list(j,2)<=retval.solution)
                fleq_eta = fleq_eta+1;
            end
        end
        fleq_eta = fleq_eta/num_samples;
        fleq_eta_list = [fleq_eta_list fleq_eta];
    end
end