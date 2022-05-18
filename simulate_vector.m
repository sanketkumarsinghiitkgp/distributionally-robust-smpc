function [z_history,e_history, u_history, radius_history, eta_history] = simulate_vector(A,B,K,a,z,v,e,mu,sigma,support_a,N,px)
    SIG = [sigma 0; 0 sigma];
    AA = [0 1; 1 0 ; 0 -1 ; -1 0];
    BB = support_a*[1; 1; 1; 1];
    w_history = []; %O(N2) method, optimize if this is causes ineffeciency
    z_history = [];
    e_history = [];
    u_history = [];
    radius_history = [];
    eta_history = [];
    for i = 1:N
        w = rmvnrnd(mu,SIG,1,AA,BB);
        w = w';
        w_history = [w_history; w'];
        e_list = [];
        sz = size(w_history);
        for j = 1:sz(1)
            e_list = [e_list; ((A+B*K)*e+w_history(j,:))'];
        end
        [radius_list, solution_list] = find_eta_vector(a , e_list , 1-px, sigma,1e-3);
        z_history = [z_history [z(1); z(2)]];
        e_history = [e_history [e(1); e(2)]];
        radius_history = [radius_history radius_list];
        eta_history = [eta_history solution_list];
        [z, e, u] = compute_next_state(A,B,z,v,e,w,K);
        u_history = [u_history u];
    end
    eta_history
    %legend_name = strcat('px=',num2str(px))
    %plot(radius_history,'DisplayName',legend_name);
    %plot(eta_history, 'DisplayName',legend_name);
    
end