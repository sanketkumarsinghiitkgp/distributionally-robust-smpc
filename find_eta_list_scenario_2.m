function eta_x_history = find_eta_list_scenario_2(A,B,K,H,N,M,w_samples)
    eta_x_history = {};
    Acl = A+B*K;
    p = size(H,1);
    for i=1:N %pass N+1 for eta_x
        eta_x = [];
        for j=1:p
            a = H(j,:)';
            eta_x = [eta_x find_eta_scenario_2(a,i,Acl,M,w_samples)];
        end
        eta_x_history{i} = eta_x;
    end
end
