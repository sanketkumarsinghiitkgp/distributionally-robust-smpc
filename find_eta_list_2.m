function eta_x_history=find_eta_list_2(A, B, K, H, N, M, e0, w_samples,sigma,theta,eps_x,w_support)
    eta_x_history = {};
    Acl = A+B*K;
    p = size(H,1);
    e_list = [];
    for i=1:M
        e_list = [e_list; e0'];
    end
    n_x = size(e0,1);
    e_support = Polyhedron('Ae',eye(n_x), 'be',e0);
    for i=1:N %pass N+1 for eta_x
        eta_x = [];
        support_C = e_support.H(:,1:size(e0,1));
        support_h = e_support.H(:,end);
        for j=1:p
            a = H(j,:)';
            retval = find_eta(a,e_list,eps_x/p,sigma,theta,support_C, support_h);
            eta_x = [eta_x retval.solution];
        end
        eta_x_history{i} = eta_x;
        e_support = Acl*e_support+w_support;
        for j=1:M
            e_list(j,:) = (Acl*e_list(j,:)'+w_samples{j}(i,:)')';
        end
    end
end