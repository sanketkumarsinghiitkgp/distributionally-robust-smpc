function f = find_v_drsmpc(P,Q,R,A,B,K,H,h,G,g,z_k,e_k,N, M_orig, w_samples, sigma, theta, eps_x, eps_u, w_support)
    eta_x_history = find_eta_list_2(A, B, K, H, N+1, M_orig, e_k, w_samples,sigma,theta,eps_x,w_support);
    eta_u_history = find_eta_list_2(A, B, K, G, N, M_orig,  e_k, w_samples,sigma,theta,eps_u,w_support);
    yalmip('clear');
    n_u = size(G,2);
    n_x = size(H,2);
    v = sdpvar(N,n_u,'full');
    z0 = sdpvar(n_x,1,'full');
    z_var = [ z0 ];
    
    for i=1:N
        z_last = z_var(:,end); 
        z_var = [z_var A*z_last+B*v(i,:)'];
    end
    objective = z_var(:,N+1)'*P*z_var(:,N+1);
    for i=1:N
        zz = z_var(:,i);
        objective = objective + (zz'*Q*zz+v(i,:)*R*v(i,:)');
    end
    z_var(:,1)
    
    constraints = [z_var(:,1) == z_k];
    for i=2:N+1 %changed
        zz = z_var(:,i);
        constraints = [constraints, H*zz<=h-eta_x_history{i}']; % eta_x_history length should be increased
    end
    for i=1:N
        constraints = [constraints, G*v(i,:)' <= g-eta_u_history{i}'];
    end
    sol = optimize(constraints,objective,sdpsettings('solver', 'mosek'));
    solution = 0; %setting this as -1 was causing problems
    
    if sol.problem == 0
        solution = value(v);
        f = {solution(1),eta_x_history{2},eta_u_history{1},eta_x_history{6},eta_u_history{6}};
    else
        f = {-1, -1};
        display('Something went wrong in find_v_drsmpc');
        sol.info;
        yalmiperror(sol.problem)
    end
end