function f = find_v_drsmpc_2(P,Q,R,A,B,K,H,h,G,g,z_k,e_k,N, M, w_samples, sigma, theta, eps_x, eps_u, w_support,c_theta)
    yalmip('clear');
    n_u = size(G,2);
    n_x = size(H,2);
    v = sdpvar(N,n_u,'full');
    z0 = sdpvar(n_x,1,'full');
    Acl = A+B*K;
    z_var = [ z0 ];
    for i=1:N
        z_last = z_var(:,end); 
        z_var = [z_var A*z_last+B*v(i,:)'];
    end
    Hh = H;
    Gg = G;
    p = size(H,1);
    q = size(G,1);
    for i=1:p
        Hh(i,:) = H(i,:)/h(i);
    end
    for i=1:q
        Gg(i,:) = G(i,:)/g(i);
    end
    
    e_list = [];
    for i=1:M
        e_list = [e_list; e_k'];
    end
    e_list_big = {};
    for i=1:N+1
        e_list_big{i} = e_list;
        for j=1:M
            e_list(j,:) = (Acl*e_list(j,:)'+w_samples{j}(i,:)')';
        end
    end
    

    constraints = [];
    Ns_h = M;
    tau_h = sdpvar(p,N+1,'full');
    lmbd_h = sdpvar(p,N+1,'full');
    s_h = sdpvar(p,Ns_h,N+1,'full');
    alpha_h = eps_x/p;
    big_theta = sdpvar(N+1,1,'full');
    for i=1:p
        for t=1:N+1 %actually 0 to N
            sm=0;
            constraints = [constraints lmbd_h(i,t)>=0];
            for j=1:Ns_h
                sm = sm + s_h(i,j,t);
                gamma_i = Hh(i,:)*(z_var(:,t)+e_list_big{t}(j,:)')-1;
                constraints = [constraints max(0,gamma_i+tau_h(i,t))<=s_h(i,j,t)];
                constraints = [constraints norm(Hh(i,:))<=lmbd_h(i,t)];
                constraints = [constraints s_h(i,j,t)];
            end
            constraints = [constraints big_theta(t)>=0];
            constraints = [constraints -alpha_h*tau_h(i,t)+theta*lmbd_h(i,t)+sm/Ns_h<=big_theta];
        end
    end


    Ns_g = M;
    tau_g = sdpvar(q,N,'full');
    lmbd_g = sdpvar(q,N,'full');
    s_g = sdpvar(q,Ns_g,N,'full');
    alpha_g = eps_u/q;
    
    for i=1:q
        for t=1:N %actually 0 to N
            sm=0;
            constraints = [constraints lmbd_g(i,t)>=0];
            for j=1:Ns_g
                sm = sm + s_g(i,j,t);
                gamma_i = Gg(i,:)*(v(t,:)+K*e_list_big{t}(j,:)')-1;
                constraints = [constraints max(0,gamma_i+tau_g(i,t))<=s_g(i,j,t)];
                constraints = [constraints norm(Gg(i,:))<=lmbd_g(i,t)];
                constraints = [constraints s_g(i,j,t)];
            end
            constraints = [constraints -alpha_g*tau_g(i,t)+theta*lmbd_g(i,t)+sm/Ns_g<=0];
        end
    end


    objective = z_var(:,N+1)'*P*z_var(:,N+1)+c_theta*norm(big_theta,inf);
    for i=1:N
        zz = z_var(:,i);
        objective = objective + (zz'*Q*zz+v(i,:)*R*v(i,:)');
    end
    constraints = [z_var(:,1) == z_k];
    sol = optimize(constraints,objective,sdpsettings('solver', 'mosek'));
    solution = 0; %setting this as -1 was causing problems
    
    if sol.problem == 0
        solution = value(v);
        f = {solution(1),-1,-1};
    else
        f = {-1, -1};
        display('Something went wrong in find_v_drsmpc');
        sol.info;
        yalmiperror(sol.problem)
    end


    
end