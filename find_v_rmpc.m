function f=find_v_rmpc(X,S,A,B,K,U,x_k,T,Q,R)
    Acl = A+B*K;
    d = size(x_k);
    d = d(1);
    z_var0 = sdpvar(d,1,'full');
    z_var = [z_var0];
    n_u = size(B,2);
    v = sdpvar(T,n_u,'full');
    for i=1:T
        z_last = z_var(:,end);
        z_var = [z_var A*z_last+B*v(i,:)'];
    end
    objective = z_var(:,T+1)'*Q*z_var(:,T+1);
    for i=1:T
        zz = z_var(:,i);
        objective = objective + (zz'*Q*zz+v(i,:)*R*v(i,:)');
    end
    constraints = [];
    for i=2:T %changed
        Z = X-S(i);
        sz = size(Z.b);
        for j=1:sz(1)
            constraints = [constraints, Z.A(j,:)*z_var(:,i)<=Z.b(j)];
        end
    end
    for i=1:T %changed
        V = U-K*S(i);
        sz = size(V.b);
        for j=1:sz(1)
            constraints = [constraints, V.A(j,:)*v(i,:)'<=V.b(j)];
        end
    end
    Z = X-S(T+1);
    %Z.plot
    sz = size(Z.b);
    for j=1:sz(1)
        constraints = [constraints, Z.A(j,:)*z_var(:,T+1)<=Z.b(j)];
    end
    constraints = [constraints, z_var(:,1)==x_k];
    sol = optimize(constraints,objective,sdpsettings('solver', 'mosek'));
    solution = 0; %setting this as -1 was causing problems
    
    if sol.problem == 0
        solution = value(v);
        f = solution(1);
    else
        
        f = [-1,-1];
        display('Something went wrong in find_v_rmpc');
        sol.info;
        constraints
        yalmiperror(sol.problem)
    end
end