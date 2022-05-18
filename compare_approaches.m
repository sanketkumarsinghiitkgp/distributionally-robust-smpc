function [x_list,u_list,z_list,e_list,v_list,X_1_ul_list,X_1_ll_list,X_2_ul_list,X_2_ll_list,U_ul_list,U_ll_list,prob_constraint_satisfied_x_list,prob_constraint_satisfied_u_list] = compare_approaches(A,B,K,G,g,H,h,x0,e0,P,Q,R,eps_x, eps_u,sigma,w_support,num_iter,M,T,method_list, beta_scenario,theta,X,W,U,S0,w_distribution,sigma_multiplier,c_theta)
    %method_list = ["precompute_eta_drsmpc" "robust_mpc" "scenario_precompute_eta_smpc" "scenario_precompute_eta_2_smpc"]
    %remember double quotes
    %TODO a way to write the code so that it gives bounds for any kind of
    %constraints
    %{
    Typical  values of X.
    X = Polyhedron('A',[1 0;-1 0;0 1;0 -1],'b',[4;4;1;1]);
    W = Polyhedron('A',[1 0;-1 0;0 1;0 -1],'b',[0.2;0.2;0.2;0.2]);
    U = Polyhedron('A',[1;-1],'b',[0.7;0.7]);
    %}
    
    
    %rng('default')  % For reproducibility
    
    n_x = size(A,2);
    n_u = size(B,2);
    z0 = x0-e0;
    x_list = containers.Map('KeyType','char','ValueType','any');
    z_list = containers.Map('KeyType','char','ValueType','any');
    u_list = containers.Map('KeyType','char','ValueType','any');
    e_list = containers.Map('KeyType','char','ValueType','any');
    v_list = containers.Map('KeyType','char','ValueType','any');
    X_1_ul_list = containers.Map('KeyType','char','ValueType','any');
    X_1_ll_list = containers.Map('KeyType','char','ValueType','any');
    X_2_ul_list = containers.Map('KeyType','char','ValueType','any');
    X_2_ll_list = containers.Map('KeyType','char','ValueType','any');
    U_ul_list = containers.Map('KeyType','char','ValueType','any');
    U_ll_list = containers.Map('KeyType','char','ValueType','any');
    prob_constraint_satisfied_x_list = containers.Map('KeyType','char','ValueType','any');
    prob_constraint_satisfied_u_list = containers.Map('KeyType','char','ValueType','any');
    %% TODO add asserts to check the size of all matrices
    %% if anything in your method requires more than T+1 samples of w, modify M here
    
    M_x = M;
    M_u = M;
    p = size(H,1);
    q = size(G,1);
    % no of cases stuff must be sorted.
    %if (any(method_list == "scenario_precompute_eta_smpc") || any(method_list == "scenario_precompute_eta_2_smpc"))
    %    M_x = ceil((2/(eps_x/p))*(log(1/beta_scenario)+n_x)); %eps_x/p or eps_x
    %    M_u = ceil((2/(eps_u/q))*(log(1/beta_scenario)+n_u));
    %    M = max(M, M_x);
    %    M = max(M, M_u);
    %end
    %% Generate w samples
    w = [];
    w_samples = {};
    num_samples = M;
    %if (1 == 1)
    if (w_distribution == "RMPC_example")
        for i=1:num_samples
            w_cur = [];
            for j=1:T+1
                w_cur = [w_cur; w_generator()];
            end
            w_samples{i} = w_cur;
        end
        for i = 1:num_iter
            w = [w; w_generator()];
        end
    elseif (w_distribution == "SMPC_reachability_truncated_example")
        mu = [0;0];
        Sigma = sigma_multiplier*[0.01 0; 0 1];
        AA = w_support.A;
        bb = w_support.b;
        for i=1:num_samples
            w_cur = [];
            for j=1:T+1
                w_cur = [w_cur; w_generator("rmvnrnd",mu,Sigma,1,AA,bb)];
            end
            w_samples{i} = w_cur;
        end
        for i = 1:M+num_iter
            w = [w; w_generator("rmvnrnd",mu,Sigma,1,AA,bb)];
        end
    elseif (w_distribution == "RMPC_gaussian_example")
        mu = [0;0];
        Sigma = sigma_multiplier*[0.01 0; 0 1];
        AA = w_support.A;
        bb = w_support.b;
        for i=1:num_samples
            w_cur = [];
            
            for j=1:T+1
                w_cur = [w_cur; w_generator("rmvnrnd",mu,Sigma,1,AA,bb)];
            end
            w_samples{i} = w_cur;
        end
        for i = 1:num_iter
            w = [w; w_generator("rmvnrnd",mu,Sigma,1,AA,bb)];
        end
    elseif (w_distribution == "Two_mass_spring")
        mu = [0;0;0;0];
        Sigma = sigma_multiplier*eye(4);
        for i=1:num_samples
            w_cur = [];
            for j=1:T+1
                w_cur = [w_cur; w_generator("mvnrnd",mu,Sigma)];
            end
            w_samples{i} = w_cur;
        end
        for i = 1:num_iter
            w = [w; w_generator("mvnrnd",mu,Sigma)];
        end
    elseif (w_distribution == "Two_mass_spring_triangle")
        pd = makedist('Triangular','A',-100*sigma_multiplier,'B',0,'C',100*sigma_multiplier);
        for i=1:num_samples
            w_cur = [];
            for j=1:T+1
                w_cur = [w_cur; random(pd,1,n_x)];
            end
            w_samples{i} = w_cur;
        end
        for i = 1:num_iter
            w = [w; random(pd,1,n_x)];
        end
    end

    w_samples = update_w_samples_cur(w_samples, w(1,:)',T,num_samples);
                    
%%%%%%%
%SUBJECT TO REPLACEMENT
%%%%%%%
    %% Simulate all methods
    for method=method_list
        wbar = waitbar(0,"Percentage of iterations completed...");
        switch method
            case "scenario_precompute_eta_2_smpc"
                x = [x0'];
                u = [];
                v = [];
                z = [z0'];
                e = [e0'];
                eta_x_history = find_eta_list_scenario_2(A, B, K, H, T+1, M_x, w_samples);
                eta_u_history = find_eta_list_scenario_2(A, B, K, G, T, M_u, w_samples);
                X_1_ll = [];
                X_1_ul = [];
                X_2_ul = [];
                X_2_ll = [];
                U_ul = [];
                U_ll = [];
                for i=1:T+1
                    temp = h-eta_x_history{i}';
                    X_1_ul = [X_1_ul temp(1)];
                    X_1_ll = [X_1_ll -temp(2)];
                    X_2_ul = [X_2_ul temp(3)];
                    X_2_ll = [X_2_ll -temp(4)];
                end
                for i=1:T
                    temp = g-eta_u_history{i}';
                    U_ul = [U_ul temp(1)];
                    U_ll = [U_ll -temp(2)];
                end
                z_cur = z0;
                x_cur = x0;
                e_cur = e0;
                for k = 1:num_iter
                    waitbar(k/num_iter, wbar, "Iterations completed for "+method);
                    w_cur = w(k,:)';
                     x_cur
                    A
                    "yes"
                    v_cur = find_v_drsmpc_precompute_eta(P,Q,R,A,B,H,h,G,g,x_cur,T,eta_x_history,eta_u_history);
                    if(size(v_cur,2)==2)
                        "Infeasibility"
                        abc=input("abc");
                        v_cur = 0;
                    end
                    u = [u; (K*e_cur+v_cur)'];
                    v = [v;(v_cur)'];
                    [new_z, new_e, ~] = compute_next_state(A,B,z_cur,v_cur,e_cur,w_cur,K); %new_u must be removed, it is useless.
                    z_cur = new_z;
                    e_cur = new_e;
                    x_cur = z_cur+e_cur;
                    x = [x; x_cur'];
                    z = [z; z_cur'];
                    e = [e; e_cur'];
                end
                prob_constraint_satisfied_x = find_prob_constraint_satisfied(H,h,x);
                prob_constraint_satisfied_u = find_prob_constraint_satisfied(G,g,u);
                X_1_ul_list(method) = X_1_ul;
                X_1_ll_list(method) = X_1_ll;
                X_2_ul_list(method) = X_2_ul;
                X_2_ll_list(method) = X_2_ll;
                U_ul_list(method) = U_ul;
                U_ll_list(method) = U_ll;
                x_list(method) = x;
                u_list(method) = u;
                z_list(method) = z;
                e_list(method) = e;
                v_list(method) = v;
                prob_constraint_satisfied_u_list(method) = prob_constraint_satisfied_u;
                prob_constraint_satisfied_x_list(method) = prob_constraint_satisfied_x;
                
                % save x, z, u, v, e, prob_c... etc for later plotting
            delete(wbar)
            %DEFUNCT METHOD
            case "scenario_precompute_eta_smpc" %DEFUNCT
                x = [x0'];
                u = [];
                v = [];
                z = [z0'];
                e = [e0'];
                %replace this hack asap
                [eta_x_history, ~] = find_eta_list_scenario(A, B, K, eps_x, eps_u, H, G, e0, T, w(1:max(M,T+1),:),M, w_support, sigma, theta);
                [~, eta_u_history] = find_eta_list_scenario(A, B, K, eps_x, eps_u, H, G, e0, T, w(1:max(M,T+1),:),M, w_support, sigma, theta);
                X_1_ll = [];
                X_1_ul = [];
                X_2_ul = [];
                X_2_ll = [];
                U_ul = [];
                U_ll = [];
                for i=1:T+1
                    temp = h-eta_x_history{i}';
                    X_1_ul = [X_1_ul temp(1)];
                    X_1_ll = [X_1_ll -temp(2)];
                    X_2_ul = [X_2_ul temp(3)];
                    X_2_ll = [X_2_ll -temp(4)];
                end
                for i=1:T
                    temp = g-eta_u_history{i}';
                    U_ul = [U_ul temp(1)];
                    U_ll = [U_ll -temp(2)];
                end
                z_cur = z0;
                x_cur = x0;
                e_cur = e0;
                for k = 1:num_iter
                    waitbar(k/num_iter, wbar, "Iterations completed for "+method);
                    w_cur = w(k,:)';
                    v_cur = find_v_drsmpc_precompute_eta(P,Q,R,A,B,H,h,G,g,x_cur,T,eta_x_history,eta_u_history);
                    if(size(v_cur,2)==2)
                        "Infeasibility"
                        abc=input("abc");
                        v_cur = 0;
                    end
                    u = [u; (K*e_cur+v_cur)'];
                    v = [v;(v_cur)'];
                    [new_z, new_e, ~] = compute_next_state(A,B,z_cur,v_cur,e_cur,w_cur,K); %new_u must be removed, it is useless.
                    z_cur = new_z;
                    e_cur = new_e;
                    x_cur = z_cur+e_cur;
                    x = [x; x_cur'];
                    z = [z; z_cur'];
                    e = [e; e_cur'];
                end
                prob_constraint_satisfied_x = find_prob_constraint_satisfied(H,h,x);
                prob_constraint_satisfied_u = find_prob_constraint_satisfied(G,g,u);
                X_1_ul_list(method) = X_1_ul;
                X_1_ll_list(method) = X_1_ll;
                X_2_ul_list(method) = X_2_ul;
                X_2_ll_list(method) = X_2_ll;
                U_ul_list(method) = U_ul;
                U_ll_list(method) = U_ll;
                x_list(method) = x;
                u_list(method) = u;
                z_list(method) = z;
                e_list(method) = e;
                v_list(method) = v;
                prob_constraint_satisfied_u_list(method) = prob_constraint_satisfied_u;
                prob_constraint_satisfied_x_list(method) = prob_constraint_satisfied_x;
                
                % save x, z, u, v, e, prob_c... etc for later plotting.
                delete(wbar)
            %DEFUNCT METHOD
            case "precompute_eta_drsmpc"
                x = [x0'];
                u = [];
                v = [];
                z = [z0'];
                e = [e0'];
                [~, ~, eta_x_history, eta_u_history] = find_eta_list(A,B,K,eps_x,eps_u,H,G,x0-x0,T,w(1:M_orig,:),M_orig,w_support,sigma,theta);
                X_1_ll = [];
                X_1_ul = [];
                X_2_ul = [];
                X_2_ll = [];
                U_ul = [];
                U_ll = [];
                for i=1:T+1
                    temp = h-eta_x_history{i};
                    X_1_ul = [X_1_ul temp(1)];
                    X_1_ll = [X_1_ll -temp(2)];
                    X_2_ul = [X_2_ul temp(3)];
                    X_2_ll = [X_2_ll -temp(4)];
                end
                for i=1:T
                    temp = g-eta_u_history{i};
                    U_ul = [U_ul temp(1)];
                    U_ll = [U_ll -temp(2)];
                end
                %store X_2_ul, X_2_ll, U_ul and U_ll for plotting later
                z_cur = z0;
                x_cur = x0;
                e_cur = e0;
                for k = 1:num_iter
                    waitbar(k/num_iter, wbar, "Iterations completed for "+method);
                    w_cur = w(k,:)';
                    
                    v_cur = find_v_drsmpc_precompute_eta(P,Q,R,A,B,H,h,G,g,x_cur,T,eta_x_history,eta_u_history);
                    if(size(v_cur,2)==2)
                        "Infeasibility"
                        abc=input("abc");
                        v_cur = 0;
                    end
                    u = [u; (K*e_cur+v_cur)'];
                    v = [v;(v_cur)'];
                    [new_z, new_e, ~] = compute_next_state(A,B,z_cur,v_cur,e_cur,w_cur,K); %new_u must be removed, it is useless.
                    z_cur = new_z;
                    e_cur = new_e;
                    x_cur = z_cur+e_cur;
                    x = [x; x_cur'];
                    z = [z; z_cur'];
                    e = [e; e_cur'];
                end
                prob_constraint_satisfied_x = find_prob_constraint_satisfied(H,h,x);
                prob_constraint_satisfied_u = find_prob_constraint_satisfied(G,g,u);
                X_1_ul_list(method) = X_1_ul;
                X_1_ll_list(method) = X_1_ll;
                X_2_ul_list(method) = X_2_ul;
                X_2_ll_list(method) = X_2_ll;
                U_ul_list(method) = U_ul;
                U_ll_list(method) = U_ll;
                x_list(method) = x;
                u_list(method) = u;
                z_list(method) = z;
                e_list(method) = e;
                v_list(method) = v;
                prob_constraint_satisfied_u_list(method) = prob_constraint_satisfied_u;
                prob_constraint_satisfied_x_list(method) = prob_constraint_satisfied_x;
                delete(wbar)
                % save x, z, u, v, e, prob_c... etc for later plotting.
            case "precompute_eta_2_drsmpc"
                x = [x0'];
                u = [];
                v = [];
                z = [z0'];
                e = [e0'];
                eta_x_history = find_eta_list_2(A, B, K, H, T+1, M, e0, w_samples,sigma,theta,eps_x,w_support); %M_orig was here
                eta_u_history = find_eta_list_2(A, B, K, G, T, M,  e0, w_samples,sigma,theta,eps_u,w_support);%M_orig was here
                X_1_ll = [];
                X_1_ul = [];
                X_2_ul = [];
                X_2_ll = [];
                U_ul = [];
                U_ll = [];
                for i=1:T+1
                    temp = h-eta_x_history{i};
                    X_1_ul = [X_1_ul temp(1)];
                    X_1_ll = [X_1_ll -temp(2)];
                    X_2_ul = [X_2_ul temp(3)];
                    X_2_ll = [X_2_ll -temp(4)];
                end
                for i=1:T
                    temp = g-eta_u_history{i};
                    U_ul = [U_ul temp(1)];
                    U_ll = [U_ll -temp(2)];
                end
                %store X_2_ul, X_2_ll, U_ul and U_ll for plotting later
                z_cur = z0;
                x_cur = x0;
                e_cur = e0;
                for k = 1:num_iter
                    waitbar(k/num_iter, wbar, "Iterations completed for "+method);
                    M
                    k
                    w_cur = w(k,:)';
                    
                    v_cur = find_v_drsmpc_precompute_eta(P,Q,R,A,B,H,h,G,g,x_cur,T,eta_x_history,eta_u_history);
                    if(size(v_cur,2)==2)
                        "Infeasibility"
                        abc=input("abc");
                        v_cur = 0;
                    end
                    u = [u; (K*e_cur+v_cur)'];
                    v = [v;(v_cur)'];
                    [new_z, new_e, ~] = compute_next_state(A,B,z_cur,v_cur,e_cur,w_cur,K); %new_u must be removed, it is useless.
                    z_cur = new_z;
                    e_cur = new_e;
                    x_cur = z_cur+e_cur;
                    x = [x; x_cur'];
                    z = [z; z_cur'];
                    e = [e; e_cur'];
                end
                prob_constraint_satisfied_x = find_prob_constraint_satisfied(H,h,x);
                prob_constraint_satisfied_u = find_prob_constraint_satisfied(G,g,u);
                X_1_ul_list(method) = X_1_ul;
                X_1_ll_list(method) = X_1_ll;
                X_2_ul_list(method) = X_2_ul;
                X_2_ll_list(method) = X_2_ll;
                U_ul_list(method) = U_ul;
                U_ll_list(method) = U_ll;
                x_list(method) = x;
                u_list(method) = u;
                z_list(method) = z;
                e_list(method) = e;
                v_list(method) = v;
                prob_constraint_satisfied_u_list(method) = prob_constraint_satisfied_u;
                prob_constraint_satisfied_x_list(method) = prob_constraint_satisfied_x;
                delete(wbar)
            case "robust_mpc"
                x = [x0'];
                u = [];
                v = [];
                z = [z0'];
                e = [e0'];
                Acl = A+B*K;
                S_cur = S0;
                S = [S_cur];
                for i=1:T
                    S_cur = Acl*S_cur+W;
                    S = [S, S_cur];
                end
                X_1_ul = [];
                X_1_ll = [];
                X_2_ul = [];
                X_2_ll = [];
                U_ul = [];
                U_ll = [];
                for i=1:T+1
                    Z = X-S(i);
                    V = U-K*S(i);
                    X_1_ul = [X_1_ul Z.b(5)]; %fishy method, relies on MPT preserving the order of inequalities.
                    X_1_ll = [X_1_ll -Z.b(6)];
                    X_2_ul = [X_2_ul Z.b(7)];
                    X_2_ll = [X_2_ll -Z.b(8)];
                    U_ul = [U_ul V.b(3)];
                    U_ll = [U_ll -V.b(4)];
                end
                %store X_2_ul, X_2_ll, U_ul and U_ll for plotting later
                z_cur = z0;
                x_cur = x0;
                e_cur = e0;
                for k = 1:num_iter
                    waitbar(k/num_iter, wbar, "Iterations completed for "+method);
                    w_cur = w(k,:)';
                    v_cur = find_v_rmpc(X,S,A,B,K,U,x_cur,T,Q,R);
                    if(size(v_cur,2)==2)
                        "Infeasibility"
                        abc=input("abc");
                        v_cur = 0;
                    end
                    u = [u; (v_cur+K*e_cur)'];
                    z_cur = A*z_cur+B*v_cur;
                    e_cur = (Acl)*e_cur+w_cur;
                    x_cur = z_cur+e_cur;
                    x = [x;x_cur'];
                    v = [v;v_cur'];
                    z = [z;z_cur'];
                    e = [e;e_cur'];
                end
                prob_constraint_satisfied_x = find_prob_constraint_satisfied(H,h,x);
                prob_constraint_satisfied_u = find_prob_constraint_satisfied(G,g,u);
                X_1_ul_list(method) = X_1_ul;
                X_1_ll_list(method) = X_1_ll;
                X_2_ul_list(method) = X_2_ul;
                X_2_ll_list(method) = X_2_ll;
                U_ul_list(method) = U_ul;
                U_ll_list(method) = U_ll;
                x_list(method) = x;
                u_list(method) = u;
                z_list(method) = z;
                e_list(method) = e;
                v_list(method) = v;
                prob_constraint_satisfied_u_list(method) = prob_constraint_satisfied_u;
                prob_constraint_satisfied_x_list(method) = prob_constraint_satisfied_x;
                % save x, z, u, v, e, prob_c... etc for later plotting.
                delete(wbar)
            case "online_drsmpc"
                x = [x0'];
                u = [];
                v = [];
                z = [z0'];
                e = [e0'];
                z_cur = z0;
                x_cur = x0;
                e_cur = e0;
                w_samples_cur = w_samples;
                
                X_1_ll = [];
                X_1_ul = [];
                X_2_ul = [];
                X_2_ll = [];
                U_ul = [];
                U_ll = [];
                for k = 1:num_iter
                    waitbar(k/num_iter, wbar, "Iterations completed for "+method);
                    %v_cur = find_v_drsmpc_precompute_eta(P,Q,R,A,B,H,h,G,g,x_cur,T,eta_x_history,eta_u_history);
                    temp = find_v_drsmpc(P,Q,R,A,B,K,H,h,G,g,z_cur,e_cur,T, M, w_samples_cur, sigma, theta, eps_x, eps_u, w_support); %M_orig was here
                    v_cur = temp{1};
                    if(size(temp,2)==2)
                        "Infeasibility"
                        abc=input("abc");
                        v_cur = 0;
                    end
                    eta_x_history_1 = temp{2};
                    eta_u_history_1 = temp{3};
                    eta_x_history_6 = temp{4};
                    eta_u_history_6 = temp{5};
                    temp = h-eta_x_history_1;
                    temp2 = h-eta_x_history_6;
                    X_1_ul = [X_1_ul [temp(1);temp2(1)]];
                    X_1_ll = [X_1_ll [-temp(2);-temp2(2)]];
                    X_2_ul = [X_2_ul [temp(3);temp2(3)]];
                    X_2_ll = [X_2_ll [-temp(4);-temp2(4)]];
                    
                    temp = g-eta_u_history_1;
                    temp2 = g-eta_u_history_6;
                    U_ul = [U_ul [temp(1);temp2(1)]];
                    U_ll = [U_ll [-temp(2);-temp2(2)]];

                    u = [u; (K*e_cur+v_cur)'];
                    v = [v;(v_cur)'];
                    w_cur = w(k,:)';
                    [new_z, new_e, ~] = compute_next_state(A,B,z_cur,v_cur,e_cur,w_cur,K); %new_u must be removed, it is useless.
                    w_samples_cur = update_w_samples_cur(w_samples_cur, w_cur,T,num_samples);
                    z_cur = new_z;
                    e_cur = new_e;
                    x_cur = z_cur+e_cur;
                    x = [x; x_cur'];
                    z = [z; z_cur'];
                    e = [e; e_cur'];
                end
                prob_constraint_satisfied_x = find_prob_constraint_satisfied(H,h,x);
                prob_constraint_satisfied_u = find_prob_constraint_satisfied(G,g,u);
                X_1_ul_list(method) = X_1_ul;
                X_1_ll_list(method) = X_1_ll;
                X_2_ul_list(method) = X_2_ul;
                X_2_ll_list(method) = X_2_ll;
                U_ul_list(method) = U_ul;
                U_ll_list(method) = U_ll;
                x_list(method) = x;
                u_list(method) = u;
                z_list(method) = z;
                e_list(method) = e;
                v_list(method) = v;
                prob_constraint_satisfied_u_list(method) = prob_constraint_satisfied_u;
                prob_constraint_satisfied_x_list(method) = prob_constraint_satisfied_x;
                delete(wbar)
            case "online_drsmpc_2"
                x = [x0'];
                u = [];
                v = [];
                z = [z0'];
                e = [e0'];
                z_cur = z0;
                x_cur = x0;
                e_cur = e0;
                w_samples_cur = w_samples;
                
                X_1_ll = [];
                X_1_ul = [];
                X_2_ul = [];
                X_2_ll = [];
                U_ul = [];
                U_ll = [];
                for k = 1:num_iter
                    waitbar(k/num_iter, wbar, "Iterations completed for "+method);
                    %v_cur = find_v_drsmpc_precompute_eta(P,Q,R,A,B,H,h,G,g,x_cur,T,eta_x_history,eta_u_history);
                    temp = find_v_drsmpc_2(P,Q,R,A,B,K,H,h,G,g,z_cur,e_cur,T, M, w_samples_cur, sigma, theta, eps_x, eps_u, w_support,c_theta); %M_orig was here
                    v_cur = temp{1};
                    if(size(temp,2)==2)
                        "Infeasibility"
                        abc=input("abc");
                        v_cur = 0;
                    end
                    eta_x_history_1 = temp{2};
                    eta_u_history_1 = temp{3};
                    temp = h-eta_x_history_1;
                    X_1_ul = [X_1_ul temp(1)];
                    X_1_ll = [X_1_ll -temp(2)];
                    X_2_ul = [X_2_ul temp(3)];
                    X_2_ll = [X_2_ll -temp(4)];
                    
                    temp = g-eta_u_history_1;
                    U_ul = [U_ul temp(1)];
                    U_ll = [U_ll -temp(2)];

                    u = [u; (K*e_cur+v_cur)'];
                    v = [v;(v_cur)'];
                    w_cur = w(k,:)';
                    [new_z, new_e, ~] = compute_next_state(A,B,z_cur,v_cur,e_cur,w_cur,K); %new_u must be removed, it is useless.
                    w_samples_cur = update_w_samples_cur(w_samples_cur, w_cur,T,num_samples);
                    z_cur = new_z;
                    e_cur = new_e;
                    x_cur = z_cur+e_cur;
                    x = [x; x_cur'];
                    z = [z; z_cur'];
                    e = [e; e_cur'];
                end
                prob_constraint_satisfied_x = find_prob_constraint_satisfied(H,h,x);
                prob_constraint_satisfied_u = find_prob_constraint_satisfied(G,g,u);
                X_1_ul_list(method) = X_1_ul;
                X_1_ll_list(method) = X_1_ll;
                X_2_ul_list(method) = X_2_ul;
                X_2_ll_list(method) = X_2_ll;
                U_ul_list(method) = U_ul;
                U_ll_list(method) = U_ll;
                x_list(method) = x;
                u_list(method) = u;
                z_list(method) = z;
                e_list(method) = e;
                v_list(method) = v;
                prob_constraint_satisfied_u_list(method) = prob_constraint_satisfied_u;
                prob_constraint_satisfied_x_list(method) = prob_constraint_satisfied_x;

            otherwise
                disp("Invalid method: "+method);
        end
    end

end