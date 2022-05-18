function f = isProblemFeasible(x0,z0,e0,A,B,K,H,h,G,g,P,Q,R,X,W,U,M_x,M_u,M,T,sigma,theta,eps_x,eps_u,w_support,S0,w_samples,method)
    f=1;
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
            v_cur = find_v_drsmpc_precompute_eta(P,Q,R,A,B,H,h,G,g,x_cur,T,eta_x_history,eta_u_history);
            if size(v_cur,2) == 2
                f = 0;
                return;
            end
            
            % save x, z, u, v, e, prob_c... etc for later plotting
        case "scenario_precompute_eta_smpc"
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
            v_cur = find_v_drsmpc_precompute_eta(P,Q,R,A,B,H,h,G,g,x_cur,T,eta_x_history,eta_u_history);    
            if size(v_cur,2) == 2
                f = 0;
                return;
            end
            % save x, z, u, v, e, prob_c... etc for later plotting.
        case "precompute_eta_drsmpc"
            x = [x0'];
            u = [];
            v = [];
            z = [z0'];
            e = [e0'];
            [~, ~, eta_x_history, eta_u_history] = find_eta_list(A,B,K,eps_x,eps_u,H,G,x0-x0,T,w(1:max(M,T+1),:),M,w_support,sigma,theta);
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
            v_cur = find_v_drsmpc_precompute_eta(P,Q,R,A,B,H,h,G,g,x_cur,T,eta_x_history,eta_u_history);
            if size(v_cur,2) == 2
                f = 0;
                return;
            end
            % save x, z, u, v, e, prob_c... etc for later plotting.
        case "precompute_eta_2_drsmpc"
            x = [x0'];
            u = [];
            v = [];
            z = [z0'];
            e = [e0'];
            eta_x_history = find_eta_list_2(A, B, K, H, T+1, max(M,T+1), e0, w_samples,sigma,theta,eps_x,w_support);
            eta_u_history = find_eta_list_2(A, B, K, G, T, max(M,T+1),  e0, w_samples,sigma,theta,eps_u,w_support);
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
            v_cur = find_v_drsmpc_precompute_eta(P,Q,R,A,B,H,h,G,g,x_cur,T,eta_x_history,eta_u_history);
            if size(v_cur,2) == 2
                f = 0;
                return;
            end
            
        case "robust_mpc"
            x = [x0'];
            u = [];
            v = [];
            z = [z0'];
            e = [e0'];
            Acl = A+B*K;
            S_cur = S0;
            S = [S_cur];
            for i=1:1+T
                S_cur = Acl*S_cur+W;
                S = [S, S_cur];
            end
            
            %store X_2_ul, X_2_ll, U_ul and U_ll for plotting later
            z_cur = z0;
            x_cur = x0;
            e_cur = e0;
            v_cur = find_v_rmpc(X,S(1:end),A,B,K,U,x_cur,T,Q,R);
            if size(v_cur,2) == 2
                f = 0;
                return;
            end
        otherwise
            disp("Invalid method: "+method);
    f=1;
end