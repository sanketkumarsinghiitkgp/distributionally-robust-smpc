function [radii_x_history, radii_u_history, eta_x_history, eta_u_history] = find_eta_list(A,B,K,eps_x, eps_u, H,G,e0,N,w_observed,M,w_support,sigma,theta)
    w_history = [];
    radii_x_history = {};
    radii_u_history = {};
    e_support = w_support;
    eta_x_history = {};
    eta_u_history = {};
    e = e0;
    sz = size(H);
    p = sz(1);
    sz = size(G);
    q = sz(1);
    Acl = A+B*K;
    for i=1:N
        
        %e_list = [e0']; COMMENTED OUT
        e_list = [];
        sz = size(w_observed);%pehle it was w_history
        %%%CONTINUE DEBUGGING FROM HERE
        for j=1:M
            e_list = [e_list; (Acl*e+w_observed(j,:)')']; %pehle it was w_history
        end
        eta_x = [];
        radii_x = [];
        sz = size(e)
        support_C = e_support.H(:,1:sz(1));
        support_h = e_support.H(:,end);
        for j=1:p
            a = H(j,:)';
            retval = find_eta(a,e_list,eps_x/p,sigma,theta,support_C, support_h);
            radii_x = [radii_x, retval.radius];
            eta_x = [eta_x retval.solution];
        end
        eta_u = [];
        radii_u = [];
        for j=1:q
            a = (G(j,:)'*K)';
            retval = find_eta(a,e_list,eps_u/q,sigma,theta,support_C, support_h);
            radii_u = [radii_u, retval.radius];
            eta_u = [eta_u retval.solution];
        end
        radii_x_history{i} = radii_x;
        eta_x_history{i} = eta_x;
        radii_u_history{i} = radii_u;
        eta_u_history{i} = eta_u;
        %if i==N
        %    break
        %end
        w = w_observed(i,:);
        w = w';
        w_history = [w_history; w'];
        e = Acl*e+w;
        e_support = Acl*e_support+w_support;
    end
    %e_list = [e0'];
    e_list = [];
    sz = size(w_observed); %pehle it was w_history
    for j=1:M
        e_list = [e_list; (Acl*e+w_observed(j,:)')']; %pehle it was w_observed
    end
    eta_x = [];
    radii_x = [];
    sz = size(e)
    support_C = e_support.H(:,1:sz(1));
    support_h = e_support.H(:,end);
    for j=1:p
        a = H(j,:)';
        retval = find_eta(a,e_list,eps_x/p,sigma,theta,support_C, support_h);
        radii_x = [radii_x, retval.radius];
        eta_x = [eta_x retval.solution];
    end
    radii_x_history{N+1} = radii_x;
    eta_x_history{N+1} = eta_x;
end