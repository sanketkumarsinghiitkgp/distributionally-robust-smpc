function f = find_eta(a,xi,alpha, sigma,theta,C,h)
    %ask sir, taken beta = alpha
    yalmip('clear');
    sz= size(xi);
    N = sz(1);
    n = sz(2);
    m=n;
    sz = size(C);
    p = sz(1);
    beta = 0.05;
    c = (1+2^0.5)*(1+3^0.5)*3^(3.5-1.0/n)*(2^7)*(sigma^3)*(n^1.5);
    if(theta == -1)
        theta = ((2*n*sigma*sigma*log(1/beta)/N)^0.5)+c*N^(-1/max(n,2));
    end
    x = sdpvar(1);
    pii = sdpvar(N,p,'full');
    Objective = x;
    t = sdpvar(1);
    lmbd = sdpvar(1);
    s = sdpvar(N,1,'full');
    sm = 0;
    for i=1:N
        sm = sm + s(i);
    end
    sm = sm / N;
    Constraints = [ lmbd*theta + sm-t*alpha<=0, lmbd >= 0 ];
    for i=1:N
        Constraints = [Constraints, max(0, -x+t+(a'-pii(i,:)*C)*(xi(i,:)')+pii(i,:)*h)-s(i) <= 0];
        Constraints = [Constraints, pii(i,:)>=0];
        Constraints = [ Constraints, norm(a'-pii(i,:)*C)-lmbd<=0 ]; 
    end
    sol = optimize(Constraints,Objective,sdpsettings('solver', 'mosek'));
    solution = 0; %CHANGED FROM -1
    if sol.problem == 0
        solution = value(x);
    else
        display('Something went wrong in find_eta');
        abc = input('Continue, ignoring the error in find_eta');
        sol.info;
        yalmiperror(sol.problem)
    end
    f.solution = solution;
    f.radius = theta;
end