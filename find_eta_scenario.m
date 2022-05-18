function f = find_eta_scenario(a,xi,alpha, sigma,theta,C,h)
    yalmip('clear');
    sz= size(xi);
    N = sz(1);
    x = sdpvar(1)
    Objective = x;
    Constraints = [];
    for i=1:N
        e_i = xi(i,:)';
        size(e_i)
        size(a')
        Constraints = [Constraints a'*e_i<=x];
    end
    sol = optimize(Constraints,Objective,sdpsettings('solver', 'mosek'));
    solution = 0; %CHANGED FROM -1
    if sol.problem == 0
        solution = value(x);
    else
        display('Something went wrong in find_eta_scenario');
        abc = input('Continue, ignoring the error in find_eta_scenario');
        sol.info;
        yalmiperror(sol.problem)
    end
    f=solution;
end