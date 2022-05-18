function f = find_eta_scenario_2(a,ii,Acl,M,w_samples)
    yalmip('clear');
    x = sdpvar(1);
    Objective = x;
    Constraints = [];
    for i=1:M
        s = 0;
        A_power = 1;
        for j = 1:ii-1
            s = s + A_power*w_samples{i}(ii-j,:)'; %can be optimized powers can be calculated beforehand
            A_power = A_power*Acl;
        end
        Constraints = [Constraints a' * s <= x];
    end
    sol = optimize(Constraints,Objective,sdpsettings('solver', 'mosek'));
    solution = 0; %CHANGED FROM -1
    if sol.problem == 0
        solution = value(x);
    else
        display('Something went wrong in find_eta_scenario_2');
        abc = input('Continue, ignoring the error in find_eta_scenario_2');
        sol.info;
        yalmiperror(sol.problem)
    end
    f=solution;
end