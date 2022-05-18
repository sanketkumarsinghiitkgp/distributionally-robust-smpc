function [radius_list, solution_list]= find_eta_vector(a,xi,alpha, sigma,theta)
    radius_list = [];
    solution_list = [];
    sz = size(xi);
    n = sz(2);
    for i=1:n
        
        f = find_eta(a,xi(:,i),alpha,sigma,theta);
        radius_list = [radius_list; f.radius];
        solution_list = [solution_list; f.solution];
    end
end