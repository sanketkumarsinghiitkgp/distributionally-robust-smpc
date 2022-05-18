function objective = find_cost(z,v,P,Q,R)
    N = size(z,1)-1;
    objective = z(N+1,:)*P*z(N+1,:)';
    for i=1:N
        zz = z(i,:)';
        objective = objective + (zz'*Q*zz+v(i,:)*R*v(i,:)');
    end
end