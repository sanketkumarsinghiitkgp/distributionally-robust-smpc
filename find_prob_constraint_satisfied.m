function p = find_prob_constraint_satisfied(H,h,x)
    p = 0;
    sz = size(x);
    for i=1:sz(1)
        if H*x(i,:)' <= h
            p = p+1;
        end
    end
    p = p/sz(1);
end