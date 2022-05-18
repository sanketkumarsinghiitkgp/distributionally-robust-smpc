function [new_z, new_e, new_u] = compute_next_state(A,B,z,v,e,w,K)
    new_z = A*z+B*v;
    new_e = (A+B*K)*e+w;
    new_u = K*e+v;
end
