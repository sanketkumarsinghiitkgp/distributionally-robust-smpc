px_list = [0.95];
hold on;
xlabel('k');
global_e_history = {};
global_z_history = {};
global_u_history = {};
global_radius_history = {};
for i=1:length(px_list)
    px = px_list(i);
    A = [ 1 1; 0 1 ];
    B = [0;1];
    K = [-0.4 -1.2];
    a = [1; 1];
    z = [0; 0];
    v = 0;
    e = [0; 0];
    mu = [0; 0];
    sigma = 1;
    support_a = 1;
    N = 100;
    [z_history,e_history, u_history, radius_history, eta_history] = simulate(A,B,K,a,z,v,e,mu,sigma,support_a,N,px);
    global_e_history{i} = e_history;
    global_z_history{i}= z_history;
    global_u_history{i} = u_history;
    global_radius_history{i} = radius_history;
end
reduced_e_history = reducer(global_e_history{1});
display(reduced_e_history)
hold on
plot(reduced_e_history,'DisplayName','e(1)+e(2)');
plot(eta_history, 'DisplayName','eta');
legend
hold off;
function f = reducer(arr)
    
    shape = size(arr)
    f=zeros(shape(2),1);
    for i=1:shape(2)
        for j=1:shape(1)
            f(i) = f(i) + arr(j,i);
        end
    end
end