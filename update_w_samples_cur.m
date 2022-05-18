function updated_w_samples = update_w_samples_cur(w_samples_cur, w_cur,T,num_samples)
    w_samples_cur_flattened = [];
    for j=2:T+1
        w_samples_cur_flattened = [w_samples_cur_flattened; w_samples_cur{1}(j,:)]
    end
    for i=2:num_samples
        w_samples_cur_flattened = [w_samples_cur_flattened; w_samples_cur{i}(T+1,:)]
    end
    w_samples_cur_flattened = [w_samples_cur_flattened; (w_cur')];
    updated_w_samples = {};
    for i=1:num_samples
        updated_w_samples{i} = w_samples_cur_flattened(i:i+T,:);
    end
    
end