function [ensemble_average_fleq_eta_list, ensemble_average_eta_history] = ensemble_average_fleq_eta_list_generator(theta,px)
    num_realizations = 20;
    [ensemble_average_fleq_eta_list, ensemble_average_eta_history] = fleq_eta_list_generator(theta,px);
    for i=1:(num_realizations-1)
        [fleq_eta_list, eta_history] = fleq_eta_list_generator(theta,px);
        ensemble_average_fleq_eta_list = ensemble_average_fleq_eta_list+fleq_eta_list;
        ensemble_average_eta_history = ensemble_average_eta_history+eta_history;
    end
    ensemble_average_fleq_eta_list = ensemble_average_fleq_eta_list/num_realizations;
    ensemble_average_eta_history = ensemble_average_eta_history/num_realizations;
end