theta_list = [1e-3];
px_list = [0.5,0.7,0.95];
figure(1)
title("Fraction <= eta")
figure(2)
title("Eta history")
for theta=theta_list
    for px=px_list
        [ensemble_average_fleq_eta_list, ensemble_average_eta_history] = ensemble_average_fleq_eta_list_generator(theta,px);
        figure(1)
        hold on
        legendName = "theta = "+theta+", px = "+px;
        plot(ensemble_average_fleq_eta_list, 'DisplayName',legendName);
        legend
        figure(2)
        hold on
        plot(ensemble_average_eta_history,'DisplayName',legendName);
        legend
    end 
end
legend
hold off
