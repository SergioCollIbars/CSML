function [] = save_results(X_f, P_f, folder)
    % Save current results
    path = "data/" + folder + "/state.mat";
    save(path, 'X_f');

    path = "data/" + folder + "/covariance.mat";
    save(path, 'P_f');
end

