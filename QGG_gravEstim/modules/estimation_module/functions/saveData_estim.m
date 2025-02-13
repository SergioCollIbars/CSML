function saveData_estim(estimObj, name, count, CS_err)
    %%                       SAVE DATA FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 04/11/2022
    %
    %   Description: This function save object data into txt file.
    %
    %   Input:
    %       estimObj: attitude object
    %       name: file name
    %       count: number of iterations
    %
    %   Output: null
    %
    % --------------------------------------------------------------------%

    % obtain complete file path
    path1 = "./data_files/" + name + "_1";
    path2 = "./data_files/" + name + "_2";
    path3 = "./data_files/" + name + "_3";
    path4 = "./data_files/" + name + "_4";
    path5 = "./data_files/" + name + "_5";
    path6 = "./data_files/" + name + "_6";
    
    % obtain variables to save
    t = estimObj.TIME;

    X0 = estimObj.Xfilter(:, 1);
    Xt = estimObj.Xfilter;
    P = estimObj.Pj0(:, :, count);
    Np = length(P);
    Nc = estimObj.Nc;
    Ns = estimObj.Ns;
    posf = estimObj.postfit;
    pref = estimObj.prefit;
    N = ones(Np, 1) * estimObj.n_max;
    sigma_t = estimObj.sigma_t;

    % add zeros due to extra values
    if(Np > (Nc + Ns))
        NE = Np - (Nc + Ns);
        cs_err = [CS_err; zeros(NE, 1)]; 
    else
        cs_err = CS_err;
    end

    % covariance value
    sigma = sqrt(diag(P));

    % postfit values
    dP11 = posf(1, :)';
    dP12 = posf(2, :)';
    dP13 = posf(3, :)';

    dP22 = posf(4, :)';
    dP23 = posf(5, :)';

    dP33 = posf(6, :)';

    % prefit values
    dp11 = pref(1, :)';
    dp12 = pref(2, :)';
    dp13 = pref(3, :)';

    dp22 = pref(4, :)';
    dp23 = pref(5, :)';

    dp33 = pref(6, :)';

    % create table
    T = table(t, dP11, dP12, dP13, dP22, dP23, ...
        dP33, ...
        'VariableNames', {'TIME','dP11', 'dP12', 'dP13',...
        'dP22', 'dP23', 'dP33'});

    % save table
    writetable(T, path1);

    % create table
    T = table(t, dp11, dp12, dp13, dp22, dp23, ...
        dp33, ...
        'VariableNames', {'TIME','dp11', 'dp12', 'dp13',...
        'dp22', 'dp23', 'dp33'});

    % save table
    writetable(T, path3);

    % create table
    T = table(X0, sigma, cs_err, N,...
        'VariableNames', {'X', 'sigma', 'CS_err', 'n_max'});

     % save table
    writetable(T, path2);

    % create table
    T = table(reshape(P, [Np^2, 1]),...
        'VariableNames', {'P'});

    % save table
    writetable(T, path4);

    % create table
    T = table(sigma_t,...
        'VariableNames', {'sigma_t'});

    % save table
    writetable(T, path5);

     % create table
    T = table(Xt,...
        'VariableNames', {'X_t'});

    % save table
    writetable(T, path6);
end