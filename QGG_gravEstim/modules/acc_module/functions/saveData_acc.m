function saveData_acc(AccObj, name)
    %%                       SAVE DATA FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 04/11/2022
    %
    %   Description: This function save object data into txt file.
    %
    %   Input:
    %       AccObj: attitude object
    %       name: file name
    %
    %   Output: null
    %
    % --------------------------------------------------------------------%

    % obtain complete file path
    path = "./data_files/" + name;

    % obtain variables to save
    t = AccObj.t';
    Nt = length(t);
    
    Ub_x = AccObj.dUb(1, :)';
    Ub_y = AccObj.dUb(2, :)';
    Ub_z = AccObj.dUb(3, :)';
    
    % tensor components
    ad_xx = zeros(Nt, 1);
    ad_xy = zeros(Nt, 1);
    ad_xz = zeros(Nt, 1);

    ad_yx = zeros(Nt, 1);
    ad_yy = zeros(Nt, 1);
    ad_yz = zeros(Nt, 1);

    ad_zx = zeros(Nt, 1);
    ad_zy = zeros(Nt, 1);
    ad_zz = zeros(Nt, 1);

    for j = 1:Nt
        up = 3*j;
        down = up - 2;
        accTk = AccObj.accT(down:up, :);

        ad_xx(j) =  accTk(1, 1);
        ad_xy(j) =  accTk(1, 2);
        ad_xz(j) =  accTk(1, 3);

        ad_yx(j) =  accTk(2, 1);
        ad_yy(j) =  accTk(2, 2);
        ad_yz(j) =  accTk(2, 3);

        ad_zx(j) =  accTk(3, 1);
        ad_zy(j) =  accTk(3, 2);
        ad_zz(j) =  accTk(3, 3);
    end


    % create table
    T = table(t, ad_xx, ad_xy, ad_xz, ad_yx, ad_yy, ad_yz, ...
        ad_zx, ad_zy, ad_zz, Ub_x, Ub_y, Ub_z, ...
        'VariableNames', {'TIME','ad_xx', 'ad_xy', 'ad_xz', 'ad_yx',...
        'ad_yy', 'ad_yz', 'ad_zx', 'ad_zy', 'ad_zz', 'Ub_x',...
        'Ub_y', 'Ub_z'});

    % save table
    writetable(T, path);
end