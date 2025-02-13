function saveData_orbit(OrbitObj, name)
    %%                       SAVE DATA FUNCTION
    % ------------------------------------------------------------------- %
    %   Author: Sergio Coll Ibars
    %
    %   Date: 04/11/2022
    %
    %   Description: This function save object data into txt file.
    %
    %   Input:
    %       OrbitObj: orbit object
    %       name: file name
    %
    %   Output: null
    %
    % --------------------------------------------------------------------%

    % obtain complete file path
    path = "./data_files/" + name;

    % obtain variables to save
    t = OrbitObj.t';

    ri_x = OrbitObj.ri(1, :)';
    ri_y = OrbitObj.ri(2, :)';
    ri_z = OrbitObj.ri(3, :)';
    
    rb_x = OrbitObj.rb(1, :)';
    rb_y = OrbitObj.rb(2, :)';
    rb_z = OrbitObj.rb(3, :)';

    vi_x = OrbitObj.vi(1, :)';
    vi_y = OrbitObj.vi(2, :)';
    vi_z = OrbitObj.vi(3, :)';

    vb_x = OrbitObj.vb(1, :)';
    vb_y = OrbitObj.vb(2, :)';
    vb_z = OrbitObj.vb(3, :)';

    e = OrbitObj.alpha(1, :)';
    h = OrbitObj.alpha(2, :)';
    a = OrbitObj.alpha(3, :)';
    rho = OrbitObj.alpha(4, :)';
    f = OrbitObj.alpha(5, :)';
    i = OrbitObj.alpha(6, :)';
    Omega = OrbitObj.alpha(7, :)';
    omega = OrbitObj.alpha(8, :)';

    % create table
    T = table(t, ri_x, ri_y, ri_z, rb_x, rb_y, rb_z, vi_x, vi_y, vi_z, ...
        vb_x, vb_y, vb_z, e, h, a, rho, f, i, Omega, omega, ...
        'VariableNames', ...
        {'TIME', 'ri_x', 'ri_y', 'ri_z', 'rb_x', 'rb_y', 'rb_z', ...
        'vi_x', 'vi_y', 'vi_z', 'vb_x', 'vb_y', 'vb_z', 'e', 'h', 'a', ...
        'rho', 'f', 'i', 'Omega', 'omega'});

    % save table
    writetable(T, path);

    % save rotation matrix
    BN = OrbitObj.BN;
    T = table(BN, 'VariableNames', {'BN'});
    writetable(T, './data_files/ECI2Body');

    SAR_N = OrbitObj.SAR_N;
    T = table(SAR_N, 'VariableNames', {'SAR_N'});
    writetable(T, './data_files/N2SAR');

    ACAF_N = OrbitObj.ACAF_N;
    T = table(ACAF_N, 'VariableNames', {'ACAF_N'});
    writetable(T, './data_files/N2ACAF');
end