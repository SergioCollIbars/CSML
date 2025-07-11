function [planetParams, poleParams, Kaula, r, coeffs, I] = loadPlanet(option)
    if(option == "Earth") % select Earth
        path = "HARMCOEFS_EARTH_1.txt";
        [Cnm, Snm, Re] = readCoeff(path);
        path = "SIGMACOEFS_EARTH_1.txt";
        GM = 3.986004418E14;
        n_max  = 100;
        normalized = 1;
        W = 2 * pi / (24*3600);     % Rotation ang. vel   [rad/s]
        W0 = 0;                     % Initial asteroid longitude
        RA = deg2rad(15);           % Right Ascension     [rad]
        DEC = 0;                    % Declination         [rad]

        % orbit radius
        r = Re + 250E3;         % [m] 

        % Kaula rule
        [Kaula] = compute_Kaula(n_max, 1E-5);

        % Inertia matrix
        I = [153,-23,-6;-23,2691,1;-6,1,2653];
    elseif(option == "Bennu")   % select Bennu
        path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
        [Cnm, Snm, Re] = readCoeff(path);
        GM = 5.2;
        n_max  = 6;
        normalized = 1;
        W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
        W0 = 0;                   % Initial asteroid longitude
        RA = deg2rad(86.6388);    % Right Ascension     [rad]
        DEC = deg2rad(-65.1086);  % Declination         [rad]
        
        % orbit radius
        r = 0.35E3;         % [m] 

        Kaula = compute_Kaula(n_max, 0.813);

        % Inertia matrix
        I = [153,-23,-6;-23,2691,1;-6,1,2653];
    end

    [Nc, Ns, ~] = count_num_coeff(n_max); 
    poleParams   = [W, W0, RA, DEC];
    planetParams = [GM, Re, n_max, normalized];
    coeffs       = mat2list(Cnm, Snm, Nc, Ns);
end

