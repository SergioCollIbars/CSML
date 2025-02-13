function [Cmat, Smat, R, normalize] = readCoeff(path)
            list = table2array(readtable(path));
            degree = list(1);
            X = list(4:end);
            R = list(2);
            normalize = list(3);

            % count number of coefficients
            Nc = 1;
            for k = 2:degree
                Nc = Nc + k + 1;
            end
            Ns = 0;
            for k = 2:degree
                Ns = Ns + k;
            end
            
            % define matrices
            Cmat = zeros(degree + 1, degree + 1);
            Smat = zeros(degree + 1, degree + 1);
            Cmat(1, 1) = 1;

            n = 2;
            m = 0;
            for j = 2:Nc
                N = n + 1;
                M = m + 1;
       
                Cmat(N, M) = X(j);
                if(m < n)
                    m = m + 1;
                else
                    m = 0;
                    n = n +1;
                end
            end
        
            n = 2;
            m = 0;
            for j = Nc + 1:Ns + Nc
                N = n + 1;
                M = m + 2;
        
                Smat(N, M) = X(j);
                if(m < n - 1)
                    m = m + 1;
                else
                    m = 0;
                    n = n + 1;
                end
           
            end
end
