function [a, f, ff] = solver_initCond(j, vj, At, STM, gamma, aj)
   scale = 1;
   aj = aj*scale;
   vj = vj*scale;
   % compute alpha matrix
   PHI_next2current = reshape(STM(j+1, 1:end), [6,6])/(reshape(STM(j, 1:end), [6,6]));
   PHI_prev2current = reshape(STM(j-1, 1:end), [6,6])/(reshape(STM(j, 1:end), [6,6]));

   alpha = PHI_next2current - PHI_prev2current;
   betta = PHI_next2current*scale + PHI_prev2current*scale - 2*eye(6,6);

   a21 = alpha(4:6, 1:3);
   a22 = alpha(4:6, 4:6);

   b21 = betta(4:6, 1:3);
   b22 = betta(4:6, 4:6);

   % compute grad accelarion
   gamma_j = reshape(gamma(:, j), [3, 3])*scale;
   gammaDot_j = (1/At).*(reshape(gamma(:, j+1), [3, 3]) - ...
       reshape(gamma(:, j-1), [3, 3]))*scale;
   
   % rotation frame
   Wtilde = [0, -1, 0;1, 0, 0;0, 0, 0];
   Wtilde  = zeros(3,3); % non rotating frame
   
   % compute A matrix
   A = a22 + 2*Wtilde.*At;
   A = b22 - gamma_j*(At^2)/4;

   % compute B matrix
   B = (gamma_j*vj).*At - a21*vj;
   B = gammaDot_j*vj*(At^2)/4 - b21*vj;
    
   % solve matrix equation: Ax = B
   a = A\B;
   
   % test equations. Linearization 1st order
   f1 = alpha * [vj;aj];
   c = gamma_j*vj;
   f2 = [aj;c]*At;
   f = f2 - f1;

   % test eqautions. Linearization 2nd order
   f1 = (betta) * [vj;aj];
   c = gammaDot_j*vj + gamma_j*aj;
   b = gamma_j * vj;
   f2 = [b;c]*At*At/4;
   ff = f2 -f1;
end

