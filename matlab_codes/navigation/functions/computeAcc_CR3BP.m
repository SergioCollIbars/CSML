function [dU] = computeAcc_CR3BP(x)
    m_1 = 5.974E24;  % [Kg]
    m_2 = 7.348E22;  % [Kg]
    mu =  m_2 / (m_1 + m_2); % mass ratio
            
    r1 = sqrt((x(1) + mu)^2 + x(2)^2 + x(3)^2);
    r2 = sqrt((x(1) + mu - 1)^2 + x(2)^2 + x(3)^2);
    
    dU(1) = 2*x(5) + x(1) - (1-mu)*(x(1)+mu)/(r1^3) - ...
        mu*(x(1)+mu-1)/(r2^3);
    dU(2) = -2*x(4) + x(2) - (1-mu)*x(2)/(r1^3) - mu*x(2)/(r2^3);
    dU(3) = -(1-mu)*x(3)/(r1^3) - mu*x(3)/(r2^3);
end

