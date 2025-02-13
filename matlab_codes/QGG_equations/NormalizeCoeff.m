clear;
clc;

%%      NORMALIZED COEFFICIENT CODE

  % DESCRIPTION: This code converts from normalized coefficients to non
  % nomalized

  % normalized coefficient
  Xtrue = [5.2; -0.017511; -0.0000023; 0.0058194; 0.00561002; 0.001547; ...
            0.0001115;  0.0026660; 0.0102498; 0.0004360; -0.0021919; -0.0010761; 0.0021356;...
            -0.0013767; 0.0005161; 0.0005464; 0.0004250; 0.0013801; 0.0004869;...
            -0.0024871; 0.0011514; 0.0007281; 0.0013361; -0.0005090; 0.0001236; 0.00071695;...
            0; -0.0000197; 0.0015368; 0.0000653; -0.0009332; 0.0018562; 0.0007749; 0.0001024;...
            0.0030684; -0.0000786; -0.0012414; 0.0002269; -0.0005038; 0.0005241;...
            -0.0000084; -0.0002936; -0.0009706; -0.0006485; -0.0005993; -0.0000500];

  n = 2;
  m = 0;
  for j = 2:26 
      if(m == 0)
          delta = 1;
      else 
          delta = 0; 
      end
      f1 = factorial(n+m);
      f2 = factorial(n-m);

      N = 1/sqrt(f2*(2*n + 1)*(2 - delta)/f1);
      Xtrue(j) = Xtrue(j)/N;
      
      if(m < n)
          m = m +1;
      else
          n = n+1;
          m = 0;
      end
  end
  n = 2;
  m = 0;
  for j = 27:46             %WARNING: There's an error here. This has to be fixed. Problem with indexing
      if(m ~= 0)
          delta = 0; 
          f1 = factorial(n+m);
          f2 = factorial(n-m);
    
          N = 1/sqrt(f2*(2*n + 1)*(2 - delta)/f1);
          Xtrue(j) = Xtrue(j)/N;
      end
      
      if(m < n)
          m = m +1;
      else
          n = n+1;
          m = 0;
      end
  end