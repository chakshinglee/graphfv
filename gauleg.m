function QuadRule = gauleg(a,b,n,tol)
% GL_RULE 1D Gauss-Legendre quadrature rule.
%
%   QUADRULE = GAULEG(A,B,N,TOL) computes the N-point Gauss-Legendre
%   quadrature rule on the interval [A,B] up to the prescribed tolerance
%   TOL. If no tolerance is prescribed GAULEG uses the machine precision
%   EPS.
%
%   Note that all quadrature rules obtained from GAULEG are of order 2*N-1.
%
%   The struct QUADRULE contains the following fields:
%    W N-by-1 matrix specifying the weights of the quadrature rule.
%    X N-by-1 matrix specifying the abscissae of the quadrature rule.
%   
%   Example:
%
%   QuadRule = gauleg(0,1,10,1e-6);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Check for right number of arguments 

  if(nargin == 3)
    tol = eps;          
  end
  
  % Initalize variables
  
  m = floor((n+1)/2);
  xm = (b+a)/2;
  xl = (b-a)/2;
   
  for i = 1:m
    
    % Initial guess of root (starting value)
    
    z = cos(pi*(i-1/4)/(n+1/2));
    
    delta = tol+1;
    
    while(delta > tol)
        
      p1 = 0;
      p2 = 1;
      
      for k = 0:(n-1)

        % Computing value of n-th Legendre polynomial at point z using the
        % recursion:
        %
        %   (j+1)*P_(j+1)(z) = (2*j+1)*z*P_(j)(z)-j*P_(j-1)(z)
          
        p3 = ((2*k+1)*z*p2-k*p1)/(k+1);
        
        % Computing value of first derivative of n-th Legendre polynomial
        % at point z using the recursion:
        %
        %   (1-z^2)*P'_(j)(z) = j*[z*P_(j)(z)-P_(j-1)(z)]
        
        dp = n*(z*p3-p2)/(z^2-1);
        p1 = p2;
        p2 = p3;
        
      end    
    
      % Performing Newton update
      
      z_old = z;
      z = z_old-p3/dp;
      
      delta = abs(z-z_old);
      
    end
    
    % Computing weights and abscissae
    
    x(i) = xm-xl*z;
    x(n+1-i) = xm+xl*z;
    w(i) = 2*xl/((1-z^2)*dp^2);
    w(n+1-i) = w(i);
    
  end
  
  % Assign output arguments
  
  QuadRule.w = w(:);
  QuadRule.x = x(:);
        
return