classdef ERANataf
   %% Nataf Transformation of random variables
   %{
---------------------------------------------------------------------------
Developed by Sebastian Geyer, Iason Papaioannou and Felipe Uribe
Engineering Risk Analysis Group
Technische Universitat Munchen
www.era.bgu.tum.de
Version 2016-06
---------------------------------------------------------------------------
* This software performs the Nataf transformation of random variables.
* It is possible to generate random numbers according to their Nataf
joint pdf and to evaluate this pdf.
* The inverse Nataf transformation is also defined as a function.
* It requires the use of Objects of the class ERADist which is also
published on the homepage of the ERA Group of TUM.
---------------------------------------------------------------------------
References:
 1. http://www.dynardo.de/fileadmin/Material_Dynardo/bibliothek/WOST_8.0/Paper_Bucher.pdf
---------------------------------------------------------------------------
   %}
   
   %% MATLAB class: definition of the 'properties' block
   properties
      Rho_X         % Correlation matrix of the vector X
      Rho_Z         % Correlation matrix of the correlated normal vector Z
      A             % Lower triangular matrix of the Cholesky decomposition of Rho_Z
      Marginals     % Contains all marginal distribution objects
   end
   
   %% MATLAB class: definition of the 'methods' block
   %{
    Definition of all member functions of the ERANataf class. Those are:
   - ERANataf   (Constructor)
   - X2U        samples from physical to standard
   - U2X        samples from standard to physical
   - random     Generates random numbers according to their joint pdf
   - pdf        Evaluaters the joint pdf of a sample vector
   %}
   methods
      % The constructor transforms the Correlation matrix Rho_x to Rho_Z.
      % Furthermore it evaluates whether Rho_Z is postdtive definite. Matlab will
      % exit the code if that is not the case.
      % A is the lower triangular matrix generated by the cholesky
      % decompostdtion such that Rho_Z = A*A'
      function Obj = ERANataf(M,Correlation)
         Obj.Marginals = M;
         n_dist        = length(M);
         Obj.Rho_X     = Correlation;
         
         % Check whether all distributions have finite moments
         for i = 1:n_dist
            if ~(isfinite(M(i).mean) && isfinite(M(i).std))
               error('The marginal distributions need to have finite mean and variance');
            end
         end
         
         % Calculation of the transformed correlation matrix. This is achieved by a
         % quadratic two-dimensional Gauss-Legendre integration
         zmax = 6;       % Integration bounds
         zmin = -zmax;   % Integration bounds
         n    = 1024;    % number of integration points along each dimension
         
         % Legendre-Gauss nodes and weights  on the interval [zmin,zmax]
         [points,w_1D] = quad_GL(n,zmin,zmax);
         
         xi    = reshape(repmat(points,1,n)',[n^2,1]);
         eta   = repmat(points,n,1);
         w_2D  = reshape(repmat(w_1D,1,n)',[n^2,1]).*repmat(w_1D,n,1);
         f_xi  = zeros(n^2,n_dist);
         f_eta = zeros(n^2,n_dist);
         
         % check is X is the identity
         Obj.Rho_Z = eye(n_dist);
         if norm(Obj.Rho_X-eye(n_dist)) >= 1e-5   % samples are uncorrelated
            
            % Transformation of parameters of the distributions for the
            % matlab cdf/icdf functions and calculations of those
            for i = 1:n_dist
               f_eta(:,i) = (M(i).icdf(normcdf(eta))-M(i).mean)/M(i).std;
               f_xi(:,i)  = (M(i).icdf(normcdf(xi))-M(i).mean)/M(i).std;
            end
            
            % This loop construction makes sure that no unnecessary
            % calculations are done -> Only off-diagonal elements of the upper part of
            % the matrix are calculated -> afterwards they are copied to the lower part
            % options = optimoptions(@fsolve,'Display','off');
            options   = optimset('Display','off');
            for i = 1:n_dist
               for j = i+1:n_dist
                  if Obj.Rho_X(i,j) == 0
                     continue;
                  end
                  if strcmp(M(i).Name,'standardnormal') && strcmp(M(j).Name,'standardnormal')
                     Obj.Rho_Z(i,j) = Obj.Rho_X(i,j);
                     Obj.Rho_Z(j,i) = Obj.Rho_Z(i,j);
                     continue;
                  elseif strcmp(M(i).Name,'normal') && strcmp(M(j).Name,'normal')
                     Obj.Rho_Z(i,j) = Obj.Rho_X(i,j);
                     Obj.Rho_Z(j,i) = Obj.Rho_Z(i,j);
                     continue;
                  elseif strcmp(M(i).Name,'normal') && strcmp(M(j).Name,'lognormal')
                     Vj             = M(j).std/M(j).mean;
                     Obj.Rho_Z(i,j) = Obj.Rho_X(i,j)*Vj/sqrt(log(1+Vj^2));
                     Obj.Rho_Z(j,i) = Obj.Rho_Z(i,j);
                     continue;
                  elseif strcmp(M(i).Name,'lognormal') && strcmp(M(j).Name,'normal')
                     Vi             = M(i).std/M(i).mean;
                     Obj.Rho_Z(i,j) = Obj.Rho_X(i,j)*Vi/sqrt(log(1+Vi^2));
                     Obj.Rho_Z(j,i) = Obj.Rho_Z(i,j);
                     continue;
                  elseif strcmp(M(i).Name,'lognormal') && strcmp(M(j).Name,'lognormal')
                     Vi             = M(i).std/M(i).mean;
                     Vj             = M(j).std/M(j).mean;
                     Obj.Rho_Z(i,j) = log(1+Obj.Rho_X(i,j)*Vi*Vj)/sqrt(log(1+Vi^2)*log(1+Vj^2));
                     Obj.Rho_Z(j,i) = Obj.Rho_Z(i,j);
                     continue;
                  end
                  % solving Nataf
                  coef = f_xi(:,j).*f_eta(:,i).*w_2D;
                  fun  = @(rho0) sum(coef.*1/(2*pi*sqrt(1-rho0^2)).*exp( -1/(2*(1-rho0^2))*...
                     (xi.^2 - 2*rho0*xi.*eta + eta.^2) ))-Obj.Rho_X(i,j);
                  [xs,~,exitflag] = fzero(fun,Obj.Rho_X(i,j),options);
                  if exitflag > 0
                     Obj.Rho_Z(i,j) = xs;
                     Obj.Rho_Z(j,i) = Obj.Rho_Z(i,j);
                  else
                     fun = @(rho0) sum(coef.*1/(2*pi*sqrt(1-rho0^2)).*exp( -1/(2*(1-rho0^2)) *...
                        (xi.^2 - 2*rho0*xi.*eta + eta.^2) ))-Obj.Rho_X(i,j);
                     [xs,~,exitflag] = fzero(fun,-Obj.Rho_X(i,j),options);
                     if exitflag > 0
                        Obj.Rho_Z(i,j) = xs;
                        Obj.Rho_Z(j,i) = Obj.Rho_Z(i,j);
                     else
                        error('fzero could not converge to a solution of the Nataf integral equation');
                     end
                  end
               end
            end
         end
         % perform Choleski decomposition
         [Obj.A,p] = chol(Obj.Rho_Z,'lower');
         if p ~= 0
            error('Transformed correlation matrix is not positive definite --> Nataf transformation is not applicable');
         end
      end
      
      %-----------------------------------------------------------------------
      % This function performs the transformation from X to U by taking
      % the inverse standard normal cdf of the cdf of every value. Then it
      % performs the transformation from Z to U. A is the lower
      % triangular matrix of the cholesky decomposition of Rho_Z and
      % U is the resulting independent standard normal vector
      % Afterwards it calculates the Jacobian of this Transformation if it is needed
      function [U,Jac] = X2U(Nataf,X,opt)
         if size(Nataf.A,1) ~= size(X,1)
            X = X';
         end
         [m,n] = size(X);
         if length(Nataf.Marginals)==1
            Z = zeros(n,m);
            for i = 1:m
               Z(i,:) = norminv(Nataf.Marginals(i).cdf(X(i,:)));
            end
            diag = zeros(m);
         else
            Z = zeros(m,n);
            for i = 1:m
               Z(i,:) = norminv(Nataf.Marginals(i).cdf(X(i,:)));
            end
            diag = zeros(m);
         end
         U = Nataf.A\Z;
         
         % Jacobian of X to U
         if nargin == 2
            Jac = [];
         elseif strcmp('Jac',opt) == 1
            for i = 1:m
               diag(i,i) = normpdf(Z(i))/Nataf.Marginals(i).pdf(X(i));
            end
            Jac = diag*Nataf.A;
         else
            error('Wrong Input');
         end
         
      end
      
      %-----------------------------------------------------------------------
      % This function performs the transformation from U to X
      function [X,Jac] = U2X(Nataf,U,opt)
         if size(Nataf.A,1) ~= size(U,1)
            U = U';
         end
         Z     = Nataf.A*U;
         [m,n] = size(U);
         if length(Nataf.Marginals)==1
            X = zeros(m,n);
            for i = 1:m
               X(i,:) = Nataf.Marginals.icdf(normcdf(Z(i,:)));
            end
            diag = zeros(m);
         else
            X = zeros(m,n);
            for i = 1:m
               X(i,:) = Nataf.Marginals(i).icdf(normcdf(Z(i,:)));
            end
            diag = zeros(m);
         end
         
         % Jacobian of U to X
         if nargin == 2
            Jac = [];
         elseif strcmp('Jac',opt) == 1
            for i = 1:m
               diag(i,i) = Nataf.Marginals(i).pdf(X(i))/normpdf(Z(i));
            end
            Jac = Nataf.A\diag;
         else
            error('Wrong Input');
         end
      end
      
      %-----------------------------------------------------------------------
      % This function generates random numbers according to their joint
      % distribution
      function jointrandom = random(Nataf,N)
         % Generate uncorrelated standard normal random variables u and
         % transform them to x --> Every column corresponds to one
         % correlated sample of each variable
         [m,n] = size(Nataf.Marginals);
         if m > n
            Size = m;
         else
            Size = n;
         end
         U           = randn(Size,N);
         Z           = Nataf.A*U;
         jointrandom = zeros(Size,N);
         for i = 1:Size
            jointrandom(i,:) = Nataf.Marginals(i).icdf(normcdf(Z(i,:)));
         end
      end
      
      %-----------------------------------------------------------------------
      % This function returns the joint PDF of the Nataf distribution
      function jointpdf = pdf(Nataf,X)
         if iscolumn(X) == 1
            X = X';
         end
         [m,n] = size(Nataf.Marginals);
         if m > n
            Size = m;
         else
            Size = n;
         end
         n    = size(X,1);
         U    = zeros(Size,n);
         phi  = zeros(Size,n);
         f    = zeros(Size,n);
         mu   = zeros(1,Size);
         for i = 1:Size
            U(i,:)   = norminv(Nataf.Marginals(i).cdf(X(:,i)));
            phi(i,:) = normpdf(U(i,:));
            f(i,:)   = Nataf.Marginals(i).pdf(X(:,i));
         end
         phi_n    = mvnpdf(U',mu,Nataf.Rho_Z);   % multivariate with mean 0 std 1 and corr Rho_Z
         jointpdf = zeros(n,1);
         for i = 1:n
            jointpdf(i) = (prod(f(:,i),1)/prod(phi(:,i),1))*phi_n(i);
            if isnan(jointpdf(i)) == 1
               jointpdf(i) = 0;
            end
         end
      end
      
      %-----------------------------------------------------------------------
      % This function returns the joint CDF of the Nataf distribution
      function jointcdf = cdf(Nataf,X)
         if iscolumn(X) == 1
            X = X';
         end
         Size = size(Nataf.Marginals,2);
         n    = size(X,1);
         U    = zeros(Size,n);
         mu   = zeros(1,Size);
         for i = 1:Size
            U(i,:) = norminv(Nataf.Marginals(i).cdf(X(:,i)),0,1);
         end
         jointcdf = mvncdf(U',mu,Nataf.Rho_Z);
      end
   end
   
end

%% nested function: Gauss-Legendre quadrature
function [x,w] = quad_GL(N,a,b)
%{
Written by Greg von Winckel - 02/25/2004
This script is for computing definite integrals using Legendre-Gauss
Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
[a,b] with truncation order N

Suppose you have a continuous function f(x) which is defined on [a,b]
which you can evaluate at any x in [a,b]. Simply evaluate it at all of
the values contained in the x vector to obtain a vector f. Then compute
the definite integral using sum(f.*w);
%}
N  = N-1;
N1 = N+1;
N2 = N+2;
xu = linspace(-1,1,N1)';
%
y  = cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);   % Initial guess
L  = zeros(N1,N2);   % Legendre-Gauss Vandermonde Matrix
Lp = zeros(N1,N2);   % Derivative of LGVM

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0 = 2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
   L(:,1)  = 1;
   L(:,2)  = y;
   Lp(:,1) = 0;
   Lp(:,2) = 1;
   for k = 2:N1
      L(:,k+1) = ( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
   end
   Lp = (N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
   y0 = y;
   y  = y0-L(:,N2)./Lp;
end

% Linear map from [-1,1] to [a,b]
x = (a*(1-y)+b*(1+y))/2;

% Compute the weights
w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

end