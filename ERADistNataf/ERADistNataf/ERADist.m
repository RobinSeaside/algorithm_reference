classdef ERADist
%% Generation of distribution objects
%{
---------------------------------------------------------------------------
Developed by Sebastian Geyer, Iason Papaioannou and Felipe Uribe
Engineering Risk Analysis Group
Technische Universitat Munchen
www.era.bgu.tum.de
Version 2016-06
---------------------------------------------------------------------------
This software generates distribution objects according to the parameters
and definitions used in the distribution table of the ERA Group of TUM.
They can be defined either by their parameters, the first and second
moment or by data, given as a vector.
---------------------------------------------------------------------------
%}
   
%% MATLAB class: definition of the 'properties' block
properties
   Name         % Name of the distribution
   Par          % Parameters of the distribution
   Dist         % matlab distribution object of the marginal distribution
end

%% MATLAB class: definition of the 'methods' block
%{
Definition of all member functions of the ERADist class. Those are:
- ERADist    (Constructor)
- mean       returns mean
- std        returns standard deviaton
- cdf        returns value of CDF
- icdf       returns value of inverse CDF
- pdf        returns value of pdf
- random     generates random numbers according to the distribution of the object
%}
methods
   function Obj = ERADist(name,opt,val)
      Obj.Name = name;
      switch opt
         %-----------------------------------------------------------------
         case 'PAR'
            if strcmp(name,'standardnormal')
               val = [0,1];
            end
            Obj.Par = val;
            switch name
               case 'binomial'
                  if (0 <= val(2)) && (val(2) <= 1) && (val(1) > 0) && (mod(val(1),1) == 0)
                     Obj.Dist = makedist(name,'N',val(1),'p',val(2));
                  else
                     error('The binomial distribution is not defined for your parameters');
                  end
                  
               case 'geometric'
                  % special case of negative binomial distribution
                  if (0 <= val) && (val <= 1)
                     Obj.Dist = makedist('negativebinomial','r',1,'p',val);
                  else
                     error('The geometric distribution is not defined for your parameter');
                  end
                  
               case 'negativebinomial'
                  if (0 <= val(2)) && (val(2) <= 1) && (val(1) > 0) && (mod(val(1),1) == 0)
                     Obj.Dist = makedist(name,'r',val(1),'p',val(2));
                  else
                     error('The negative binomial distribution is not defined for your parameters');
                  end
                  
               case 'poisson'
                  if (val(1) > 0) && (val(2) > 0)
                     Obj.Dist = makedist(name,'lambda',val(1)*val(2));
                  else
                     error('The poisson distribution is not defined for your parameters');
                  end
                  
               case 'uniform'
                  Obj.Dist = makedist(name,'lower',val(1),'upper',val(2));
                  
               case 'standardnormal'
                  % special case of normal distribution
                  Obj.Dist = makedist('normal');
                  
               case 'normal'
                  if val(2) > 0
                     Obj.Dist = makedist(name,'mu',val(1),'sigma',val(2));
                  else
                     error('The normal distribution is not defined for your parameters');
                  end
                  
               case 'lognormal'
                  if val(2) > 0
                     Obj.Dist = makedist(name,'mu',val(1),'sigma',val(2));
                  else
                     error('The lognormal distribution is not defined for your parameters');
                  end
                  
               case 'exponential'
                  if val > 0
                     Obj.Dist = makedist(name,'mu',1/val);
                  else
                     error('The exponential distribution is not defined for your parameter');
                  end
                  
               case 'gamma'
                  if (val(1) > 0) && (val(2) > 0)
                     Obj.Dist = makedist(name,'a',val(1),'b',1/val(2));
                  else
                     error('The gamma distribution is not defined for your parameters');
                  end
                  
               case 'beta'
                  if (val(1) > 1) && (val(2) > 1) && (val(3) < val(4))
                     Obj.Dist = makedist(name,'a',val(1),'b',val(2));
                  else
                     error('The beta distribution is not defined for your parameters');
                  end
                  
              case 'gumbel'  % this distribution can be used to model minima
                  if val(2) > 0
                     Obj.Dist = makedist('extremevalue','mu',val(1),'sigma',val(2));
                  else
                     error('The gumbel distribution is not defined for your parameters');
                  end 
                  
               case 'gumbelMax'  % mirror image of this distribution can be used to model maxima
                  if val(2) > 0
                     Obj.Dist = makedist('extremevalue','mu',val(1),'sigma',val(2));
                  else
                     error('The gumbel distribution is not defined for your parameters');
                  end
                  
               case 'frechet'
                  if (val(1) > 0) && (val(2) > 0)
                     Obj.Dist = makedist('GeneralizedExtremeValue','k',1/val(2),'sigma',val(1)/val(2),'mu',val(1));
                  else
                     error('The frechet distribution is not defined for your parameters');
                  end
                  
               case 'weibull'
                  if (val(1) > 0) && (val(2) > 0)
                     Obj.Dist = makedist(name,'a',val(1),'b',val(2));
                  else
                     error('The weibull distribution is not defined for your parameters');
                  end
                  
               case 'GEV'
                  if val(2) > 0
                     Obj.Dist = makedist('GeneralizedExtremeValue','k',val(1),'sigma',val(2),'mu',val(3));
                  else
                     error('The Generalized extreme value distribution is not defined for your parameters');
                  end
                  
               case 'GEVMin'
                  if val(2) > 0
                     Obj.Dist = makedist('GeneralizedExtremeValue','k',val(1),'sigma',val(2),'mu',val(3));
                  else
                     error('The Generalized extreme value distribution is not defined for your parameters');
                  end
                  
               case 'pareto'
                  if (val(1) > 0) && (val(2) > 0)
                     Obj.Dist = makedist('Generalizedpareto','k',1/val(2),'sigma',val(1)/val(2),'theta',val(1));
                  else
                     error('The pareto distribution is not defined for your parameters');
                  end
                  
               case 'rayleigh'
                  if val > 0
                     Obj.Dist = makedist(name,'b',val);
                  else
                     error('The rayleigh distribution is not defined for your parameters');
                  end
                  
               case 'chisquare'
                  % special case of gamma distribution
                  if (val>0) && (mod(val,1)==0)
                     Obj.Dist = makedist('gamma','a',val/2,'b',2);
                  else
                     error('The chisquare distribution is not defined for your parameters');
                  end
                  
               otherwise
                  disp('Not available distribution');
            end
         %-----------------------------------------------------------------   
         case 'MOM'
            if val(2) < 0
               error('The standard deviation must not be negative');
            else
               switch name
                  case 'binomial'
                     % Solve System of two equations for the parameters of the distribution
                     Obj.Par(2) = (1-val(2)^2)/val(1);
                     Obj.Par(1) = val(1)/Obj.Par(2);
                     % Evaluate if distribution can be defined on the parameters
                     if mod(Obj.Par(1),1) <= 1e-4
                        Obj.Par = round(Obj.Par,0);
                     elseif (0 <= Obj.Par(2)) && (Obj.Par(2) <= 1)
                     else
                        error('Please select other moments');
                     end
                     % Define distribution
                     Obj.Dist = makedist(name,'N',Obj.Par(1),'p',Obj.Par(2));
                     
                  case 'geometric'                     
                     % Solve Equation for the parameter of the distribution based on the first moment
                     Obj.Par = 1/val(1);                     
                     % Evaluate if distribution can be defined on the paramater and if the moments are well defined
                     if (0 <= Obj.Par) && (Obj.Par <= 1)
                     else
                        error('Please select other moments')
                     end                     
                     % Define distribution
                     Obj.Dist = makedist('negativebinomial','r',1,'p',Obj.Par);
                                          
                  case 'negativebinomial'
                     % Solve System of two equations for the parameters of the distribution
                     Obj.Par(2) = val(1)/(val(1)+val(2)^2);
                     Obj.Par(1) = Obj.Par(2)*val(1);
                     % Evaluate if distribution can be defined on the parameters
                     if mod(Obj.Par(1),1) <= 1e-4
                        Obj.Par = round(Obj.Par,0);
                     elseif (0 <= Obj.Par(2)) && (Obj.Par(2) <= 1)
                     else
                        error('Please select other moments');
                     end
                     % Define distribution
                     Obj.Dist = makedist(name,'r',Obj.Par(1),'p',Obj.Par(2));
                     
                  case 'poisson'
                     Obj.Par = val(1);
                     % Evaluate if moments match
                     if 0 < Obj.Par
                     else
                        error('Please select other moments');
                     end
                     % Define distribution
                     Obj.Dist = makedist(name,'lambda',val(1));
                     
                  case 'uniform'
                     syms a b;
                     sol        = solve((a+b)/2==val(1),((b-a)^2)/12==val(2)^2);
                     Obj.Par(1) = double(sol.a);
                     Obj.Par(2) = double(sol.b);   
                     % Define distribution
                     Obj.Dist = makedist(name,'lower',Obj.Par(1),'upper',Obj.Par(2));
                     
                  case 'standardnormal'
                     % special case of normal distribution
                     Obj.Par  = [0,1];
                     Obj.Dist = makedist('normal');
                     
                  case 'normal'
                     Obj.Par  = val;
                     Obj.Dist = makedist(name,'mu',val(1),'sigma',val(2));
                     
                  case 'lognormal'
                     % Solve two equations for the parameters of the distribution
                     Obj.Par(1) = log(val(1)) - log(sqrt(1+(val(2)/val(1))^2));  % mean normal
                     Obj.Par(2) = sqrt(log(1+(val(2)/val(1))^2));   % sigma normal
                     Obj.Dist   = makedist(name,'mu',Obj.Par(1),'sigma',Obj.Par(2));
                     
                  case 'exponential'
                     % Solve Equation for the parameter of the distribution based on the first moment
                     Obj.Par = 1/val(1);
                     % Evaluate if distribution can be defined on the paramater and if the moments are well defined
                     if (0 <= Obj.Par) && (Obj.Par <= Inf)
                     else
                        error('Please select other moments');
                     end
                     % Define distribution
                     Obj.Dist = makedist(name,'mu',val(1));
                     
                  case 'gamma'
                     % Solve System of two equations for the parameters of the distribution
                     syms a b;
                     sol        = solve(a/b==val(1),a/b^2==val(2)^2);
                     Obj.Par(1) = double(sol.a);
                     Obj.Par(2) = double(sol.b);   
                     % Evaluate if distribution can be defined on the parameters
                     if (0 < Obj.Par(1)) && (0 < Obj.Par(2))
                     else
                        error('Please select other moments');
                     end
                     % Define distribution
                     Obj.Dist = makedist(name,'a',Obj.Par(1),'b',1/Obj.Par(2));                     
                     
                  case 'beta'                     
                     if val(4) <= val(3)
                        error('Upper bound of Beta distribution is smaller or equal to lower bound');
                     end
                     % Solve System of two equations for the parameters of the distribution
                     Obj.Par(1) = ((val(4)-val(1))*(val(1)-val(3))/val(2)^2-1)*...
                                  (val(1)-val(3))/(val(4)-val(3));
                     Obj.Par(2) = Obj.Par(1)*(val(4)-val(1))/(val(1)-val(3));
                     Obj.Par(3) = val(3);
                     Obj.Par(4) = val(4);                     
                     % Evaluate if distribution can be defined on the parameters
                     if (Obj.Par(1)>0) && (Obj.Par(2)>0)
                     else
                        error('Please select other moments');
                     end                     
                     % Define distribution
                     Obj.Dist = makedist(name,'a',Obj.Par(1),'b',Obj.Par(2));                     
                     
                  case 'gumbel'   % gumbel can be used to model minima
                     % Solve two equations for the parameters of the distribution
                     Obj.Par(2) = val(2)*sqrt(6)/pi;                      % scale parameter
                     Obj.Par(1) = val(1) + double(eulergamma)*Obj.Par(2); % location parameter
                     % Define distribution
                     Obj.Dist = makedist('extremevalue','mu',Obj.Par(1),'sigma',Obj.Par(2));  
                     
                  case 'gumbelMax'   % mirror image of gumbel can be used to model maxima
                     % Solve two equations for the parameters of the distribution
                     Obj.Par(2) = val(2)*sqrt(6)/pi;                     % scale parameter
                     Obj.Par(1) = -(val(1) + double(eulergamma)*Obj.Par(2)); % location parameter
                     % Define distribution
                     Obj.Dist = makedist('extremevalue','mu',-Obj.Par(1),'sigma',Obj.Par(2));                     
                     
                  case 'frechet'
                     % Solve equation for the parameters of the distribution
                     options = optimset('Display','off');
                     par0    = [2.001,1.0e3];
                     fun     = @(par) sqrt(gamma(1-2/par)-(gamma(1-1/par)).^2)./...
                                      gamma(1-1/par)-val(2)/val(1);
                     [xs,~,exitflag] = fzero(fun,par0,options);
                     if exitflag > 0
                        Obj.Par(2) = xs;
                        Obj.Par(1) = val(1)/gamma(1-1/Obj.Par(2));
                     else
                        error('fzero could not converge to a solution for determining the parameters of the frechet distribution');
                     end
                     % Evaluate if distribution can be defined on the parameters
                     if (Obj.Par(1) > 0) && (Obj.Par(2) > 0)
                     else
                        error('Please select other moments');
                     end
                     % Define distribution
                     Obj.Dist = makedist('GeneralizedExtremeValue','k',1/Obj.Par(2),'sigma',Obj.Par(1)/Obj.Par(2),'mu',Obj.Par(1));
                     
                  case 'weibull'
                     % Solve equation for the parameters of the distribution
                     options = optimset('Display','off');
                     par0    = [0.02,1.0e3];
                     fun     = @(par) sqrt(gamma(1+2/par)-(gamma(1+1/par)).^2)./gamma(1+1/par)-val(2)/val(1);
                     [xs,~,exitflag] = fzero(fun,par0,options);
                     if exitflag > 0
                        Obj.Par(2) = xs;
                        Obj.Par(1) = val(1)/gamma(1+1/Obj.Par(2));
                     else
                        error('fzero could not converge to a solution for determining the parameters of the weibull distribution')
                     end
                     % Evaluate if distribution can be defined on the parameters
                     if (Obj.Par(1) > 0) && (Obj.Par(2) > 0)
                     else
                        error('Please select other moments');
                     end
                     % Define distribution
                     Obj.Dist = makedist(name,'a',Obj.Par(1),'b',Obj.Par(2));
                     
                  case 'GEV'
                     % Solve System of two equations for the parameters of the distribution
                     if val(1) == val(3)
                        Obj.Par(1) = -1;
                        Obj.Par(2) = val(2);
                        Obj.Par(3) = val(3);
                     else
                        options = optimset('Display','off');
                        if val(1) > val(3)
                           par0 = 0.3;
                        else
                           par0 = -1.5;
                        end
                        fun = @(par) (gamma(1-2*par)-gamma(1-par).^2)./(gamma(1-par)-1).^2-...
                                     (val(2)/(val(1)-val(3)))^2;
                        [xs,~,exitflag] = fzero(fun,par0,options);
                        if exitflag > 0
                           Obj.Par(1) = xs;
                           Obj.Par(2) = (val(1)-val(3))*Obj.Par(1)/(gamma(1-Obj.Par(1))-1);
                           Obj.Par(3) = val(3);
                        else
                           error('fzero could not converge to a solution for determining the parameters of the GEV distribution')
                        end
                     end
                     % Evaluate if distribution can be defined on the parameters
                     if Obj.Par(2) > 0
                     else
                        error('Please select other moments');
                     end
                     % Define distribution
                     Obj.Dist = makedist('GeneralizedExtremeValue','k',Obj.Par(1),'sigma',Obj.Par(2),'mu',Obj.Par(3));
                     
                  case 'GEVMin' % mirror image of gumbel can be used to model maxima
                     % Solve System of two equations for the parameters of the distribution
                     if val(1) == val(3)
                        Obj.Par(1) = -1;
                        Obj.Par(2) = val(2);
                        Obj.Par(3) = val(3);
                     else
                        options = optimset('Display','off');
                        if val(1) > val(3)
                           par0 = 0.3;
                        else
                           par0 = -1.5;
                        end
                        fun = @(par) (gamma(1-2*par)-gamma(1-par).^2)./(gamma(1-par)-1).^2-...
                                     (val(2)/(val(1)-val(3)))^2;
                        [xs,~,exitflag] = fzero(fun,par0,options);
                        if exitflag > 0
                           Obj.Par(1) = xs;
                           Obj.Par(2) = (val(1)-val(3))*Obj.Par(1)/(gamma(1-Obj.Par(1))-1);
                           Obj.Par(3) = val(3);
                        else
                           error('fzero could not converge to a solution for determining the parameters of the GEV distribution')
                        end
                     end
                     % Evaluate if distribution can be defined on the parameters
                     if Obj.Par(2) > 0
                     else
                        error('Please select other moments');
                     end
                     % Define distribution
                     Obj.Dist = makedist('GeneralizedExtremeValue','k',Obj.Par(1),'sigma',Obj.Par(2),'mu',Obj.Par(3));
                     
                  case 'pareto'
                     % Solve System of two equations for the parameters of the distribution
                     Obj.Par(2) = 1 + sqrt(1+(val(1)/val(2))^2);
                     Obj.Par(1) = val(1)*(Obj.Par(2)-1)/Obj.Par(2);
                     % Evaluate if distribution can be defined on the parameters
                     if (Obj.Par(1) > 0) && (Obj.Par(2 )> 0)
                     else
                        error('Please select other moments');
                     end
                     % Define distribution
                     Obj.Dist = makedist('Generalizedpareto','k',1/Obj.Par(2),'sigma',Obj.Par(1)/Obj.Par(2),'theta',Obj.Par(1));
                     
                  case 'rayleigh'                     
                     % Solve Equation for the parameter of the distribution based on the first moment
                     Obj.Par = val(1)/sqrt(pi/2);                     
                     % Evaluate if distribution can be defined on the paramater and if the moments are well defined
                     if (0 < Obj.Par) && (Obj.Par <= Inf)
                     else
                        error('Please select other moments');
                     end
                     % Define distribution
                     Obj.Dist = makedist(name,'b',Obj.Par);                     
                     
                  case 'chisquare'                     
                     % Solve Equation for the parameter of the distribution based on the first moment
                     Obj.Par = val(1);                     
                     % Evaluate if distribution can be defined on the paramater and if the moments are well defined
                     if mod(Obj.Par,1) <= 1e-4
                        Obj.Par = round(Obj.Par,0);
                     elseif (0 < Obj.Par) && (Obj.Par <= Inf)
                     else
                        error('Please select other moments');
                     end
                     % Define distribution
                     Obj.Dist = makedist('gamma','a',Obj.Par/2,'b',2);    % chi-squared special case of the gamma                
                     
                  otherwise
                     disp('Not available distribution');
               end
            end
         %-----------------------------------------------------------------   
         case 'DATA'
            switch name
               case 'binomial' % Error occurs by using fitdist
                  Obj.Dist = fitdist(val,name);
                  Obj.Par  = (Obj.Dist.ParameterValues);
                  
               case 'geometric' % see negativebinomial
                  error('The geometric distribution is not supported in DATA');
                  
               case 'negativebinomial' % Error occurs by using fitdist
                  Obj.Dist = fitdist(val,name);
                  Obj.Par  = (Obj.Dist.ParameterValues);
                  
               case 'poisson'
                  Obj.Dist = fitdist(val,name);
                  Obj.Par  = (Obj.Dist.ParameterValues);
                  
               case 'normal'
                  Obj.Dist = fitdist(val,name);
                  Obj.Par  = (Obj.Dist.ParameterValues);
                  
               case 'lognormal'
                  Obj.Dist = fitdist(val,name);
                  Obj.Par  = (Obj.Dist.ParameterValues);
                  
               case 'exponential'
                  Obj.Dist = fitdist(val,name);
                  Obj.Par  = (1/Obj.Dist.ParameterValues);
                  
               case 'gamma'
                  Obj.Dist = fitdist(val,name);
                  Obj.Par  = [Obj.Dist.ParameterValues(1) 1/Obj.Dist.ParameterValues(2)];
                  
               case 'beta' % Upper and lower bound have to be implemented
                  error('The beta distribution is not supported in DATA.');
                  
               case 'gumbel'   % to model the minimum value
                  Obj.Dist = fitdist(val,'extremevalue');
                  Obj.Par  = [Obj.Dist.ParameterValues(1) Obj.Dist.ParameterValues(2)];
                  
               case 'gumbelMax' % to model the maximum value, use the negative of the original values
                  Obj.Dist = fitdist(-val,'extremevalue');
                  Obj.Par  = [Obj.Dist.ParameterValues(1) Obj.Dist.ParameterValues(2)];
                  
               case 'frechet'
                  error('The frechet distribution is not supported in DATA.');
                  
               case 'weibull'
                  Obj.Dist = fitdist(val,name);
                  Obj.Par  = Obj.Dist.ParameterValues;
                  
               case 'GEV'
                  Obj.Dist = fitdist(val,'GeneralizedExtremeValue');
                  Obj.Par  = Obj.Dist.ParameterValues;
                  
               case 'GEVMin'
                  Obj.Dist = fitdist(-val,'GeneralizedExtremeValue');
                  Obj.Par  = Obj.Dist.ParameterValues;
                  
               case 'pareto'
                  error('The pareto distribution is not supported in DATA.')
                  
               case 'rayleigh'
                  Obj.Dist = fitdist(val,name);
                  Obj.Par  = Obj.Dist.ParameterValues;
                  
               case 'chisquare'
                  error('The chi-square distribution is not supported in DATA.')
                  
               otherwise
                  disp('Not available distribution');
            end
            
         otherwise
            disp('Not available distribution');
      end
   end
   
   %-----------------------------------------------------------------------
   function Mean = mean(Obj)
      switch Obj.Name
         case 'geometric'
            Mean = Obj.Dist.mean+1;            
         case 'negativebinomial'
            Mean = Obj.Dist.mean+Obj.Dist.R;            
         case 'beta'
            Mean = (Obj.Par(2)*Obj.Par(3)+Obj.Par(1)*Obj.Par(4))/(Obj.Par(1)+Obj.Par(2));   
         case 'gumbelMax'
            Mean = -Obj.Dist.mean;
         case 'GEVMin'
            Mean = -Obj.Dist.mean;
         case 'frechet'
            Mean = Obj.Par(1)*gamma(1-1/Obj.Par(2));
         otherwise
            Mean = Obj.Dist.mean;
      end
   end
   
   %-----------------------------------------------------------------------
   function Standarddeviation = std(Obj)
      switch Obj.Name
         case 'beta'
            Standarddeviation = Obj.Dist.std*(Obj.Par(4)-Obj.Par(3));
         otherwise
            Standarddeviation = Obj.Dist.std;
      end
   end
   
   %-----------------------------------------------------------------------
   function CDF = cdf(Obj,x)
      switch Obj.Name
         case 'geometric'
            CDF = Obj.Dist.cdf(x-1);
         case 'negativebinomial'
            CDF = Obj.Dist.cdf(x-Obj.Par(1));
         case 'beta'
            CDF = Obj.Dist.cdf((x-Obj.Par(3))/(Obj.Par(4)-Obj.Par(3)));
         case 'gumbelMax'
            CDF = Obj.Dist.cdf(-x);
         case 'GEVMin'
            CDF = Obj.Dist.cdf(-x);
         otherwise
            CDF = Obj.Dist.cdf(x);
      end
   end
   
   %-----------------------------------------------------------------------
   function InverseCDF = icdf(Obj,y)
      switch Obj.Name
         case 'geometric'
            InverseCDF = Obj.Dist.icdf(y)+1;
         case 'negativebinomial'
            InverseCDF = Obj.Dist.icdf(y)+Obj.Par(1);
         case 'beta'
            InverseCDF = Obj.Dist.icdf(y)*(Obj.Par(4)-Obj.Par(3))+Obj.Par(3);
         case 'gumbelMax'
            InverseCDF = -Obj.Dist.icdf(y);
         case 'GEVMin'
            InverseCDF = -Obj.Dist.icdf(y);
         otherwise
            InverseCDF = Obj.Dist.icdf(y);
      end
   end
   
   %-----------------------------------------------------------------------
   function PDF = pdf(Obj,x)
      switch Obj.Name
         case 'geometric'
            PDF = Obj.Dist.pdf(x-1);
         case 'negativebinomial'
            PDF = Obj.Dist.pdf(x-Obj.Par(1));
         case 'beta'
            PDF = Obj.Dist.pdf((x-Obj.Par(3))/(Obj.Par(4)-Obj.Par(3)));
         case 'gumbelMax'
            PDF = Obj.Dist.pdf(-x);
         case 'GEVMin'
            PDF = Obj.Dist.pdf(-x);
         otherwise
            PDF = Obj.Dist.pdf(x);
      end
   end
   
   %-----------------------------------------------------------------------
   function Random = random(Obj,m,n)
      if nargin == 2
         switch Obj.Name
            case 'geometric'
               Random = random(Obj.Dist,m)+1;
            case 'negativebinomial'
               Random = random(Obj.Dist,m)+Obj.Par(1);
            case 'beta'
               Random = random(Obj.Dist,m)*(Obj.Par(4)-Obj.Par(3))+Obj.Par(3);
            case 'gumbelMax'
               Random = -random(Obj.Dist,m);
            case 'GEVMin'
               Random = -random(Obj.Dist,m);
            otherwise
               Random = random(Obj.Dist,m);
         end
      elseif nargin == 3
         switch Obj.Name
            case 'geometric'
               Random = random(Obj.Dist,m,n)+1;
            case 'negativebinomial'
               Random = random(Obj.Dist,m,n)+Obj.Par(1);
            case 'beta'
               Random = random(Obj.Dist,m,n)*(Obj.Par(4)-Obj.Par(3))+Obj.Par(3);
            case 'gumbelMax'
               Random = -random(Obj.Dist,m,n);
            case 'GEVMin'
               Random = -random(Obj.Dist,m,n);
            otherwise
               Random = random(Obj.Dist,m,n);
         end
      else
         error('Wrong Input');
      end
   end
end

end