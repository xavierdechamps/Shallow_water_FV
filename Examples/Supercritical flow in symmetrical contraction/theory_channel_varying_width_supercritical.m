clear,clc,close all
format long e

p.B1 = 0.6096;
p.B2 = 0.3048;
y(1) = 0.03048 ;
p.g = 9.81;
p.Q = 0.041 ;
V(1) = p.Q / ( p.B1 * y(1)) ;
Fr(1)= V(1) / sqrt(p.g*y(1)) ;
p.Q = V(1) * p.B1 * y(1);
p.n = 0.015 ; %0.0107;
p.S0= 0.05664 ;
p.s = 0;
p.eps = 1.e-6;

L1 = 1.          * 0.3048; % Length upstream
L2 = 4.757182227 * 0.3048; % Length contraction
L3 = 7.5         * 0.3048; % Length downstream

q = p.Q / p.B1 ;
yc(1) = ( q^2 / p.g )^(1/3) ;

% Downstream section
z(1) = 0;
numsteps = 101;
x     = linspace (0.,L1,numsteps) ;
B     = linspace (p.B1,p.B1,numsteps);
gradB = linspace (0,0,numsteps);
gradB(end) = 0.5*(p.B2 - p.B1)/L2 ; % Average between gradient and parallel
ind(1) = numsteps ; 

% Contraction section
numsteps = 501;
x(end:end+numsteps-1)     = linspace (L1,L1+L2,numsteps) ;
B(end:end+numsteps-1)     = linspace (p.B1,p.B2,numsteps) ;
gradB(end+1:end+numsteps-1) =  (p.B2 - p.B1)/L2 ;
gradB(end) = 0.5*(p.B2 - p.B1)/L2 ; % Average between gradient and parallel
ind(2) = ind(1) + numsteps - 1; 

% Downstream section
numsteps = 501;
x(end:end+numsteps-1)     = linspace (L1+L2,L1+L2+L3,numsteps) ;
B(end:end+numsteps-1)     = linspace (p.B2,p.B2,numsteps) ;
gradB(end+1:end+numsteps-1) = 0.;
ind(3) = ind(2) + numsteps - 1; 

%
A(1)   = B(1)*y(1) + p.s*y(1)^2 ; % Area
P(1)   = B(1) + 2*y(1)*sqrt(1+p.s^2); % Wetted perimeter
R(1)   = A(1)/P(1) ;
Sf(1)  = p.n^2 * V(1)^2 / ( R(1)^(4/3) ) ;
f(1)   = ( p.S0 - Sf(1) + p.Q^2 * y(1)* gradB(1) / (p.g*A(1).^3) ) / ( 1 - p.Q^2 * B(1) / (p.g*A(1).^3) ) ;

for i = 2:ind(3)
  go=1;
  [ynew,z(i),A(i),P(i),R(i),V(i),Sf(i),f(i),go] = iterative_solver( p,x(i-1),x(i),y(i-1),z(i-1),A(i-1),P(i-1),R(i-1),V(i-1),Sf(i-1),f(i-1),B(i),gradB(i),go );
  if (go==0)
    disp("------- Break out")
    break;
  endif
  y(i) = ynew ;
  Fr(i)  = V(i) / sqrt( p.g * y(i) );
  q = p.Q / B(i) ;
  yc(i) = ( q^2 / p.g )^(1/3) ;
end
##
##disp("Computing the normal depth")
##for i=1:ind(5)
##  yn(i) = 1.1 * yc(i) ; % Initial value
##  fun = @(x)(p.Q * p.n - sqrt(p.S0) * B(i) * x * (B(i)*x/(B(i)+2*x))^(2/3)) ;
##  yn(i) = fzero(fun,yn(i));
##endfor
##
##tosave = [x' y' yc' yn' Fr'];
len = length(y) ;

figure(1)
plot(x(1:len) / y(1), y / y(1)    ,'b'); grid on; hold on ;
%plot(x(1:ind(2)) , ydiv  , 'r') ;
##plot(x(1:ind(1)) , hconj , 'r--') ;
plot(x(1:len) / y(1), yc / y(1)    , 'k--') ;
##plot(x(1:ind(1)) , yn    , 'r--') ;


figure(2)
plot(x(1:len) / y(1), Fr     ,'b'); grid on;

tosave = [ x(1:len)' , y' , Fr' ];