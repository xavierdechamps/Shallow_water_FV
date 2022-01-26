clear,clc,close all
format long e

B = 2;
g = 9.81;
Q = 8;
n = 0.00708;
S0= 0.002;
s = 0;
eps = 1.e-6;

y(1) = 3.4;
x    = linspace (1000,000,1001) ;
z(1) = 0;

A(1)   = B*y(1) + s*y(1)^2 ; % Area
P(1)   = B + 2*y(1)*sqrt(1+s^2); % Wetted perimeter
R(1)   = A(1)/P(1) ;
V(1)   = Q / A(1) ;       % Q / A
kin(1) = 0.5*V(1)^2 /g ; % V^2 / 2g
H(1)   = kin(1) + y(1) + z(1) ;
Sf(1)  = n^2 * V(1)^2 / ( R(1)^(4/3) ) ;
f(1)   = ( S0 - Sf(1) ) / ( 1-Q^2 * B / (g*A(1).^3) ) ;

%%%% 

for i = 2:length(x)
  
 % Predictor
 y(i) = y(i-1) + f(i-1) * (x(i)-x(i-1)) ;
 
 z(i) = z(i-1) - S0*(x(i)-x(i-1)) ;
 A(i)   = B*y(i) + s*y(i)^2 ; % Area
 P(i)   = B + 2*y(i)*sqrt(1+s^2); % Wetted perimeter
 R(i)   = A(i)/P(i) ;
 V(i)   = Q / A(i) ;       % Q / A
 kin(i) = 0.5*V(i)^2 /g ; % V^2 / 2g
 H(i)   = kin(i) + y(i) + z(i) ;
 Sf(i)  = n^2 * V(i)^2 / ( R(i)^(4/3) ) ;
 f(i)   = ( S0 - Sf(i) ) / ( 1-Q^2 * B / (g*A(i).^3) ) ;

 ynew = y(i-1) + (f(i-1)+f(i))*0.5 * (x(i)-x(i-1)) ;
 while ( abs( y(i) - ynew ) > eps )
  % Corrector
   y(i) = ynew ;
   z(i) = z(i-1) - S0*(x(i)-x(i-1)) ;
   A(i)   = B*y(i) + s*y(i)^2 ; % Area
   P(i)   = B + 2*y(i)*sqrt(1+s^2); % Wetted perimeter
   R(i)   = A(i)/P(i) ;
   V(i)   = Q / A(i) ;       % Q / A
   kin(i) = 0.5*V(i)^2 /g ; % V^2 / 2g
   H(i)   = kin(i) + y(i) + z(i) ;
   Sf(i)  = n^2 * V(i)^2 / ( R(i)^(4/3) ) ;
   f(i)   = ( S0 - Sf(i) ) / ( 1-Q^2 * B / (g*A(i).^3) ) ;

   ynew = y(i-1) + (f(i-1)+f(i))*0.5 * (x(i)-x(i-1)) ;
   
 endwhile
 
 y(i) = ynew ;
 
end

ii = imag(y)==0 ; % select only physical solution where water depth > 0
xx = x(ii)' ;
yy = y(ii)' ;
zz = z(ii)' ;

plot(xx,yy+zz,'b'); hold on ; plot(x,z,'k'); grid on
%axis( [ 0 1000 0 3.25 ] )

step = -100;
%[  xx(1:step:end) yy(1:step:end) ]
[  xx(end:step:1) yy(end:step:1) ]