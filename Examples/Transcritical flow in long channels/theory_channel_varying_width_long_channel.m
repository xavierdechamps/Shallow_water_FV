clear,clc,close all
format long e

p.B1 = 40;
p.B2 = 26.5;
p.g = 9.81;
p.Q = 500;
p.n = 0.0389;
p.S0= 0.002;
p.s = 0;
p.eps = 1.e-6;

Lu = 1000; % Length upstream
Ld = 500; % Length downstream
Lb = 10; % Length of narrowest section
Le = 20 ; % Length of expansion / contraction
lengL = Lu+Ld;

q = p.Q / p.B1 ;
yc(1) = ( q^2 / p.g )^(1/3) ;

%    1 ---- Ld ---- xx --- Le --- xx --- Lb --- xx --- Le --- xx --- Lu --- Last Point
%ind(1)           ind(4)        ind(3)        ind(2)        ind(1)
% Downstream section
y(1) = yc(1) * 1.01;
z(1) = 0;
numsteps = 501;
x     = linspace (lengL,lengL-Ld+Le+Lb*0.5,numsteps) ;
B     = linspace (p.B1,p.B1,numsteps);
gradB = linspace (0,0,numsteps);
gradB(end) = 0.5*(p.B1 - p.B2)/Le ; % Average between gradient and parallel
ind(1) = numsteps ; 

% Divergent section
numsteps = 301;
x(end:end+numsteps-1)     = linspace (Lu+Le+Lb*0.5,Lu+Lb*0.5,numsteps) ;
B(end:end+numsteps-1)     = linspace (p.B1,p.B2,numsteps) ;
gradB(end+1:end+numsteps-1) =  (p.B1 - p.B2)/Le ;
gradB(end) = 0.5*(p.B1 - p.B2)/Le ; % Average between gradient and parallel
ind(2) = ind(1) + numsteps - 1; 

% Narrowest section
numsteps = 101;
x(end:end+numsteps-1)     = linspace (Lu+Lb*0.5,Lu-Lb*0.5,numsteps) ;
B(end:end+numsteps-1)     = linspace (p.B2,p.B2,numsteps) ;
gradB(end+1:end+numsteps-1) = 0.;
gradB(end) = 0.5*(p.B2 - p.B1)/Le ; % Average between gradient and parallel
ind(3) = ind(2) + numsteps - 1; 

% Convergent section
numsteps = 301;
x(end:end+numsteps-1)     = linspace (Lu-Lb*0.5,Lu-Lb*0.5-Le,numsteps) ;
B(end:end+numsteps-1)     = linspace (p.B2,p.B1,numsteps) ;
gradB(end+1:end+numsteps-1) = (p.B2 - p.B1)/Le;
gradB(end) = 0.5*(p.B2 - p.B1)/Le ; % Average between gradient and parallel
ind(4) = ind(3) + numsteps - 1; 

% Upstream section
numsteps = 501;
x(end:end+numsteps-1)     = linspace (Lu-Lb*0.5-Le,0.,numsteps) ;
B(end:end+numsteps-1)     = linspace (p.B1,p.B1,numsteps) ;
gradB(end+1:end+numsteps-1) = 0.;
ind(5) = ind(4) + numsteps - 1; 


%
A(1)   = B(1)*y(1) + p.s*y(1)^2 ; % Area
P(1)   = B(1) + 2*y(1)*sqrt(1+p.s^2); % Wetted perimeter
R(1)   = A(1)/P(1) ;
V(1)   = p.Q / A(1) ;       % Q / A
Sf(1)  = p.n^2 * V(1)^2 / ( R(1)^(4/3) ) ;
f(1)   = ( p.S0 - Sf(1) + p.Q^2 * y(1)* gradB(1) / (p.g*A(1).^3) ) / ( 1 - p.Q^2 * B(1) / (p.g*A(1).^3) ) ;
Fr(1)  = V(1) / sqrt( p.g * y(1) );

%%%%  Downstream section
disp("Downstream section")
for i = 2:ind(1)
  go=1;
  [ynew,z(i),A(i),P(i),R(i),V(i),Sf(i),f(i),go] = iterative_solver( p,x(i-1),x(i),y(i-1),z(i-1),A(i-1),P(i-1),R(i-1),V(i-1),Sf(i-1),f(i-1),B(i),gradB(i),go );
  y(i) = ynew ;
  Fr(i)  = V(i) / sqrt( p.g * y(i) );
  q = p.Q / p.B1 ;
  yc(i) = ( q^2 / p.g )^(1/3) ;
end

%%%% Divergent section
disp("Divergent section")
for i = ind(1)+1:ind(2)
  go=1;
  [ynew,z(i),A(i),P(i),R(i),V(i),Sf(i),f(i),go] = iterative_solver( p,x(i-1),x(i),y(i-1),z(i-1),A(i-1),P(i-1),R(i-1),V(i-1),Sf(i-1),f(i-1),B(i),gradB(i),go );
  if (go==0)
    break;
  endif
  y(i) = ynew ;
  Fr(i)  = V(i) / sqrt( p.g * y(i) );
  q = p.Q / B(i) ;
  yc(i) = ( q^2 / p.g )^(1/3) ;
end

ii = i-1;
if (ii<ind(2) && ii>ind(1))
  disp("Trying to compute the shock... ")
  
  q = p.Q / B(ind(2)) ;
  yc(ind(2)) = ( q^2 / p.g )^(1/3) ;
  ydiv(ind(2))    = yc(ind(2)) *0.995 ; % Slightly less than critical depth, otherwise difficulty to find the correct value
  z(ind(2))  = z(1) - p.S0*(x(ind(2))-x(1)) ;
  A(ind(2))  = B(ind(2))*ydiv(1) + p.s*ydiv(1)^2 ; % Area
  P(ind(2))  = B(ind(2)) + 2*ydiv(1)*sqrt(1+p.s^2); % Wetted perimeter
  R(ind(2))  = A(ind(2))/P(ind(2)) ;
  V(ind(2))  = p.Q / A(ind(2)) ;       % Q / A
  Sf(ind(2)) = p.n^2 * V(ind(2))^2 / ( R(ind(2))^(4/3) ) ;
  f(ind(2))  = 0. ;
  
  Fr(ind(2))  = V(ind(2)) / sqrt( p.g * ydiv(1) );
  hconj(ind(2)) = yc(ind(2));
  
  % Start from yc at throat downstream
  for i = ind(2)-1:-1:1
    go=1;
    [ynew,z(i),A(i),P(i),R(i),V(i),Sf(i),f(i),go] = iterative_solver( p,x(i+1),x(i),ydiv(i+1),z(i+1),A(i+1),P(i+1),R(i+1),V(i+1),Sf(i+1),f(i+1),B(i),gradB(i),go );
    if (go==0)
      break;
    endif
    ydiv(i) = ynew ;
    
    tmp  = V(i) / sqrt( p.g * ydiv(i) );
    hconj(i) = 0.5*ydiv(i)*(sqrt(1+8*tmp*tmp)-1);
    q = p.Q / B(i) ;
    yc(i) = ( q^2 / p.g )^(1/3) ;
  endfor
  ydiv(i:-1:1) = yc(i);
  hconj(i:-1:1) = yc(i);
endif

for i = 1:length(y)
  if (y(i)<hconj(i))
    disp("Switching to the conjugate height...")
    break;
  endif
endfor
y(i:ind(2)) = ydiv(i:ind(2));
Fr(i:ind(2))  = V(i:ind(2)) ./ sqrt( p.g * ydiv(i:ind(2)) );

  y(ind(2))    = yc(ind(2)) *1.001 ; % Slightly more than critical depth, otherwise difficulty to find the correct value
  z(ind(2))  = z(1) - p.S0*(x(ind(2))-x(1)) ;
  A(ind(2))  = B(ind(2))*y(1) + p.s*y(1)^2 ; % Area
  P(ind(2))  = B(ind(2)) + 2*y(1)*sqrt(1+p.s^2); % Wetted perimeter
  R(ind(2))  = A(ind(2))/P(ind(2)) ;
  V(ind(2))  = p.Q / A(ind(2)) ;       % Q / A
  Sf(ind(2)) = p.n^2 * V(ind(2))^2 / ( R(ind(2))^(4/3) ) ;
  f(ind(2))  = 0. ;
  
  Fr(ind(2))  = 1. ;
  hconj(ind(2)) = yc(ind(2));

%%%% Narrowest section upstreamwards
disp("From narrowest section to inlet")
for i = ind(2)+1:ind(5)
  go=1;
  [ynew,z(i),A(i),P(i),R(i),V(i),Sf(i),f(i),go] = iterative_solver( p,x(i-1),x(i),y(i-1),z(i-1),A(i-1),P(i-1),R(i-1),V(i-1),Sf(i-1),f(i-1),B(i),gradB(i),go );
  if (go==0)
    break;
  endif
  y(i) = ynew ;
  Fr(i)  = V(i) / sqrt( p.g * y(i) );
  q = p.Q / B(i) ;
  yc(i) = ( q^2 / p.g )^(1/3) ;
  
%  hconj(i) = 0.5*y(i)*(sqrt(1+8*Fr(i)*Fr(i))-1);
end

disp("Computing the normal depth")
for i=1:ind(5)
  yn(i) = 1.1 * yc(i) ; % Initial value
  fun = @(x)(p.Q * p.n - sqrt(p.S0) * B(i) * x * (B(i)*x/(B(i)+2*x))^(2/3)) ;
  yn(i) = fzero(fun,yn(i));
endfor

tosave = [x' y' yc' yn' Fr'];

figure(1)
plot(x(1:ind(5)) , y     ,'b'); grid on; hold on ;
%plot(x(1:ind(2)) , ydiv  , 'r') ;
plot(x(1:ind(2)) , hconj , 'r--') ;
plot(x(1:ind(5)) , yc    , 'k--') ;
plot(x(1:ind(5)) , yn    , 'r--') ;

axis([Lu-125 Lu+125 min(y) max(y)])

figure(2)
plot(x(1:ind(5)) , Fr     ,'b'); grid on;
axis([Lu-125 Lu+125 min(Fr) max(Fr)])