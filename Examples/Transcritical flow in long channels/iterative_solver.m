function [ynew,z2,A2,P2,R2,V2,Sf2,f2,go] = iterative_solver(p,x1,x2,y1,z1,A1,P1,R1,V1,Sf1,f1,B,gradB,go)
  % Predictor
  y2 = y1 + f1 * (x2-x1) ;
  
  z2   = z1 - p.S0*(x2-x1) ;
  A2   = B*abs(y2) + p.s*y2^2 ; % Area
  P2   = B + 2*abs(y2)*sqrt(1+p.s^2); % Wetted perimeter
  R2   = A2/P2 ;
  V2   = p.Q / A2 ;       % Q / A
  Sf2  = p.n^2 * V2^2 / ( R2^(4/3) ) ;
  f2   = ( p.S0 - Sf2 + p.Q^2 * y2 * gradB / (p.g*A2^3) ) / ( 1-p.Q^2 * B / (p.g*A2^3) ) ;

  ynew = y1 + (f1+f2)*0.5 * (x2-x1) ;
  ncurrent = 1;
  nmax = 100;
  while ( abs( y2 - ynew ) > p.eps )
   % Corrector
    y2 = ynew ;
    z2 = z1 - p.S0*(x2-x1) ;
    A2   = B*abs(y2) + p.s*y2^2 ; % Area
    P2   = B + 2*abs(y2)*sqrt(1+p.s^2); % Wetted perimeter
    R2   = A2/P2 ;
    V2   = p.Q / A2 ;       % Q / A
    Sf2  = p.n^2 * V2^2 / ( R2^(4/3) ) ;
    f2   = ( p.S0 - Sf2 + p.Q^2 * y2 * gradB / (p.g*A2^3) ) / ( 1-p.Q^2 * B / (p.g*A2^3) ) ;
        
    ynew = y1 + (f1+f2)*0.5 * (x2-x1) ;
    ncurrent = ncurrent + 1 ;
    if (ncurrent > nmax)
      disp(strcat("Didn't converge at position ",num2str(x2)))
      go=0;
      break
    endif
  endwhile
 
endfunction