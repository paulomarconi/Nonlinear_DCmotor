function [sys, x0] = relay_slab(t,x,u,flag)
%
% Relay nonlinear model, configured a la SFUNC(T,X,U,FLAG,P1,P2, ....)
% for use with SimuLink ODE integrators
% 
% JH Taylor, 21 February 1995
%
% If no right hand arguments then return the number of states, outputs 
% and inputs.
if nargin==0,  sys=[2,0,0,0,0,0]; x0 = zeros(2,1); return, end
%
if abs(flag) == 1 
    sys(1) = x(2);
    sys(2) = -sign(x(1));
elseif abs(flag) == 0
    sys=[2,0,0,0,0,0]; x0 = zeros(2,1);
else
    sys = [];                   % 'roots' and 'feed-thru' flag (ignored).
end

