function[like]=Likelihood(theta,beta,V,data)
%% Log-Likelihood function
% PURPOSE: Compute the Likelihood Function
%--------------------------------------------------------------------
% USAGE: [like] = Likelihood(theta,beta,V, data)
% where: theta  = 2 x 1 parameter vector to be estimated 
%        beta   = discount rate (=0 for myopic case, 0.95 for NFXP)
%        V      = 1 x 5 vector of estimated V(s) value funtions from the
%                 NFXP Algorithm in the Main File "ECON890_HW4_Ackerman.m"
%        data   = 2400 x 4 data matrix
%--------------------------------------------------------------------
% RETURNS: like = the negative of the likelihood function
%--------------------------------------------------------------------
% REFERENCE:
% Rust, John (1987)., Optimal Replacement of GMC Bus Engines: An Empirical
% Model of Harold Zucher., Econometrica., Vol. 55., No. 5 (Sep., 1987). 
% pp 999-1033
% -------------------------------------------------------------------
% Written by Robert Ackerman and Anonymous 

like=0;

for i = 1:100; % bus loop
    for t = 1:24; % time loop
        if data((t + (i-1)*24),3) == 5 % for s=5 case, we don't want to add in the value function since it is never the case that s=6
        l = (1-data((t + (i-1)*24),4)) * log(exp(- theta(1,1) *data((t + (i-1)*24),3) + (beta * V(data((t + (i-1)*24),3)))) ... 
        /(exp(- theta(1,1) *data((t + (i-1)*24),3) + (beta * V(data((t + (i-1)*24),3)))) + exp(-theta(2,1) + (beta*V(1,1))))) ...
        + data((t + (i-1)*24),4) * log(exp(-theta(2,1) + (beta*V(1,1))) ... 
        /(exp(- theta(1,1) *data((t + (i-1)*24),3) + (beta * V(data((t + (i-1)*24),3)))) + exp(-theta(2,1) + (beta*V(1,1)))));
        else % for the s=1-4 case we want to add since the bus ages, since we want the correct value function to account for the aging process
        l = (1-data((t + (i-1)*24),4)) * log(exp(- theta(1,1) *data((t + (i-1)*24),3) + (beta * V(data((t + 1 + (i-1)*24),3)))) ... 
        /(exp(- theta(1,1) *data((t + (i-1)*24),3) + (beta * V(data((t + 1 + (i-1)*24),3)))) + exp(-theta(2,1) + (beta*V(1,1))))) ...
        + data((t + (i-1)*24),4) * log(exp(-theta(2,1) + (beta*V(1,1))) ... 
        /(exp(- theta(1,1) *data((t + (i-1)*24),3) + (beta * V(data((t + 1 + (i-1)*24),3)))) + exp(-theta(2,1) + (beta*V(1,1)))));
        end
    like = like + l;
    end
end
like = -like; % to min rather than max using fminunc    