% Robert Ackerman and Anonymous

% Homework 4
% March 21, 2014
% Rust, John (1987)., Optimal Replacement of GMC Bus Engines: An Empirical
% Model of Harold Zucher., Econometrica., Vol. 55., No. 5 (Sep., 1987). 
% pp 999-1033
%% Step 1: Prelminary Settings
clear; clc;

% Load the data
data = csvread('bus_data_reshape_hw4_nolabel.csv');
% col 1 = bus, col 2 = time, col 3 = state, col 4 = action
% dimensions of data: 2400 x 4 

% myopic model
% beta = 0;
% forward looking model
beta = 0.95;

theta = [1 1]'; % 2 x 1 intitial parameter value guess
V = zeros(1,5); % 1 x 5 vector of initial value function guesses for the value function iteration process
Vplus = zeros(1,5); % 1 x 5 vector of value functions that updates in the value function iteration process
eta = 0.0000000000001; % tolerance for the value function loops and likelihood loops
%% Step 2: Summary Statistics

% mean of a for all bus and periods
mean(data(:,4))

% mean of s for all bus periods
mean(data(:,3))

% get mean of s for a=1
B = sortrows(data,4);
mean(B(1509:2400,3))

% get mean of s for a=0
mean(B(1:1508,3))

% get mean of a for each s
B = sortrows(data,3);

% get percent replaced replaced for all s (prob a=1 | s)
mean(B(1:876,4)) % percent replaced if s = 1
mean(B(877:1597,4)) % percent replaced if s = 2
mean(B(1598:2078,4)) % percent replaced if s = 3
mean(B(2079:2306,4)) % percent replaced if s = 4
mean(B(2307:2400,4)) % percent replaced if s = 5

%% Step 3: Value Loop and Likelihood Estimation
for i =1:100 % maximum number of likelihood iterations
    for r=1:2000 % maximum number of value function iterations
        for s = 5:-1:1 % for each s 
            if s==5
                Vplus(1,s) = log(exp(- theta(1,1) * s + beta * V(1,s)) + exp(-theta(2,1) + beta * V(1,1)));
            else
                Vplus(1,s) = log(exp(-theta(1,1) * s + beta * V(1,s+1)) + exp(-theta(2,1) + beta * V(1,1)));
            end
        end
        if abs((Vplus(1,1)- V(1,1))) <= eta ... % check if change is below tolerance
            && abs((Vplus(1,2)- V(1,2))) <= eta ...
            && abs((Vplus(1,3)- V(1,3))) <= eta ...
            && abs((Vplus(1,4)- V(1,4))) <= eta ...
            && abs((Vplus(1,4)- V(1,4))) <= eta ...
        ,break 
        else
            V(1,1)=Vplus(1,1);
            V(1,2)=Vplus(1,2);
            V(1,3)=Vplus(1,3);
            V(1,4)=Vplus(1,4);
            V(1,5)=Vplus(1,5);
        end
    end


% maximum likelihood part 
% calls fminunc to minimize the user written function  "Likelihood.m"
% note: the @(x) part creates an "anonymous function" that allows for
% fminunc to call "Likelihood.m", which has multiple inputs, but minimize
% it with respect to just theta
options = optimoptions('fminunc', 'TolFun', 0.0000000001, 'TolX', 0.0000000001);
options.Display ='iter'; % this sets it so matlab reports each iteration while minimizing the likelihood function
[mintheta, fval, exitflag, output,hessian] = fminunc(@(x) Likelihood(x,beta,V,data),theta,options);
if (abs((theta(1,1)-mintheta(1,1))/(theta(1,1))) <=eta) && (abs((theta(1,1)-mintheta(1,1))/(theta(1,1))) <= eta)
   break
else
end
theta=mintheta;   
end
%% Step 4: Report the Estimated Probabilities 
for s = 5:-1:1
    if s==5
        p0(1,s) = exp(-theta(1,1) * s + (beta * V(1,s))) / (exp(-theta(1,1) * s + (beta * V(1,s))) + exp(-theta(2,1) + (beta * V(1,1))));
        p1(1,s) = exp(-theta(2,1) + (beta * V(1,1))) / (exp(-theta(1,1) * s + (beta * V(1,s))) + exp(-theta(2,1) + (beta * V(1,1))));
    else
        p0(1,s) = exp(-theta(1,1) * s + (beta * V(1,s+1))) / (exp(-theta(1,1) * s + (beta * V(1,s+1))) + exp(-theta(2,1) + (beta * V(1,1))));
        p1(1,s) = exp(-theta(2,1) + (beta * V(1,1))) / (exp(-theta(1,1) * s + (beta * V(1,s+1))) + exp(-theta(2,1) + (beta * V(1,1))));
    end
end



