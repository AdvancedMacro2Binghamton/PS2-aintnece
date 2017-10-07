clear all
close all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
pi=[0.977, 0.023; 0.926, 0.074];
A=[0.678;1.1];
%%%% Set up discretized state space
k_min = 0;
k_max = 1.1*(alpha*A(2)/(1/beta-1+delta))^(alpha/(1-alpha))+(1-delta)*(alpha*1.1/(1/beta-1+delta))^(1/(1-alpha));;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k+1);
k(1)=[];
k_mat = repmat(k', [1 num_k]); % this will be useful in a bit
sd=1;
while sd(1)>0.018 | sd(2)>0.018
 A=0.5*(A+1);
%Simulate A
consl = A(1)*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 
consh = A(2)*k_mat .^ alpha + (1 - delta) * k_mat - k_mat';
retl = consl .^ (1 - sigma) / (1 - sigma);
reth = consh .^ (1 - sigma) / (1 - sigma);
retl(consl < 0) = -Inf;
reth(consh < 0) = -Inf;
%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(2, num_k);
while dis > tol
    % compute the utility value for all possible combinations of k and k':
    value_matl = retl + beta *(pi(2,1)* repmat(v_guess(1,:), [num_k 1])+pi(2,2)*repmat(v_guess(2,:), [num_k 1]));
    value_math = reth + beta *(pi(1,2)* repmat(v_guess(1,:), [num_k 1])+pi(1,1)*repmat(v_guess(2,:), [num_k 1]));
    
    % find the optimal k' for every k:
    [vfnl, pol_indxl] = max(value_matl, [], 2);
    vfnl = vfnl';
    [vfnh, pol_indxh] = max(value_math, [], 2);
    vfnh = vfnh';
    % what is the distance between current guess and value function
    dis = [max(abs(vfnl - v_guess(1,:))) ; max(abs(vfnh - v_guess(2,:)))] ;
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess= [vfnl;vfnh];
end

gl = k(pol_indxl); 
gh = k(pol_indxh); % policy function
[c indexl] = min(abs(k-gl));
[c indexh] = min(abs(k-gh));
k_ss=[k(indexl);k(indexh)];


A_bar=zeros(2,5000);
A_bar(1,1)=A(1);
A_bar(2,1)=A(2);
for j = 1:2
    for i=2:5000
    if A_bar(j,i-1)==A(1);
        R=binornd(1,pi(2,1));
        if R==1
            A_bar(j,i)=A(1);
        else
            A_bar(j,i)=A(2);
        end
    else
        R=binornd(1,pi(1,1));
        if R==1
            A_bar(j,i)=A(2);
        else
            A_bar(j,i)=A(1);
        end
    end
    end
end
Y=zeros(2,5000);
for i = 1:2
    K=k_ss(i);
    for j = 1:5000
        [c index] = min(abs(K-k));
        if A_bar(i,j)==A(1)
            Y(i,j)=A(1)*(k(index))^(alpha);
            K=gl(index);
        else
            Y(i,j)=A(2)*(k(index))^(alpha);
            K=gh(index);
        end
    end
end
sd=std(Y,0,2)
end