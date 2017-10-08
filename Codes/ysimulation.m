close all
clear all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
A=[0.678,1.1];
Q=[0.977,0.023;0.074,0.926];

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]);

psd=100;
while psd(1)>1.8||psd(2)>1.8
A=0.5*(A+1);

cons(:,:,1) = k_mat .^ alpha.*A(1) + (1 - delta) * k_mat - k_mat'; 
cons(:,:,2) = k_mat .^ alpha.*A(2) + (1 - delta) * k_mat - k_mat'; 

ret = cons .^ (1 - sigma) ./ (1 - sigma); % return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret(cons < 0) = -Inf;

%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(1, num_k);
for i=1:2
while dis > tol
    % compute the utility value for all possible combinations of k and k':
    value_mat = ret(:,:,i) + beta * repmat(v_guess, [num_k 1]);
    
    % find the optimal k' for every k:
    [vfn(i,:), I] = max(value_mat, [], 2);
    
    % what is the distance between current guess and value function
    dis = max(abs(vfn(i,:) - v_guess));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn(i,:);
end
pol_indx(:,i)=I;
dis=1;
end

gl = k(pol_indx(:,1));% policy function
gh = k(pol_indx(:,2));

k_star=((1+delta*beta-beta)./(beta*alpha.*A)).^(1/(alpha-1));
[k_a,indexl]=min(abs(k_star(1)-k));
[k_a,indexh]=min(abs(k_star(2)-k));
k_a=[gl(indexl),gh(indexh)];
At(1,1)=A(1);
At(2,1)=A(2);
for j = 1:2
    for i=2:5000
    if At(j,i-1)==A(1)
        z=binornd(1,Q(2,1));
        if z==1
            At(j,i)=A(1);
        else
            At(j,i)=A(2);
        end
    else
        z=binornd(1,Q(1,1));
        if z==1
            At(j,i)=A(2);
        else
            At(j,i)=A(1);
        end
    end
    end
end
Y=zeros(2,5000);
for i = 1:2
    K=k_a(i);
    for j = 1:5000
        [c index] = min(abs(K-k));
        if At(i,j)==A(1)
            Y(i,j)=A(1)*(k(index))^(alpha);
            K=gl(index);
        else
            Y(i,j)=A(2)*(k(index))^(alpha);
            K=gh(index);
        end
    end
end
psd=100.*std(Y,0,2)./mean(Y,2);
end