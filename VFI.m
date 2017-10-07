close all
clear all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
A=[1.1,0.678];
Q=[0.023,0.977;0.074,0.926];

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
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

gh = k(pol_indx(:,1));% policy function
gl = k(pol_indx(:,2));

plot(k,vfn)
legend('High State','Low State','Location','southeast')
figure
plot(k,[gh;gl])
legend('High State','Low State','Location','southeast')

P=Q^1000; % invariant distribution
z=random(makedist('Binomial','N',1,'P',P(1,1)),1000,1); % At generator
for i=1:1000
    if z(i)==1
         At(i)=A(1);
    else At(i)=A(2);
    end
end
k_star=((1+delta*beta-beta)/(beta*alpha*A(1)))^(1/(alpha-1));
[k_a,index]=min(abs(k_star-k));
z=1;
for i=2:1000
    if z==1
        k_a(i)=gh(index(i-1));
    else
        k_a(i)=gl(index(i-1));
    end
    index(i)=find(k==k_a(i));
    z=random(makedist('Binomial','N',1,'P',Q(z+1,1)),1,1);
end
y=k_a.^alpha.*At;
100*std(y)/mean(y)  %% I've tried to calibrate Ah, but even if I decrease Ah to 1.001, the standard deviation is 2.9%.