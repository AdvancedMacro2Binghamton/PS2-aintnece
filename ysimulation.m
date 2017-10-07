close all
clear all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
A=[1.1,0.678];
Q=[0.977,0.023;0.074,0.926];

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]);

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