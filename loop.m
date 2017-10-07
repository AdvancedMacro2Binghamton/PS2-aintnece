close all
clear all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
A=[1.1,0.678];

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k+1);
k = k(2:1001);
%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v=zeros(num_k,2,2);
v(:,2,:)=10;
for i=1:2
    p=1;
    v(:,end,i)=10;
while norm(v(:,end,i)-v(:,end-1,i))>=tol
    p=p+1;
    m=zeros(num_k,num_k);
    for j=1:num_k
        for h=1:num_k
            a=A(i)*k(j)^alpha+(1-delta)*k(j)-k(h);
            if a>=0
                m(j,h)=a^(1-sigma)/(1-sigma)+beta*v(h,p-1,i);
            else m(j,h)=-1000;
            end
        end
    end
    [v(:,p,i),Ir]=max(m,[],2);
end
pol_indx(:,i)=Ir;
V(:,i)=v(:,p,i);
end
gh = k(pol_indx(:,1));% policy function
gl = k(pol_indx(:,2));

plot(k,V)
legend('High State','Low State','Location','southeast')
figure
plot(k,[gh;gl])
legend('High State','Low State','Location','southeast')