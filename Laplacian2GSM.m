% Gaussian Mixture Expectation Maximization

clc, clear all, close all;

Ndim = 4;

mx_0 = zeros(Ndim,1);
Pbar_0 = diag([1e-2*ones(1,Ndim)]);

% Gaussian Mixture Expectation Maximization

N_training = 10000;

x_initial = mx_0(1:Ndim) + chol(Pbar_0(1:Ndim,1:Ndim)/2)*randl(Ndim,N_training);

K = 20; % Number of Gaussians in the GMM

mu_k = randn(Ndim,K);
P_k = repmat(eye(Ndim),[1,1,K]);
w_k = 1/K*ones(K,1);

Initial_LLF_inner = zeros(N_training,1);
for n = 1:N_training
    for k = 1:K
        Initial_LLF_inner(n) = Initial_LLF_inner(n) + w_k(k)*mvnpdf(x_initial(:,n),mu_k(:,k),P_k(:,:,k));
    end
end
Initial_LLF = sum(log(Initial_LLF_inner));
Prev_LLF = Initial_LLF;
counter = 0;
while(true)
    counter = counter + 1
    % E-step
    for k = 1:K
        for n = 1:N_training
            gamma_prior(k,n) = w_k(k)*mvnpdf(x_initial(:,n),mu_k(:,k),P_k(:,:,k));
        end
    end
    gamma = gamma_prior./sum(gamma_prior(:,:),1);
    
    % M-step
    mu_k = zeros(Ndim,K);
    P_k = zeros(Ndim,Ndim,K);
    w_k = zeros(K,1);
    for k = 1:K
        for n = 1:N_training
            mu_k(:,k) = mu_k(:,k) + gamma(k,n)*x_initial(:,n);
            P_k(:,:,k) = P_k(:,:,k) + gamma(k,n)*(x_initial(:,n) - mu_k(:,k))*(x_initial(:,n) - mu_k(:,k)).';
            w_k(k) = w_k(k) + gamma(k,n);
        end
        mu_k(:,k) = mu_k(:,k)/sum(gamma(k,:),2);
        P_k(:,:,k) = P_k(:,:,k)/sum(gamma(k,:),2);
        w_k(k) = w_k(k)/N_training;
    end
    
    Current_LLF_inner = zeros(N_training,1);
    for n = 1:N_training
        for k = 1:K
            Current_LLF_inner(n) = Current_LLF_inner(n) + w_k(k)*mvnpdf(x_initial(:,n),mu_k(:,k),P_k(:,:,k));
        end
    end
    Current_LLF = sum(log(Current_LLF_inner))
    if abs(Current_LLF - Prev_LLF)/abs(Prev_LLF) <= 0.01
        break;
    end
    Prev_LLF = Current_LLF;
    
end