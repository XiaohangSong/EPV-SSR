function Data=compressed_sensing_EM_VMP_EPV_modeling(M,b,SNR,damp,sim_frames_str)
%  compressed_sensing_EM_VMP_EPV_modeling carries a simluation on a sparse
%  bayesian learning problem with unknown priors using the EM-VMP-EPV
%  (Algorithm4) in paper
%
%  "Unifying Message Passing Algorithms Under the Framework of Constrained Bethe Free Energy Minimization"
%   Authors: Dan Zhang, Xiaohang Song, Wenjin Wang, Gerhard Fettweis, Xiqi Gao
%   https://arxiv.org/abs/1703.10932
%
%   The input variables are set as sequences for easy deployments on a compute
%   farm. The calling syntax is:
%   Data=compressed_sensing_EM_VMP_EPV_modeling('M','b','SNR_dB','damp','sim_frames_str');
%
%   Data contains the simulation results, including
%        Data.sim_frames: Number of frames being successfully simulated.
%               Data.mse: Mean squre error (MSE) w.r.t. the number of iterations.
%         Data.mse_final: The final MSE;
%              Data.iter: The mean number of iterations applied;
%
%   'M': Number of unknown variables, and the number of observations is set
%        to half of this value. This ratio can be changed within the implementation;
%   'b': Sparsity, i.e., the ratio between the number of non-zero parameters
%        of insterest and its total number;
%   'SNR': Signal-to-noise ratio in dB;
%   'damp': Damping factor;
%   'sim_frames_str': The desired maximum number of frames per simulation.
%        For earlier simulation termination, the actual number of frames can be smaller.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Copyright (C) 2019~2021 Xiaohang Song
%  Vodafone Chair, Technische Universitaet Dresden.
%  Last updated on Jan. 17, 2021
%
%  This code is distributed for the users to reproduce and extend the results. 
%  When reporting, reproducing or extending the results, please cite the paper above.
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

foldername = ['./sim_results/2021-01-17/UnknownModel/EM-VMP-EPV/'];
if isdir(foldername) == 0
    mkdir(foldername);
end

rng(12342,'twister');

M         = str2double(M);
dim_ratio = 0.5;
N         = ceil(M*dim_ratio);

b       = str2double(b);
SNR     = str2double(SNR); %dB;
mu_x0   = zeros(M,1);
tau_x0  = ones(M,1);

schedule='Pal'; %Sceduling operation for further algorithm extention;
rho     = '0'; %Correlation bettwen elements for further algorithm extention;

filename = [foldername 'CS-BernolliGauss--Rho-' num2str(b) '-correlation-' rho '-SNR-' num2str(SNR) 'dB-M' num2str(M) '-N' num2str(N) '-schedule-' schedule '-Damp-' damp '.mat'];
damp    = str2double(damp);

Ex   = b*(abs(mu_x0).^2+tau_x0);
N0   = (sum(Ex)/N)/power(10,SNR/10);
Nr_iter = 500;
mse    = zeros(1,Nr_iter);

mse_final  = 0;
sim_frames = str2double(sim_frames_str);

convergence_analysis = zeros(1,sim_frames);
for n_frame = 1:sim_frames
    
    flag = rand(1,M);
    x    = sqrt(0.5*tau_x0).*complex(randn(M,1),randn(M,1)) + mu_x0;
    %x    = sqrt(tau_x0).*randn(M,1) + mu_x0;%complex(randn(M,1),randn(M,1)) + mu_x0;
    x(flag < 1-b) = 0;
    w    = sqrt(N0/2)*complex(randn(N,1),randn(N,1));
    
    A     = sqrt(0.5/N)*complex(randn(N,M),randn(N,M));
    
    y     = A*x + w;
    S     = abs(A)'.^2; %simplify the computionation later
    beta  = zeros(N,1);
    
    lambda  = 100/var(y);
    theta   = 1/M*ones(M,1);
    mean_x = zeros(M,1);
    est_x  = zeros(M,1);
    var_x  = theta;
    epsilon = 1.5;
    eta     = 1;
    e_old  = 0;
    support_x1 = M;
    
    tau_t = zeros(M,N+1);
    tau_0 = zeros(M,1);
    tmp_mse = zeros(1,Nr_iter);
    
    for iter = 1:Nr_iter
        
        tau_t(:,1+N) = 1./var_x;                                                                %Alg2-Step 3.1
        
        %tic
        tau_t(:,1+N)=tau_t(:,1+N)-tau_0;                                                        %Alg2-Step 3.2
        tau_n_mtx = max(1e-8,repmat(sum(tau_t,2),1,N) - tau_t(:,1:N));                          %Alg2-Step 4
        gamma_n   = sum(S./tau_n_mtx(:,1:N),1)';                                                %Alg2-Step 5
        mu_z_n    = A*mean_x - beta.*gamma_n;                                                   %Alg2-Step 6
        
        var_z_n  = gamma_n./(1 + lambda*gamma_n);                                               %Alg4-Step 6.0
        beta      = (1-damp)*beta + damp*(y - mu_z_n)./(1/lambda + gamma_n);                    %Alg2-Step 7
        tau_t(:,1:N) = transpose(S'./(repmat(gamma_n,1,M) + 1/lambda-S'./transpose(tau_n_mtx)));%Alg2-Step 8
        %toc
        
        tau_0 = max(1e-8,sum(tau_t,2) - tau_t(:,N+1));                                          %Alg2-Step 9
        mu_x  = mean_x + (A'*beta)./tau_0;                                                      %Alg2-Step 10
        var_x = 1./tau_0;        
        mean_x  = (1-damp)*mean_x+ damp*((mu_x.*theta)./(var_x+theta));                         %Alg2-Step 11
        var_x   = (theta.*var_x)./(var_x+theta);                                                %Alg2-Step 11
 
        tmp_mse(iter:end) = sum(abs(mean_x-x).^2)/sum(Ex);
        nse = sum(abs(mean_x-x).^2)/sum(Ex);
 
        e_new = (y-A*mean_x)'*(y-A*mean_x)/N;                                                     %Alg4-Step 5
        lambda = 1/(e_new + mean(var_z_n));                                                       %Alg4-Step 6
        theta  = min(max(0,(epsilon-2)/(2*eta) + sqrt((epsilon-2)^2+4*eta*(abs(mean_x).^2+var_x))/(2*eta)),1e8);
        %Alg4-Step 7
        
        support_x = sum(abs(mean_x) > 1e-8);
        
        if abs(e_new - e_old) < 1e-6 && support_x >= support_x1
            epsilon = 0.95*epsilon;                                                               %Alg4-Step 8
        end
        support_x1 = support_x;
        e_old = e_new;                                                                            %Alg4-Step 9
        
        if ~isdeployed
            if max(abs(est_x-mean_x)) < 1e-5
                break
            end
            est_x = mean_x;
        end
        
    end
    mse        = mse + tmp_mse;
    mse_final = mse_final + sum(abs(mean_x-x).^2)/sum(Ex);
    convergence_analysis(n_frame)=iter;
    if ~isdeployed
        if mod(n_frame,50)==0
            fprintf('EM_VMP_EPV:rho=%.2f, Nr_frames %i, MSE: %.2f %.2f %.2f, Iteration Num %.2f\n'...
                ,b, n_frame,10*log10(mse(end)/n_frame),10*log10(mse_final/n_frame),10*log10(nse), mean(convergence_analysis(1:n_frame)));
        end
    end
    
    if mod(n_frame,100) == 0
        Data.sim_frames = n_frame;
        Data.mse  = mse/n_frame;
        Data.mse_final = mse_final/n_frame;
        Data.iter = mean(convergence_analysis(1:n_frame));
        save(filename,'Data');
    end
end

end
