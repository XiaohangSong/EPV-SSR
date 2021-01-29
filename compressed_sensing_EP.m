function [Data]=compressed_sensing_EP(M,b,SNR, damp, sim_frames_str)
%  compressed_sensing_EP carries a simluation on a sparse signal
%  recovery problem using the Expectation Propagation (Algorithm1)
%  in paper
%
%  "Unifying Message Passing AlgorithmsUnder the Framework ofConstrained Bethe Free Energy Minimization"
%   Authors: Dan Zhang, Xiaohang Song, Wenjin Wang, Gerhard Fettweis, Xiqi Gao
%   https://arxiv.org/abs/1703.10932
%
%   The input variables are set as sequences for easy deployments on a compute
%   farm. The calling syntax is:
%   Data=compressed_sensing_EP('M','b','SNR_dB','damp','sim_frames_str');
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

damp_tmp    = str2double(damp);
if  damp_tmp==1
    foldername = ['./sim_results/2021-01-17/Alg1-EP/'];
else
    foldername = ['./sim_results/2021-01-17/Alg1-EP-damping/'];
end

if isdir(foldername) == 0
    mkdir(foldername);
end

rng(12342,'twister');

M         = str2double(M);  %BS anntenna number
dim_ratio = 0.5;
N         = ceil(M*dim_ratio); %UE number

b       = str2double(b);     %channel sparsity
SNR     = str2double(SNR); %dB;
mu_x0   = zeros(M,1);
tau_x0  = ones(M,1);
schedule='Pal'; %Sceduling operation for further algorithm extention;
rho     = '0'; %Correlation bettwen elements for further algorithm extention;

filename = [foldername 'CS-BernolliGauss--Rho-' num2str(b) '-correlation-' rho '-SNR-' num2str(SNR) 'dB-M' num2str(M) '-N' num2str(N) '-schedule-' schedule '-Damp-' damp '.mat'];
damp    = str2double(damp);

Ex   = b*(abs(mu_x0).^2+tau_x0);
N0   = (sum(Ex)/N)/power(10,SNR/10);
Nr_iter = 500;%Initial 500;
mse    = zeros(1,Nr_iter);

mse_final  = 0;
sim_frames = str2double(sim_frames_str);
convergence_analysis = zeros(1,sim_frames);

for n_frame = 1:sim_frames
    
    flag = rand(1,M);
    x    = sqrt(0.5*tau_x0).*complex(randn(M,1),randn(M,1)) + mu_x0;
    
    x(flag < 1-b) = 0;
    w    = sqrt(N0/2)*complex(randn(N,1),randn(N,1));
    
    A     = sqrt(0.5/N)*complex(randn(N,M),randn(N,M));
    
    y     = A*x + w;                                                                         %linear Model
    S     = abs(A)'.^2;
    beta  = zeros(N,1);                                                                      %Alg1-Step 1
    
    mean_x = zeros(M,1);
    est_x  = zeros(M,1);
    var_x  = b*tau_x0;
    
    tau_t = zeros(M,N+1);
    alpha_t=zeros(M,N+1);
    tau_0 = zeros(M,1);
    alpha_0=zeros(M,1);
    alpha_n_mtx=zeros(M,N);
    tau_n_mtx=zeros(M,N);
    tmp_mse = zeros(1,Nr_iter);
    
    for iter = 1:Nr_iter                                                                %Alg1-Step 2
        
        alpha_t(:,1+N) = mean_x./var_x;                                                 %Alg1-Step 3.1.1
        tau_t(:,1+N) = 1./var_x;                                                        %Alg1-Step 3.2.1
        
        %tic
        alpha_t(:,1+N)=alpha_t(:,1+N)-alpha_0;                                          %Alg1-Step 3.1.2
        tau_t(:,1+N)=tau_t(:,1+N)-tau_0;                                                %Alg1-Step 3.2.2
        
        if  damp==1
            alpha_n_mtx = repmat(sum(alpha_t,2),1,N) - alpha_t(:,1:N);                  %Alg1-Step 4.1
            tau_n_mtx = max(1e-8,repmat(sum(tau_t,2),1,N) - tau_t(:,1:N));              %Alg1-Step 4.2
        else
            alpha_n_mtx = (1-damp)*alpha_n_mtx +...
                damp*(repmat(sum(alpha_t,2),1,N) - alpha_t(:,1:N));                     %Alg1-Step 4.1
            tau_n_mtx = (1-damp)*tau_n_mtx + ...
                damp*max(1e-8,repmat(sum(tau_t,2),1,N) - tau_t(:,1:N));                 %Alg1-Step 4.2
        end
        gamma_n_m = S./tau_n_mtx;                                                       %Alg1-Step 5.1.0
        gamma_n   = sum(gamma_n_m,1)';                                                  %Alg1-Step 5.1
        mu_z_n_m    = A.*transpose(alpha_n_mtx./tau_n_mtx);                             %Alg1-Step 5.2.0
        mu_z_n    = sum(mu_z_n_m,2);                                                    %Alg1-Step 5.2
        
        denominator=repmat(gamma_n,1,M) + N0-transpose(gamma_n_m);                      %Alg1-Step 6.0
        alpha_t(:,1:N) = transpose(conj(A).*(repmat((y - mu_z_n),1,M)+...
            mu_z_n_m )./denominator);                                                   %Alg1-Step 6.1
        tau_t(:,1:N) = transpose(S'./denominator);                                      %Alg1-Step 6.2
        
        %toc
        
        alpha_0 = sum(alpha_t,2) - alpha_t(:,N+1);                                          %Alg1-Step 7.1
        tau_0 = max(1e-8,sum(tau_t,2) - tau_t(:,N+1));                                      %Alg1-Step 7.2
        mu_x  = (alpha_0)./tau_0;                                                           %Alg1-Step 8.0.1
        var_x = 1./tau_0;                                                                   %Alg1-Step 8.0.2
        
        tmp_mu  = (mu_x.*tau_x0+mu_x0.*var_x)./(var_x+tau_x0);                              %Alg1-Step 8.1, calculate mean
        tmp_var = (tau_x0.*var_x)./(var_x+tau_x0);                                          %Alg1-Step 8.2, calculate variance
        mean_x  = tmp_mu./(1+...
            (1-b)/b*tau_x0./tmp_var.*exp(-abs(tmp_mu).^2./tmp_var+abs(mu_x0).^2./tau_x0));
        var_x   = (abs(tmp_mu).^2+tmp_var)./(1+...
            (1-b)/b*tau_x0./tmp_var.*exp(-abs(tmp_mu).^2./tmp_var+abs(mu_x0).^2./tau_x0));
        var_x   = max(1e-8,var_x - abs(mean_x).^2);
        
        tmp_mse(iter:end) = sum(abs(mean_x-x).^2)/sum(Ex);                                  %mse per x realization
        nse = sum(abs(mean_x-x).^2)/sum(Ex);
        if ~isdeployed
            %if iter > 100 && max(abs(est_x-mean_x)) < 1e-5
            if  max(abs(est_x-mean_x)) < 1e-5                                               %Alg1-Step 9, if a station point is found
                break
            end
            est_x = mean_x;
        end
        
    end
    mse        = mse + tmp_mse;
    mse_final = mse_final + sum(abs(mean_x-x).^2)/sum(Ex);
    convergence_analysis(n_frame)=iter;
    if ~isdeployed
        if mod(n_frame,100)==0
            fprintf('EP:rho=%.2f, Nr_frames %i, MSE: %.2f %.2f %.2f, Iteration Num %.2f\n',b, n_frame,10*log10(mse(end)/n_frame),10*log10(mse_final/n_frame),10*log10(nse), mean(convergence_analysis(1:n_frame)));
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
