%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sample script for message passing algorithms in the paper 
% "Unifying Message Passing AlgorithmsUnder the Framework of Constrained Bethe Free Energy Minimization"
% Authors: Dan Zhang, Xiaohang Song, Wenjin Wang, Gerhard Fettweis, Xiqi Gao
%  https://arxiv.org/abs/1703.10932
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Copyright (C) 2019~2021 Xiaohang Song
%  Vodafone Chair, Technische Universitaet Dresden.
%  Last updated on Jan. 17, 2021
% 
%  This code is distributed for the users to reproduce and extend the results. 
%  When reporting, reproducing or extending the results, please cite the paper above.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
sim_frames=1e2; % For results shown in the paper, the simulation is carried out at a compute farm with 10^5 frames per sparsity ratio
                % per algorithm until the reserved computing time is reached. In addition, the number of iterations per realizaiton is 
                % set to 500, while the early termination is disabled. 
% Note: With current implementation, the output data is only created after
% 100 realizations, i.e., sim_frames >=100.
rho_range=0.1:0.05:0.4;  %Sparsity
MSE=zeros(size(rho_range));
Iteration_Nr=zeros(size(rho_range));
algorithms_to_be_test=[1 1.1 2 3 4];
for algorithm_under_test=algorithms_to_be_test
    counter=1;
    for rho=rho_range
        if algorithm_under_test==1
            fprintf('EP started\n');
            Data=compressed_sensing_EP('500',num2str(rho),'30','1',num2str(sim_frames));
            color_opt='-ks';
        elseif algorithm_under_test==1.1
            fprintf('EP with damping started\n');
            Data=compressed_sensing_EP('500',num2str(rho),'30','0.9',num2str(sim_frames));
            color_opt='-g*';             
        elseif algorithm_under_test==2
            fprintf('EPV started\n');
            Data=compressed_sensing_EP_Variant('500',num2str(rho),'30','1',num2str(sim_frames));
            color_opt='-ro';
        elseif algorithm_under_test==3
            fprintf('AMP started\n');
            Data=compressed_sensing_amp('500',num2str(rho),'30','1',num2str(sim_frames));
            color_opt='-bx';
        elseif algorithm_under_test==4
            fprintf('Hybrid EM-VMP-EPVariant started\n'); 
            Data=compressed_sensing_EM_VMP_EPV_modeling('500',num2str(rho),'30','1',num2str(sim_frames)); 
            color_opt='-bo';
        end
        
        MSE(counter)=Data.mse_final;
        Iteration_Nr(counter)=Data.iter;
        counter=counter+1;
    end
    
    figure(1);
    plot(rho_range,10*log10(MSE),color_opt);
    xlim([0.1 0.4]);
    ylim([-40 0]);
    grid on;
    hold on;
    xlabel('Sparsity ratio \rho');
    ylabel('NMSE');
    
   % figure(2);
   % plot(rho_range,Iteration_Nr,color_opt);
   % xlim([0.1 0.4]);
   % grid on;
   % hold on;
   % xlabel('Sparsity ratio \rho');
   % ylabel('The Mean of the Number of Iterations until Convergence');
end
