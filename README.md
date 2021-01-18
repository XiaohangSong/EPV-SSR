# EPV-SSR
This is the companion code for the message passing algorithms reported in the paper "Unifying Message Passing Algorithms Under the Framework of Constrained Bethe Free Energy Minimization" authored by Dan Zhang, Xiaohang Song, Wenjin Wang, Gerhard Fettweis, Xiqi Gao. For more details, one can refer to the paper at https://arxiv.org/abs/1703.10932. This code is distributed for the readers to reproduce and extend the results. When reporting, reproducing or extending the results, please cite the paper above.

# Purpose of the project
This code is developed for and published as an application example of the framework proposed in [https://arxiv.org/abs/1703.10932]. It will neither be maintained nor monitored in any way.

# License
EPV-SSR is open-sourced under the AGPL-3.0 license. See the LICENSE file for details.

# Getting Start
One can run the Sample_script_main.m file to obtain some quick results. For the results shown in the paper, the simulation is carried out at a compute farm with 10^5 frames per sparsity ratio per algorithm until running out of the permitted resources, and the early termination should be disabled for convergence analysis. 

Here, algorithms one to four are provided for a linear sparse signal recovery (SSR) problem, namely the expectation propagation (EP) algorithm, the expectation propagation variant (EPV)  algorithm, the approximate message passing (AMP) algorithm, and the hybrid 'expectation maximization'-'variational message passing'-'expectation propagation variant' (EM-VMP-EPV)  algorithm. The first three algorithms require prior information, while the last one can deal with SSR problems without a perfect statistical model.
