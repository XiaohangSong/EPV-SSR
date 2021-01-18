# EPV-SSR
This is the companion code for the message passing algorithms reported in the paper "Unifying Message Passing Algorithms Under the Framework of Constrained Bethe Free Energy Minimization" authored by Dan Zhang, Xiaohang Song, Wenjin Wang, Gerhard Fettweis, Xiqi Gao. For more details, one can refer to the paper at https://arxiv.org/abs/1703.10932. This code is distributed for the readers to reproduce and extend the results. When reporting, reproducing or extending the results, please cite the paper above.

# Purpose of the project
This code is a developed for and published as an application example of the framework proposed in [https://arxiv.org/abs/1703.10932]. It will neither be maintained nor monitored in any way.

# License
EPV-SSR is open-sourced under the AGPL-3.0 license. See the LICENSE file for details.

# Getting Start
One can run the Sample_script_main.m file to see how to obtain some quick results.

To fully reproduce the results in the paper, we suggest the users to use batch processing. For results shown in the paper, the simulation is carried out at a compute farm with 10^5 frames per sparsity ratio per algorithm until running out of the permitted resources. In addition, the number of iterations per realizaiton is set to 500, while the early termination is disabled. 
