# Computation of star spin foam amplitude

**This code can be used simply interacting with the file "configs_to_compute.jl", specifying the configurations that you want to compute and the M-H parameters for the computation.**


The usage, as well as the meaning of the various flags, is explained in this document (see sections below). 

***The code builds a Markov chain for each julia process and the latter are parallelized on the available cores.***

After choosing the configurations to compute, in order to execute the code (on a single machine with the synthax below) just run the following command:

```
julia -p [N-1] execute_me.jl
```

where [N-1] is the number of Markov chains that you want to assemble, and the julia executable must be linked. If you want to use a single chain just omit the "-p [N-1]".

**The code can run on multiple machines.** In this case the above synthax must be replaced, for example (on CC clusters) with:

```
srun hostname -s > hostfile
sleep 5
julia --machine-file ./hostfile ./execute_me.jl
```

*It is advisable for the performance to use a number of Markov chains equal to or less than the physical number of cores present on the system.*

At the end of the computation of the draws and operators requested by the user, the code assembles together the data computed by every chain involved and writes it in csv format in the "final_data" folder, where the various computations are separated based on the date and time in which they ended.

A full list of the packages used in this code can be found in the file "/inc/pkgs.jl". *Before executing the code, all packages must be installed.*

**Finally, notice that the julia's Just-in-Time compiler is such that the first execution of the random walk and operators' computation is considerably slower that following ones, and it also allocates much more memory**. To avoid this, you can use the [DaemonMode package](https://github.com/dmolina/DaemonMode.jl).

![Alt text](https://github.com/PietropaoloFrisoni/Star-spinfoam-model/blob/master/RW_benchmark.pdf?raw=true "Random walk benchmark")

## Usage:

```
[j, model to compute, number of iterations, number of burnin iterations, sigma of gaussian, *[random walk flags]*, *[angle operators flags]*, *[volume operators flags]*, *[entropy operator flags]*]

```



## Meaning of parameters:


- "j" ---> Common spin of the boundary links of the star spin foam. It must be a number between 0 and 10 (10 is the highest number I'm interested in, so the vertex ampls were not pre-computed for higher values). For halfintegers values, please use the float notation, for example "2.5" etc.
  
- "model to compute" ---> It must be set to "BF" or "EPRL"

- "number of iterations" ---> Number of M-H iterations *for each chain*

- "number of burnin iterations" ---> Number of M-H burnin iterations *for each chain*

- "sigma of gaussian" ---> Standard deviation of M-H Gaussian proposal density


*[random walk flags]* = [random_walk, add_chains]

- "random_walk" --->  If "true" perform the random walk for such configuration, if "false" it doesn't (in this case the draws must have been been previously computed). The code automatically searches for previously computed operators in every chain. *The stored values are used if "random_walk" is false, otherwise they are computed again*

- "add_chains"  --->  If this parameter is true, then *all new computed draws and operators will be stored in addition to those previously found in "data_folder_path", and nothing is overwritten*
      
                                                       
*[angle operator flags]* = [compute_angles, compute_angles_correlations]

- "compute_angles"              --->  If "true" computes the angles for such configuration

- "compute_angles_correlations" --->  If "true" computes all the angles correlations C(angle_i,angle_j) for such configuration


*[volume operator flags]* = [compute_volumes, compute_volumes_correlations, node_1, node_2]

- "computes_volumes"             ---> If "true" computes the volumes for such configuration. With respect to angles, the code computes only the volumes in the node 1
                                      
- "compute_volumes_correlations" ---> If "true" computes volume correlations C(vol_node_1, vol_node_2) for such configuration
                                      
- "node_1"                       ---> It must be an integer between 1 and 20. It refers to the first node for the computation of volumes correlations

- "node_2"                       ---> It must be an integer between 1 and 20. It refers to the second node for the computation of volumes correlations


*[entropy operator flags]* = [compute_entropy, [nodes_in_subsystem]]

- "compute_entropy"              ---> If "true" computes the entropy of the specified subsystem for such configuration
                                      
- "[nodes_in_subsystem]"         ---> It must be a vector of integers between 1 and 20, with maximum length equal to 20. It specifies the nodes in the subsystem
                                      


## Example of usage:

See "configs_to_compute".



#### Current limitations:

- The contraction of vertex amplitudes on the GPU is not currently implemented. This would improve significantly the performance for large spins

- The algorithm is written in such a way as to perform the random walk and compute the observables in the same run, using the same number of resources. The two phases could be totally separated, in such a way as to parallelize the computation of observables over an arbitrary number of tasks (which already happens) by also exploiting an arbitrary number of CPUs for each task, that is, for each Markov chain. So far, this has only been done for the density matrix computation, where each chain can use an arbitrary number of threads, making it possible to compute such a matrix for a subsystem with many nodes

- If "add_chain" is true, even if in some cases it is possible to assemble more chains than those chosen by the user (this depends on the operators previously stored and on those that the user wants to compute for each configuration), the code only assembles a number of chains corresponding to that chosen by the user



