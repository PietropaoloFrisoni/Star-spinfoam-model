function compute_density_matrix_function!(j, D, vertex, draws, number_of_draws, ampls, density_matrix, N, b, density_matrices_folder, chain_id, subsystem, number_of_nodes_in_subsystem, density_matrix_linear_dim, total_density_matrices_already_stored=0)

  # containers for contraction
  v1 = zeros(Float64, D)
  v2 = zeros(Float64, D)
  v3 = zeros(Float64, D)
  v4 = zeros(Float64, D)
  v5 = zeros(Float64, D)
 
  indices = zeros(Int64,20) 
  vector_unflatted_indices = zeros(Int, number_of_nodes_in_subsystem)
  comparison_vector = zeros(Int, number_of_nodes_in_subsystem)
 
  for flat_index_row = 1:density_matrix_linear_dim
  
      for n = 1:number_of_draws 
      
          from_flat_index_to_vector_indices!(D, vector_unflatted_indices, number_of_nodes_in_subsystem, flat_index_row)
          
          for k=1:number_of_nodes_in_subsystem
          @inbounds comparison_vector[k] = abs(vector_unflatted_indices[k] - draws[subsystem[k], n])
          end
          
          if maximum(comparison_vector) > 0
          continue
          end
          
          for i=1:20
          @inbounds indices[i] = draws[i, n]
          end

              for flat_index_column = flat_index_row:density_matrix_linear_dim

                  from_flat_index_to_vector_indices!(D, vector_unflatted_indices, number_of_nodes_in_subsystem, flat_index_column)

                  for k=1:number_of_nodes_in_subsystem
                  @inbounds indices[subsystem[k]] = vector_unflatted_indices[k] 
                  end
              
                    if (flat_index_column == flat_index_row) 
                    @inbounds density_matrix[flat_index_row, flat_index_column] += 1 
                    continue
                    end   
              
                  @inbounds density_matrix[flat_index_row, flat_index_column] += star_amplitude(D, vertex, v1, v2, v3, v4, v5, indices)/ampls[n]
              
              end
               
      end # draw cycle
      
  end # flat_index_row cycle
 
  density_matrix = (density_matrix + transpose(density_matrix) - Diagonal(density_matrix))/sum(Diagonal(density_matrix))
  
  @save "$(density_matrices_folder)/density_matrix_chain=$(chain_id + total_density_matrices_already_stored).jld2" density_matrix 

end





function entropy_assemble(conf::Configuration, chains_to_assemble::Int64, density_matrix_linear_dim)

  density_matrix_all_chains = zeros(Float64, density_matrix_linear_dim, density_matrix_linear_dim, chains_to_assemble) 

  for id_chain=1:chains_to_assemble
    
    @load "$(conf.density_matrices_folder)/density_matrix_chain=$(chain_id).jld2" density_matrix 
      
    density_matrix_all_chains[:, :, id_chain] = density_matrix[:, :]

  end # cycle in id_chain
      
  density_matrix_all_chains = sum(density_matrix_all_chains, dims = 3)
  density_matrix_all_chains[:] ./= chains_to_assemble
  
  # trick to make disappear the third fictitious dimension (N×N×1 Array{Float64, 3} ---> N×N Matrix{Float64})
  density_matrix_all_chains = density_matrix_all_chains[:,:]
  
  for i=1:density_matrix_linear_dim, j=1:density_matrix_linear_dim  
  @inbounds density_matrix_all_chains[j,i] = round(density_matrix_all_chains[j,i], digits = 5)
  end      
      
  decomp = eigen(Symmetric(density_matrix_all_chains))
    
  density_matrix_eigenvalues = decomp.values
    
  entropy = 0.0
  
  for i=1:density_matrix_linear_dim
  entropy -= density_matrix_eigenvalues[i]*log(density_matrix_eigenvalues[i])
  end 

  entropy_dataframe = DataFrame(to_rename = entropy)
  column_name = "j=$(conf.j)"
  rename!(entropy_dataframe, :to_rename => column_name) # julia is weird 
 
  entropy_table_name = "/entropy_$(chains_to_assemble)_chains_combined.csv"
      
  entropy_table_full_path = conf.tables_folder*entropy_table_name

  CSV.write(entropy_table_full_path, entropy_dataframe)
  
  return entropy_dataframe
  
end





