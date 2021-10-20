mutable struct Configuration
   
    j::Float64       
    M::String
    N::Int64 
    b::Int64 
    σ::Float64 
    j_half_int::Half{Int64}
    D::Int64
    d::Int64
    
    random_walk::Bool 
    add_chains::Bool   
    draws_stored_for_this_chain::Bool
    total_draws_already_stored::Int64
  
    compute_angles::Bool   
    compute_angles_correlations::Bool
    angles_already_stored_for_this_chain::Bool 
    angles_pseudo_correlations_already_stored_for_this_chain::Bool 
    total_angles_already_stored::Int64
    total_angles_pseudo_correlations_already_stored::Int64
    
    compute_volumes::Bool
    compute_volumes_correlations::Bool
    volumes_correlations_node_1::Int64
    volumes_correlations_node_2::Int64
    volumes_already_stored_for_this_chain::Bool 
    volumes_pseudo_correlations_already_stored_for_this_chain::Bool
    total_volumes_already_stored::Int64 
    total_volumes_pseudo_correlations_already_stored::Int64
    
    compute_entropy::Bool
    subsystem::Vector{Int64}
    density_matrix_already_stored_for_this_chain::Bool
    total_density_matrices_already_stored::Int64
    
    draws_folder::String
    ampls_folder::String 
    operators_folder::String
    tables_folder::String
    
    angles_folder::String  
    angles_sq_folder::String
    angles_spread_folder::String
    angles_pseudo_correlations_folder::String 
    angles_correlations_folder::String 
    
    volumes_folder::String 
    volumes_sq_folder::String
    volumes_spread_folder::String
    volumes_pseudo_correlations_folder::String
    volumes_correlations_folder::String
    
    density_matrices_folder::String

end


function init_config(user_conf::Vector{Any}, data_folder_path::String, draws_stored_for_this_chain::Bool=false, total_draws_already_stored::Int64=0, angles_already_stored_for_this_chain::Bool=false, total_angles_already_stored::Int64=0, angles_pseudo_correlations_already_stored_for_this_chain::Bool=false, total_angles_pseudo_correlations_already_stored::Int64=0, volumes_already_stored_for_this_chain=false, total_volumes_already_stored::Int64=0, volumes_pseudo_correlations_already_stored_for_this_chain=false, total_volumes_pseudo_correlations_already_stored::Int64=0, density_matrix_already_stored_for_this_chain=false, total_density_matrices_already_stored=0)

    j_half_int = half(2*user_conf[1])
    j = user_conf[1]
    M = user_conf[2]
    N = user_conf[3]
    b = user_conf[4]
    σ = user_conf[5]    
    
    D = Dimension(user_conf[1])
    d = dimension(user_conf[1])

    volumes_correlations_node_1 = user_conf[8][3]
    volumes_correlations_node_2 = user_conf[8][4]
    
    nodes_in_subsystem = user_conf[9][2]

    intermediate_path_data_folder = "/data_star_model/$(M)/j=$(j)/N=$(N)_b=$(b)_σ=$(σ)"
    draws_folder = data_folder_path*intermediate_path_data_folder*"/draws"
    ampls_folder = data_folder_path*intermediate_path_data_folder*"/ampls"
    operators_folder = data_folder_path*intermediate_path_data_folder*"/operators"
    tables_folder = data_folder_path*intermediate_path_data_folder*"/tables"
    angles_folder = operators_folder*"/angles" 
    angles_sq_folder = operators_folder*"/angles_sq"
    angles_spread_folder = operators_folder*"/angles_spread"
    angles_pseudo_correlations_folder = operators_folder*"/angles_pseudo_correlations"
    angles_correlations_folder = operators_folder*"/angles_correlations"
    volumes_folder = operators_folder*"/volumes" 
    volumes_sq_folder = operators_folder*"/volumes_sq"
    volumes_spread_folder = operators_folder*"/volumes_spread"
    volumes_pseudo_correlations_folder = operators_folder*"/volumes_pseudo_correlations/node_1=$(volumes_correlations_node_1)_node_2=$(volumes_correlations_node_2)"
    volumes_correlations_folder = operators_folder*"/volumes_correlations"
    density_matrices_folder = operators_folder*"/density_matrices/subsystem_$(nodes_in_subsystem)"
  
    conf = Configuration(user_conf[1], user_conf[2], user_conf[3], user_conf[4], user_conf[5], j_half_int, D, d, user_conf[6][1], user_conf[6][2], draws_stored_for_this_chain, total_draws_already_stored, user_conf[7][1], user_conf[7][2], angles_already_stored_for_this_chain, angles_pseudo_correlations_already_stored_for_this_chain, total_angles_already_stored, total_angles_pseudo_correlations_already_stored, user_conf[8][1], user_conf[8][2], user_conf[8][3], user_conf[8][4], volumes_already_stored_for_this_chain, volumes_pseudo_correlations_already_stored_for_this_chain, total_volumes_already_stored, total_volumes_pseudo_correlations_already_stored, user_conf[9][1], user_conf[9][2], density_matrix_already_stored_for_this_chain, total_density_matrices_already_stored, draws_folder, ampls_folder, operators_folder, tables_folder, angles_folder, angles_sq_folder, angles_spread_folder, angles_pseudo_correlations_folder, angles_correlations_folder, volumes_folder, volumes_sq_folder, volumes_spread_folder, volumes_pseudo_correlations_folder, volumes_correlations_folder, density_matrices_folder)
    
    return conf
    
end 

function make_folders(conf::Configuration)

  mkpath(conf.draws_folder)
  mkpath(conf.ampls_folder) 
  mkpath(conf.operators_folder) 
  mkpath(conf.tables_folder)
   
  mkpath(conf.angles_folder)  
  mkpath(conf.angles_sq_folder)    
  mkpath(conf.angles_spread_folder)
   
  mkpath(conf.angles_pseudo_correlations_folder)
  mkpath(conf.angles_correlations_folder)
   
  mkpath(conf.volumes_folder)
  mkpath(conf.volumes_sq_folder) 
  mkpath(conf.volumes_spread_folder)    
 
  mkpath(conf.volumes_pseudo_correlations_folder)
  mkpath(conf.volumes_correlations_folder)
  
  mkpath(conf.density_matrices_folder)
    
end
