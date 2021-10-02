
# Deprecated slower versions
#=
function star_amplitude_old_1(D::Int64, A::Array{Float64, 5}, indices::Vector{Int64})
    
      amp = 0.0   
    
      for i=1:1: D
      I = A[i,indices[1],indices[2],indices[3],indices[4]]   
      if (I == 0) continue   end    
        
      for j=1:1: D
      J = A[j,indices[5],indices[6],indices[7],indices[8]]  
      if (J == 0) continue   end     
            
      for k=1:1: D
      K = A[k,indices[9],indices[10],indices[11],indices[12]]  
      if (K == 0) continue   end      
                
      for l=1:1: D
      L = A[l,indices[13],indices[14],indices[15],indices[16]]
      if (L == 0) continue   end   
                
      for m=1:1: D
      M = A[m,indices[17],indices[18],indices[19],indices[20]]
      if (M == 0) continue   end    
                        
      bulk = A[m,l,k,j,i]                   
                                                
      @inbounds amp += I*J*K*L*M*bulk                    
        
                  end # m cycle
               end # l cycle
            end # k cycle
         end # j cycle
      end  # i cycle   
    
      return amp
      
end  


# Deprecated
function star_amplitude_old_2(A::Array{Float64, 5}, indices::Vector{Int64})

  @tullio amp := A[q,indices[1],indices[2],indices[3],indices[4]]*A[l,indices[5],indices[6],indices[7],indices[8]]*A[k,indices[9],indices[10],indices[11],indices[12]]*
                 A[j,indices[13],indices[14],indices[15],indices[16]]*A[m,indices[17],indices[18],indices[19],indices[20]]*A[m,j,k,l,q]   
    
  return amp        
    
end 
=#



# computes star amplitude for given dimensions, vertex amplitude, containers for contraction and indices
function star_amplitude(D,A,v1,v2,v3,v4,v5,indices)  

  @turbo for i in 1:D    
        
  v1[i] = A[i,indices[1],indices[2],indices[3],indices[4]]       
  v2[i] = A[i,indices[5],indices[6],indices[7],indices[8]]     
  v3[i] = A[i,indices[9],indices[10],indices[11],indices[12]]
  v4[i] = A[i,indices[13],indices[14],indices[15],indices[16]]
  v5[i] = A[i,indices[17],indices[18],indices[19],indices[20]]        
        
  end   
    
  s = 0.0
  
  @turbo for q in 1:D, l in 1:D, k in 1:D, j in 1:D, m in 1:D
         s += A[m,j,k,l,q]*v1[q]*v2[l]*v3[k]*v4[j]*v5[m]      
         end
  return s 
    
end




