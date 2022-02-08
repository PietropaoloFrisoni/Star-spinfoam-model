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




