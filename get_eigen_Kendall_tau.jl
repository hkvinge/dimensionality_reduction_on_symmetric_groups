using SnFFT

n = 5

# Generate permutations
perms = permutations(1:n) |> collect
record_transform = zeros(factorial(n),factorial(n))

for i = 1:factorial(n)
    
    # Print permutation
    println("                  ")
    println("                  ")
    println("The decomposition for permutation ")
    println(perms[i])
    println("is")
    println("******************")

    # Set proper coefficient equal to 1
    perm_coef = zeros(factorial(n))
    perm_coef[i] = 1

    # Generate group algebra element
    x = snf(n,perms,perm_coef)

    # Get information related to Young's orthogonal form
    RA, PT = yor(n)

    # Get Fourier transform
    x_trans = sn_fft(n,x,RA,PT)

    println(x_trans)
 
    # Record the transform in one matrix. This will be used to 
    # find which eigenvalues of the MDS operator correspond to
    # which representations   
    count = 1
    transform_size = size(x_trans)
    # Iterate through different terms in product
    for j = 1:transform_size[1]
        # Identify size of representation and record
        # its entries.
        rep_size = size(x_trans[j])
        rep = x_trans[j]
        for k = 1:rep_size[1]
            for l = 1:rep_size[1]
                record_transform[i,count] = rep[k,l]
                count = count + 1
        
            end
        end
    end
end

# Get coordinates of identity element to center embedding


distance_matrix = zeros(factorial(n),factorial(n))

# Compute the pairwise Kendall tau distance between all points
for i = 1:factorial(n)
    println(i)
    for j = 1:factorial(n)
        p1 = preferencematrix(perms[i])
        p2 = preferencematrix(perms[j])
        distance_matrix[i,j] = kendalldistance(p1,p2)
    end
end

# Perform MDS
D = map(x->.5*x^2,distance_matrix)
H = eye(factorial(n),factorial(n)) - (1/factorial(n))*ones(factorial(n),factorial(n))
M = H*D*H

U,S,V = svd(M)
U = map(x->round(x,3),U)
S = map(x->round(x,6),S)
println(S)
for i = 1:factorial(n)
    println("The coordinates of")
    println(perms[i])
    println("is")
    println(U[i,:])
    println("----------------")
end

# Solving linear systems to determine which eigenvalue corresponds
# to each representation

for i = 1:factorial(n)
    x = record_transform\U[:,i]
    println("This eigenvalue corresponds to")
    println(i)
    println(S[i])
    println(round(x,3))
    println("------------------------------")
end
