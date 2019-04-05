using SnFFT

# Choose the size of symmetric group to compute eigenvalues for
n = 5

# If print_perm is equal to 1, then all permutations and their FFT
# are printed out (if n > 6, do not do this).
print_perm = 0

# If print_embed is equal to 1, then he MDS coordinates of each
# permutation are printed out (if n > 6, do not do this).
print_embed = 0

# Generate permutations
perms = permutations(1:n) |> collect

# Set array to hold matrix coefficients for FFT
record_transform = zeros(factorial(n),factorial(n))

# Loop through all permutations of S_n, computing FFT of each
for i = 1:factorial(n)

    # Set proper coefficient equal to 1
    perm_coef = zeros(factorial(n))
    perm_coef[i] = 1

    # Generate group algebra element
    x = snf(n,perms,perm_coef)

    # Get information related to Young's orthogonal form
    RA, PT = yor(n)

    # Get Fourier transform
    x_trans = sn_fft(n,x,RA,PT)

    # If print_perm is equal to 1, print each permutation and its
    # FFT
    if (print_perm == 1)
        println("                  ")
        println("                  ")
        println("The decomposition for permutation ")
        println(perms[i])
        println("is")
        println("******************")
        println(x_trans)
    end

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

# Create matrix to store Kendall tau distances
distance_matrix = zeros(factorial(n),factorial(n))

# Compute the pairwise Kendall tau distance between all points
# and store these in the distance matrix
for i = 1:factorial(n)
    for j = 1:factorial(n)
        p1 = preferencematrix(perms[i])
        p2 = preferencematrix(perms[j])
        distance_matrix[i,j] = kendalldistance(p1,p2)
    end
end

# Perform MDS on the distance matrix that has been found
D = map(x->.5*x^2,distance_matrix)
H = eye(factorial(n),factorial(n)) - (1/factorial(n))*ones(factorial(n),factorial(n))
M = H*D*H

U,S,V = svd(M)
U = map(x->round(x,3),U)
S = map(x->round(x,6),S)

# Print eigenvalues for MDS operator
println("The full list of MDS eigenvalues is")
println(S)

# Print MDS coordinates for permutations
if print_embed == 1
    for i = 1:factorial(n)
        println("The coordinates of")
        println(perms[i])
        println("is")
        println(U[i,:])
        println("----------------")
    end
end

# Solve the linear system to determine which eigenvalue corresponds
# to each representation
for i = 1:factorial(n)
    # Only print this data when the eigenvalue is not zero
    if S[i] != 0.0
        x = record_transform\U[:,i]
        println("The eigenvalue")
        println(S[i])
        println("corresponds to coefficients:")
        println(round(x,3))
        println("------------------------------")
    end
end
