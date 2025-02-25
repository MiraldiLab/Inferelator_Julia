
# Method 1: Matrix Inversion Method: 
#           this methods computes partial correlation from precision-matrix (the inverse of variance-covariance matrix)
using LinearAlgebra 
function partialCorrelationMat(X::Matrix{Float64}; epsilon = 1e-7,  first_vs_all = false)
    # Mean centering of the columns (mean subtraction)
    X_centered = X .- mean(X, dims=1)

    # Compute the covariance matrix
    sigma = cov(X_centered)
    # regularizing the covariance matrix to avoid ill-conditining and singularity
    sigma = sigma + epsilon * I

    # Precision matrix (inverse of the covariance matrix)
    theta = inv(sigma)
    p = size(X, 2)  # number of features/predictors/explanatory variables

    if first_vs_all
        # Compute partial correlation for the first variable vs all others
        P = ones(1, p) # Initialize the partial correlation matrix
        for j in 2:p
            P[j] = -theta[1, j] / sqrt(theta[1, 1] * theta[j, j])
        end 
    else
        # Compute the full partial correlation matrix
        P = ones(p, p) # Initialize the partial correlation matrix
        for i in 1:(p-1)
            for j in (i+1):p
                P[i , j] = -theta[i, j] / sqrt(theta[i, i] * theta[j, j])
                P[j, i] = P[i, j]  # Symmetry
            end
        end
    end

    return P
end


#= # Method 2: Using Regression
using GLM, Statistics, DataFrames

function partialCorrReg(X::Matrix{Float64}; first_vs_all = false)
    p = size(X, 2)  
    P = ones(p, p)  

    # Convert the matrix to a DataFrame for easier variable handling
    df = DataFrame(X, :auto);
    col_names = names(df)
    
    if first_vs_all
        P = ones(1, p)  
        for j in 2:p 
            keep_indices = setdiff(2:p, j)
            covariates =  col_names[keep_indices] # All columns except i and j
            
            # First, regress the first covariate on all other features except 
            model_1 = lm(term.(col_names[1]) ~ sum(term.(covariates)), df)
            res_1 = residuals(model_1)  # Get residuals for 1

            # Then, regress j on all others except i
            model_j = lm(term.(col_names[j]) ~ sum(term.(covariates)), df)
            res_j = residuals(model_j)  # Get residuals for j

            # Compute c1rrelation between the residuals
            r = cor(res_1, res_j)
            P[j] = r
        end
    else
        P = ones(p, p)  
        for i in 1:p
            for j in (i+1):p 

                keep_indices = setdiff(1:p, [i ,j])
                covariates =  col_names[keep_indices] # All columns except i and j

                # reference: https://discourse.julialang.org/t/using-all-independent-variables-with-formula-in-a-multiple-linear-model/43691/4

                # First, regress i on all other features except j (i is the response)
                model_i = lm(term.(col_names[i]) ~ sum(term.(covariates)), df)
                res_i = residuals(model_i)  # Get residuals for i

                # Then, regress j on all others except i
                model_j = lm(term.(col_names[j]) ~ sum(term.(covariates)), df)
                res_j = residuals(model_j)  # Get residuals for j

                # Compute correlation between the residuals
                r = cor(res_i, res_j)
                P[i, j] = r
                P[j, i] = r  # Symmetric matrix
            end
        end
    end
    return P
end


 =#