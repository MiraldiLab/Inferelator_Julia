function updatePenaltyMatrix(penalty_matrix, targGenes, allPredictors, priorFile_penalties, lambdaBias)
    # Create a dictionary to store the minimum lambda for each interaction
    min_lambda_dict = Dict{Tuple{String,String}, Float64}()
    
    # Iterate through each prior file and its associated lambda
    for (file_path, lambda) in zip(priorFile_penalties, lambdaBias)
        # Read the prior file
        prior_data = readdlm(file_path)
        
        # Extract gene names and TF names from the prior file
        prior_genes = prior_data[2:end, 1]
        # Get TF names, filtering out any empty or missing entries
        prior_tfs = filter(tf -> !ismissing(tf) && !isempty(string(tf)), prior_data[1, :])
        
        # Create indices mapping for faster lookup
        gene_to_idx = Dict(gene => i for (i, gene) in enumerate(targGenes))
        tf_to_idx = Dict(tf => i for (i, tf) in enumerate(allPredictors))
        
        # Process the interactions
        for (i, gene) in enumerate(prior_genes)
            for (j, tf) in enumerate(prior_tfs)
                if prior_data[i+1, j+1] != 0 && haskey(gene_to_idx, gene) && haskey(tf_to_idx, tf)
                    interaction = (gene, tf)
                    if !haskey(min_lambda_dict, interaction) || lambda < min_lambda_dict[interaction]
                        min_lambda_dict[interaction] = lambda
                    end
                end
            end
        end
    end
    
    # Apply the penalties using the minimum lambda values
    for ((gene, tf), min_lambda) in min_lambda_dict
        gene_idx = findfirst(==(gene), targGenes)
        tf_idx = findfirst(==(tf), allPredictors)
        penalty_matrix[gene_idx, tf_idx] = min_lambda
    end
    
    return penalty_matrix
end
