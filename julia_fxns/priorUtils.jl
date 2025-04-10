using DataFrames
using CSV
using LinearAlgebra
using TickTock


# COnvert Wide/long data to long/wide
function convertToLong(data)
  #dfs = [convertToLong(df) for df in dfs]
    if ncol(data) > 3
        return stack(data, Not(1) )
    else
        return data
    end
end


function convertToWide(Data; indices::Union{Nothing, NTuple{3, Int}}=nothing)
    """
    Converts a 3‑column long-format DataFrame to wide‑format using the unstack function.
    If the input DataFrame has exactly 3 columns, the conversion is performed.
    
    - Data columns  is a wide matrix with columns as TF and rows as target genes or 
            a long data with columns in the order TF, Gene, Weights.
    indices::Union{Nothing, NTuple{3, Int}} = nothing:
            A tuple specifying the column indices in the order (pivot, key, value).
                - pivot: The column that provides the row identifier.
                - key: The column whose unique values will become new column names.
                - value: The column from which the cell values are taken.

    If no indices are provided, the function defaults to (1, 2, 3).
    If the DataFrame has more than 3 columns, the original DataFrame is returned.
    If the DataFrame has less than 3 columns, an error is thrown.

    # USAGE
    dfs = [convertToWide(df) for df in dfs]
    dfs = [convertToWide(df; indices = (1,2,3)) for df in dfs]
    """

    ncols = ncol(Data)

    if ncols < 3
        error("DataFrame has less than 3 columns. A 3‑column DataFrame or a wide-matrix is required.")
    elseif ncols > 3
        # More than 3 columns: return data unchanged.
        return Data
    else # Exactly 3 columns
        # Use provided indices, or default to (2, 1, 3)
        inds = isnothing(indices) ? (2,1,3) : indices
        # Convert the specified columns to Symbols for unstack.
        idSym    = Symbol(names(Data)[inds[1]])
        keySym   = Symbol(names(Data)[inds[2]])
        valueSym = Symbol(names(Data)[inds[3]])

        return unstack(Data, idSym, keySym, valueSym)
    end
end
    


function frobeniusNormalize1(df::DataFrame, dims::Symbol =:column)

    """
    frobeniusNormalize(df::DataFrame; dims::Symbol = :row)

    Normalize a DataFrame based on the Frobenius (L2) norm.

    • If dims = :row (default), it normalizes each row (ignoring the first column).
    • If dims = :column, it normalizes each column (ignoring the first column).

    The first column is assumed to contain non‐numeric data (e.g., row identifiers)
    and is left unchanged.
    """
    # Create a copy so that the original DataFrame remains unchanged.
    dfNorm = deepcopy(df)
    dfMat = Matrix(dfNorm[!, 2:end])  # Extract numerical part
    dfMat = convert(Matrix{Float64}, dfMat)  # Ensure it's Float64

    if dims == :row
        # Normalize each row.
        for i in 1:size(dfMat, 1)
            nrm = norm(dfMat[i, :], 2)  # Compute the L2 norm for the row.
            if nrm ≠ 0
                dfMat[i, :] ./= nrm
            end
        end
    elseif dims == :column
        # Normalize each column.
        for j in 1:size(dfMat, 2)
            nrm = norm(dfMat[:, j], 2)  # Compute the L2 norm for the column.
            if nrm ≠ 0
                dfMat[:, j] ./= nrm
            end
        end
    else
        throw(ArgumentError("dims must be :row or :column"))
    end
    
    # Update the DataFrame with the normalized numeric values.
    dfNorm[!, 2:end] .= dfMat
    return dfNorm
end



function completeDF(df::DataFrame, id::Symbol, idsALL, allCols)
    """
    completeDF(df::DataFrame, id::Symbol, idsALL, allCols)

    Aligns a DataFrame to a common set of row identifiers and columns.

    ### Arguments:
    - `df::DataFrame`: The input DataFrame to align.
    - `id::Symbol`: The identifier column (e.g., gene or sample ID).
    - `idsALL`: A vector of all unique row identifiers across multiple DataFrames.
    - `allCols`: A vector of all unique column names across multiple DataFrames.

    ### Returns:
    - A new DataFrame with:
    - Rows corresponding to `idsALL` (missing rows filled with `missing`).
    - Columns matching `allCols` (missing columns filled with `missing`).
    - The original values retained where available.

    """
    dfNew = DataFrame()
    dfNew[!, id] = idsALL  # Ensure all IDs are included

    for col in allCols
        if col in Symbol.(names(df))
            mapping = Dict(row[id] => row[col] for row in eachrow(df))  # Store existing values
            dfNew[!, col] = [get(mapping, rid, missing) for rid in idsALL]  # Align data
        else
            dfNew[!, col] = fill(missing, length(idsALL))  # Fill missing columns
        end
    end
    return dfNew
end


function mergeDFs(dfs::Vector{DataFrame}, id::Symbol = :Gene, option::String = "sum")
    """
    "Merges a vector of DataFrames by summing or averaging them cell-wise.
    
    Arguments:
    - dfs::Vector{DataFrame}: A vector of DataFrames to merge.
    - id::Symbol: The identifier column (common to every DataFrame).
    - option::String: The operation to perform: either "sum" or "avg".
    
    Returns:
    A DataFrame where the first column is the identifier and the remaining columns are the cell-wise
    sum or average (depending on the option) of the numeric values from all data frames (with columns sorted in alphabetical order).
    """

    # 1. Compute the full set of row identifiers
    idsALL = reduce(union, [unique(df[!, id]) for df in dfs])
    sort!(idsALL)

    # 2. Compute the full set of numeric columns
    allCols = reduce(union, [setdiff(names(df), [string(id)]) for df in dfs])
    allCols = sort(Symbol.(allCols))

    # 3. Complete all DataFrames to have the same structure
    completed_dfs = [completeDF(df, id, idsALL, allCols) for df in dfs]

    # 4. Initialize matrices for sum and count tracking
    nRows, nCols = length(idsALL), length(allCols)
    sumMat = zeros(nRows, nCols)
    countMat = zeros(Int, nRows, nCols)

    # 5. Compute sum and count matrices
    for df in completed_dfs
        mat = Matrix(df[!, allCols])
        mat = coalesce.(mat, 0.0)
        sumMat .+= mat
        countMat .+= (mat .!= 0)  # Count nonzero interactions
    end

    # 6. Compute sum or average
    resultMat = option == "avg" ? sumMat ./ max.(countMat, 1) : sumMat

    # 7. Create the merged DataFrame
    merged = DataFrame(id => idsALL)
    for (j, col) in enumerate(allCols)
        merged[!, col] = resultMat[:, j]
    end

    return merged
end


function check_column_norms(df::DataFrame)

    """
    Checks if the columns are of length 1 (L2 Norm). Returns the number of true and false cases
    """
    details = Dict{String,Bool}()
    count_true = 0
    count_false = 0

    for col in names(df)[2:end]
        # Check if the column is numeric by testing its element type.
        if !(eltype(df[!, col]) <: Number)
            println("Column: ", col, " is non-numeric, skipping...")
            continue
        end
    # Convert the column values to Float64.
    values = Float64.(df[!, col])
    # Compute the L₂ norm.
    col_norm = norm(values, 2)
    # Check whether the norm is approximately 1.
    is_norm_one = isapprox(col_norm, 1.0; atol=1e-8)
    
    details[string(col)] = is_norm_one
    if is_norm_one
        count_true += 1
    else
        count_false += 1
    end
    # println("Column: ", col, " Norm: ", col_norm, " ~ 1? ", is_norm_one)
    end

    println("Total numeric columns normalized (True): ", count_true)
    println("Total numeric columns not normalized (False): ", count_false)
    return details, count_true, count_false
end

function writeTSVWithEmptyFirstHeader(df::DataFrame, filepath::String; delim::Char='\t')
    open(filepath, "w") do io
    # Create custom header:
    # Make the first header cell an empty string.
    # The remaining header cells are the remaining column names.
    hdr = [""; names(df)[2:end]]
    # Write the header row using the delimiter.
    println(io, join(hdr, delim))

       # Write each row. Each row is converted into strings.
       for row in eachrow(df)
        # Convert every element of the row to a string.
        row_values = [string(x) for x in collect(row)]
        println(io, join(row_values, delim))
        end
    end
end



#=
# Test workflow with small data

Example DataFrames with row names (represented as the first column in the DataFrame)
df1 = DataFrame(RowName = ["W", "X"], A = [1, 2], B = [3, 4], C = [5, 6])
df2 = DataFrame(Gene = ["Y", "X"], B = [7, 8], C = [9, 10], D = [11, 12])
df3 = DataFrame(Blue = ["P", "Q"], A = [13, 14], D = [15, 16], E = [17,18])

# # # # List of all DataFrames you want to combine
dfs = [df1, df2, df3]
dfs = [convertToLong(df) for df in dfs]
dfs = [convertToWide(df; indices = (1,2,3)) for df in dfs]

tick()
dfNorm = [frobeniusNormalize1(df, :row) for df in dfs]
tock()

commonID = :Gene
dfNorm = [rename!(df, names(df)[1] => commonID) for df in dfNorm]

tick()
mergeDFs(dfNorm, :Gene, "sum")
tock()

=#



#Compute Frobenius norm for each row (excluding the first column)
# norm_dfs = [DataFrame(A = df[!, 1], FrobeniusNorm = [norm(row[2:end], 2) for row in eachrow(df)]) for df in dfs]   # row
# norm_dfs = [DataFrame(Column = names(df)[2:end], FrobeniusNorm = [ norm(col, 2) for col in eachcol(df[!, 2:end]) ]) for df in dfs]  # col
    