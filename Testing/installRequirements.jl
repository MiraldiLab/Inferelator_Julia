# Pkg module for managing packages
using Pkg
# Define the path to the file containing required package names
package_file = "/data/miraldiNB/Katko/Projects/Julia2/Testing/Packages.txt"

# Read package names from the file
required_packages = []
if isfile(package_file)
    # Open the file and read package names
    open(package_file) do file
        for line in eachline(file)
            line = strip(line) 
            if !isempty(line)
                push!(required_packages, line)
            end
        end
    end
else
    println("Error: Package requirements file '$package_file' not found.")
    exit(1)
end

# Check and install required packages
for pkg in required_packages
    try
        eval(Meta.parse("using $pkg"))    # check if package already exist by loading
    catch
        println("Installing $pkg...") #install package if not installed
        Pkg.add(pkg)
        eval(Meta.parse("using $pkg"))       # Load the package after installation
    end
end
