
# This script will allow one to run the notebooks from the REPL
# IJulia.jl is required
# First, navigate to the Examples directory in Multitaper.jl and issue 
#
# (v1.0) pkg> activate .
#
# (Examples) pkg> instantiate
#  
# Then run this script from the REPL

using IJulia

notebook(dir=".")

