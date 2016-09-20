output="BetaCore"
g++ -I ../EigenLib -fopenmp BetaCore.cpp -o $output
echo "Compile done, output is $output"
