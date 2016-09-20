output="Paral2"
g++-4.8 -std=c++11 -I ../EigenLib -fopenmp Paral.2.cpp -o $output
echo "Compile done, output is $output"
