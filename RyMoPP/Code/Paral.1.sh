output="Beta40SP1007"
g++ -std=c++0x -I ../EigenLib -fopenmp Beta40SP.10.07.cpp -o $output
echo "Compile done, output is $output"
