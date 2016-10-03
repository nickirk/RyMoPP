output="RyMoPP"
g++ -std=c++11 -I ../EigenLib -fopenmp Kore.cpp -o $output
if [ $? -eq 0 ]; then
	echo "Compile done, output is $output"
else
	echo "Compile failed!!"
fi
