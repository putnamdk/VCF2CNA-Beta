bash genutil_build.sh
g++ -std=c++0x -c consprep.cpp 
g++ -std=c++0x -c snvcounts.cpp 
g++ -std=c++0x -o consprep consprep.o genutil.o 
g++ -std=c++0x -o snvcounts snvcounts.o genutil.o
