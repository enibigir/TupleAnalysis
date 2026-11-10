g++ -std=c++11 -fPIC `root-config --cflags` -I inc -c src/RegionMassPlot.cc 
g++ -std=c++11 -fPIC `root-config --cflags` -I inc -c src/MassTools.cc 
g++ -std=c++11 -fPIC `root-config --cflags` -I inc -c src/CPlots.cc 
g++ -fPIC -shared CPlots.o MassTools.o RegionMassPlot.o `root-config --glibs` -o libTools.so
rm *.o
