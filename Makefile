
test: config.cpp config.h config_FNAL.cpp 
	g++ -ggdb --std=c++11 -o test config.cpp config_FNAL.cpp -lboost_system -lboost_program_options `root-config --cflags --glibs`
