#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <cmath>
#include "function.hpp"

using namespace std;

void execute_julia(string jl_source){
    stringstream command_jl, delete_jl ;
    //command_jl << "julia -J TWPA_sysimg.so " << jl_source;
    command_jl << "julia " << jl_source;
    delete_jl << "rm -rf " << jl_source;
    if(system(command_jl.str().c_str()) == -1){ //Juliaの実行が失敗した場合
        cout << "error:julia was not executed correctly." << endl;

    
    }
    system(delete_jl.str().c_str());

}

