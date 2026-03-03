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
#include "function.hpp"


using namespace std;

 double gain_th = 20;
 int generation_num = 100;
 int generation_time = 200;

int main(int argc, const char *argv[]) {

    vector<ele_unit> ele;
    vector<string> jl_source;
    vector<string> arg_arr(argv, argv+argc);  //コマンドライン引数が格納されている動的配列  コマンドライン引数が格納されている静的配列(argv)の要素を動的配列(arg_arr)に格


    cout << " Julia source file: " << arg_arr[1] << endl;

    int menu;

    cout << " Please select targets of optimization" << endl << endl;

    cout << " 1. C, Ip, wp" << endl;
    cout << " 2. C, Ip" << endl;
    cout << " 3. C" << endl;
    cout << "  Selected: ";

    cin >> menu;
    cout << endl;

    cout << " gain threshold (gain_th): " ;
    cin >> gain_th;
    cout << endl << " number of generation: ";
    cin >> generation_num;
    cout << endl << " genetic times: ";
    cin >> generation_time;
    cout << endl;

    read_jl(ele, jl_source, arg_arr[1]);
    calculation(ele, jl_source);

    switch(menu){
        case 1: run_nsga2_pagmo(generation_num, generation_time, ele, jl_source, out_value(ele, "Lj"), 1e-17, 1e-13, 1e-17, 1e-13, 5e-6, 6e-6, 8e9, 9e9);  //C, Ip, wp
        case 2: run_nsga2_pagmo_Ip(generation_num, generation_time, ele, jl_source, out_value(ele, "Lj"), 1e-17, 1e-13, 1e-17, 1e-13, 5e-6, 6e-6);  //C, Ip
        case 3: run_nsga2_pagmo_without(generation_num, generation_time, ele, jl_source, out_value(ele, "Lj"), 1e-17, 1e-13, 1e-17, 1e-13);  //C
        default: return 0;
    }


    return 0;
}
