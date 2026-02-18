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
    cout << " gain threshold (gain_th): " ;
    cin >> gain_th;
    cout << endl << " number of generation: ";
    cin >> generation_num;
    cout << endl << " genetic times: ";
    cin >> generation_time;

    read_jl(ele, jl_source, arg_arr[1]);
    calculation(ele, jl_source);
    run_nsga2_pagmo(generation_num, generation_time, ele, jl_source, out_value(ele, "Lj"), 1e-17, 1e-13, 1e-17, 1e-13, 5e-6, 6e-6, 8e9, 9e9);

    return 0;
}
