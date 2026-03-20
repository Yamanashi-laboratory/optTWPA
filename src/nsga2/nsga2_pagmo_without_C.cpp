#include "function.hpp"
#include <random>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <fstream>

#if __has_include(<pagmo/pagmo.hpp>)
#include <pagmo/pagmo.hpp>
#include <pagmo/algorithms/nsga2.hpp>
#include <pagmo/population.hpp>
#include <pagmo/types.hpp>
#define PAGMO_AVAILABLE 1
#else
#pragma message("pagmo library not found, run_nsga2_pagmo will be disabled")
#define PAGMO_AVAILABLE 0
#endif

#if PAGMO_AVAILABLE
using pagmo::vector_double;

using pagmo::problem;
using pagmo::algorithm;
using pagmo::population;
using pagmo::nsga2;


//======================================
// UDP: JosephsonCircuits を呼び出して (gain, bandwidth, ripple) を得る
//======================================
struct josephson_problem_without_C {
    double Lj;
    std::vector<ele_unit> ele;
    std::vector<std::string> jl_source;
    double Ip_min, Ip_max, wp_min, wp_max;

    josephson_problem_without_C() = default;
    josephson_problem_without_C(double Lj_, const std::vector<ele_unit>& ele_,
                      const std::vector<std::string>& jl_src_)
        : Lj(Lj_), ele(ele_), jl_source(jl_src_) {}

    result eval(double Ip, double wp) const {
        auto tmp = ele;
        change_param(tmp, "Ip", Ip);
        change_param(tmp, "wp", wp);
        return calculation(tmp, jl_source);
    }

    vector_double fitness(const vector_double &x) const {
        result r = eval(x[0], x[1]);
        double f1 = -r.gain;
        double f2 = -r.bandwidth;
        if (r.ripple > 0.1) {
            f1 = 1e6 + r.ripple;
            f2 = 1e6 + r.ripple;
        }
        return {f1, f2};
    }

    std::pair<vector_double, vector_double> get_bounds() const {
        return {{Ip_min, wp_min}, {Ip_max, wp_max}};
    }

    std::size_t get_nobj() const { return 2u; }
    std::size_t get_nic() const { return 0u; }
    std::size_t get_nec() const { return 0u; }

    void set_bounds(double _Ip_min, double _Ip_max,
                    double _wp_min, double _wp_max) {
        Ip_min = _Ip_min;  Ip_max = _Ip_max;
        wp_min = _wp_min;  wp_max = _wp_max;
    }

    std::string get_name() const { return "josephson_problem"; }
};

//======================================
// run_nsga2_pagmo ：Pagmo を使った NSGA‐II 実行
//======================================
void run_nsga2_pagmo_without_C(int pop_size,
                     int generations,
                     const std::vector<ele_unit>& ele,
                     const std::vector<std::string>& jl_source,
                     double Lj,
                     double Ip_min, double Ip_max,
                     double wp_min, double wp_max) {
    josephson_problem_without_C prob_udp(Lj, ele, jl_source);
    prob_udp.set_bounds(Ip_min, Ip_max, wp_min, wp_max);
    pagmo::problem prob{prob_udp};

    pagmo::algorithm algo{ pagmo::nsga2(
        generations,
        0.9,
        20.0,
        1.0/2.0,
        20.0
    ) };

    pagmo::population pop{prob, static_cast<unsigned int>(pop_size)};
    pop = algo.evolve(pop);

    std::ofstream ofs("nsga2_result_without_C.csv");
    ofs << "Ip,wp,gain,bandwidth,ripple\n";

    std::cout << "# Ip\twp\tgain\tbandwidth\tripple\n";
    for (std::size_t i = 0; i < pop.size(); ++i) {
        const auto xv = pop.get_x()[i];
        double Ip = xv[0];
        double wp = xv[1];
        result r = prob_udp.eval(Ip, wp);
        if (r.ripple <= 0.1) {
            std::cout << Ip << "\t"
                      << wp << "\t"
                      << r.gain << "\t"
                      << r.bandwidth << "\t"
                      << r.ripple << "\n";
            ofs << Ip << ','
                << wp << ','
                << r.gain << ','
                << r.bandwidth << ','
                << r.ripple << '\n';
        }
    }
    ofs.close();
}
#else
void run_nsga2_pagmo(int pop_size,
                     int generations,
                     const std::vector<ele_unit>& ele,
                     const std::vector<std::string>& jl_source,
                     double Lj,
                     double Cg_min, double Cg_max,
                     double Cc_min, double Cc_max,
                     double Ip_min, double Ip_max,
                     double wp_min, double wp_max) {
    (void)pop_size; (void)generations; (void)ele; (void)jl_source;
    (void)Lj; (void)Cg_min; (void)Cg_max; (void)Cc_min; (void)Cc_max;
    (void)Ip_min; (void)Ip_max; (void)wp_min; (void)wp_max;
    std::cerr << "run_nsga2_pagmo is disabled because pagmo library is not available." << std::endl;
}
#endif
