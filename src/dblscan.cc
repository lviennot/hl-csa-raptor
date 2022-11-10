#include <iostream>

#include "timetable.hh"
#include "double_scan.hh"
#include "connection_scan.hh"
#include "logging.hh"
#include "file_util.hh"

void usage_exit (char **argv) {
    auto paragraph = [](std::string s, int width=80) -> std::string {
        std::string acc;
        while (s.size() > 0) {
            int pos = s.size();
            if (pos > width) pos = s.rfind(' ', width);
            std::string line = s.substr(0, pos);
            acc += line + "\n";
            s = s.substr(pos);
        }
        return acc;
    };
    
    std::cerr <<"Usage: "<< argv[0]
              <<" [options] gtfs_dir\n"
              <<"\n";
    exit(1);
}

std::vector<std::string> get_args(int argc, char **argv) {
    std::vector<std::string> a;
    for (int i = 0; i < argc; ++i) {
        std::string s(argv[i]);
        if (s[0] == '-') continue; // option
        a.push_back(s);
    }
    return a;
}

bool has_opt(int argc, char **argv, std::string opt) {
    assert(opt[0] == '-');
    for (int i = 0; i < argc; ++i) {
        std::string s(argv[i]);
        if (s == opt) return true;
    }
    return false;
}

std::string get_opt(int argc, char **argv,
                    std::string opt_prefix, std::string dft) {
    int len = opt_prefix.size();
    assert(opt_prefix[0] == '-'
           && (opt_prefix[len-1] == '=' || opt_prefix[len-1] == ':'));
    for (int i = 0; i < argc; ++i) {
        std::string s(argv[i]);
        if (s.size() >= len && s.substr(0, len) == opt_prefix) {
            return s.substr(len);
        }
    }
    return dft;
}

int get_int_opt(int argc, char **argv,
                std::string opt_prefix, int dft) {
    return std::stoi(get_opt(argc, argv, opt_prefix, std::to_string(dft)));
}

typedef pareto_rev<int> pset;

int main (int argc, char **argv) {
    logging main_log("--");

    // ------------------------ usage -------------------------
    std::vector<std::string> args = get_args(argc, argv);
    if (args.size() != 2) {
        usage_exit(argv);
    }
    std::string dir{args[1]};
    std::cerr <<"--------------- "<< dir <<" ---------------------\n";
    dir += "/";

    // ------------------------ time -------------------------
    main_log.cerr() << "start\n";
    double t = main_log.lap();

    //std::cerr <<"first rands: "<< rand() <<" "<< rand() <<" "<< rand() <<"\n";


    // ------------------------- load timetable ----------------------
    /*
    timetable ttbl{args[1], args[2],
            dir+"calendar.txt", dir+"calendar_dates.txt",
            dir+"trips.txt", dir+"stop_times.txt", dir+"transfers.txt", true};
    */
    timetable ttbl{dir+"stop_times.csv.gz",
            dir+"walk_and_transfer_inhubs.gr.gz",
            dir+"walk_and_transfer_outhubs.gr.gz",
            dir+"transfers.csv.gz", true};
    //dir+"walking_and_transfers.gr", t_from, t_to};
    std::cerr << ttbl.n_r <<" routes, "<< ttbl.n_st <<" sations, "
              << ttbl.n_s <<" stops\n";
    main_log.cerr(t) << "timetable\n";
    t = main_log.lap();
    //exit(0);

    // --------------- earliest arrival time through CSA and DS ---------
    double_scan<cost_eat> ds(ttbl);
    main_log.cerr(t) << "double_scan initialized\n";
    t = main_log.lap();

    connection_scan csa(ttbl);
    main_log.cerr(t) << "csa initialized\n";
    t = main_log.lap();

    const bool hub=true, trf=false;
    const int chg=std::stoi(get_opt(argc, argv, "-min-change-time=", "60")),
        km=48;

    uint64_t sum = 0, n_ok = 0;


    // ------------------ read queries -------------------
    int n_q = 0;
    std::vector<std::tuple<int, int, int> > queries;
    if (get_opt(argc, argv, "-query-file=", "queries-unif.csv") != "") {
        auto rows = read_csv
            (dir + get_opt(argc, argv, "-query-file=", "queries-unif.csv"), 3,
             "source", "destination", "departure_time");
             //"source", "target", "time");
        int n_q_max = std::stoi(get_opt(argc, argv, "-nq=", "10000"));
        for (auto r : rows) {
            if (n_q >= n_q_max) break;
            if (has_opt(argc, argv, "-1") && n_q >= 1) break;
            if (has_opt(argc, argv, "-10") && n_q >= 10) break;
            if (has_opt(argc, argv, "-100") && n_q >= 100) break;
            int src = ttbl.id_to_station[r[0]];
            int dst = ttbl.id_to_station[r[1]];
            int t = std::stoi(r[2]);
            queries.push_back(std::make_tuple(src, dst, t));
            ++n_q;
        }
        main_log.cerr(t) << n_q << " queries\n";
        t = main_log.lap();
    }
    // */


    // make andom successful queries
    if (get_opt(argc, argv, "-random-queries=", "") != "") {
        n_q = std::stoi(get_opt(argc, argv, "-random-queries=", "1000"));
        int max_delay = std::stoi(get_opt(argc, argv, "-max-delay=",
                                          std::to_string(ttbl.t_max)));
        int t_beg = 0*3600, t_end = 24*3600;
        int n_try = 0, n_err = 0;
        while (queries.size() < n_q) {
            ++n_try;
            int src = rand() % ttbl.n_st;
            int dst = rand() % ttbl.n_st;
            int t = t_beg + rand() % (t_end - t_beg);
            int arr = csa.earliest_arrival_time(src, dst, t, false, true,
                                                 chg);
            if (arr <= t + max_delay) {
                ++n_ok;
                queries.push_back(std::make_tuple(src, dst, t));
            }
        }
        main_log.cerr(t) << n_q <<" Random queries, success rate : "
                         << (n_q*100/n_try) <<"%\n";
        
        t = main_log.lap();
    }
    // */



    //* Arrival times
    if (has_opt(argc, argv, "-arrival-times")) {
        std::cout <<"src,dst,tdep,eat,eat_Unrestricted_Walking,eat_Walk_Only\n";
        n_ok = 0;
        for (auto q : queries) {
            int src = std::get<0>(q);
            int dst = std::get<1>(q);
            int t = std::get<2>(q);
            int arr1 = csa.earliest_arrival_time(src, dst, t, false, true, chg);
            int arr2 = csa.earliest_arrival_time_opt(src, dst, t, false, true, chg);
            assert(arr1 == arr2); // can fail if chg == 0
            std::cout << ttbl.station_id[src] <<","<< ttbl.station_id[dst]
                      <<","<< t <<","<< arr1
                      <<","<< (t + ttbl.walking_time(src, dst)) <<"\n";
            std::cout.flush();
            ++n_ok;
        }
        main_log.cerr(t) << n_ok << " arrival times\n";
        t = main_log.lap();
    }
    // */


    // ------------------------ end -------------------------
    main_log.cerr() << "end "<< dir <<"\n";

}
