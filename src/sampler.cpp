#include "common.h"
#include "nonsymmetric.h"
#include "symmetric.h"
#include "readwrite.h"

void write_dags(const std::vector<std::vector<int>>& dags) {
    for(const std::vector<int>& dag : dags) {
        int size = dag.size();

        for(int i = 0; i < size; ++i) {
            if(i) {
                std::cout << ", ";
            }

            std::cout << i << " <- {";

            bool first = true;
            for(int j = 0; j < size; ++j) {
                if(dag[i] & (1 << j)) {
                    if(!first) {
                        std::cout << ", ";
                    }
                    first = false;
                    std::cout << j;
                }
            }
            std::cout << "}";
        }
        std::cout << "\n";
    }
}

template <class Sampler>
void run_sampler(int number_of_dags, typename Sampler::WeightT weights) {
    std::cerr << "Sampling " << number_of_dags << " DAGs\n";

    clock_t begin = clock();

    Sampler sampler(std::move(weights));

    clock_t mid = clock();
    
    std::vector<std::vector<int>> dags;
    for (int i = 0; i < number_of_dags; ++i) {
        dags.push_back(sampler.sample());
    }

    clock_t end = clock();

    write_dags(dags);
    double pre_elapsed_secs = double(mid - begin) / CLOCKS_PER_SEC;
    double samp_elapsed_secs = double(end - mid) / CLOCKS_PER_SEC;
    std::cerr << "Precomputation: " << pre_elapsed_secs << "s\n";
    std::cerr << "Per DAG: " << samp_elapsed_secs / number_of_dags << "s\n";
}

void usage() {
    std::cerr << "Usage:\n";
    std::cerr << "    ./sampler symmetric uniform <number_of_nodes> <number_of_dags>\n";
    std::cerr << "    ./sampler symmetric input <input_file> <number_of_dags>\n";
    std::cerr << "    ./sampler nonsymmetric <input_file> <number_of_dags>\n";
}

int main(int argc, char* argv[]) {
    int argi = 1;
    auto getArg = [&]() {
        if(argi >= argc) {
            std::cerr << "Too few command line arguments\n";
            usage();
            exit(1);
        }
        return argv[argi++];
    };
    auto argsDone = [&]() {
        if(argi != argc) {
            std::cerr << "Extra command line arguments\n";
            usage();
            exit(1);
        }
    };

    std::string symmetry_type = getArg();

    if(symmetry_type == "symmetric") {
        std::string weight_arg = getArg();

        if(weight_arg == "uniform") {
            int size = std::stoi(getArg());
            int n_dags = std::stoi(getArg());
            argsDone();

            std::vector<Lognum> weights(size, Lognum::one());
            run_sampler<SymmetricSampler<Lognum>>(n_dags, std::move(weights));
        } else if (weight_arg == "input") {
            std::string input = getArg();
            int n_dags = std::stoi(getArg());
            argsDone();

            std::vector<Lognum> weights = read_symmetric_weights<Lognum>(input);
            run_sampler<SymmetricSampler<Lognum>>(n_dags, std::move(weights));
        } else {
            std::cerr << "Unknown weight type " << weight_arg << "\n";
            usage();
            return 1;
        }
    } else if (symmetry_type == "nonsymmetric") {
        std::string input = getArg();
        int n_dags = std::stoi(getArg());
        argsDone();
        
        std::vector<std::vector<Lognum>> weights = read_nonsymmetric_weights<Lognum>(input);
        run_sampler<NonSymmetricSampler<Lognum>>(n_dags, std::move(weights));
    } else {
        std::cerr << "Unknown symmetry type " << symmetry_type << "\n";
        usage();
        return 1;
    }

    return 0;
}
