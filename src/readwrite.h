#pragma once

#include "common.h"

template <typename T>
std::vector<T> read_symmetric_weights(const std::string& filename) {
    std::ifstream file;
    file.exceptions(file.failbit | file.badbit);
    file.open(filename);

    int size;
    file >> size;

    if(size <= 0) {
        std::cerr << "Invalid symmetric weight file\n";
        exit(1);
    }

    std::vector<T> weights(size);
    for(int i = 0; i < size; ++i) {
        double log_value;
        file >> log_value;
        weights[i] = T::from_log(log_value);
    }

    return weights;
}

template <typename T>
std::vector<std::vector<T>> read_nonsymmetric_weights(const std::string& filename) {
    std::ifstream file;
    file.exceptions(file.failbit | file.badbit);
    file.open(filename);

    int size;
    file >> size;

    if(size <= 0) {
        std::cerr << "Invalid nonsymmetric weight file\n";
        exit(1);
    }
    if(size >= 31) {
        std::cerr << "Too many nodes in nonsymmetric weight file\n";
        exit(1);
    }

    std::map<std::string, int> name_to_idx;
    for(int i = 0; i < size; ++i) {
        std::string name;
        int score_count;
        file >> name >> score_count;

        if(name_to_idx.count(name)) {
            std::cerr << "Invalid nonsymmetric weight file\n";
            exit(1);
        }
        name_to_idx[name] = i;

        for(int j = 0; j < score_count; ++j) {
            double log_score;
            int parent_count;
            file >> log_score >> parent_count;
            if(parent_count < 0) {
                std::cerr << "Invalid nonsymmetric weight file\n";
                exit(1);
            }

            for(int k = 0; k < parent_count; ++k) {
                std::string parent;
                file >> parent;
            }
        }
    }

    file.seekg(0);

    int ignore;
    file >> ignore;

    std::vector<std::vector<T>> weights(size);
    for (int i = 0; i < size; ++i) {
        weights[i].resize(1 << size, T::zero());

        std::string name;
        int score_count;
        file >> name >> score_count;
        
        for(int j = 0; j < score_count; ++j) {
            double log_score;
            int parent_count;
            file >> log_score >> parent_count;

            std::bitset<32> parents;

            for(int k = 0; k < parent_count; ++k) {
                std::string parent;
                file >> parent;

                auto it = name_to_idx.find(parent);
                if(it == name_to_idx.end() || it->second == i) {
                    std::cerr << "Invalid nonsymmetric weight file\n";
                    exit(1);
                }
                parents[it->second] = 1;
            }

            weights[i][parents.to_ulong()] = T::from_log(log_score);
        }
    }
    
    return weights;
}
