#pragma once

#include "common.h"
#include "lognum.h"
#include "subtable.h"

namespace nonsymmetric_ {

template <class T>
std::vector<SubTable<T>> calculate_hat_weights(int size, const std::vector<std::vector<T>>& weights) {
	/*
		Section 3.1 in the article
	*/

    std::vector<std::vector<T>> hat_weights_1;
    std::vector<SubTable<T>> hat_weights_2;

    for (int i = 0; i < size+1; ++i)
    {
        std::vector<T> vectA1(((size_t)1 << size), T::zero());
        SubTable<T> vectB1(size);
        hat_weights_1.push_back(vectA1);
        hat_weights_2.push_back(vectB1);
    }


    for (int i = 0; i < size; ++i)
    {
        std::bitset<32> V(((size_t)1 << size)-1);
        V[i] = 0;
        int V_sub_i = (int)V.to_ulong();

        hat_weights_1[i][0] = weights[i][0];
        hat_weights_2[i](0, 0) = weights[i][0];

        //Go through all nonempty subsets of V\{i};
        for (int t = 0; (t=(t-V_sub_i)&V_sub_i);)
        {
            T sum_1 = weights[i][0];
            for (int S1 = 0; (S1=(S1-t)&t);)
            {
                sum_1 = sum_1 + weights[i][S1];
            }
            hat_weights_1[i][t] = sum_1;

            /*SINGLETON CASE*/
            for(int p = 0; p < size; p++) {
                int node = (size_t)1 << p;
                if((t&node) != 0) {
                    T sum_2 = T::zero();
                    for(int S = 0; (S=(S-t)&t);) {
                        if((S&node) != 0) {
                            sum_2 = sum_2 + weights[i][S];
                        }
                    }
                    hat_weights_2[i](node, t) = sum_2;
                }
            }

            for(int R = 0; (R=(R-t)&t);) {
                T sum = T::zero();
                std::bitset<32> R_bits(R);
                int previous = 0;
                for(int k = 0; k < size; k++) {
                    if(R_bits[k] == 1) {
                        int k_num = (size_t)1 << k;
                        int set = t;

                        set = set&(~previous);
                        previous = previous|k_num;
                        sum = sum + hat_weights_2[i](k_num, set);
                    }
                }
                hat_weights_2[i](R, t) = sum;
            }

            hat_weights_2[i](0,t) = hat_weights_1[i][t];
        }
    }

    return hat_weights_2;
}

template <class T>
SubTable<T> monotone_calculate_fs(int size, const std::vector<SubTable<T>> &hws) {
	/*
	Section 3.1.1
	MONOTONE VERSION.
	*/
    SubTable<T> fs(size);

    fs(0, 0) = T::one();

    std::bitset<32> V(((size_t)1 << size)-1);
    int V_sub = (int) V.to_ulong();

    for (int U=0; (U=(U-V_sub)&V_sub);) {
        for (int S_0 = 0; (S_0=(S_0-U)&U);) {
        int upmask = U&(~S_0);
        if(S_0 == U) {
            fs(S_0,U) = T::one();
        } else {
            T sum1 = T::zero();
            for (int S_1 = 0; (S_1=(S_1-upmask)&upmask);) {

                T product = T::one();
                std::bitset<32> S_1_bits(S_1);

                for (int i = 0; i < size; i++) {
                    if(S_1_bits[i] == 1) {
                        product = product*hws[i](S_0, V_sub&(~upmask));
                    }
                }
                product = product*fs(S_1, U&(~S_0));
                sum1 = sum1 + product;
            }
            fs(S_0, U) = sum1;
            }
        }
    }

    return fs;
}

template <class T>
std::vector<int> sample_layering(int size, const std::vector<SubTable<T>>& hws, SubTable<T> fs) {
	/*
	Section 3.2.
	*/
    std::vector<int> layering;
    layering.push_back(0);
    
    int partition_count = 0;
    int j = 1;
    int previous_rs = 0;
    int U;
    int V = ((size_t)1 << size)-1;


    while(partition_count < size) {
        U = V&(~previous_rs);
        std::bitset<32> U_bits(U);
        T upper_bound = T::zero();

        std::vector<T> bounds(U+1, T::zero());
        for(int R = 0; (R=(R-U)&U);) {
            T product = T::one();
            T inner_product = T::one();
            int r_k_union = 0;

            for (int k = 1; k <= j-1; k++) {
                int r_mk = layering[k-1];
                r_k_union = r_k_union|r_mk;

                std::bitset<32> r_k_bits(layering[k]);
                for (int i = 0; i < size; i++) {
                    if(r_k_bits[i] == 1) {
                        T prev = hws[i](r_mk, r_k_union);
                        inner_product = inner_product*prev;
                    }
                }
            }

            product = inner_product;

            int r_mj = layering[j-1];
            std::bitset<32> r_bits(R);

            for (int i = 0; i < size; i++) {
                if(r_bits[i] == 1) {
                    T prev = hws[i](r_mj, V&(~U));
                    product = product*prev;
                }
            }
            product = product*fs(R, U);
            upper_bound = upper_bound + product;

            bounds[R] = upper_bound;
        }

        T random_number = T::uniform_rand(upper_bound);
        for(int R = 0; (R=(R-U)&U);) {
            if(bounds[R] > random_number) {
                layering.push_back(R);
                previous_rs = previous_rs|R;
                break;
            }
        }

        std::bitset<32> added(layering[j]);
        for (int i = 0; i < size; i++) {
            if(added[i] == 1) {
                partition_count = partition_count +1;
            }
        }
        j = j + 1;
    }

    return layering;
}

template <class T>
std::vector<int> sample_parents_ns(int size, const std::vector<int>& layering,
	const std::vector<std::vector<T>>& weights) {
	/*
	Section 3.2.
	*/


    std::vector<int> dag(size, 0);

    int U = 0;
    U = U|layering[1];

    int previous_partition = U;

    for(int j = 2; j < (int) layering.size(); j ++) {
        std::bitset<32> layer_bits(layering[j]);
        for(int p = 0; p < size; p++) {
            if(layer_bits[p] == 0) {
                continue;
            }
            int node = p;
            T upper_bound = T::zero();
            for(int G = 0; (G=(G-U)&U);) {
                if((G&previous_partition) != 0) {
                    upper_bound = upper_bound + weights[node][G];
                }
            }

            T random = T::uniform_rand(upper_bound);
            T cumulative = T::zero();
            for(int G = 0; (G=(G-U)&U);) {
                if((G&previous_partition) != 0) {
                    cumulative = cumulative + weights[node][G];
                    if(cumulative > random) {
                        std::bitset<32> parent_bits(G);
                        for (int i = 0; i < size; i++) {
                            if(parent_bits[i] == 1) {
                                dag[node] |= 1 << i;
                            }
                        }
                        break;
                    }
                }
            }
        }
        previous_partition = layering[j];
        U = U|previous_partition;
    }

    return dag;
}

}

template <class T>
class NonSymmetricSampler {
public:
    typedef std::vector<std::vector<T>> WeightT;

    NonSymmetricSampler(WeightT weights) : weights(std::move(weights)) {
        preprocess();
    }

    std::vector<int> sample() const {
        using namespace nonsymmetric_;

        std::vector<int> layering = sample_layering<T>(weights.size(), h, non_symmetric_fs2);
        return sample_parents_ns<T>(weights.size(), layering, weights);
    }

private:
    int size;
    WeightT weights;
    std::vector<SubTable<T>> h;
    SubTable<T> non_symmetric_fs2;

    void preprocess() {
        using namespace nonsymmetric_;

        h = calculate_hat_weights<T>(weights.size(), weights);
        non_symmetric_fs2 = monotone_calculate_fs<T>(weights.size(), h);
    }
};
