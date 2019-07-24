#pragma once

#include "common.h"
#include "lognum.h"

namespace symmetric_ {

template <class T>
std::vector<std::vector<T>> calculate_hat_weights(int size, int l_bound, const std::vector<T>& weights) {
    /*
    Computes the \hat{w}(r, t)-values following the description in section 4.2 of the article.
    */
    std::vector<std::vector<T>> auxiliary(size+1);
    std::vector<std::vector<T>> hw(size+1);

    for (int i = 0; i < size+1; ++i)
    {
        std::vector<T> vect(size+1);
        std::vector<T> aux(size+1);
        hw[i] = vect;
        auxiliary[i] = aux;
    }


    hw[0][0] = T::one();

    for(int t = 1; t < size; t++) 
    {
        T sum;
        for(int j = 0; j <= t; j++) 
        {
            sum = sum+(T::binomial(t, j)*weights[j]);
        }
        hw[0][t] = sum;
    }


    for(int t = 1; t < size; t++) 
    {
        T sum;
        for(int j = 1; j <= t; j++) 
        {
            sum = sum+(T::binomial(t-1, j-1)*weights[j]);
        }
        hw[1][t] = sum;
    }

    for(int t = 1; t < size; t++) 
    {

        int bnd = (int) std::min(l_bound, t);
        
        for(int r = 2; r <= bnd; r++) 
        {
            T sum;
            sum = sum + hw[1][t];
            sum = sum + hw[r-1][t-1];
            hw[r][t] = sum;
        }
    }

    return hw;
}

template <class T>
std::vector<std::vector<T>> calc_ru_recursively(int size, int l_bound, std::vector<std::vector<T>> &hw) {
    /*
    Calculates the values f(r, u) (as they appear in the section 4.2. of the paper).
    */

    std::vector<std::vector<T>> rus;
    for (int i = 0; i < size+1; ++i)
    {
        std::vector<T> vect(size+1);
        rus.push_back(vect);
    }

    int bnd = (int) std::min(size, l_bound);

    for (int i = 0; i < bnd+1; i++)
    {
        rus[i][i] = T::one();
    }

    for (int i = 1; i < size+1; i++)
    {
        rus[0][i] = T::zero();
    }

    for (int u = 1; u < size+1; u++) 
    {
        int bnd1 = std::min(u-1, l_bound);

        for (int r = 1; r <= bnd1; r++) 
        {
            T sum;
            int bnd2 = std::min(u-r, l_bound);
            for (int r_prime = 1; r_prime <= bnd2; r_prime++) 
            {
                T temp;
                temp = hw[r][size-u+r].powi(r_prime);

                temp = temp*T::binomial(u-r, r_prime);
                temp = temp*rus[r_prime][u-r];
                sum = sum+temp;
            }
            rus[r][u] = sum;
        }
    }

    return rus;
}

template <class T>
T number_of_compatible_dags(int size, int u, int r, int previous_size,
    const std::vector<std::vector<T>>& rus, const std::vector<std::vector<T>>& hw) {
    /*
        See section 4.3 in the article.
    */
    
    T weightnumber = hw[previous_size][size-u];
    weightnumber = weightnumber.powi(r);
    weightnumber = weightnumber*(rus[r][u]*T::binomial(u, r));
    return weightnumber;
}

template <class T>
std::vector<int> sample_partition(int size, int l_bound, const std::vector<std::vector<T>>& rus, const std::vector<std::vector<T>>& hw) {
    /*
    Section 4.3 in the article.
    */
    std::vector<int> partition(size);

    int partition_count = 0;
    int previous_size = 0;

    int r, u;
    int partition_length = 0;

    while(partition_count < size) 
    {

        r = 1;
        u = size-partition_count;

        T upper_bound;

        if(previous_size == 0) {
            T up;
            int bnd = (int) std::min(u, l_bound);
            for (int i = 1; i <= bnd; ++i)
            {
                up = up + number_of_compatible_dags<T>(size, u, i, 0, rus, hw);
            }
            upper_bound = up;
        } else {
            upper_bound = rus[previous_size][previous_size + u];
        }

        T random_number = T::uniform_rand(upper_bound);

        bool notfound = true;

        T sum;

        while(notfound) 
        {

            sum = sum + number_of_compatible_dags<T>(size, u, r, previous_size, rus, hw);

            if(sum > random_number) {
                break;
            }

            r++;
        }

        partition[partition_length] = r;
        partition_length = partition_length + 1;
        partition_count += r;
        previous_size = r;
    }

    partition.resize(partition_length);


    return partition;
}


template <class T>
std::vector<int> sample_parents(int size, const std::vector<T>& weights, const std::vector<std::vector<T>>& hws, const std::vector<int>& partition) {
    /*
    Section 4.3 in the article.
    */
 
    std::vector<int> dag(size, 0);

    int partition_length = (int) partition.size();

    std::vector<std::vector<int>> layering(partition_length+1);
    std::vector<int> nodenames;
    
    for (int i = 0; i < size; ++i)
    {
        nodenames.push_back(i);
    }

    shuffle(nodenames.begin(), nodenames.end(), std::random_device());

    int index = 0;
    for(int i = 0; i < (int) partition.size(); i++) {
        
        std::vector<int> layer_vector;
        
        for(int j = 0; j < partition[i]; j++) {
            layer_vector.push_back(nodenames[index]);
            index++;
        }

        layering[i] = layer_vector;
    }

    std::vector<int> ancestors = {};
    for (int j = 1; j < partition_length; ++j)
    {
        std::vector<int> parent_layer = layering[j-1];
        std::vector<int> current_layer = layering[j];

        int parent_layer_size = (int) parent_layer.size();
        int current_layer_size = (int) current_layer.size();

        std::vector<T> upper_bounds_x(parent_layer_size);
        upper_bounds_x[0] = T::zero();
        T total = T::zero();
        std::vector<std::vector<int>> pxs;

        for(int xi = 0; xi < parent_layer_size; xi++) {
            std::vector<int> px = ancestors;
            for(int yi = 0; yi < parent_layer_size; yi++) {
                if(parent_layer[yi] > parent_layer[xi]) {
                    px.push_back(parent_layer[yi]);
                }
            }
            upper_bounds_x[xi] = hws[1][(int) px.size()+1];
            total = total + upper_bounds_x[xi];
            pxs.push_back(px);
        }

        for(int c = 0; c < current_layer_size; c++) {
            int current_node = current_layer[c];

            T random_number1 = T::uniform_rand(total);

            int xi = 0;

            T total2 = T::zero();

            for(int i = 0; i < parent_layer_size; i++) {
                total2 = total2 + upper_bounds_x[i];
                if(total2 > random_number1) {
                    xi = i;
                    break;
                }
            }

            int size_of_gi = 1;
            std::vector<int> px = pxs[xi];
            int size_of_px = (int) px.size() + 1;

            std::vector<T> upper_bounds_size(size_of_px + 1);
            upper_bounds_size[0] = T::zero();
            for(int gi = 1; gi <= size_of_px; gi++) {
                upper_bounds_size[gi] = upper_bounds_size[gi-1] + (T::binomial(size_of_px-1, gi-1)*weights[gi]);
            }

            T random_number2 = T::uniform_rand(upper_bounds_size[size_of_px]);

            for(int gi = 1; gi <= size_of_px; gi++) {
                if(upper_bounds_size[gi] > random_number2) {
                    size_of_gi = gi;
                    break;
                }
            }

            shuffle(px.begin(), px.end(), std::random_device());

            for (int i = 0; i < size_of_gi-1; ++i)
            {
                dag[current_node] |= 1 << px[i];
            }

            dag[current_node] |= 1 << parent_layer[xi];

        }

        for (int p = 0; p < parent_layer_size; ++p)
        {
            ancestors.push_back(parent_layer[p]);
        }

    }
    return dag;
}

}

template <class T>
class SymmetricSampler {
public:
    typedef std::vector<T> WeightT;
    
    SymmetricSampler(WeightT weights) : weights(std::move(weights)) {
        preprocess();
    }

    std::vector<int> sample() const {
        using namespace symmetric_;
        
        std::vector<int> partition = sample_partition<T>(weights.size(), weights.size(), rus, hw);
        return sample_parents<T>(weights.size(), weights, hw, partition);
    }
private:
    WeightT weights;
    std::vector<std::vector<T>> hw;
    std::vector<std::vector<T>> rus;

    void preprocess() {
        using namespace symmetric_;
        
        hw = calculate_hat_weights<T>(weights.size(), weights.size(), weights);
        rus = calc_ru_recursively<T>(weights.size(), weights.size(), hw);
    }
};
