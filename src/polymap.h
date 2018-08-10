//
// Created by Yue Hao on 4/4/17.
//
#include <vector>

#ifndef PYTPSA_POLYMAP_H
#define PYTPSA_POLYMAP_H
class CPolyMap{
private:
    int dim;
    int max_order;
    std::vector<int> decomposite(const int& i);
public:
    CPolyMap(const int& dim, const int& order);
    std::vector<std::vector<int> > map;
    std::vector<int> getindexmap(const int& i) {
        return map[i];
    }
    void setindexmap();
    
    
};
#endif //PYTPSA_POLYMAP_H
