#ifndef CFD_BASICS_CORE_H_
#define CFD_BASICS_CORE_H_

#include <Eigen/Core>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <vector>

namespace cfd_basics
{
    using std::function;
    using std::size_t;

    template<typename T>
    inline T maxNumLimit()
    {
        return std::numeric_limits<T>::max();
    }

    template<typename T>
    using List = std::list<T>;

    template<typename Key, typename T>
    using Map = std::map<Key, T>;

    using string = std::string;

    template<typename T>
    using VectorSTL = std::vector<T>;

    template<typename Scalar, int SizeAtCompileTime = Eigen::Dynamic>
    using Vector = Eigen::Matrix<Scalar, SizeAtCompileTime, 1>;

    template<typename Scalar,
            int RowsAtCompileTime = Eigen::Dynamic,
            int ColsAtCompileTime = Eigen::Dynamic,
            int Options = 0,
            int MaxRowsAtCompileTime = RowsAtCompileTime,
            int MaxColsAtCompileTime = ColsAtCompileTime>
    using Matrix = Eigen::Matrix<Scalar,
                                 RowsAtCompileTime,
                                 ColsAtCompileTime,
                                 Options,
                                 MaxRowsAtCompileTime,
                                 MaxColsAtCompileTime>;

} // cfd_basics

#endif /* CFD_BASICS_CORE_H_ */
