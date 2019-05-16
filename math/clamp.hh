#include <algorithm>

namespace sfm {
namespace math {

template <class T>
inline T clamp(const T x, const T lower, const T upper) {
    return std::max(std::min(x, upper), lower);
}

}
}
