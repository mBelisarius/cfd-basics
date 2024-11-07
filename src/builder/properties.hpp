#ifndef CFD_BASICS_BUILDER_PROPERTIES_HPP_
#define CFD_BASICS_BUILDER_PROPERTIES_HPP_

#include <map>

namespace cfd_basics {

enum class Property {
  kCondutivity,
  kDensity,
  kHeatSource,
  kThermalDiffusivity
};

template<typename Scalar>
class PropertiesList {
 public:
  PropertiesList() = default;

  explicit PropertiesList(std::map<Property, Scalar> props) : props_(props) {}

  std::map<Property, Scalar> Props() { return props_; }

  Scalar operator[](Property prop) { return props_[prop]; }

 private:
  std::map<Property, Scalar> props_;
};

}

#endif // CFD_BASICS_BUILDER_PROPERTIES_HPP_
