#ifndef CFD_BASICS_MESH_CELL_HPP_
#define CFD_BASICS_MESH_CELL_HPP_

#include <unordered_set>

#include "nuenv/core"

#include "mesh/face.hpp"

namespace cfd_basics {

class Cell {
 public:
  Cell() = default;

  explicit Cell(nuenv::Index id)
      : id_(id), facesId_() {}

  Cell(nuenv::Index id, std::unordered_set<nuenv::Index> facesId)
      : id_(id), facesId_(facesId) {}

  nuenv::Index Id() const { return id_; }

  void SetId(nuenv::Index id) { id_ = id; }

  std::unordered_set<nuenv::Index> FacesId() const { return facesId_; }

  nuenv::Index NFaces() { return facesId_.size(); }

  void AddFace(nuenv::Index faceId) { facesId_.insert(faceId); }

 private:
  nuenv::Index id_;
  std::unordered_set<nuenv::Index> facesId_;
};

} // cfd_basics

#endif // CFD_BASICS_MESH_CELL_HPP_
