#ifndef CFD_BASICS_MESH_CELL_HPP_
#define CFD_BASICS_MESH_CELL_HPP_

#include <unordered_set>

#include "nuenv/core"

#include "mesh/face.hpp"

namespace cfd_basics {

class Cell {
 public:
  Cell() = default;

  Cell(nuenv::Index id, std::unordered_set<nuenv::Index> facesId)
      : id_(id), facesId_(facesId) {}

  nuenv::Index Id() { return id_; }

  std::unordered_set<nuenv::Index> FacesId() const { return facesId_; }

  nuenv::Index NFaces() { return facesId_.size(); }

  void AddFace(nuenv::Index faceId) { facesId_.insert(faceId); }

 private:
  nuenv::Index id_;
  std::unordered_set<nuenv::Index> facesId_;
};

} // cfd_basics

#endif // CFD_BASICS_MESH_CELL_HPP_
