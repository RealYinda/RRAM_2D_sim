//
// 文件名:      Nedelec2D.h
// 描述  :      2D 三角形一阶 Nedelec 棱边元.
//              basis() 为 3 参数签名, 与 TetQuad::quadcalculateKij 模板兼容.
//              edge_order 方向处理内置于 basis() 中.
//

#ifndef included_appu_Nedelec2D
#define included_appu_Nedelec2D

#include "Pointer.h"
#include "Array.h"
#include "Patch.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "EdgeData.h"
#include "CellData.h"
#include "NodeData.h"

namespace JAUMIN {
namespace appu {

class Nedelec2D {
public:
  Nedelec2D(const hier::Patch<NDIM>& patch,
            tbox::Pointer<pdat::EdgeData<NDIM, bool> > edge_order,
            tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian);

  int nbas() const;   // 3 (三角形 3 条边)
  int dim() const;    // 2

  /// 3 参数签名, 与 TetQuad/TriQuad::quadcalculateKij 模板兼容.
  void basis(const int cell, const double* lambda, double* values) const;

  /// 2D 旋度 (标量).
  void curl(const int cell, const double* lambda, double* values) const;

private:
  const hier::Patch<NDIM>& d_patch;
  tbox::Pointer<pdat::EdgeData<NDIM, bool> > d_edge_order;
  tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_jacobian;
  tbox::Pointer<pdat::NodeData<NDIM, double> > d_node_coord;
  tbox::Array<int> d_cell_edge_ext;
  tbox::Array<int> d_cell_edge_idx;
  tbox::Array<int> d_edge_node_ext;
  tbox::Array<int> d_edge_node_idx;
  tbox::Array<int> d_cell_node_ext;
  tbox::Array<int> d_cell_node_idx;
};

}  // namespace appu
}  // namespace JAUMIN

#endif  // included_appu_Nedelec2D
