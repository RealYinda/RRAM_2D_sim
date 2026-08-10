//
// 文件名:      Nedelec2D.C
// 描述  :      2D 三角形一阶 Nedelec 棱边元实现.
//              basis() 从 d_cell_jacobian 自取 nabla,
//              签名与 3D Nedelec 一致 (cell, lambda, values).
//

#include "Nedelec2D.h"

namespace JAUMIN {
namespace appu {

Nedelec2D::Nedelec2D(const hier::Patch<NDIM>& patch,
                     tbox::Pointer<pdat::EdgeData<NDIM, bool> > edge_order,
                     tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian)
    : d_patch(patch),
      d_edge_order(edge_order),
      d_cell_jacobian(cell_jacobian) {
  patch.getPatchTopology()->getCellAdjacencyEdges(d_cell_edge_ext,
                                                   d_cell_edge_idx);
  patch.getPatchTopology()->getEdgeAdjacencyNodes(d_edge_node_ext,
                                                   d_edge_node_idx);
  patch.getPatchTopology()->getCellAdjacencyNodes(d_cell_node_ext,
                                                   d_cell_node_idx);
  d_node_coord = patch.getPatchGeometry()->getNodeCoordinates();
}

int Nedelec2D::nbas() const { return 3; }

int Nedelec2D::dim() const { return 2; }

void Nedelec2D::basis(const int cell, const double* lambda,
                       double* values) const {
  // φ_ij = λ_i ∇λ_j − λ_j ∇λ_i
  double(*phi)[2] = (double(*)[2])values;
  double(*nabla)[NDIM + 1] =
      (double(*)[NDIM + 1])(&((*d_cell_jacobian)(0, cell)));

  int nedge = 3;
  for (int ie = 0; ie < nedge; ie++) {
    int edge = d_cell_edge_idx[d_cell_edge_ext[cell] + ie];
    int v0 = d_edge_node_idx[d_edge_node_ext[edge] + 0];
    int v1 = d_edge_node_idx[d_edge_node_ext[edge] + 1];

    // 与 3D Nedelec::basis 一致: 若边反转则交换顶点
    if (!(*d_edge_order)(0, edge)) {
      int t = v0; v0 = v1; v1 = t;
    }

    // 全局顶点 → 局部顶点编号
    int lv0 = -1, lv1 = -1;
    for (int n = 0; n < 3; n++) {
      int nidx = d_cell_node_idx[d_cell_node_ext[cell] + n];
      if (nidx == v0) lv0 = n;
      if (nidx == v1) lv1 = n;
    }

    phi[ie][0] = lambda[lv0] * nabla[lv1][0] - lambda[lv1] * nabla[lv0][0];
    phi[ie][1] = lambda[lv0] * nabla[lv1][1] - lambda[lv1] * nabla[lv0][1];
  }
}

void Nedelec2D::curl(const int cell, const double* lambda,
                      double* values) const {
  // curl φ_ij = 2 ∇λ_i × ∇λ_j  (2D 标量)
  (void)lambda;
  double(*curlphi)[1] = (double(*)[1])values;
  double(*nabla)[NDIM + 1] =
      (double(*)[NDIM + 1])(&((*d_cell_jacobian)(0, cell)));

  int nedge = 3;
  for (int ie = 0; ie < nedge; ie++) {
    int edge = d_cell_edge_idx[d_cell_edge_ext[cell] + ie];
    int v0 = d_edge_node_idx[d_edge_node_ext[edge] + 0];
    int v1 = d_edge_node_idx[d_edge_node_ext[edge] + 1];

    if (!(*d_edge_order)(0, edge)) {
      int t = v0; v0 = v1; v1 = t;
    }

    int lv0 = -1, lv1 = -1;
    for (int n = 0; n < 3; n++) {
      int nidx = d_cell_node_idx[d_cell_node_ext[cell] + n];
      if (nidx == v0) lv0 = n;
      if (nidx == v1) lv1 = n;
    }

    double a = nabla[lv0][0], b = nabla[lv0][1];
    double d = nabla[lv1][0], e = nabla[lv1][1];
    curlphi[ie][0] = 2.0 * (a * e - d * b);
  }
}

}  // namespace appu
}  // namespace JAUMIN
