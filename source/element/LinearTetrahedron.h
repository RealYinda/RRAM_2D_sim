//
// 文件名:     LinearTetrahedron.h
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     单元计算类
// 类别　:     %Internal File% ( Don't delete this line )
//
#ifndef included_LinearTetrahedron
#define included_LinearTetrahedron
#include "BaseElement.h"
#include "BaseIntegrator.h"
#include "BaseShapeFunction.h"
#include "DoubleVector.h"
#include "UserConstant.h"

#include "IntegratorManager.h"
#include "ShapeFunctionManager.h"
#include "ElementManager.h"

using namespace JAUMIN;
class LinearTetrahedron : public BaseElement<NDIM> {
public:
  /**
   * @brief 构造函数.
   *
   * @param name 输入参数, 单元名字.
   *
   */
  LinearTetrahedron(const string& name);

  /**
   * @brief 析构函数.
   *
   */
  ~LinearTetrahedron();

  /**
   * @brief 计算单元质量矩阵.
   *
   * @param ele_mat        输出参数, 二维数组, 单元质量矩阵.
   * @param integrator     输入参数, 指针, 指向积分器.
   * @param shape_func     输入参数, 指针, 指向形函数.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   */
  virtual void buildStiffElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > coord, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat);

  /// 电场刚度矩阵 (电流连续性): ∫σ(T,nD) ∇N_i·∇N_j dΩ
  virtual void buildEStiffElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Array<double> T, tbox::Array<double> nD,
      int cell_flag, tbox::Pointer<tbox::Matrix<double> > ele_mat);

  /**
   * @brief 计算单元右端项.
   *
   * @param ele_rhs        输出参数, 数组, 单元右端项.
   * @param integrator     输入参数, 指针, 指向积分器.
   * @param shape_func     输入参数, 指针, 指向形函数.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   */
  virtual void buildElementRHS(tbox::Array<hier::DoubleVector<NDIM> > coord,
                               tbox::Array<double> u_val, const double dt,
                               const double time,
                               tbox::Pointer<tbox::Vector<double> > ele_vec);

  /// 电场右端项 (电流连续性): ∫f·N_i dΩ, 无体源时 f=0
  virtual void buildEElementRHS(tbox::Array<hier::DoubleVector<NDIM> > coord,
                               const double dt, const double time,
                               int cell_flag,
                               tbox::Pointer<tbox::Vector<double> > ele_vec);

  /// 电场计算: E = -∇φ
  virtual void calculateElementEx(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Array<double> phi,
      tbox::Pointer<tbox::Vector<double> > ele_Ex);

  /// 温度场刚度矩阵 (热传导): ∫K(T,nD) ∇N_i·∇N_j dΩ
  virtual void buildTStiffElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Array<double> T, tbox::Array<double> nD,
      int cell_flag, tbox::Pointer<tbox::Matrix<double> > ele_mat);

  /// 温度场右端项: 焦耳热源 ∫(|E|²σ)·N_i dΩ
  virtual void buildTElementRHS(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, double u_val, tbox::Array<double> phi,
      tbox::Array<double> E_cell, tbox::Array<double> T, tbox::Array<double> nD,
      int cell_flag, tbox::Pointer<tbox::Vector<double> > ele_rhs);

private:
  /*!@brief 对象名.  */
  string d_object_name;
};
#endif
