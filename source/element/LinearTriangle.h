//
// 文件名:     LinearTriangle.h
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     单元计算类
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "BaseElement.h"
#include "BaseIntegrator.h"
#include "BaseShapeFunction.h"
#include "DoubleVector.h"
#include "Vector.h"
#include "Matrix.h"
#include "UserConstant.h"

#include "IntegratorManager.h"
#include "ShapeFunctionManager.h"
#include "ElementManager.h"

using namespace JAUMIN;
class LinearTriangle : public BaseElement<NDIM> {
public:
  /**
   * @brief 构造函数.
   *
   * @param name 输入参数, 单元名字.
   *
   */
  LinearTriangle(const string& name);

  /**
   * @brief 析构函数.
   *
   */
  ~LinearTriangle();

  /**
   * @brief 计算单元刚度矩阵.
   *
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_mat        输出参数, 指针, 指向矩阵.
   */
  virtual void buildStiffElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > coord, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat);
  /// 电场刚度矩阵 (电流连续性): ∫σ(T,nD) ∇N_i·∇N_j dΩ
  virtual void buildEStiffElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Array<double> T, tbox::Array<double> nD,
      int cell_flag, tbox::Pointer<tbox::Matrix<double> > ele_mat);
  virtual void thermal_buildStiffElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > coord, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,
	   double density,double K,double Cp);


  ///////////////////////////////////////////
  virtual void calculateElementEx(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Array<double> phi,
      tbox::Pointer<tbox::Vector<double> > ele_Ex);
  /////////////////////////////////////////////////////////////
  virtual void calculateElementSigma(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Array<double> T, tbox::Pointer<double> ele_Ex);
  //////////////////////////////////////////////////////////////////////////////////
  /**
   * @brief 计算单元右端项.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_vec        输出参数, 指针, 指向单元向量.
   */
  virtual void buildElementRHS(tbox::Array<hier::DoubleVector<NDIM> > coord,
                               const double dt, const double time,
                               tbox::Pointer<tbox::Vector<double> > ele_vec);
  virtual void thermal_buildElementRHS(tbox::Array<hier::DoubleVector<NDIM> > coord,
                               const double dt, const double time,
                               double u_val,
							   tbox::Array<double> T_val,
							   tbox::Pointer<tbox::Vector<double> > ele_vec);

  /// 漂移扩散矩阵: ele_mat = M + ele_K * dt
  virtual void buildDDStiffElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_K,
      tbox::Array<double> phi, tbox::Array<double> T,
      tbox::Pointer<tbox::Matrix<double> > ele_mat);

  /// 漂移扩散右端项: ele_rhs = M * nD
  virtual void buildDDElementRHS(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Array<double> nD,
      tbox::Pointer<tbox::Vector<double> > ele_rhs);

  /// 电场右端项 (电流连续性): ∫f·N_i dΩ, 无体源时 f=0
  virtual void buildEElementRHS(tbox::Array<hier::DoubleVector<NDIM> > coord,
                               const double dt, const double time,
                               int cell_flag,
                               tbox::Pointer<tbox::Vector<double> > ele_vec);

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
};
