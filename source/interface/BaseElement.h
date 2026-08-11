//
// 文件名:     BaseElement.h
// 软件包:     2D FEM application
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     有限元单元类基类
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_BaseElement
#define included_fem_BaseElement

#include <vector>
#include "Array.h"
#include "Pointer.h"
#include "Matrix.h"
#include "Vector.h"

using namespace std;
using namespace JAUMIN;

namespace JAUMIN {

namespace hier {
template <int DIM>
class DoubleVector;
}
}
template <int DIM>
class BaseIntegrator;
template <int DIM>
class BaseShapeFunction;

template <int DIM>
class BaseElement {
public:
  /**
   * @brief 构造函数.
   *
   * @param name 输入参数, 字符串, 对象名称.
   */
  BaseElement(const string& name);

  /**
   * @brief 析构函数.
   */
  virtual ~BaseElement();

  /**
   * @brief 计算单元刚度矩阵.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_mat        输出参数, 指针, 指向单元刚度矩阵.
   *
   */
  virtual void buildStiffElementMatrix(
      tbox::Array<hier::DoubleVector<DIM> > ele_info, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,double sigma);
  virtual void buildEStiffElementMatrix(
          tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
          const double time,  tbox::Array<double> T,tbox::Array<double> nD,int cell_flag,tbox::Pointer<tbox::Matrix<double> > ele_mat);
  virtual void buildTStiffElementMatrix(
          tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
          const double time,  tbox::Array<double> T, tbox::Array<double> nD,int cell_flag,tbox::Pointer<tbox::Matrix<double> > ele_mat);
  virtual void buildTElementRHS(
          tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
          const double time, double u_val,tbox::Array<double> phi,tbox::Array<double> E_cell,
          tbox::Array<double> T,tbox::Array<double> nD,int cell_flag,tbox::Pointer<tbox::Vector<double> > ele_rhs);
  virtual void thermal_buildStiffElementMatrix(
          tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
          const double time,  tbox::Array<double> T, tbox::Array<double> nD,int cell_flag,tbox::Pointer<tbox::Matrix<double> > ele_mat);
  virtual void thermal_buildStiffElementMatrix(
          tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
          const double time, double u_val,tbox::Array<double> phi,tbox::Array<double> T,tbox::Array<double> nD,int cell_flag,tbox::Pointer<tbox::Vector<double> > ele_rhs);
//zgd 20171225
  virtual void buildDDStiffElementMatrix(
      tbox::Array<hier::DoubleVector<DIM> > ele_info, const double dt,
      const double time,tbox::Pointer<tbox::Matrix<double> > ele_K,tbox::Array<double> phi,tbox::Array<double> T, tbox::Pointer<tbox::Matrix<double> > ele_mat);


  virtual void calculateElementEx(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Array<double> phi, tbox::Pointer<tbox::Vector<double> > ele_Ex);
  virtual void calculateElementSigma(
          tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
          const double time, tbox::Array<double> T, tbox::Array<double> nD,double* ele_sigma);
  virtual void calculateElementK(
          tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
          const double time, tbox::Array<double> T, tbox::Array<double> nD,double *ele_K);
  /**
   * @brief 计算单元质量矩阵.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_mat        输出参数, 指针, 指向单元质量矩阵.
   */
  virtual void buildMassElementMatrix(
      tbox::Array<hier::DoubleVector<DIM> > ele_info, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat);

  /**
   * @brief 计算单元阻尼矩阵.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_mat        输出参数, 指针, 指向单元阻尼矩阵.
   */
  virtual void buildDumpElementMatrix(
      tbox::Array<hier::DoubleVector<DIM> > ele_info, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat);

  /**
   * @brief 计算单元右端项.
   *
   * @param ele_info       输入参数, 指针, 指向单元信息对象.
   * @param dt             输入参数, 双精度, 时间步长.
   * @param time           输入参数, 双精度, 当前时刻.
   * @param ele_vec        输出参数, 指针, 指向单元右端项.
   */
  virtual void buildElementRHS(tbox::Array<hier::DoubleVector<DIM> > ele_info,
                               const double dt, const double time,
                               tbox::Pointer<tbox::Vector<double> > ele_vec);
  virtual void buildEElementRHS(tbox::Array<hier::DoubleVector<DIM> > ele_info,
                               const double dt, const double time,int cell_flag,
                               tbox::Pointer<tbox::Vector<double> > ele_vec);
  virtual void thermal_buildElementRHS(tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
                                        const double time, double u_val,tbox::Array<double> phi,tbox::Array<double> E_cell,tbox::Array<double> T,tbox::Array<double> nD,int cell_flag,tbox::Pointer<tbox::Vector<double> > ele_rhs);
  virtual void buildDDElementRHS(tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
                                  const double time,tbox::Array<double> nD, tbox::Pointer<tbox::Vector<double> > ele_rhs);


  /**
   * @brief 判断给定的名字是否为该对象的名字.
   *
   * @param name 输入参数, 字符串, 约束名字.
   *
   * @return 布尔型, 是否一致.
   */
  virtual bool isMatching(const string& name);

  /**
   * @brief 获取对象的名字.
   *
   * @return 字符串, 对象的名字.
   */
  virtual const string& getName();

  /**
   * @brief 获取单元求解问题维数.
   *
   * @return 整型, 单元求解问题的维数.
   */
  virtual const int getProblemDim();

private:
  string d_object_name;
};
#endif
