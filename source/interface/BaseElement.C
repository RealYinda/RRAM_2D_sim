
//
// 文件名:     BaseElement.C
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     有限元单元基类
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "JAUMIN_config.h"
#include "Utilities.h"
#include "BaseIntegrator.h"
#include "DoubleVector.h"
#include "BaseShapeFunction.h"
#include "BaseElement.h"

#include <string>
#include <list>

#ifndef NULL
#define NULL (0)
#endif

/************************************************************************
 * 构造函数
 *************************************************************************/
template <int DIM>
BaseElement<DIM>::BaseElement(const string& name) {
  d_object_name = name;
}

/************************************************************************
 * 析构函数
 *************************************************************************/
template <int DIM>
BaseElement<DIM>::~BaseElement() {}

/************************************************************************
 * 建立单元刚度矩阵
 *************************************************************************/
//电
template <int DIM>
void BaseElement<DIM>::buildStiffElementMatrix(
    tbox::Array<hier::DoubleVector<DIM> > ele_info, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,double sigma) {
  TBOX_WARNING(
      "Function buidStiffElementMatrix is not implemented, please make sure "
      "your code is right! ");
}
template <int DIM>
 void BaseElement<DIM>::calculateElementEx(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Array<double> phi, tbox::Pointer<tbox::Vector<double> > ele_Ex){
	  TBOX_WARNING(
	      "Function LinearTriangle::calculateElementEx is not implemented, please make sure "
	      "your code is right! ");
	}
 template <int DIM>
 void BaseElement<DIM>::buildEStiffElementMatrix(
         tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
         const double time,  tbox::Array<double> T,tbox::Array<double> nD,int cell_flag,tbox::Pointer<tbox::Matrix<double> > ele_mat) {
   TBOX_WARNING(
       "Function buidStiffElementMatrix is not implemented, please make sure "
       "your code is right! ");
 }
//漂移扩散
 template <int DIM>
 void BaseElement<DIM>::buildDDStiffElementMatrix(
     tbox::Array<hier::DoubleVector<DIM> > ele_info, const double dt,
     const double time,tbox::Pointer<tbox::Matrix<double> > ele_K,tbox::Array<double> phi,tbox::Array<double> T, tbox::Pointer<tbox::Matrix<double> > ele_mat) {
   TBOX_WARNING(
       "Function buidStiffElementMatrix is not implemented, please make sure "
       "your code is right! ");
 }
//热
template <int DIM>
void BaseElement<DIM>::buildTStiffElementMatrix(
        tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
        const double time,  tbox::Array<double> T, tbox::Array<double> nD,int cell_flag,tbox::Pointer<tbox::Matrix<double> > ele_mat) {
  TBOX_WARNING(
      "Function buildTStiffElementMatrix is not implemented, please make sure "
      "your code is right! ");
}
template <int DIM>
void BaseElement<DIM>::buildTElementRHS(
        tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
        const double time, double u_val,tbox::Array<double> phi,tbox::Array<double> E_cell,
        tbox::Array<double> T,tbox::Array<double> nD,int cell_flag,tbox::Pointer<tbox::Vector<double> > ele_rhs) {
  TBOX_WARNING(
      "Function buildTElementRHS is not implemented, please make sure "
      "your code is right! ");
}
template <int DIM>
void BaseElement<DIM>::thermal_buildStiffElementMatrix(
        tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
        const double time,  tbox::Array<double> T, tbox::Array<double> nD,int cell_flag,tbox::Pointer<tbox::Matrix<double> > ele_mat) {
  TBOX_WARNING(
      "Function buidStiffElementMatrix is not implemented, please make sure "
      "your code is right! ");
}
template <int DIM>
void BaseElement<DIM>::thermal_buildStiffElementMatrix(
        tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
        const double time, double u_val,tbox::Array<double> phi,tbox::Array<double> T,tbox::Array<double> nD,int cell_flag,tbox::Pointer<tbox::Vector<double> > ele_rhs) {
  TBOX_WARNING(
      "Function buidStiffElementMatrix is not implemented, please make sure "
      "your code is right! ");
}

/************************************************************************
 * 建立单元质量矩阵
 *************************************************************************/
template <int DIM>
void BaseElement<DIM>::buildMassElementMatrix(
    tbox::Array<hier::DoubleVector<DIM> > ele_info, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat) {
  TBOX_WARNING(
      "Function buidMassElementMatrix is not implemented, please make sure "
      "your code is right! ");
}

/************************************************************************
 * 建立单元阻尼矩阵
 *************************************************************************/
template <int DIM>
void BaseElement<DIM>::buildDumpElementMatrix(
    tbox::Array<hier::DoubleVector<DIM> > ele_info, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat) {
  TBOX_WARNING(
      "Function buidDumpElementMatrix is not implemented, please make sure "
      "your code is right! ");
}

/************************************************************************
 * 建立单元右端项
 *************************************************************************/
template <int DIM>
void BaseElement<DIM>::buildElementRHS(
    tbox::Array<hier::DoubleVector<DIM> > ele_info, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec) {
  TBOX_WARNING(
      "Function buidElementRHS is not implemented, please make sure your code "
      "is right! ");
}
template <int DIM>
void BaseElement<DIM>::buildEElementRHS(
    tbox::Array<hier::DoubleVector<DIM> > ele_info, const double dt,
    const double time, int cell_flag,tbox::Pointer<tbox::Vector<double> > ele_vec) {
  TBOX_WARNING(
      "Function buidElementRHS is not implemented, please make sure your code "
      "is right! ");
}
template <int DIM>
void BaseElement<DIM>::thermal_buildElementRHS(
        tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
        const double time, double u_val,tbox::Array<double> phi,tbox::Array<double> E_cell,
        tbox::Array<double> T,tbox::Array<double> nD,
        int cell_flag,tbox::Pointer<tbox::Vector<double> > ele_rhs) {
  TBOX_WARNING(
      "Function buidElementRHS is not implemented, please make sure your code "
      "is right! ");

}
template <int DIM>
void BaseElement<DIM>::buildDDElementRHS(
            tbox::Array<hier::DoubleVector<NDIM> > real_vertex,
            const double dt,
            const double time,
            tbox::Array<double> nD,
            tbox::Pointer<tbox::Vector<double> > ele_rhs) {
  TBOX_WARNING(
      "Function buidElementRHS is not implemented, please make sure your code "
      "is right! ");
}
/************************************************************************
 *
 *
 *
 * 计算单元的材料参数
 *
 *
 *
 *
 *************************************************************************/
/*****************************************************************************
 * 计算电导率参数.
 *****************************************************************************/
template <int DIM>
void BaseElement<DIM>::calculateElementSigma(
        tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
        const double time, tbox::Array<double> T, tbox::Array<double> nD,double* ele_sigma){
	  TBOX_WARNING(
	      "Function calculateElementSigma is not implemented, please make sure your code "
	      "is right! ");
	}
template <int DIM>
void BaseElement<DIM>::calculateElementK(
        tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
        const double time, tbox::Array<double> T, tbox::Array<double> nD,double*ele_K){
	  TBOX_WARNING(
	      "Function calculateElementSigma is not implemented, please make sure your code "
	      "is right! ");
	}

/*****************************************************************************
 * 判断给定名字是否以本对象名字匹配.
 *****************************************************************************/
template <int DIM>
bool BaseElement<DIM>::isMatching(const string& name) {
  return (d_object_name == name);
}

/*****************************************************************************
 * 获取对象名字.
 *****************************************************************************/
template <int DIM>
const string& BaseElement<DIM>::getName() {
  return d_object_name;
}

template <int DIM>
const int BaseElement<DIM>::getProblemDim() {
  TBOX_WARNING(
      "Function getProblemDim is not implemented, please make sure your code "
      "is right! ");
  return 0;
}

template class BaseElement<NDIM>;
