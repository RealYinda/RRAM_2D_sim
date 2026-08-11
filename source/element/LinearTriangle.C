//
// 文件名:     LinearTriangle.C
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     单元计算类实现
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "LinearTriangle.h"
using namespace JAUMIN;
LinearTriangle::LinearTriangle(const string& name) : BaseElement<NDIM>(name) {}
LinearTriangle::~LinearTriangle() {}
/*************************************************************************
 *
 *  电处理
 *
 *
 ************************************************************************/
void LinearTriangle::buildStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat) {
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Triangle");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Triangle");

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, quad_pnt);

  /// 计算单元刚度矩阵.
  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < n_dof; ++j) {
      for (int l = 0; l < num_quad_pnts; ++l) {
        double JxW = volume * jac * weight[l];

        (*ele_mat)(i, j) += JxW * (bas_grad[l][i][0] * bas_grad[l][j][0] +
                                   bas_grad[l][i][1] * bas_grad[l][j][1]);
      }
    }
  }
}


void LinearTriangle::buildElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_rhs) {
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Triangle");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Triangle");

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  /// 计算单元右端项.
  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxW = volume * jac * weight[l];
      //const double PI = 4.0 * atan(1.0);
     // double f_val = 104.0 * PI * PI * sin(2.0 * PI * quad_pnt[l][0]) *
                     //sin(10.0 * PI * quad_pnt[l][1]);
      double f_val=0;
      (*ele_rhs)[i] += JxW * f_val * bas_val[l][i];
    }
  }
}
/*************************************************************************
 *计算单元的Ex
 ************************************************************************/
void LinearTriangle::calculateElementEx(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time,  tbox::Array<double> phi,
    tbox::Pointer<tbox::Vector<double> > ele_Ex)
{
	/// 取出积分器对象.
	  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
	      IntegratorManager<NDIM>::getManager();
	  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
	      integrator_manager->getIntegrator("Triangle");

	  /// 取出形函数对象.
	  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
	      ShapeFunctionManager<NDIM>::getManager();
	  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
	      shape_manager->getShapeFunction("Triangle");

	  /// 取出单元上自由度数目.
	  int n_dof = shape_func->getNumberOfDof();
	  /// 取出积分点数目.
	  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
	  /// 取出积分点.
	  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
	      integrator->getQuadraturePoints(real_vertex);

	  /// 取出基函数在积分点梯度值 (线性单元梯度为常数).
	  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
	      shape_func->gradient(real_vertex, quad_pnt);

	  /// E = -∇φ = -Σ φ_i ∇N_i, 与 3D LinearTetrahedron::calculateElementEx 同构.
	  double Ex = 0, Ey = 0;
	  for (int i = 0; i < n_dof; ++i) {
	    Ex += -phi[i] * bas_grad[0][i][0];
	    Ey += -phi[i] * bas_grad[0][i][1];
	  }
	  (*ele_Ex)[0] = Ex;
	  (*ele_Ex)[1] = Ey;
	};
/*************************************************************************
 *
 *  热处理
 *
 *
 ************************************************************************/
/*************************************************************************
 *  计算热处理矩阵
 ************************************************************************/
void LinearTriangle::thermal_buildStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,double density,double K,double Cp) {
	 /// 取出积分器对象.
	  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
	      IntegratorManager<NDIM>::getManager();
	  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
	      integrator_manager->getIntegrator("Triangle");

	  /// 取出形函数对象.
	  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
	      ShapeFunctionManager<NDIM>::getManager();
	  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
	      shape_manager->getShapeFunction("Triangle");

	  /// 取出单元上自由度数目.
	  int n_dof = shape_func->getNumberOfDof();
	  /// 取出积分点数目.
	  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
	  /// 取出模板单元的面积.
	  double volume = integrator->getElementVolume();
	  /// 取出积分点.
	  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
	      integrator->getQuadraturePoints(real_vertex);
	  /// 取出jacobian矩阵行列式.
	  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
	  /// 取出积分点的积分权重.
	  tbox::Array<double> weight = integrator->getQuadratureWeights();

	  /// 取出基函数在积分点的值和梯度值.
	  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
	      shape_func->gradient(real_vertex, quad_pnt);
	  tbox::Array<tbox::Array<double> > bas_val =
	      shape_func->value(real_vertex, quad_pnt);
	 // cout << dt << endl;
	  /// 计算单元刚度矩阵.
	  for (int i = 0; i < n_dof; ++i) {
	    for (int j = 0; j < n_dof; ++j) {
	      for (int l = 0; l < num_quad_pnts; ++l) {
	        double JxW = volume * jac * weight[l];
	        (*ele_mat)(i, j) += JxW * (bas_val[l][i] * bas_val[l][j] +
	                                   dt * bas_grad[l][i][0] * bas_grad[l][j][0] +
	                                   dt * bas_grad[l][i][1] * bas_grad[l][j][1]);//在这里热传导系数k/(密度*等压比热容)=1
	      }
	     // cout << (*ele_mat)(i, j) << "  "<<i<<"   "<<j<<endl;
	    }
	  }
	}
void LinearTriangle::thermal_buildElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, double u_val, tbox::Array<double> T_val, tbox::Pointer<tbox::Vector<double> > ele_rhs) {
	  /// 取出积分器对象.
	  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
	      IntegratorManager<NDIM>::getManager();
	  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
	      integrator_manager->getIntegrator("Triangle");
	  //cout << u_val<<" sdfsdf" << endl;
	  //cout << T_val[0]<<" sdfsdf" << endl;
	  /// 取出形函数对象.
	  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
	      ShapeFunctionManager<NDIM>::getManager();
	  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
	      shape_manager->getShapeFunction("Triangle");

	  /// 取出单元上自由度数目.
	  int n_dof = shape_func->getNumberOfDof();
	  /// 取出积分点数目.
	  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
	  /// 取出模板单元的面积.
	  double volume = integrator->getElementVolume();
	  /// 取出积分点.
	  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
	      integrator->getQuadraturePoints(real_vertex);
	  /// 取出jacobian矩阵行列式.
	  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
	  /// 取出积分点的积分权重.
	  tbox::Array<double> weight = integrator->getQuadratureWeights();

	  /// 取出基函数在积分点的值和梯度值.
	  tbox::Array<tbox::Array<double> > bas_val =
	      shape_func->value(real_vertex, quad_pnt);
	  double f;
	  /// 计算单元右端项.
	  for (int i = 0; i < n_dof; ++i) {
	    for (int l = 0; l < num_quad_pnts; ++l) {
	      double u_value = 0.0;
	      for (int i1 = 0; i1 < n_dof; ++i1) {
	        double tmp_val = T_val[i1];
	        u_value += tmp_val * bas_val[l][i1];
	      }
	      double JxW = volume * jac * weight[l];
	     // double expt = exp((quad_pnt[l][0] + quad_pnt[l][1]) * (time + dt));
	      double f_val =u_val*dt+u_value;
	      //double f_val=1;
	      f=f_val;
	      (*ele_rhs)[i] += JxW * f_val * bas_val[l][i];
	    }
	    //cout << (*ele_rhs)[i]<<"(*ele_rhs)["<<i<<"]"<<f<<endl;
	  }
}
/*************************************************************************
 *
 *  材料参数处理
 *
 *
 ************************************************************************/
/*************************************************************************
 *计算单元的sigma
 ************************************************************************/
void LinearTriangle::calculateElementSigma(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Array<double> T, tbox::Pointer<double> ele_sigma){

	/// 取出积分器对象.
	  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
	      IntegratorManager<NDIM>::getManager();
	  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
	      integrator_manager->getIntegrator("Triangle");

	  /// 取出形函数对象.
	  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
	      ShapeFunctionManager<NDIM>::getManager();
	  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
	      shape_manager->getShapeFunction("Triangle");

	  /// 取出单元上自由度数目.
	  int n_dof = shape_func->getNumberOfDof();
	  /// 取出积分点数目.
	  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
	  /// 取出模板单元的面积.
//	  double volume = integrator->getElementVolume();
	  /// 取出积分点.
	  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
	      integrator->getQuadraturePoints(real_vertex);
	  /// 取出jacobian矩阵行列式.
//	  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
	  /// 取出积分点的积分权重.
//	  tbox::Array<double> weight = integrator->getQuadratureWeights();

	  /// 取出基函数在积分点梯度值.
	  tbox::Array<tbox::Array<double> > bas_val =
	      shape_func->value(real_vertex, quad_pnt);

	  double ele_T=0;

	  for (int i = 0; i < n_dof; ++i) {
	    for (int l = 0; l < num_quad_pnts; ++l) {
	      ele_T+= T[i]* bas_val[l][i];
	    }
	  }
	  (*ele_sigma)= (*ele_sigma)-0.1*ele_T;
}

void LinearTriangle::buildDDStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_K,
    tbox::Array<double> phi, tbox::Array<double> T,
    tbox::Pointer<tbox::Matrix<double> > ele_mat) {
  (void)time; (void)phi; (void)T;

  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Triangle");
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Triangle");

  int n_dof = shape_func->getNumberOfDof();
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  double volume = integrator->getElementVolume();
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  tbox::Array<double> weight = integrator->getQuadratureWeights();
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  tbox::Pointer<tbox::Matrix<double> > M = new tbox::Matrix<double>();
  M->resize(n_dof, n_dof);
  for (int i = 0; i < n_dof; ++i)
    for (int j = 0; j < n_dof; ++j)
      (*M)(i, j) = 0.0;

  for (int i = 0; i < n_dof; ++i)
    for (int j = 0; j < n_dof; ++j)
      for (int l = 0; l < num_quad_pnts; ++l)
        (*M)(i, j) += volume * jac * weight[l] * bas_val[l][i] * bas_val[l][j];

  for (int i = 0; i < n_dof; ++i)
    for (int j = 0; j < n_dof; ++j)
      (*ele_mat)(i, j) = (*M)(i, j) + (*ele_K)(i, j) * dt;
}

void LinearTriangle::buildDDElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Array<double> nD,
    tbox::Pointer<tbox::Vector<double> > ele_rhs) {
  (void)dt; (void)time;

  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Triangle");
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Triangle");

  int n_dof = shape_func->getNumberOfDof();
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  double volume = integrator->getElementVolume();
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  tbox::Array<double> weight = integrator->getQuadratureWeights();
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  tbox::Pointer<tbox::Matrix<double> > M = new tbox::Matrix<double>();
  M->resize(n_dof, n_dof);
  for (int i = 0; i < n_dof; ++i)
    for (int j = 0; j < n_dof; ++j)
      (*M)(i, j) = 0.0;

  for (int i = 0; i < n_dof; ++i)
    for (int j = 0; j < n_dof; ++j)
      for (int l = 0; l < num_quad_pnts; ++l)
        (*M)(i, j) += volume * jac * weight[l] * bas_val[l][i] * bas_val[l][j];

  for (int i = 0; i < n_dof; ++i) {
    (*ele_rhs)[i] = 0.0;
    for (int j = 0; j < n_dof; ++j)
      (*ele_rhs)[i] += (*M)(i, j) * nD[j];
  }
}

void LinearTriangle::buildEStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Array<double> T, tbox::Array<double> nD,
    int cell_flag, tbox::Pointer<tbox::Matrix<double> > ele_mat) {
  (void)dt; (void)time;
  /// 电场刚度矩阵 (电流连续性方程 ∇·(σ∇V) = 0).
  /// 与 3D LinearTetrahedron 同构, 积分器换为 Triangle.

  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Triangle");

  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Triangle");

  int n_dof = shape_func->getNumberOfDof();
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  double volume = integrator->getElementVolume();
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, quad_pnt);

  for (int i = 0; i < n_dof; ++i)
    for (int j = 0; j < n_dof; ++j)
      (*ele_mat)(i, j) = 0.0;

  double *pnts_T = new double[num_quad_pnts];
  double *pnts_nD = new double[num_quad_pnts];
  double *sigma = new double[num_quad_pnts];
  for (int l = 0; l < num_quad_pnts; ++l) {
    pnts_T[l] = 0;
    pnts_nD[l] = 0;
    sigma[l] = 0;
  }

  for (int i = 0; i < n_dof; ++i)
    for (int l = 0; l < num_quad_pnts; ++l) {
      pnts_T[l] += T[i] * bas_val[l][i];
      pnts_nD[l] += nD[i] * bas_val[l][i];
    }

  for (int l = 0; l < num_quad_pnts; ++l) {
    if (cell_flag == 1 || cell_flag == 2) {
      sigma[l] = Sigma_0(pnts_nD[l]) * exp(-Eac_var(pnts_nD[l]) / (K_b * pnts_T[l]));
    } else if (cell_flag == 3) {
      sigma[l] = Sigmate;
    } else if (cell_flag == 4) {
      sigma[l] = Sdiodeon;
    } else if (cell_flag == 5 || cell_flag == 7) {
      sigma[l] = Sdiodeoff;
    } else if (cell_flag == 6) {
      sigma[l] = Sigmame;
    }
  }

  /// 组装: K_ij = Σ_l JxW * σ * ∇N_i · ∇N_j
  for (int i = 0; i < n_dof; ++i)
    for (int j = 0; j < n_dof; ++j)
      for (int l = 0; l < num_quad_pnts; ++l) {
        double JxW = volume * jac * weight[l];
        for (int k = 0; k < NDIM; ++k)
          (*ele_mat)(i, j) += JxW * sigma[l] *
                              bas_grad[l][i][k] * bas_grad[l][j][k];
      }

  delete[] pnts_T; pnts_T = NULL;
  delete[] pnts_nD; pnts_nD = NULL;
  delete[] sigma; sigma = NULL;
}

void LinearTriangle::buildEElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, int cell_flag, tbox::Pointer<tbox::Vector<double> > ele_rhs) {
  (void)dt; (void)time; (void)cell_flag;
  /// 电场右端项组装: ∫f·N_i dΩ, 电流连续性方程无体源 (f=0), 框架保留.
  /// 与 3D LinearTetrahedron::buildVElementRHS 同构, 积分器换为 Triangle.

  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Triangle");

  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Triangle");

  int n_dof = shape_func->getNumberOfDof();
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  double volume = integrator->getElementVolume();
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxW = volume * jac * weight[l];
      double f_val = 0;  // 体源项 (电流连续性无源)
      (*ele_rhs)[i] += JxW * f_val * bas_val[l][i];
    }
  }
}

void LinearTriangle::buildTStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Array<double> T, tbox::Array<double> nD,
    int cell_flag, tbox::Pointer<tbox::Matrix<double> > ele_mat) {
  (void)dt; (void)time;
  /// 热传导刚度矩阵: ∫K(T,nD) ∇N_i·∇N_j dΩ.
  /// 与 3D 同构, 积分器换为 Triangle.

  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Triangle");

  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Triangle");

  int n_dof = shape_func->getNumberOfDof();
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  double volume = integrator->getElementVolume();
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, quad_pnt);
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  double *pnts_nD = new double[num_quad_pnts];
  double *K = new double[num_quad_pnts];
  for (int l = 0; l < num_quad_pnts; ++l) {
    pnts_nD[l] = 0;
    K[l] = 0;
  }
  for (int i = 0; i < n_dof; ++i)
    for (int l = 0; l < num_quad_pnts; ++l)
      pnts_nD[l] += nD[i] * bas_val[l][i];

  for (int l = 0; l < num_quad_pnts; ++l) {
    if (cell_flag == 1 || cell_flag == 2) {
      K[l] = Kth_var(pnts_nD[l]);
    } else if (cell_flag == 3) {
      K[l] = Kte;
    } else if (cell_flag == 4) {
      K[l] = Kdiodeon;
    } else if (cell_flag == 5 || cell_flag == 7) {
      K[l] = Kdiodeoff;
    } else if (cell_flag == 6) {
      K[l] = Kme;
    }
  }

  for (int i = 0; i < n_dof; ++i)
    for (int j = 0; j < n_dof; ++j)
      for (int l = 0; l < num_quad_pnts; ++l) {
        double JxW = volume * jac * weight[l];
        for (int k = 0; k < NDIM; ++k)
          (*ele_mat)(i, j) += JxW * K[l] *
                              bas_grad[l][i][k] * bas_grad[l][j][k];
      }

  delete[] pnts_nD; pnts_nD = NULL;
  delete[] K; K = NULL;
}

void LinearTriangle::buildTElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, double u_val, tbox::Array<double> phi,
    tbox::Array<double> E_cell, tbox::Array<double> T, tbox::Array<double> nD,
    int cell_flag, tbox::Pointer<tbox::Vector<double> > ele_rhs) {
  (void)dt; (void)time;
  /// 温度场右端项: 焦耳热源 ∫(|E|²σ)·N_i dΩ.
  /// 与 3D 同构, 积分器换为 Triangle.

  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Triangle");

  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Triangle");

  int n_dof = shape_func->getNumberOfDof();
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  double volume = integrator->getElementVolume();
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  double *pnts_T = new double[num_quad_pnts];
  double *pnts_nD = new double[num_quad_pnts];
  double *sigma = new double[num_quad_pnts];
  for (int l = 0; l < num_quad_pnts; ++l) {
    pnts_T[l] = 0;
    pnts_nD[l] = 0;
    sigma[l] = 0;
  }
  for (int i = 0; i < n_dof; ++i)
    for (int l = 0; l < num_quad_pnts; ++l) {
      pnts_T[l] += T[i] * bas_val[l][i];
      pnts_nD[l] += nD[i] * bas_val[l][i];
    }

  for (int l = 0; l < num_quad_pnts; ++l) {
    if (cell_flag == 1 || cell_flag == 2) {
      sigma[l] = Sigma_0(pnts_nD[l]) * exp(-Eac_var(pnts_nD[l]) / (K_b * pnts_T[l]));
    } else if (cell_flag == 3) {
      sigma[l] = Sigmate;
    } else if (cell_flag == 4) {
      sigma[l] = Sdiodeon;
    } else if (cell_flag == 5 || cell_flag == 7) {
      sigma[l] = Sdiodeoff;
    } else if (cell_flag == 6) {
      sigma[l] = Sigmame;
    }
  }

  /// 焦耳热源
  double f_val = 0;
  for (int k = 0; k < NDIM; ++k)
    f_val += E_cell[k] * E_cell[k];
  f_val *= sigma[0];

  for (int i = 0; i < n_dof; ++i)
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxW = volume * jac * weight[l];
      (*ele_rhs)[i] += JxW * f_val * bas_val[l][i];
    }

  delete[] pnts_T; pnts_T = NULL;
  delete[] pnts_nD; pnts_nD = NULL;
  delete[] sigma; sigma = NULL;
}
