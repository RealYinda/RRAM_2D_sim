//
// 文件名:     LinearTetrahedron.C
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     单元计算类实现
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "LinearTetrahedron.h"

LinearTetrahedron::LinearTetrahedron(const string& name)
    : BaseElement<NDIM>(name) {}

LinearTetrahedron::~LinearTetrahedron() {}

void LinearTetrahedron::buildStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double time,
    const double dt, tbox::Pointer<tbox::Matrix<double> > ele_mat) {
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Tetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Tetrahedron");

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

  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < n_dof; ++j) {
      (*ele_mat)(i, j) = 0.0;
    }
  }

  /// 计算单元刚度矩阵.
  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < n_dof; ++j) {
      for (int l = 0; l < num_quad_pnts; ++l) {
        double JxW = volume * jac * weight[l];

        (*ele_mat)(i, j) += JxW * (bas_val[l][i] * bas_val[l][j] +
                                   dt * bas_grad[l][i][0] * bas_grad[l][j][0] +
                                   dt * bas_grad[l][i][1] * bas_grad[l][j][1] +
                                   dt * bas_grad[l][i][2] * bas_grad[l][j][2]);
      }
    }
  }
}

void LinearTetrahedron::buildElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex,
    tbox::Array<double> u_val, const double time, const double dt,
    tbox::Pointer<tbox::Vector<double> > ele_rhs) {
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Tetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Tetrahedron");

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
      double u_value = 0.0;
      for (int i1 = 0; i1 < n_dof; ++i1) {
        u_value += u_val[i1] * bas_val[l][i1];
      }

      double JxW = volume * jac * weight[l];
      double val =
          (quad_pnt[l][0] + quad_pnt[l][1] + quad_pnt[l][2]) * (time + dt);
      double expt = exp(val);

      double f_val = dt * ((quad_pnt[l][0] + quad_pnt[l][1] + quad_pnt[l][2]) -
                           3.0 * (time + dt) * (time + dt)) *
                         expt +
                     u_value;
      (*ele_rhs)[i] += JxW * f_val * bas_val[l][i];
    }
  }
}

void LinearTetrahedron::buildEStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Array<double> T, tbox::Array<double> nD,
    int cell_flag, tbox::Pointer<tbox::Matrix<double> > ele_mat) {
  (void)dt; (void)time;
  /// 电场刚度矩阵 (电流连续性方程 ∇·(σ∇V) = 0).
  /// 电导率 σ(T,nD) 按材料区域 cell_flag 取常数或依赖温度/浓度.

  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Tetrahedron");

  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Tetrahedron");

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

  // 高斯积分点的温度 / 浓度 / 电导率
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
        (*ele_mat)(i, j) += JxW * sigma[l] *
                            (bas_grad[l][i][0] * bas_grad[l][j][0] +
                             bas_grad[l][i][1] * bas_grad[l][j][1] +
                             bas_grad[l][i][2] * bas_grad[l][j][2]);
      }

  delete[] pnts_T; pnts_T = NULL;
  delete[] pnts_nD; pnts_nD = NULL;
  delete[] sigma; sigma = NULL;
}

void LinearTetrahedron::buildEElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, int cell_flag, tbox::Pointer<tbox::Vector<double> > ele_rhs) {
  (void)dt; (void)time; (void)cell_flag;
  /// 电场右端项组装: ∫f·N_i dΩ, 电流连续性方程无体源 (f=0), 框架保留.
  /// 与 3D LinearTetrahedron::buildVElementRHS 同构.

  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Tetrahedron");

  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Tetrahedron");

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

void LinearTetrahedron::buildTStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Array<double> T, tbox::Array<double> nD,
    int cell_flag, tbox::Pointer<tbox::Matrix<double> > ele_mat) {
  (void)dt; (void)time;
  /// 热传导刚度矩阵: ∫K(T,nD) ∇N_i·∇N_j dΩ.
  /// 与 3D LinearTetrahedron::thermal_buildStiffElementMatrix 同构.

  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Tetrahedron");

  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Tetrahedron");

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

  // 高斯积分点的浓度 / 热导率
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
        (*ele_mat)(i, j) += JxW * K[l] *
                            (bas_grad[l][i][0] * bas_grad[l][j][0] +
                             bas_grad[l][i][1] * bas_grad[l][j][1] +
                             bas_grad[l][i][2] * bas_grad[l][j][2]);
      }

  delete[] pnts_nD; pnts_nD = NULL;
  delete[] K; K = NULL;
}

void LinearTetrahedron::buildTElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, double u_val, tbox::Array<double> phi,
    tbox::Array<double> E_cell, tbox::Array<double> T, tbox::Array<double> nD,
    int cell_flag, tbox::Pointer<tbox::Vector<double> > ele_rhs) {
  (void)dt; (void)time;
  /// 温度场右端项: 焦耳热源 ∫(|E|²σ)·N_i dΩ.
  /// 与 3D LinearTetrahedron::thermal_buildElementRHS 同构.

  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Tetrahedron");

  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Tetrahedron");

  int n_dof = shape_func->getNumberOfDof();
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  double volume = integrator->getElementVolume();
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  // 高斯积分点的温度 / 浓度 / 电导率
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
  double f_val = (E_cell[0] * E_cell[0] + E_cell[1] * E_cell[1] +
                  E_cell[2] * E_cell[2]) * sigma[0];

  for (int i = 0; i < n_dof; ++i)
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxW = volume * jac * weight[l];
      (*ele_rhs)[i] += JxW * f_val * bas_val[l][i];
    }

  delete[] pnts_T; pnts_T = NULL;
  delete[] pnts_nD; pnts_nD = NULL;
  delete[] sigma; sigma = NULL;
}

void LinearTetrahedron::calculateElementEx(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Array<double> phi,
    tbox::Pointer<tbox::Vector<double> > ele_Ex) {
  (void)dt; (void)time;
  /// 电场计算: E = -∇φ = -Σ φ_i ∇N_i.
  /// 与 3D LinearTetrahedron::calculateElementEx 同构.

  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("Tetrahedron");

  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("Tetrahedron");

  int n_dof = shape_func->getNumberOfDof();
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);

  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, quad_pnt);

  double Ex = 0, Ey = 0, Ez = 0;
  for (int i = 0; i < n_dof; ++i) {
    Ex += -phi[i] * bas_grad[0][i][0];
    Ey += -phi[i] * bas_grad[0][i][1];
    Ez += -phi[i] * bas_grad[0][i][2];
  }
  (*ele_Ex)[0] = Ex;
  (*ele_Ex)[1] = Ey;
  (*ele_Ex)[2] = Ez;
}
