//
// 文件名:      TetQuad.h
// 软件包:
// 版权  :      (c) 2004-2015 北京应用物理与计算数学研究所
//              (c) 2013-2015 中物院高性能数值模拟软件中心
// 版本号:      $Revision$
// 修改  :      $Date$
// 描述  :
//

#ifndef included_appu_TriQuad
#define included_appu_TriQuad

#include <algorithm>
#include <strings.h>
#include "Pointer.h"
#include "Array.h"
#include "Patch.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "EdgeData.h"
#include "CellData.h"
#include "NodeData.h"
#include "GridInfo.h"

namespace JAUMIN {
namespace appu {
class TriQuad
{
public:
    TriQuad(const hier::Patch<NDIM>& patch,
            tbox::Pointer<pdat::CellData<NDIM, double> > cell_volume,
        tbox::Pointer<pdat::CellData<NDIM, double> > cell_jacobian);
    template<class TPYE>
    void nodeBasMultinodeBas(int order,double area,double *value)
    {
        const Quad *quad=d_face_quad_table[order];
        double *phi;
        int nbas=3;
        double(*cc)[nbas]=(double(*)[nbas])value;//基函数Ni*Nj在一个单元的积分值
        int i_order[3]={0,1,2};
        if(area<0)
            {
            i_order[2]=1;
            i_order[1]=2;
        }
        //初始化
        bzero(&(cc[0][0]), nbas * nbas * sizeof(*value));
        for(int n = 0;n < quad->npoints;n++)
        {
            phi=quad->points+n*3;
            for(int i=0;i<nbas;i++)
            {
                for(int j=0;j<nbas;j++)
                {
                    double v=0;
                    v=phi[i_order[i]]*phi[i_order[j]];
                    v*=quad->weights[n];
                    cc[i][j]+=v;
                }
            }

        }
        for(int i=0;i<nbas;i++)
        {
            for (int j=0; j<nbas; j++)
                {
                cc[i][j]*=abs(area);
            }
        }
    }
    template<class TPYE>
    void nodegradBasDotnodegradBas(int order,double area,tbox::Array<tbox::Array<double> >nodegradBas,double *value)
    {
        const Quad *quad=d_face_quad_table[order];
        //double *phi;
        int nbas=3;
        double(*cc)[nbas]=(double(*)[nbas])value;//基函数gradNi*gradNj在一个单元的积分值
        ///初始化
         bzero(&(cc[0][0]), nbas * nbas * sizeof(*value));
//         cout<<cc[0][0]<<" "<<cc[0][1]<<"// "
//             <<cc[1][0]<<" "<<cc[1][1]<<"// "
//               <<cc[2][0]<<" "<<cc[2][1]<<"// " <<endl;
//                 cout<<nodegradBas[0][0]<<" "<<nodegradBas[0][1]<<"// "
//                     <<nodegradBas[1][0]<<" "<<nodegradBas[1][1]<<"// "
//                       <<nodegradBas[2][0]<<" "<<nodegradBas[2][1]<<"// " <<endl;
         int i_order[3]={0,1,2};
         if(area<0)
             {
             i_order[2]=1;
             i_order[1]=2;
         }
        for(int n = 0;n < quad->npoints;n++)
        {
            //phi=quad->points+n*3;
            for(int i=0;i<nbas;i++)
            {
                for(int j=0;j<nbas;j++)
                {
                    double v=0;
                    for(int d=0;d<2;d++)
                    {
                        v+=nodegradBas[i_order[i]][d]*nodegradBas[i_order[j]][d];
                    }
                    v*=quad->weights[n];
                //    cout<<quad->weights[n]<<endl;
                    cc[i][j]+=v;
                }
            }

        }
//        cout<<nodegradBas[0][0]<<" "<<nodegradBas[0][1]<<"// "
//            <<nodegradBas[1][0]<<" "<<nodegradBas[1][1]<<"// "
//              <<nodegradBas[2][0]<<" "<<nodegradBas[2][1]<<"// " <<endl;
//        cout<<cc[0][0]<<" "<<cc[0][1]<<"// "
//            <<cc[1][0]<<" "<<cc[1][1]<<"// "
//              <<cc[2][0]<<" "<<cc[2][1]<<"// " <<endl;
        for(int i=0;i<nbas;i++)
        {
            for (int j=0; j<nbas; j++)
                {
                cc[i][j]*=abs(area);
            }
        }
    }
  ///////////////////////////////////////边基函数的积分/////////////////////////////////////////
    template <class ShapeFunc>
    void quadedgeBasDotedgeBas(const int face,const int order,double area, const ShapeFunc *shapefunc,
                        double *value,tbox::Array<tbox::Array<double> >nodeGradBas) const {
      check(order);
      int nbas = shapefunc->nbas();
      int dim = shapefunc->dim();
      double(*bb)[nbas] = (double(*)[nbas])value;
      const Quad *quad = d_face_quad_table[order];
      ///初始化bb
       bzero(&(bb[0][0]), nbas * nbas * sizeof(*value));
       int faceorder=0;//表示面的节点是否是逆时针编号
       if(area>0)
       {
           faceorder=1;
       }

      for (int n = 0; n < quad->npoints; n++) {
        double phi[nbas][dim];
        shapefunc->basis(face, quad->points + n * 3, &(phi[0][0]),nodeGradBas,faceorder);
//        cout<<phi[0][0]<<" "<<phi[0][1]<<" || "
//            <<phi[1][0]<<" "<<phi[1][1]<<" || "
//            <<phi[2][0]<<" "<<phi[2][1]<<" || "
//                                         <<face<<endl;
        for (int i = 0; i < nbas; i++) {
          for (int j = 0; j < nbas; j++) {
            double v = 0.;
            for (int d = 0; d < dim; d++) v += phi[i][d]*phi[j][d];
            v *= quad->weights[n];
            bb[i][j] += v;
          }
        }
      }
      for (int i = 0; i < nbas; i++)
      {
            for (int j = 0; j < nbas; j++)
            {
                bb[i][j] *= abs(area);
            }
       }
    }
    template <class ShapeFunc>
    void quadcurledgeBasDotcurledgeBas(const int face,const int order,double area, const ShapeFunc *shapefunc,
                        double *value,tbox::Array<tbox::Array<double> >nodeGradBas) const {
      check(order);
      int nbas = shapefunc->nbas();
      int dim = 1;
      double(*bb)[nbas] = (double(*)[nbas])value;
      const Quad *quad = d_face_quad_table[order];
      ///初始化bb
       bzero(&(bb[0][0]), nbas * nbas * sizeof(*value));
       int faceorder=0;//表示面的节点是否是逆时针编号
       if(area>0)
       {
           faceorder=1;
       }

      for (int n = 0; n < quad->npoints; n++) {
        double curlphi[nbas][1];
        shapefunc->curl(face, quad->points + n * 3, &(curlphi[0][0]),nodeGradBas,faceorder);
        for (int i = 0; i < nbas; i++) {
          for (int j = 0; j < nbas; j++) {
            double v = 0.;
            for (int d = 0; d < dim; d++) v += curlphi[i][d]*curlphi[j][d];
            v *= quad->weights[n];
            bb[i][j] += v;
          }
        }
      }
      for (int i = 0; i < nbas; i++)
      {
            for (int j = 0; j < nbas; j++)
            {
                bb[i][j] *= abs(area);
            }
       }
    }
    ///////////////////////////////边基函数和节点基函数/////////////////////////////////
    template <class ShapeFunc>
    void quadedgeBasDotnodegraBas(const int face,const int order,double area, const ShapeFunc *shapefunc,
                        double *value,tbox::Array<tbox::Array<double> >nodeGradBas) const {
      check(order);
      int nbas = shapefunc->nbas();
      int dim = shapefunc->dim();
      double(*bb)[nbas] = (double(*)[nbas])value;
      const Quad *quad = d_face_quad_table[order];
      ///初始化bb
       bzero(&(bb[0][0]), nbas * nbas * sizeof(*value));
       int faceorder=0;//表示面的节点是否是逆时针编号
       int i_order[3]={0,2,1};
       if(area>0)
       {
           faceorder=1;
           i_order[2]=2;
           i_order[1]=1;
       }

      for (int n = 0; n < quad->npoints; n++) {
        double phi[nbas][dim];
        shapefunc->basis(face, quad->points + n * 3, &(phi[0][0]),nodeGradBas,faceorder);
        for (int i = 0; i < nbas; i++) {
          for (int j = 0; j < nbas; j++) {
            double v = 0.;
            for (int d = 0; d < dim; d++) v += phi[i][d]*nodeGradBas[i_order[j]][d];
            v *= quad->weights[n];
            bb[i][j] += v;
          }
        }
      }
      for (int i = 0; i < nbas; i++)
      {
            for (int j = 0; j < nbas; j++)
            {
                bb[i][j] *= abs(area);
            }
       }
    }

    template <class ShapeFunc>
    void quadnodegradBasDotedgeBas(const int face,const int order,double area, const ShapeFunc *shapefunc,
                        double *value,tbox::Array<tbox::Array<double> >nodeGradBas) const {
      check(order);
      int nbas = shapefunc->nbas();
      int dim = shapefunc->dim();
      double(*bb)[nbas] = (double(*)[nbas])value;
      const Quad *quad = d_face_quad_table[order];
      ///初始化bb
       bzero(&(bb[0][0]), nbas * nbas * sizeof(*value));
       int faceorder=0;//表示面的节点是否是逆时针编号
       int i_order[3]={0,2,1};
       if(area>0)
       {
           faceorder=1;
           i_order[2]=2;
           i_order[1]=1;
       }

      for (int n = 0; n < quad->npoints; n++) {
        double phi[nbas][dim];
        shapefunc->basis(face, quad->points + n * 3, &(phi[0][0]),nodeGradBas,faceorder);
        for (int i = 0; i < nbas; i++) {
          for (int j = 0; j < nbas; j++) {
            double v = 0.;
            for (int d = 0; d < dim; d++) v += nodeGradBas[i_order[i]][d]*phi[j][d];
            v *= quad->weights[n];
            bb[i][j] += v;
          }
        }
      }
      for (int i = 0; i < nbas; i++)
      {
            for (int j = 0; j < nbas; j++)
            {
                bb[i][j] *= abs(area);
            }
       }
    }
/////////////////////////////////////////模式功率（E叉乘H）/////////////////////////////////
    template <class ShapeFunc,class TYPE, class TYPE2,class TYPE3>
    void quadfaceEsquare(const int face,const int order,double area,
                           const ShapeFunc *shapefunc,
                           const TYPE *E_edge,////边上插值系数
                           const TYPE *E_node,///节点上插值系数
                           TYPE2 *value,
                           tbox::Array<tbox::Array<double> >nodeGradBas,
                           TYPE3 beta,
                           double omega,
                           double mu) const {
         check(order);
         int nbas = shapefunc->nbas();
         int dim = shapefunc->dim();
        // TYPE(*bb) = (double(*)[1])value;
         TYPE2 (*bb)=value;
             bzero(&(bb[0]), 1* sizeof(*value));
         const Quad *quad = d_face_quad_table[order];
         ///初始化bb
              bzero(&(bb[0]), 1* sizeof(*value));
          int faceorder=0;//表示面的节点是否是逆时针编号
          int i_order[3]={0,2,1};
          if(area>0)
          {
              faceorder=1;
              i_order[2]=2;
              i_order[1]=1;
          }

         for (int n = 0; n < quad->npoints; n++) {
           double phi[nbas][dim];
           shapefunc->basis(face, quad->points + n * 3, &(phi[0][0]),nodeGradBas,faceorder);

           TYPE et[2]={0,0};//横向场//et=beta*Et
           for(int i=0;i<3;i++){//一个面上有三个基函数
               et[0]+=E_edge[i]*phi[i][0];
               et[1]+=E_edge[i]*phi[i][1];
           }
//           cout<<"phi "<<phi[0][0]<<" "<<phi[0][1]<<endl;
//           cout<<"phi "<<phi[1][0]<<" "<<phi[1][1]<<endl;
//           cout<<"phi "<<phi[2][0]<<" "<<phi[2][1]<<endl;
//           cout<<"Et "<<et[0]<<" "<<et[1]<<"   E_edge "<<E_edge[0]<<" "<<E_edge[1]<<" "<<E_edge[2]<<endl;
           TYPE gradez[2]={0,0};//纵向场ez=-j*Ez
           for(int i=0; i<3; i++){
              gradez[0]+=E_node[i]*nodeGradBas[i_order[i]][0];///节点2，3在area为负值时交换了顺序，故求出梯度顺序是0，2，1不同于面单元所规定的顺序0，1，2
              gradez[1]+=E_node[i]*nodeGradBas[i_order[i]][1];
           }
//           cout<<E_node[0]<<" "<<E_node[1]<<E_node[2]<<endl;
            TYPE2 v = 0.;
            TYPE2 Ht[2]={0,0};//Hx=(dEz/dy+j*beta*Ety)/(j*omega*mu),Hy=(-dEz/dx-j*beta*Etx)/(j*omega*mu)
            Ht[0]=-(gradez[1]*dcomplex(0,1)+dcomplex(0,1)*et[1])/(dcomplex(0,1)*omega*mu);
            Ht[1]=-(-gradez[0]*dcomplex(0,1)-dcomplex(0,1)*et[0])/(dcomplex(0,1)*omega*mu);
//            cout<<Ht[0]<<" "<<Ht[1]<<endl;
            TYPE2 Et[2]={0,0};
            Et[0]=et[0]/beta;
            Et[1]=et[1]/beta;
            //cout<<"Et "<<Et[0]<<" "<<Et[1]<<endl;
            TYPE2 S_flow=0;//E叉乘H，入射功率
            S_flow=Et[0]*Ht[1]-Et[1]*Ht[0];
           // cout<<"yita"<<Et[0]/Ht[1]<<"  "<<-Et[1]/Ht[0]<<endl;
            v += S_flow;
           // cout<<"Enorm "<<sqrt(v)<<endl;
            v *= quad->weights[n];
            (*bb)=(*bb)+ v;
         }
         (*bb)=(*bb)*abs(area);
       }

    /// 2D 三角形单元 Kij 计算.
    /// 与 3D TetQuad::quadcalculateKij 同构: 调 shapefunc->basis() + 质心单点积分.
    /// edge_order 方向处理由 shapefunc->basis() 内部完成.
    template<class ShapeFunc, class TYPE>
    void quadcalculateKij(const int cell,
                           const ShapeFunc* shapefunc,
                           const int order,
                           const int localedge,
                           const int ii,
                           double edgeJ,
                           TYPE *value) const
    {
        (void)order;
        double area = (*d_cell_volume)(0, cell);
        double(*nabla)[NDIM + 1] =
            (double(*)[NDIM + 1])(&((*d_cell_jacobian)(0, cell)));

        bzero(value, 1 * sizeof(*value));

        // 质心 barycentric 坐标 (order 0, 线性基下精确)
        double lambda[NDIM + 1];
        for (int k = 0; k <= NDIM; k++) lambda[k] = 1.0 / (NDIM + 1);

        int nbas = shapefunc->nbas();
        int dim  = shapefunc->dim();
        double phi[nbas][dim];
        shapefunc->basis(cell, lambda, &(phi[0][0]));

        double v = edgeJ * dotProduct(NDIM, phi[localedge], nabla[ii]);
        *value = v * area;
    }

public:
    struct Quad {
    const char *name;
    int dim;
    int order;
    int npoints;
    double *points;
    double *weights;
    int id;
    };

private:
    void check(int order) const {
    const int max_order = 3;
    if(order > max_order)
        TBOX_ERROR("Quadrature rules of order "
               << order << " for tetrahedra not implememted.\n");
    }

    const hier::Patch<NDIM>& d_patch;
    tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_volume;
    tbox::Pointer<pdat::CellData<NDIM, double> > d_cell_jacobian;

    tbox::Array<int> can_ext, can_idx;
#if (NDIM > 2)
    tbox::Array<int> fan_ext, fan_idx;
    tbox::Array<int> caf_ext, caf_idx;
#endif
    tbox::Array<int> ean_ext, ean_idx;
    tbox::Array<int> cae_ext, cae_idx;

    Quad** d_face_quad_table;
    Quad** d_edge_quad_table;
};
} // namespace appu
} // namespace JAUMIN

#endif // included_appu_TetQuad
