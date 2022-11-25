#include <iostream>
#include <cstddef>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Eigen/Dense"
#include "Eigen/Sparse"

using std::setprecision;
using std::setw;
using std::min;
using std::max;
using std::min;
using std::floor;
using std::sqrt;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::ArrayXi;

int main()
{
    size_t nelx = 80;
    size_t nely = 40;
    double  volfrac = 0.5,
            penal = 3,
            rmin =  2;
    MatrixXd x = MatrixXd::Constant(nely, nelx, volfrac);
    size_t loop = 0;
    double change = 1.0;
    while (change > 0.01){
        ++loop;
        MatrixXd x_old = x;
        // ELEMENT STIFFNESS MATRIX
        double E = 1.0;
        double nu = 0.3;
        // 这里向量的维度之所以设置为 9 是因为下面 KE 的那一大块是直接复制过来的
        // 但是因为 VectorXd 的下标是从 0 开始的，所以需要在 k(0) 处占用一个位置
        VectorXd k {{0, 0.5-nu/6, 0.125+nu/8, -0.25-nu/12, -0.125+3*nu/8, -0.25+nu/12, -0.125-nu/8, nu/6, 0.125-3*nu/8}};
        MatrixXd KE {
            {k(1), k(2), k(3), k(4), k(5), k(6), k(7), k(8)}
            {k(2), k(1), k(8), k(7), k(6), k(5), k(4), k(3)}
            {k(3), k(8), k(1), k(6), k(7), k(4), k(5), k(2)}
            {k(4), k(7), k(6), k(1), k(8), k(3), k(2), k(5)}
            {k(5), k(6), k(7), k(8), k(1), k(2), k(3), k(4)}
            {k(6), k(5), k(4), k(3), k(2), k(1), k(8), k(7)}
            {k(7), k(4), k(5), k(2), k(3), k(8), k(1), k(6)}
            {k(8), k(3), k(2), k(5), k(4), k(7), k(6), k(1)}
        };
        KE *= E/(1-nu*nu);
        // FE-ANALYSIS
        MatrixXd K = MatrixXd::Zero(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
        VectorXd F = VectorXd::Zero(2*(nelx+1)*(nely+1)),
                 U = VectorXd::Zero(2*(nely+1)*(nelx+1));
        for (size_t elx = 1; elx <= nelx; ++elx) {
            for (size_t ely = 1; ely <= nely; ++ely) {
                int n1 = (nely+1)*(elx-1) + ely;
                int n2 = (nely+1)*elx + ely;
                Eigen::Array<size_t, 8, 1> edof {2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2};
                // MATLAB 里的矩阵和向量是以 1 为起始下标的，而 Eigen 里的矩阵都是以 0 为起始下标的
                // 所以这里的 K 和 x 的索引相比 MATLAB 里的都要少 1
                edof -= 1;
                size_t ely_index = ely - 1;
                size_t elx_index = elx - 1;
                K(edof, edof) += pow(x(ely_index, elx_index),penalization_exponent)*KE;
            }
        }
        // DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
        // 索引相比 MATLAB 都要少 1
        F(2*(nely+1)*nelx+nely+1,0) = -1;
        // 用作 U 的索引
        std::vector<size_t> fixed_dofs;
        for (size_t i = 0; i < 2*(nely+1); i++)
        {
            fixed_dofs.push_back(i);
        }
        std::vector<size_t> all_dofs;
        for (size_t i = 0; i < 2*(nely+1)*(nelx+1); ++i)
        {
            all_dofs.push_back(i)
        }
        // 用作 U, K, F 的索引
        std::vector<size_t> free_dofs;
        for (size_t value = 2*(nely+1)+1; value < 2*(nely+1)*(nelx+1); value++)
        {
            free_dofs.push_back(value);
        }
        U(free_dofs) = K(free_dofs, free_dofs).colPivHouseholderQr().solve(F(free_dofs));
        U(fixed_dofs) = 0;

        double c = 0;
        MatrixXd dc = MatrixXd::Zero(nely, nelx);
        for (size_t elx = 1; elx <= nelx; ++elx) {
            for (size_t ely = 1; ely <= nely; ++ely) {
                size_t n1 = (nely+1)*(elx-1) + ely;
                size_t n2 = (nely+1)*elx + ely;
                Eigen::Array<size_t, 8, 1> indexes {2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2};
                // MATLAB 里的矩阵和向量是以 1 为起始下标的，而 Eigen 里的矩阵都是以 0 为起始下标的
                // 所以这里的 U 的索引相比 MATLAB 里的都要少 1
                indexes -= 1;
                VectorXi Ue(8) = U(indexes);
                size_t  ely_index = ely-1,
                        elx_index = elx-1;
                c += pow(x(index_ely,index_elx), penal)*Ue.transpose()*KE*Ue;
                dc(ely_index, elx_index) = -penal*pow(x(ely_index, elx_index), (penal-1))*Ue.transpose()*KE*Ue
            }
        }

        // FILTERING OF SENSITIVITIES
        MatrixXd dcn = MatrixXd::Zero(nely, nelx);
        for (size_t i = 1; i <= nelx; i++)
        {
            for (size_t j = 1; j <= nely; j++)
            {
                double sum = 0.0;
                size_t  j_index = j-1,
                        i_index = i-1;
                for (size_t k = max(i-floor(rmin), 1); k <= min(i+floor(rmin), nelx); k++)
                {
                    for (size_t l = max(j-floor(rmin), 1); l <= min(j+floor(rmin), nely); l++)
                    {
                        double fac = rmin-sqrt((i-k)*(i-k)+(j-l)*(j-l));
                        sum += max(0, fac);
                        size_t  l_index = l-1,
                                k_index = k-1;
                        dcn(j_index, i_index) += max(0, fac)*x(l_index, k_index)*dc(l_index, k_index);
                    }
                }
                dcn(j_index, i_index) /= x(j_index, i_index) * sum;
            }
        }
        
        // DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
        double  l1 = 0,
                l2 = 100000,
                move = 0.2;
        MatrixXd xnew(nely, nelx);
        while (l2-l1 > 1e-4)
        {
            double lmid = 0.5*(l2+l1);
            xnew = max(0.001, max(x-move, min(1.0, min(x+move, x.cwiseProduct(sqrt(-dc/lmid))))));
            if (xnew.sum() - volfrac*nelx*nely > 0){
                l1 = lmid;
            }else {
                l2 = lmid;
            }
        }
        x = xnew;
        // TODO(whoslcy@foxmail.com): 下面这个能行得通吗 cwiseAbs() 的返回值是一个 矩阵吗
        change = (x - x_old).cwiseAbs().maxCoeff();
        std::cout
            << "It: " << setw(4) << loop << " "
            << "Obj.: " << setw(10) << setprecision(6) << c << " "  // 原 matlab 代码中打印精度为小数点后4位，我目前只发现了 C++ 提供了有效数字的精度指定方式，于是指定了 6 位有效数字，达到和原 matlab 代码一样的打印效果
            << setprecision(3)
            << "Vol.: " << setw(6) << x.sum()/(nelx*nely)
            << "ch.: " << setw(6) << change
            << std::endl;
        // TODO(whoslcy@foxmail.com): 输出灰度图
    };
    return 0;
}