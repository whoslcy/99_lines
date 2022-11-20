#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Eigen/Dense"
#include "Eigen/Sparse"

using std::setprecision;
using std::setw;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

int main()
{
    int nelx = 80;
    int nely = 40;
    double volfrac = 0.5;
    double penal = 3;
    double rmin =  2;
    // TODO(whoslcy@foxmail.com): 和 Matrix<double, nely, nelx> 相比哪个更好
    MatrixXd x(nely, nelx);
    // MatrixXd 默认是列优先存储的，按照列的顺序赋值会更快一些
    // TODO(whoslcy@foxmail.com): 引用矩阵的元素要求用什么类型的值？
    for (size_t i = 0; i < nelx; i++)
    {
        for (size_t j = 0; j < nely; j++)
        {
            x(j, i) = volfrac;
        }
    }
    // TODO(whoslcy@foxmail.com): size_t 跟 string::size_type 是一回事吗
    size_t loop = 0;
    double change = 1.0;
    while (change > 0.01){
        ++loop;
        MatrixXd xold(nely, nelx) = x;

        // FE-ANALYSIS

        // ELEMENT STIFFNESS MATRIX
        double E = 1.0;
        double nu = 0.3;
        VectorXd k(8) = {0/2-nu/6, 1/8+nu/8, -1/4-nu/12, -1/8+3*nu/8, -1/4+nu/12, -1/8-nu/8, nu/6, 1/8-3*nu/8};
        MatrixXd KE(8, 8) = E/(1-nu*nu) * {
            {k(1), k(2), k(3), k(4), k(5), k(6), k(7), k(8)},
            {k(2), k(1), k(8), k(7), k(6), k(5), k(4), k(3)},
            {k(3), k(8), k(1), k(6), k(7), k(4), k(5), k(2)},
            {k(4), k(7), k(6), k(1), k(8), k(3), k(2), k(5)},
            {k(5), k(6), k(7), k(8), k(1), k(2), k(3), k(4)},
            {k(6), k(5), k(4), k(3), k(2), k(1), k(8), k(7)},
            {k(7), k(4), k(5), k(2), k(3), k(8), k(1), k(6)},
            {k(8), k(3), k(2), k(5), k(4), k(7), k(6), k(1)},
        };

        MatrixXd K(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1)),
                 F(2*(nelx+1)*(nely+1)),
                 U(2*(nely+1)*(nelx+1));
        for (int elx = 1; elx <= nelx; ++elx)
        {
            for (int ely = 1; ely <= nely; ++ely)
            {
                int n1 = (nely+1)*(elx-1) + ely;
                int n2 = (nely+1)*elx + ely;
                std::vector<int> edof {2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2};
                // MATLAB 里的矩阵和向量是以 1 为起始下标的，而 Eigen 里的矩阵都是以 0 为起始下标的，所以这里的 edof 相比 MATLAB 里的都要少 1
                for (auto &value : edof)
                {
                    value--;
                }
                K(edof, edof) += pow(x(ely, elx),penalization_exponent)*KE;
            }
        }
        // DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
        F(2*(nely+1)*nelx+nely+2,1) = -1;
        // 用作 U 的索引
        std::vector<int> fixed_dofs;
        for (int i = 0; i < 2*(nely+1); i++)
        {
            fixed_dofs.push_back(i);
        }
        std::vector<int> all_dofs;
        for (int i = 0; i < 2*(nely+1)*(nelx+1); ++i)
        {
            all_dofs.push_back(i)
        }
        // 用作 U, K, F 的索引
        VectorXi free_dofs(2*(nely+1)*nelx);
        for (int i = 0; i < 2*(nely+1)*nelx; i++)
        {
            free_dofs(i) = 2*(nely+1)+1+i;
        }
        U(free_dofs) = K(free_dofs, free_dofs).colPivHouseholderQr().solve(F(free_dofs));
        U(fixed_dofs) = 0;
        double c = 0;
        for (int elx = 1; elx <= nelx; ++elx)
        {
            for (int ely = 1; ely <= nely; ++ely)
            {
                int n1 = (nely+1)*(elx-1) + ely;
                int n2 = (nely+1)*elx + ely;
                std::vector<int> indexes {2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2};
                VectorXi Ue(8) = U(indexes);
                // MATLAB 里的矩阵和向量是以 1 为起始下标的，而 Eigen 里的矩阵都是以 0 为起始下标的，所以这里的 edof 相比 temp 里的都要少 1
                c = c + Ue.transpose()*KE*Ue*pow(x(ely-1, elx-1), penal);

            }
        }
        // TODO(whoslcy@foxmail.com): 下面这个能行得通吗 cwiseAbs() 的返回值是一个 矩阵吗
        change = (x - x_old).cwiseAbs().maxCoeff();
        std::cout
            << "It: " << setw(4) << loop << " "
            << "Obj.: " << setw(10) << setprecision(6) << c << " "  // 原 matlab 代码中打印精度为小数点后4位，我目前只发现了 C++ 提供了有效数字的精度指定方式，于是指定了 6 位有效数字，达到和原 matlab 代码一样的打印效果
            << setprecision(3)
            << "Vol.: " << setw(6) << x.sum()/(nelx*nely)  // TODO(whoslcy@foxmail.com): sum() 的返回值的类型是什么，是 double 吗，能直接除吗
            << "ch.: " << setw(6) << change
            << std::endl;
        // TODO(whoslcy@foxmail.com): 输出灰度图
    };
    return 0;
}