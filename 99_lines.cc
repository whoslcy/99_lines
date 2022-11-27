#include <iostream>
#include <cstddef>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Eigen/Eigen"

using std::setprecision;
using std::setw;
using std::floor;
using std::size_t;
using std::min;
using std::max;
using std::cout;
using std::endl;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::ArrayXi;

typedef Eigen::Triplet<double> T;

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
            {k(1), k(2), k(3), k(4), k(5), k(6), k(7), k(8)},
            {k(2), k(1), k(8), k(7), k(6), k(5), k(4), k(3)},
            {k(3), k(8), k(1), k(6), k(7), k(4), k(5), k(2)},
            {k(4), k(7), k(6), k(1), k(8), k(3), k(2), k(5)},
            {k(5), k(6), k(7), k(8), k(1), k(2), k(3), k(4)},
            {k(6), k(5), k(4), k(3), k(2), k(1), k(8), k(7)},
            {k(7), k(4), k(5), k(2), k(3), k(8), k(1), k(6)},
            {k(8), k(3), k(2), k(5), k(4), k(7), k(6), k(1)}
        };
        KE *= E/(1-nu*nu);
        
        // FE-ANALYSIS
        // https://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title3
        std::vector<T> triplet_list;                                    
        for (size_t elx = 1; elx <= nelx; ++elx) {
            for (size_t ely = 1; ely <= nely; ++ely) {
                size_t n1 = (nely+1)*(elx-1) + ely;
                size_t n2 = (nely+1)*elx + ely;
                Eigen::Array<size_t, 8, 1> edof {{2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2}};
                // MATLAB 里的矩阵和向量是以 1 为起始下标的，而 Eigen 里的矩阵都是以 0 为起始下标的
                // 所以这里的 K 和 x 的索引相比 MATLAB 里的都要少 1
                edof -= 1;
                size_t ely_index = ely - 1;
                size_t elx_index = elx - 1;
                for (size_t i = 0; i < 8; i++)
                {
                    for (size_t j = 0; j < 8; j++)
                    {
                        triplet_list.push_back(T(edof(i), edof(j), pow(x(ely_index, elx_index), penal)*KE(i, j)));
                    }
                }
            }
        }
        Eigen::SparseMatrix<double> K(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
        K.setFromTriplets(triplet_list.begin(), triplet_list.end());
        // DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
        Eigen::SparseVector<double> F(2*(nelx+1)*(nely+1)-2*(nely+1));
        F.insert(2*(nelx+1)*(nely+1)-2*(nely+1)-41) = -1;
        // 用作 U 的索引
        std::vector<size_t> fixed_dofs;
        for (size_t i = 0; i < 2*(nely+1); i++)
        {
            fixed_dofs.push_back(i);
        }
        std::vector<size_t> free_dofs;
        for (size_t i = 2*(nely+1); i < 2*(nely+1)*(nelx+1); i++)
        {
            free_dofs.push_back(i);
        }
        Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(K.bottomRightCorner(free_dofs.size(), free_dofs.size()));
        VectorXd U(2*(nelx+1)*(nely+1));
        U(free_dofs) = chol.solve(F);
        U(fixed_dofs).array() = 0;

        double c = 0;
        MatrixXd dc = MatrixXd::Zero(nely, nelx);
        for (size_t elx = 1; elx <= nelx; ++elx) {
            for (size_t ely = 1; ely <= nely; ++ely) {
                size_t n1 = (nely+1)*(elx-1) + ely;
                size_t n2 = (nely+1)*elx + ely;
                Eigen::Array<size_t, 8, 1> indexes {{2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2}};
                // MATLAB 里的矩阵和向量是以 1 为起始下标的，而 Eigen 里的矩阵都是以 0 为起始下标的
                // 所以这里的 U 的索引相比 MATLAB 里的都要少 1
                indexes -= 1;
                VectorXd Ue = U(indexes);
                size_t  ely_index = ely-1,
                        elx_index = elx-1;
                double temp = Ue.transpose()*KE*Ue;
                c += pow(x(ely_index, elx_index), penal)*temp;
                dc(ely_index, elx_index) = -penal*pow(x(ely_index, elx_index), (penal-1))*temp;
            }
        }

        // FILTERING OF SENSITIVITIES
        MatrixXd dc_new = MatrixXd::Zero(nely, nelx);
        for (size_t i = 1; i <= nelx; i++)
        {
            for (size_t j = 1; j <= nely; j++)
            {
                double sum = 0.0;
                size_t  j_index = j-1,
                        i_index = i-1;
                for (size_t k = max(i-floor(rmin), 1.0); k <= min(i+floor(rmin), (double)nelx); k++)
                {
                    for (size_t l = max(j-floor(rmin), 1.0); l <= min(j+floor(rmin), (double)nely); l++)
                    {
                        double fac = rmin-std::sqrt((i-k)*(i-k)+(j-l)*(j-l));
                        sum += max(0.0, fac);
                        size_t  l_index = l-1,
                                k_index = k-1;
                        dc_new(j_index, i_index) += max(0.0, fac)*x(l_index, k_index)*dc(l_index, k_index);
                    }
                }
                dc_new(j_index, i_index) /= x(j_index, i_index) * sum;
            }
        }
        dc = dc_new;

        // DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
        double  l1 = 0,
                l2 = 100000,
                move = 0.2;
        MatrixXd x_new(nely, nelx);
        while (l2-l1 > 1e-4) {
            double lmid = 0.5*(l2+l1);
            x_new = x.cwiseProduct((-dc/lmid).cwiseSqrt()).cwiseMin(
                            (x.array()+move).matrix()).cwiseMin(
                                                  1.0).cwiseMax(
                            (x.array()-move).matrix()).cwiseMax(0.001);
            if (x_new.sum() - volfrac*nelx*nely > 0){
                l1 = lmid;}
            else {
                l2 = lmid;}}
        x = x_new;
        change = (x - x_old).cwiseAbs().maxCoeff();
        std::cout
            << " It.: " << setw(4) << loop
            << " Obj.: " << setw(10) << setprecision(6) << c  // 原 matlab 代码中打印精度为小数点后4位，我目前只发现了 C++ 提供了有效数字的精度指定方式，于是指定了 6 位有效数字，达到和原 matlab 代码一样的打印效果
            << " Vol.: " << setw(6) << setprecision(3) << x.sum()/(nelx*nely)
            << " ch.: " << setw(6) << setprecision(2) << change
            << std::endl;
    };
    return 0;
}