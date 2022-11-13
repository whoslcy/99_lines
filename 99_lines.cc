#include <iostream>
#include <iomanip>
#include "Eigen/Dense"
#include "Eigen/Sparse"

using std::setprecision;
using std::setw;

using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::SparseMatrix;
using Eigen::SparseVector;
using Eigen::VectorXi;

typedef Matrix<double, 8, 8> Matrix8d;

void top(int x_elements_number, int y_elements_number, double volume_fraction, int penalization_exponent, int minium_radius, double convergence_criterion){
    // TODO(whoslcy@foxmail.com): 和 Matrix<double, y_elements_number, x_elements_number> 相比哪个更好
    // TODO(whoslcy@foxmail.com): 设计变量 x 重命名
    MatrixXd x(y_elements_number, x_elements_number);
    // TODO(whoslcy@foxmail.com): size_t 跟 string::size_type 是一回事吗
    size_t loop_counter = 0;
    double difference_between_two_iteration = 1.0;
    while (difference_between_two_iteration > convergence_criterion){
        ++loop_counter;
        // TODO(whoslcy@foxmail.com): x_old 重命名
        MatrixXd x_old = x;
        ?? U = finite_element_analyze(x_elements_number, y_elements_number, x, penalization_exponent);
        Matrix8d KE = lk();
        int c = 0;
        for (size_t i = 0; i < count; i++)
        {
            /* code */
        }
        // TODO(whoslcy@foxmail.com): 下面这个能行得通吗 cwiseAbs() 的返回值是一个 矩阵吗
        difference_between_two_iteration = (x - x_old).cwiseAbs().maxCoeff();
        std::cout
            << "It: " << setw(4) << loop_counter << " "
            << "Obj.: " << setw(10) << setprecision(6) << c << " "  // 原 matlab 代码中打印精度为小数点后4位，我目前只发现了 C++ 提供了有效数字的精度指定方式，于是指定了 6 位有效数字，达到和原 matlab 代码一样的打印效果
            << setprecision(3)
            << "Vol.: " << setw(6) << x.sum()/(x_elements_number*y_elements_number)  // TODO(whoslcy@foxmail.com): sum() 的返回值的类型是什么，是 double 吗，能直接除吗
            << "ch.: " << setw(6) << difference_between_two_iteration
            << std::endl;
        // TODO(whoslcy@foxmail.com): 输出灰度图
    }
}


// TODO(whoslcy@foxmail.com): 返回值类型是稀疏的吗？
VectorXi finite_element_analyze(int x_elements_number, int y_elements_number, MatrixXd &x, int penalization_exponent){
    Matrix8d KE = lk();
    // TODO(whoslcy@foxmail.com): 想一个常量名字，把 2*(y_elements_number+1)*(x_elements_number+1) 存起来
    SparseMatrix<double, RowMajor> K(2*(x_elements_number+1)*(y_elements_number+1), 2*(x_elements_number+1)*(y_elements_number+1));
    SparseVector<double>    F(2*(x_elements_number+1)*(y_elements_number+1)),
                            U(2*(y_elements_number+1)*(x_elements_number+1));
    // dof means degree of freedom.
    for (int element_x = 1; element_x <= x_elements_number; ++element_x)
    {
        for (int element_y = 1; element_y <= y_elements_number; ++element_y)
        {
            int node1 = (y_elements_number+1)*(element_x-1) + element_y;
            int node2 = (y_elements_number+1)*element_x + element_y;
            // MATLAB 里的矩阵和向量是以 1 为起始下标的，而 Eigen 里的矩阵都是以 0 为起始下标的，所以这里的 element_dof 相比 MATLAB 里的都要少 1
            std::vector<int> element_dof {2*node1-2, 2*node1-1, 2*node2-2, 2*node2-1, 2*node2, 2*node2+1, 2*node1, 2*node1+1};
            K(element_dof, element_dof) += pow(x(element_y, element_x),penalization_exponent)*KE;
        }
        
    }
    F(2,1) = -1;
    VectorXi fixed_dofs(y_elements_number+2);
    for (int i = 0; i < y_elements_number + 1; i++)
    {
        fixed_dofs(i) = 2*i+1;
    }
    fixed_dofs(y_elements_number+1) = 2*(nelx+1)*(y_elements_number+1+);
    VectorXi all_dofs(2*(y_elements_number+1)*(nelx+1));
    for (int i = 0; i < 2*(y_elements_number+1)*(nelx+1); ++i){
        all_dofs(i) = i+1;
    }
    // TODO(whoslcy@foxmail.com): free_dofs uncompleted
    VectorXi free_dofs = 
}


// 返回单元刚度矩阵
// TODO(whoslcy@foxmail.com): lk 改名
Matrix8d lk(double E, double nu){
    // TODO(whoslcy@foxmail.com): 普通 10 位精度 的 double 的范围会不会不够用，我不知道接下来 E 和 nu 会到什么程度
    // 8 维行向量
    Matrix<double, 1, 8> k {0/2-nu/6, 1/8+nu/8, -1/4-nu/12, -1/8+3*nu/8, -1/4+nu/12, -1/8-nu/8, nu/6, 1/8-3*nu/8};
    // 8 阶方阵
    // TODO(whoslcy@foxmail.com): 我目前不知道这个矩阵为什么这么安排，我只知道这是个对称矩阵，这样纯手打矩阵效率很低，容易出错，修改也不方便，有没有更优雅的方式
    // TODO(whoslcy@foxmail.com): 这个初始化方式 OK 吗？
    // TODO(whoslcy@foxmail.com): KE 改名
    Matrix8d KE = E/(1-nu*nu) *  {
        {k(1), k(2), k(3), k(4), k(5), k(6), k(7), k(8)},
        {k(2), k(1), k(8), k(7), k(6), k(5), k(4), k(3)},
        {k(3), k(8), k(1), k(6), k(7), k(4), k(5), k(2)},
        {k(4), k(7), k(6), k(1), k(8), k(3), k(2), k(5)},
        {k(5), k(6), k(7), k(8), k(1), k(2), k(3), k(4)},
        {k(6), k(5), k(4), k(3), k(2), k(1), k(8), k(7)},
        {k(7), k(4), k(5), k(2), k(3), k(8), k(1), k(6)},
        {k(8), k(3), k(2), k(5), k(4), k(7), k(6), k(1)},
    };
    return KE;
}