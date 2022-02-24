#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <iterator>
#include "matrix.h"
using namespace std;

constexpr int N = 2;
constexpr double initial_x = 1.0;
constexpr double initial_y = 1.0;
constexpr double threshold = 0.001;

struct grad{
    double fx;
    double fy;
};

inline auto func(double x, double y){
    return x * x + y * y + x * y - 2 * x;
}

inline auto fx(double x, double y){
    return 2 * x + y - 2;
}

inline auto fy(double x, double y){
    return 2 * y + x;
}

inline auto calc_grad(double x, double y)-> grad{
    return {2 * x + y - 2, 2 * y + x};
}

inline auto calc_grad_norm(grad grad){
    return sqrt(grad.fx * grad.fx + grad.fy * grad.fy);
}

inline auto Armijo(double x, double y, vector<double> d){
    double alpha = 1.0; //ステップ幅の初期値
    double raito = 0.4; //raito:係数αの減衰比率
    double coeff = 0.001; //c1 傾きの減少量

    for(int i = 0; i < 50; i++){
        if(func(x + d.at(0) * alpha, y + d.at(1) * alpha) <= func(x, y) + coeff * alpha * (fx(x, y) * d.at(0)+ fy(x, y) * d.at(1))) break;
        else alpha = alpha * raito;
    }
    return alpha;
}


inline auto Bupdate(vector<vector<double>> &B, vector<vector<double>> s, vector<vector<double>> dgrad){
    vector<vector<double>> st = make_tenchi(s); //sの転置行列
    vector<vector<double>> dgradt = make_tenchi(dgrad); //yの転置行列

//Bの更新式右辺第二項を計算
    vector<vector<double>> st_B = calc_matrix_product(st, B);
    vector<vector<double>> st_B_s = calc_matrix_product(st_B, s); //Bの更新式の右辺第二項の分母stBsを計算
    vector<vector<double>> Bs = calc_matrix_product(B, s);
    vector<vector<double>> Bs_t = make_tenchi(calc_matrix_product(B, s));
    vector<vector<double>> Bs_Bst = calc_matrix_product(Bs, Bs_t);//Bの更新式の右辺第二項の分子Bs*(Bs)tを計算
    if(fabs(st_B_s.at(0).at(0)) < 0.00001){
        st_B_s.at(0).at(0) = 0.00001;
    }
    vector<vector<double>> second_term = matrix_times_const(1 / st_B_s.at(0).at(0), Bs_Bst); //second term(右辺第二項)の結果を計算

//Bの更新式右辺第三項を計算
    vector<vector<double>> st_y = calc_matrix_product(st, dgrad); //Bの更新式の右辺第三項の分母を計算
    vector<vector<double>> y_yt = calc_matrix_product(dgrad, dgradt); //Bの更新式の右辺第三項の分子を計算
    if(fabs(st_y.at(0).at(0)) < 0.00001){
        st_y.at(0).at(0) = 0.00001;
    }
    vector<vector<double>> third_term = matrix_times_const(1 / st_y.at(0).at(0), y_yt); //third term(右辺第三項)の結果を計算

    vector<vector<double>> differ = calc_matrix_difference(B, second_term);
    B = calc_matrix_sum(differ, third_term);
}

//行列をコピーする関数
inline auto copy (vector<vector<double>> A, vector<vector<double>> &B){
    for(unsigned int i = 0; i < A.size(); i++){
        for(unsigned int j = 0; j < A.at(0).size(); j++){
            B.at(i).at(j) = A.at(i).at(j);
        }
    }
}

//行列がすでに対角行列になっているか判定(これがないとなぜかガウスの消去法でバグる)
inline auto is_diagonal_matrix(int N, vector<vector<double>> A){
    bool flag = false;
    for(int i = 0; i < N - 1; i++){ //列を変更するループ
        for(int j = 1; j < N; j++){ //行を変更するループ
            if(A.at(j).at(i) != 0){
                flag = true;
                break;
            }
        }
    }
    return flag;
}

int main(int ac, char *av[]){
    double x = initial_x;
    double y = initial_y;
    double x_new, y_new;
    double alpha;
    vector<vector<double>> B = {
        {1.0, 0.0},
        {0.0, 1.0}
    }; //ヘッセ行列の近似
    vector<vector<double>> B_for_gauss = {
        {1.0, 0.0},
        {0.0, 1.0}
    }; //ガウスの掃き出し法用のヘッセ行列の近似
    vector<double> gradient(N);
    vector<double> d(N);
    vector<vector<double>> dgrad(N, vector<double>(1)); // = yn
    vector<vector<double>> s(N, vector<double>(1));
    grad grad;
    bool loop_flag = true;

    ofstream fout(av[1]);
    if(!fout){
        cout << "Can't open the file" << endl;
        exit(8);
    }
    fout << "initial x = " << x << " initial y = " << y << " f(x, y) = " << func(x, y) << endl;

    while(loop_flag){
        grad = calc_grad(x, y);
        if(calc_grad_norm(grad) < threshold){
            loop_flag = false;
        }

        gradient.at(0) = -grad.fx;
        gradient.at(1) = -grad.fy;
        if(is_diagonal_matrix(N, B)){
            gauss(N, B_for_gauss, gradient); //Bd = -∇f(x)を計算して，探索方向dを計算
        }
        calc_ans(N, B_for_gauss, gradient, d);//Bd = -∇f(x)を計算して，探索方向dを計算

        alpha = Armijo(x, y, d); //アルミホ条件でステップ幅を更新
        x_new = x + alpha * d.at(0); //xをアップデート
        y_new = y + alpha * d.at(1); //yをアップデート
        if(isnan(x_new) || isnan(y_new)) break;
    
//行列Bをアップデート
        s.at(0).at(0) = x_new - x;
        s.at(1).at(0) = y_new - y;
        dgrad.at(0).at(0) = fx(x_new, y_new) - fx(x, y);//ynのx成分
        dgrad.at(1).at(0) = fy(x_new, y_new) - fy(x, y);//ynのy成分
        Bupdate(B, s, dgrad);
// 行列Bをアップデートここまで

        x = x_new;
        y = y_new;
        copy(B, B_for_gauss); //ガウスの掃き出し法のための行列に行列Bをコピー

        fout << " x = " << x << " y = " << y << " f(x, y) = " << func(x, y) << endl;
    }
    cout << "x = " << x << " y = " << y << " Extream value f(x, y) = " << func(x, y) << endl;

    fout.close();
    return 0;
}