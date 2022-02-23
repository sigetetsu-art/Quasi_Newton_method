#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

//Nは配列のサイズを示す

//係数Aの中身を表示する関数
void show_A(int N, vector<vector<double>> A){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            cout << A.at(i).at(j) << " ";
        }
        cout << endl;
    }
}

//定数Bの中身を表示する関数
void show_B(int N, vector<double> B){
    for(int i = 0; i < N; i++){
        cout << B.at(i) << endl;
    }
}

//解Xの中身を表示する関数
void show_X(int N, vector<double> X){
    for(int i = 0; i < N; i++){
        cout << "X" << i << " = " << X.at(i) << endl;
    }
}

//ピボット(行列の要素)の同じ列の中で最も大きい値を持つ場所の添え字を返す.引数のkは操作する列の添え字を意味する
inline int pivot_index(int N, vector<vector<double>> &A, int k){
    vector<double> D(N);
    int index = 0;

    for(int l = k; l < N; l++){
        D.at(l) = A.at(l).at(k) * A.at(l).at(k);
    }
    sort(D.begin(), D.end());
    double max = D[N - 1];//最大値を保存
    if(max < 0.00000001){
        cout << "The solution is not uniquely determined" << endl;
        exit(0);
    }

    for(int i = 0; i < N; i++){
        if(max == D.at(i)) index = i;
    }
    return index;//最大要素がある行のインデックスを返す
}

//選択したピボット(行列の要素)の行をk行目の列と総入れ替えする関数
inline void change_line(int N, vector<vector<double>> &A, vector<double> &B, int k, int index){
    double S;
    double T;

    T = B.at(k);
    B.at(k) = B.at(index);
    B.at(index) = T;

    for(int i = 0; i < N; i++){//k,indexで行，iで列を操作して，index行目をk行目と総入れ替え
        S = A.at(k).at(i);
        A.at(k).at(i) = A.at(index).at(i);
        A.at(index).at(i) = S;
    }
}

inline void remove_element(int N, vector<vector<double>> &A, vector<double> &B, int k, int i){
    if(fabs(A.at(k).at(k)) < 0.00001){
        A.at(k).at(k) = 0.00001;
    }
    double r = A.at(i).at(k) / A.at(k).at(k);

    for(int l = k; l < N; l++){
        A.at(i).at(l) = A.at(i).at(l) - r * A.at(k).at(l);
    }

    B.at(i) = B.at(i) - r * B.at(k);
}

//ガウスの消去法により行列を対角化
void gauss(int N, vector<vector<double>> &A, vector<double> &B){
    int p = 0;
    for(int k = 0; k < N - 1; k++){
        p = pivot_index(N, A, k);
        change_line(N, A, B, k, p);
        for(int i = k + 1; i < N; i++){
            remove_element(N, A, B, k, i);
        }
    }
}

void calc_ans(int N, vector<vector<double>> &A, vector<double> &B, vector<double> &X){
    if(fabs(A.at(N - 1).at(N - 1)) < 0.00001){
        A.at(N - 1).at(N - 1) = 0.00001;
    }
    X.at(N - 1) = B.at(N - 1) / A.at(N - 1).at(N - 1);
    for(int i = N - 2; i >= 0; i--){
        double S = 0;
        for(int j = i + 1; j < N; j++){
            S = S + A.at(i).at(j) * X.at(j);
            X.at(i) = (B.at(i) - S) / A.at(i).at(i);
        }
    }
}

//転置行列を作成する関数(2次元配列)
inline vector<vector<double>> make_tenchi(vector<vector<double>> A){
    int row = A.size(); //行列の行数
    int column = A.at(0).size(); //行列の列数
    vector<vector<double>> At(column, vector<double>(row));
    for(int i = 0; i < column; i++){
        for(int j = 0; j < row; j++){
            At.at(i).at(j) = A.at(j).at(i);
        }
    }
    return At;
}

//転置行列を作成する関数(1次元配列)
inline vector<vector<double>> make_tenchi(vector<double> A){
    int column = A.size(); //行列の列数
    vector<vector<double>> At(column, vector<double>(0));
    for(int i = 0; i < column; i++){
        At.at(i).at(0) = A.at(i);
    }
    return At;
}

//行列の和A+Bを計算する関数
inline vector<vector<double>> calc_matrix_sum(vector<vector<double>> A, vector<vector<double>> B){
    int row_A = A.size(); //行列Aの行数
    int column_A = A.at(0).size(); //行列Aの列数
    int row_B = B.size(); //行列Bの行数
    int column_B = B.at(0).size(); //行列Bの列数

    if(row_A != row_B || column_A != column_B){
        cout << "Can't calculate sum of matrix" << endl;
        exit(8);
    }

    vector<vector<double>> matrix_C(row_A, vector<double>(column_A));

    for(int i = 0; i < row_A; i++){
        for(int j = 0; j < column_A; j++){
            matrix_C.at(i).at(j) = A.at(i).at(j) + B.at(i).at(j);
        }
    }
    return matrix_C;
}

//行列の差A-Bを計算する関数
inline vector<vector<double>> calc_matrix_difference(vector<vector<double>> A, vector<vector<double>> B){
    int row_A = A.size(); //行列Aの行数
    int column_A = A.at(0).size(); //行列Aの列数
    int row_B = B.size(); //行列Bの行数
    int column_B = B.at(0).size(); //行列Bの列数

    if(row_A != row_B || column_A != column_B){
        cout << "Can't calculate difference of matrix" << endl;
        exit(8);
    }

    vector<vector<double>> matrix_C(row_A, vector<double>(column_A));

    for(int i = 0; i < row_A; i++){
        for(int j = 0; j < column_A; j++){
            matrix_C.at(i).at(j) = A.at(i).at(j) - B.at(i).at(j);
        }
    }
    return matrix_C;
}


//行列の積ABを計算する関数
inline vector<vector<double>> calc_matrix_product(vector<vector<double>> A, vector<vector<double>> B){
    int row_A = A.size(); //行列Aの行数
    int column_A = A.at(0).size(); //行列Aの列数
    int row_B = B.size(); //行列Bの行数
    int column_B = B.at(0).size(); //行列Bの列数 
    vector<vector<double>> matrix_C(row_A, vector<double>(column_B)); //AとBの積の結果
    vector<double> vec_A(column_A); //column_A = row_B
    vector<double> vec_B(row_B);
    int inner_product; //内積を格納する

    if(column_A != row_B){
        cout << "Can't calculate product of matrix" << endl;
        exit(8);
    }

    for(int i = 0; i < row_A; i++){
        for(int j = 0; j < column_B; j++){
            for(int k = 0; k < column_A; k++){
                vec_A.at(k) = A.at(i).at(k); //行列Aのi行目の成分だけのベクトルを用意
                vec_B.at(k) = B.at(k).at(j); //行列Bのj列目の成分だけのベクトルを用意
            }
            inner_product = 0.0;
            for(int k = 0; k < column_A; k++){
                inner_product = inner_product + vec_A.at(k) * vec_B.at(k);
            }
                matrix_C.at(i).at(j) = inner_product;
        }
    }
    return matrix_C;
}

//行列Aを定数倍する関数(二次元配列)
inline vector<vector<double>> matrix_times_const(vector<vector<double>> coeff, vector<vector<double>> A){
    int row_A = A.size(); //行列Aの行数
    int column_A = A.at(0).size(); //行列Aの列数
    for(int i = 0; i < row_A; i++){
        for(int j = 0; j < column_A; j++){
            A.at(i).at(j) = coeff.at(0).at(0) * A.at(i).at(j);
        }
    }
    return A;
}

//行列Aを定数倍する関数(二次元配列)
inline vector<vector<double>> matrix_times_const(double coeff, vector<vector<double>> A){
    int row_A = A.size(); //行列Aの行数
    int column_A = A.at(0).size(); //行列Aの列数
    for(int i = 0; i < row_A; i++){
        for(int j = 0; j < column_A; j++){
            A.at(i).at(j) = coeff * A.at(i).at(j);
        }
    }
    return A;
}
