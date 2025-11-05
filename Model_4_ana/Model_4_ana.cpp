#include "model4_ana_parameter.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

double B, NM, N;
double theta, epsilon, r, p;
double lamH, lamL, musH, musL, mufH, mufL;
double MAX_K, MAX_I, Q_MAX;

double Le_T, Lq_T, Lb_T, W_T, Wq_T, Wb_T, Pb_T, Th_T;
double Ls_T, Lf_T, Lv_T, Ws_T, Wf_T, Rs_T, Rf_T, Ths_T, Thf_T, Lvs_T, Lvf_T, Nb_T;

double Le_H, Lq_H, Lb_H, W_H, Wq_H, Wb_H, Pb_H, Th_H;
double Ls_H, Lf_H, Lv_H, Ws_H, Wf_H, Rs_H, Rf_H, Ths_H, Thf_H, Lvs_H, Lvf_H, Nb_H;

double Le_L, Lq_L, Lb_L, W_L, Wq_L, Wb_L, Pb_L, Th_L;
double Ls_L, Lf_L, Lv_L, Ws_L, Wf_L, Rs_L, Rf_L, Ths_L, Thf_L, Lvs_L, Lvf_L, Nb_L;

int iter;
double sum, convergence, normalize, error;

int outcheck, incheck;

vector<vector<vector<vector<vector<vector<vector<double>>>>>>> pi;
vector<vector<vector<vector<vector<vector<vector<double>>>>>>> old_pi;
vector<vector<vector<vector<vector<double>>>>> check;

class csv1 {
public:
    // Constructor
    csv1(string filename) : m_filename(filename) {}

    string m_filename;
    fstream m_fout; // file out

    // 清空檔案（用於新檔案）
    void clear_file() {
        m_fout.open(m_filename, ios::out); // 清空檔案內容
        m_fout.close();
    }

    // 打印標題（根據 case 動態設定標題）
    void print_title(int case_id, const vector<string>& columns) {
        m_fout.open(m_filename, ios::out | ios::app); // 追加模式
        m_fout << "Case " << case_id << " Title:" << endl;
        for (size_t i = 0; i < columns.size(); i++) {
            m_fout << columns[i];
            if (i != columns.size() - 1)
                m_fout << ", ";
        }
        m_fout << endl;
        m_fout.close();
    }

    // 打印數據
    void print_data(const vector<double>& data) {
        m_fout.open(m_filename, ios::out | ios::app);
        for (size_t j = 0; j < data.size(); j++) {
            m_fout << data[j];
            if (j != data.size() - 1)
                m_fout << ", ";
        }
        m_fout << endl;
        m_fout.close();
    }
};

void update_input(int pointer, int aa) {
    // 基準值與步長設定
    double base_value = (pointer == 1) ? DEFAULT_lamH
        : (pointer == 2) ? DEFAULT_lamL
        : (pointer == 3) ? DEFAULT_musH
        : (pointer == 4) ? DEFAULT_musL
        : (pointer == 5) ? DEFAULT_mufH
        : (pointer == 6) ? DEFAULT_mufL
        : (pointer == 7) ? DEFAULT_theta
        : (pointer == 8) ? DEFAULT_epsilon
        : (pointer == 9) ? DEFAULT_r
        : (pointer == 10) ? DEFAULT_p
        : (pointer == 11) ? DEFAULT_B
        : (pointer == 12) ? DEFAULT_NM
        : -1;

    double step = 1;
    if (base_value < 1) step = 0.1;

    // 計算新的值
    double updated_value = base_value + step * (aa - 2);

    // 使用一個數組或類型參照對應參數更新
    switch (pointer) {
    case 1:
        lamH = updated_value;
        break;
    case 2:
        lamL = updated_value;
        break;
    case 3:
        musH = updated_value;
        break;
    case 4:
        musL = updated_value;
        break;
    case 5:
        mufH = updated_value;
        break;
    case 6:
        mufL = updated_value;
        break;
    case 7:
        theta = updated_value;
        break;
    case 8:
        epsilon = updated_value;
        break;
    case 9:
        r = updated_value;
        break;
    case 10:
        p = updated_value;
        break;
    case 11:
        B = updated_value;
        break;
    case 12:
        NM = updated_value;
        break;
    }

    N = 3 * NM + 1;
    MAX_K = 2 * NM + 1;
    MAX_I = NM + 1;
    Q_MAX = 2 * B - 1;
}

void Initialization(int pointer, int aa) {
    // rate
    B = DEFAULT_B;
    NM = DEFAULT_NM;
    lamH = DEFAULT_lamH;
    lamL = DEFAULT_lamL;
    musH = DEFAULT_musH;
    musL = DEFAULT_musL;
    mufH = DEFAULT_mufH;
    mufL = DEFAULT_mufL;
    theta = DEFAULT_theta;
    epsilon = DEFAULT_epsilon;
    r = DEFAULT_r;
    p = DEFAULT_p;

    update_input(pointer, aa);

    pi = vector<vector<vector<vector<vector<vector<vector<double>>>>>>>(
        2 * B,
        vector<vector<vector<vector<vector<vector<double>>>>>>(
            2 * B,
            vector<vector<vector<vector<vector<double>>>>>(
                B + 1,
                vector<vector<vector<vector<double>>>>(
                    B + 1,
                    vector<vector<vector<double>>>(
                        2 * NM + 2, vector<vector<double>>(NM + 2, vector<double>(NM + 2, 0.0)))))));
    old_pi = vector<vector<vector<vector<vector<vector<vector<double>>>>>>>(
        2 * B,
        vector<vector<vector<vector<vector<vector<double>>>>>>(
            2 * B,
            vector<vector<vector<vector<vector<double>>>>>(
                B + 1,
                vector<vector<vector<vector<double>>>>(
                    B + 1,
                    vector<vector<vector<double>>>(
                        2 * NM + 2, vector<vector<double>>(NM + 2, vector<double>(NM + 2, 0.0)))))));

    // Initialize
    int valid_state_count = 0;

    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int yH = 0; yH <= B; yH = yH + B) {
                for (int yL = 0; yL <= B; yL = yL + B) {

                    for (int k = 0; k <= MAX_K; k++) {
                        for (int i = 0; i <= NM + 1; i++) {
                            for (int j = 0; j <= NM + 1 - i; j++) {

                                //int x = (xH + xL);
                                //int y = (yH + yL);

                                if (k == MAX_K && i + j == NM + 1) continue;
                                else if ((yH + yL) == 0 && !(k == 0 && i == 0 && j == 0)) continue;
                                else if ((xH) >= B && (xH) <= Q_MAX && (yH + yL) == 0 && k == 0 && i == 0 && j == 0) continue;
                                else if ((xL) >= B && (xL) <= Q_MAX && (yH + yL) == 0 && k == 0 && i == 0 && j == 0) continue;
                                else if ((xH + xL) > Q_MAX) continue;
                                else if (yH != 0 && yL != 0) continue;

                                else valid_state_count++;
                            }
                        }
                    }
                }
            }
        }
    }
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int yH = 0; yH <= B; yH = yH + B) {
                for (int yL = 0; yL <= B; yL = yL + B) {

                    for (int k = 0; k <= MAX_K; k++) {
                        for (int i = 0; i <= NM + 1; i++) {
                            for (int j = 0; j <= NM + 1 - i; j++) {

                                if (k == MAX_K && i + j == NM + 1) {
                                    old_pi[xH][xL][yH][yL][k][i][j] = 0.0;
                                    pi[xH][xL][yH][yL][k][i][j] = 0.0;
                                }
                                else if ((yH + yL) == 0 && !(k == 0 && i == 0 && j == 0)) {
                                    old_pi[xH][xL][yH][yL][k][i][j] = 0.0;
                                    pi[xH][xL][yH][yL][k][i][j] = 0.0;
                                }
                                else if ((xH) >= B && (xH) <= Q_MAX && (yH + yL) == 0 && k == 0 && i == 0 && j == 0) {
                                    old_pi[xH][xL][yH][yL][k][i][j] = 0.0;
                                    pi[xH][xL][yH][yL][k][i][j] = 0.0;
                                }
                                else if ((xL) >= B && (xL) <= Q_MAX && (yH + yL) == 0 && k == 0 && i == 0 && j == 0) {
                                    old_pi[xH][xL][yH][yL][k][i][j] = 0.0;
                                    pi[xH][xL][yH][yL][k][i][j] = 0.0;
                                }
                                else if ((xH + xL) > Q_MAX) {
                                    old_pi[xH][xL][yH][yL][k][i][j] = 0.0;
                                    pi[xH][xL][yH][yL][k][i][j] = 0.0;
                                }
                                else if (yH != 0 && yL != 0) {
                                    old_pi[xH][xL][yH][yL][k][i][j] = 0.0;
                                    pi[xH][xL][yH][yL][k][i][j] = 0.0;
                                }

                                else {
                                    old_pi[xH][xL][yH][yL][k][i][j] = 1.0 / valid_state_count;
                                    pi[xH][xL][yH][yL][k][i][j] = 1.0 / valid_state_count;
                                }

                            }
                        }

                    }
                }
            }
        }
    }
    double sum2 = 0.0;
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int yH = 0; yH <= B; yH = yH + B) {
                for (int yL = 0; yL <= B; yL = yL + B) {

                    for (int k = 0; k <= MAX_K; k++) {
                        for (int i = 0; i <= NM + 1; i++) {
                            for (int j = 0; j <= NM + 1 - i; j++) {
                                sum2 += pi[xH][xL][yH][yL][k][i][j];
                            }
                        }
                    }
                }
            }
        }
    }
    printf("%f\n", sum2);


    // Normalization
    convergence = 1.0;

    // Iteration
    iter = 0;

    // Performance measure
    Le_T = 0, Lq_T = 0, Lb_T = 0, W_T = 0, Wq_T = 0, Wb_T = 0, Pb_T = 0, Th_T = 0;
    Ls_T = 0, Lf_T = 0, Ws_T = 0, Wf_T = 0, Ths_T = 0, Thf_T = 0, Rs_T = 0, Nb_T = 0;
    Lv_T = 0, Lvs_T = 0, Lvf_T = 0;

    Le_H = 0, Lq_H = 0, Lb_H = 0, W_H = 0, Wq_H = 0, Wb_H = 0, Pb_H = 0, Th_H = 0;
    Ls_H = 0, Lf_H = 0, Ws_H = 0, Wf_H = 0, Ths_H = 0, Thf_H = 0, Rs_H = 0, Nb_H = 0;
    Lv_H = 0, Lvs_H = 0, Lvf_H = 0;

    Le_L = 0, Lq_L = 0, Lb_L = 0, W_L = 0, Wq_L = 0, Wb_L = 0, Pb_L = 0, Th_L = 0;
    Ls_L = 0, Lf_L = 0, Ws_L = 0, Wf_L = 0, Ths_L = 0, Thf_L = 0, Rs_L = 0, Nb_L = 0;
    Lv_L = 0, Lvs_L = 0, Lvf_L = 0;

    outcheck = 0, incheck = 0;
    return;
}

void Iteration() {
    iter++;
    // cout << "Iteration " << iteration_time << " times." << endl;
    cout << ".";

    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int yH = 0; yH <= B; yH = yH + B) {
                for (int yL = 0; yL <= B; yL = yL + B) {

                    for (int k = 0; k <= MAX_K; k++) {
                        for (int i = 0; i <= NM + 1; i++) {
                            for (int j = 0; j <= NM + 1 - i; j++) {

                                if (k == MAX_K && i + j == NM + 1) continue;
                                else if ((yH + yL) == 0 && !(k == 0 && i == 0 && j == 0)) continue;
                                else if ((xH) >= B && (xH) <= Q_MAX && (yH + yL) == 0 && k == 0 && i == 0 && j == 0) continue;
                                else if ((xL) >= B && (xL) <= Q_MAX && (yH + yL) == 0 && k == 0 && i == 0 && j == 0) continue;
                                else if ((xH + xL) > Q_MAX) continue;
                                else if (yH != 0 && yL != 0) continue;

                                else if (xH == 0) {
                                    if (xL == 0) {
                                        if (yH == 0) {
                                            if (yL == 0) {
                                                // 1
                                                double sum_musH_1 = 0.0;
                                                for (int i11 = 0; i11 <= NM; i11++) {
                                                    for (int j11 = 0; j11 <= NM - i11; j11++) {
                                                        sum_musH_1 += pi[xH][xL][B][0][MAX_K][i11][j11];
                                                    }
                                                }
                                                double sum_mufH_1 = 0.0;
                                                for (int k11 = 0; k11 <= MAX_K - 1; k11++) {
                                                    for (int i12 = 0; i12 <= NM + 1; i12++) {
                                                        sum_mufH_1 += pi[xH][xL][B][0][k11][i12][NM + 1 - i12];
                                                    }
                                                }
                                                double sum_musL_1 = 0.0;
                                                for (int i13 = 0; i13 <= NM; i13++) {
                                                    for (int j12 = 0; j12 <= NM - i13; j12++) {
                                                        sum_musL_1 += pi[xH][xL][0][B][MAX_K][i13][j12];
                                                    }
                                                }
                                                double sum_mufL_1 = 0.0;
                                                for (int k12 = 0; k12 <= MAX_K - 1; k12++) {
                                                    for (int i14 = 0; i14 <= NM + 1; i14++) {
                                                        sum_mufL_1 += pi[xH][xL][0][B][k12][i14][NM + 1 - i14];
                                                    }
                                                }
                                                pi[xH][xL][yH][yL][k][i][j] = ( musH * sum_musH_1 +
                                                                                mufH * sum_mufH_1 +
                                                                                musL * sum_musL_1 +
                                                                                mufL * sum_mufL_1) /
                                                                                (lamH + lamL);
                                                // 1 end
                                            }
                                            else if (yL == B) {
                                                // 2
                                                // k = 0
                                                if (k == 0 && i == 0 && j == 0) {
                                                    double sum_musH_2L = 0.0;
                                                    for (int i21L = 0; i21L <= NM; i21L++) {
                                                        for (int j21L = 0; j21L <= NM - i21L; j21L++) {
                                                            sum_musH_2L += pi[xH][xL + (B - 1)][B][0][MAX_K][i21L][j21L];
                                                        }
                                                    }
                                                    double sum_mufH_2L = 0.0;
                                                    for (int k21L = 0; k21L <= MAX_K - 1; k21L++) {
                                                        for (int i22L = 0; i22L <= NM + 1; i22L++) {
                                                            sum_mufH_2L += pi[xH][xL + (B - 1)][B][0][k21L][i22L][NM + 1 - i22L];
                                                        }
                                                    }
                                                    double sum_musL_2L = 0.0;
                                                    for (int i23L = 0; i23L <= NM; i23L++) {
                                                        for (int j22L = 0; j22L <= NM - i23L; j22L++) {
                                                            sum_musL_2L += pi[xH][xL + (B - 1)][0][B][MAX_K][i23L][j22L];
                                                        }
                                                    }
                                                    double sum_mufL_2L = 0.0;
                                                    for (int k22L = 0; k22L <= MAX_K - 1; k22L++) {
                                                        for (int i24L = 0; i24L <= NM + 1; i24L++) {
                                                            sum_mufL_2L += pi[xH][xL + (B - 1)][0][B][k22L][i24L][NM + 1 - i24L];
                                                        }
                                                    }
                                                    pi[xH][xL][yH][yL][k][i][j] = ( musH * sum_musH_2L +
                                                                                    mufH * sum_mufH_2L +
                                                                                    musL * sum_musL_2L +
                                                                                    mufL * sum_mufL_2L +
                                                                                    lamL * pi[xH][xL + (B - 1)][yH][yL - B][k][i][j] +
                                                                                    (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                                                    (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufL);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufL);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + musL);
                                                }
                                                // 2 end
                                            }
                                        }
                                        else if (yH == B) {
                                            if (yL == 0) {
                                                // 3
                                                // k = 0
                                                if (k == 0 && i == 0 && j == 0) {
                                                    double sum_musH_2L = 0.0;
                                                    for (int i21L = 0; i21L <= NM; i21L++) {
                                                        for (int j21L = 0; j21L <= NM - i21L; j21L++) {
                                                            sum_musH_2L += pi[xH + B][xL][B][0][MAX_K][i21L][j21L];
                                                        }
                                                    }
                                                    double sum_mufH_2L = 0.0;
                                                    for (int k21L = 0; k21L <= MAX_K - 1; k21L++) {
                                                        for (int i22L = 0; i22L <= NM + 1; i22L++) {
                                                            sum_mufH_2L += pi[xH + B][xL][B][0][k21L][i22L][NM + 1 - i22L];
                                                        }
                                                    }
                                                    double sum_musL_2L = 0.0;
                                                    for (int i23L = 0; i23L <= NM; i23L++) {
                                                        for (int j22L = 0; j22L <= NM - i23L; j22L++) {
                                                            sum_musL_2L += pi[xH + B][xL][0][B][MAX_K][i23L][j22L];
                                                        }
                                                    }
                                                    double sum_mufL_2L = 0.0;
                                                    for (int k22L = 0; k22L <= MAX_K - 1; k22L++) {
                                                        for (int i24L = 0; i24L <= NM + 1; i24L++) {
                                                            sum_mufL_2L += pi[xH + B][xL][0][B][k22L][i24L][NM + 1 - i24L];
                                                        }
                                                    }
                                                    pi[xH][xL][yH][yL][k][i][j] = ( musH * sum_musH_2L +
                                                                                    mufH * sum_mufH_2L +
                                                                                    musL * sum_musL_2L +
                                                                                    mufL * sum_mufL_2L +
                                                                                    lamH * pi[xH + (B - 1)][xL][yH - B][yL][k][i][j] +
                                                                                    (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                                                    (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufH);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufH);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + musH);
                                                }
                                                // 3 emd
                                            }
                                        }
                                    }
                                    else if (xL >= 1 && xL <= B - 1) {
                                        if (yH == 0) {
                                            if (yL == 0) {
                                                // 4
                                                double sum_musH_1 = 0.0;
                                                for (int i11 = 0; i11 <= NM; i11++) {
                                                    for (int j11 = 0; j11 <= NM - i11; j11++) {
                                                        sum_musH_1 += pi[xH][xL][B][0][MAX_K][i11][j11];
                                                    }
                                                }
                                                double sum_mufH_1 = 0.0;
                                                for (int k11 = 0; k11 <= MAX_K - 1; k11++) {
                                                    for (int i12 = 0; i12 <= NM + 1; i12++) {
                                                        sum_mufH_1 += pi[xH][xL][B][0][k11][i12][NM + 1 - i12];
                                                    }
                                                }
                                                double sum_musL_1 = 0.0;
                                                for (int i13 = 0; i13 <= NM; i13++) {
                                                    for (int j12 = 0; j12 <= NM - i13; j12++) {
                                                        sum_musL_1 += pi[xH][xL][0][B][MAX_K][i13][j12];
                                                    }
                                                }
                                                double sum_mufL_1 = 0.0;
                                                for (int k12 = 0; k12 <= MAX_K - 1; k12++) {
                                                    for (int i14 = 0; i14 <= NM + 1; i14++) {
                                                        sum_mufL_1 += pi[xH][xL][0][B][k12][i14][NM + 1 - i14];
                                                    }
                                                }
                                                pi[xH][xL][yH][yL][k][i][j] = ( musH * sum_musH_1 +
                                                                                mufH * sum_mufH_1 +
                                                                                musL * sum_musL_1 +
                                                                                mufL * sum_mufL_1 +
                                                                                lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                                                (lamH + lamL);
                                                // 4 end
                                            }
                                            else if (yL == B) {
                                                // 5
                                                // k = 0
                                                if (k == 0 && i == 0 && j == 0) {
                                                    double sum_musH_2L = 0.0;
                                                    for (int i21L = 0; i21L <= NM; i21L++) {
                                                        for (int j21L = 0; j21L <= NM - i21L; j21L++) {
                                                            sum_musH_2L += pi[xH][xL + B][B][0][MAX_K][i21L][j21L];
                                                        }
                                                    }
                                                    double sum_mufH_2L = 0.0;
                                                    for (int k21L = 0; k21L <= MAX_K - 1; k21L++) {
                                                        for (int i22L = 0; i22L <= NM + 1; i22L++) {
                                                            sum_mufH_2L += pi[xH][xL + B][B][0][k21L][i22L][NM + 1 - i22L];
                                                        }
                                                    }
                                                    double sum_musL_2L = 0.0;
                                                    for (int i23L = 0; i23L <= NM; i23L++) {
                                                        for (int j22L = 0; j22L <= NM - i23L; j22L++) {
                                                            sum_musL_2L += pi[xH][xL + B][0][B][MAX_K][i23L][j22L];
                                                        }
                                                    }
                                                    double sum_mufL_2L = 0.0;
                                                    for (int k22L = 0; k22L <= MAX_K - 1; k22L++) {
                                                        for (int i24L = 0; i24L <= NM + 1; i24L++) {
                                                            sum_mufL_2L += pi[xH][xL + B][0][B][k22L][i24L][NM + 1 - i24L];
                                                        }
                                                    }
                                                    pi[xH][xL][yH][yL][k][i][j] = ( musH * sum_musH_2L +
                                                                                    mufH * sum_mufH_2L +
                                                                                    musL * sum_musL_2L +
                                                                                    mufL * sum_mufL_2L +
                                                                                    lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                                                    (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                                                    (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + mufL);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + mufL);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + musL);
                                                }
                                                // 5 end
                                            }
                                        }
                                        else if (yH == B) {
                                            if (yL == 0) {
                                                // 6
                                                // k = 0
                                                if (k == 0 && i == 0 && j == 0) {
                                                    double sum_musH_2L = 0.0;
                                                    for (int i21L = 0; i21L <= NM; i21L++) {
                                                        for (int j21L = 0; j21L <= NM - i21L; j21L++) {
                                                            sum_musH_2L += pi[xH + B][xL][B][0][MAX_K][i21L][j21L];
                                                        }
                                                    }
                                                    double sum_mufH_2L = 0.0;
                                                    for (int k21L = 0; k21L <= MAX_K - 1; k21L++) {
                                                        for (int i22L = 0; i22L <= NM + 1; i22L++) {
                                                            sum_mufH_2L += pi[xH + B][xL][B][0][k21L][i22L][NM + 1 - i22L];
                                                        }
                                                    }
                                                    double sum_musL_2L = 0.0;
                                                    for (int i23L = 0; i23L <= NM; i23L++) {
                                                        for (int j22L = 0; j22L <= NM - i23L; j22L++) {
                                                            sum_musL_2L += pi[xH + B][xL][0][B][MAX_K][i23L][j22L];
                                                        }
                                                    }
                                                    double sum_mufL_2L = 0.0;
                                                    for (int k22L = 0; k22L <= MAX_K - 1; k22L++) {
                                                        for (int i24L = 0; i24L <= NM + 1; i24L++) {
                                                            sum_mufL_2L += pi[xH + B][xL][0][B][k22L][i24L][NM + 1 - i24L];
                                                        }
                                                    }
                                                    pi[xH][xL][yH][yL][k][i][j] = ( musH * sum_musH_2L +
                                                                                    mufH * sum_mufH_2L +
                                                                                    musL * sum_musL_2L +
                                                                                    mufL * sum_mufL_2L +
                                                                                    lamH * pi[xH + (B - 1)][xL][yH - B][yL][k][i][j] +
                                                                                    lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                                                    (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                                                    (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + mufH);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + mufH);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = ((N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j] +
                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                        (lamH + lamL + musH);
                                                }
                                                // 6 emd
                                            }
                                        }
                                    }
                                    else if (xL >= B && xL <= Q_MAX - 1) {
                                        if (yH == 0) {
                                            if (yL == B) {
                                                // 7
                                                if (k == 0 && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufL);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufL);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + musL);
                                                }
                                                // 7 end
                                            }
                                        }
                                        else if (yH == B) {
                                            if (yL == 0) {
                                                // 8
                                                // K = 0
                                                if (k == 0 && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufH);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufH);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + musH);
                                                }
                                                // 8 emd
                                            }
                                        }
                                    }
                                    else if (xL == Q_MAX) {
                                        if (yH == 0) {
                                            if (yL == B) {
                                                // 9
                                                // K = 0
                                                if (k == 0 && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        ((j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (mufL);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        ((j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (mufL);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (musL);
                                                }
                                                // 9 end
                                            }
                                        }
                                        else if (yH == B) {
                                            if (yL == 0) {
                                                // 10
                                                // K = 0
                                                if (k == 0 && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        ((j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (mufH);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        ((j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (mufH);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (musH);
                                                }
                                                // 10 emd
                                            }
                                        }
                                    }
                                }
                                else if (xH >= 1 && xH <= B - 1) {
                                    if (xL == 0) {
                                        if (yH == 0) {
                                            if (yL == 0) {
                                                // 11
                                                double sum_musH_1 = 0.0;
                                                for (int i11 = 0; i11 <= NM; i11++) {
                                                    for (int j11 = 0; j11 <= NM - i11; j11++) {
                                                        sum_musH_1 += pi[xH][xL][B][0][MAX_K][i11][j11];
                                                    }
                                                }
                                                double sum_mufH_1 = 0.0;
                                                for (int k11 = 0; k11 <= MAX_K - 1; k11++) {
                                                    for (int i12 = 0; i12 <= NM + 1; i12++) {
                                                        sum_mufH_1 += pi[xH][xL][B][0][k11][i12][NM + 1 - i12];
                                                    }
                                                }
                                                double sum_musL_1 = 0.0;
                                                for (int i13 = 0; i13 <= NM; i13++) {
                                                    for (int j12 = 0; j12 <= NM - i13; j12++) {
                                                        sum_musL_1 += pi[xH][xL][0][B][MAX_K][i13][j12];
                                                    }
                                                }
                                                double sum_mufL_1 = 0.0;
                                                for (int k12 = 0; k12 <= MAX_K - 1; k12++) {
                                                    for (int i14 = 0; i14 <= NM + 1; i14++) {
                                                        sum_mufL_1 += pi[xH][xL][0][B][k12][i14][NM + 1 - i14];
                                                    }
                                                }
                                                pi[xH][xL][yH][yL][k][i][j] = ( musH * sum_musH_1 +
                                                                                mufH * sum_mufH_1 +
                                                                                musL * sum_musL_1 +
                                                                                mufL * sum_mufL_1 +
                                                                                lamH * pi[xH - 1][xL][yH][yL][k][i][j]) /
                                                                                (lamH + lamL);
                                                // 11 end
                                            }
                                            else if (yL == B) {
                                                // 12
                                                // k = 0
                                                if (k == 0 && i == 0 && j == 0) {
                                                    double sum_musH_10 = 0.0;
                                                    for (int i101 = 0; i101 <= NM; i101++) {
                                                        for (int j101 = 0; j101 <= NM - i101; j101++) {
                                                            sum_musH_10 += pi[xH][xL + B][B][0][MAX_K][i101][j101];
                                                        }
                                                    }
                                                    double sum_mufH_10 = 0.0;
                                                    for (int k101 = 0; k101 <= MAX_K - 1; k101++) {
                                                        for (int i102 = 0; i102 <= NM + 1; i102++) {
                                                            sum_mufH_10 += pi[xH][xL + B][B][0][k101][i102][NM + 1 - i102];
                                                        }
                                                    }
                                                    double sum_musL_10 = 0.0;
                                                    for (int i103 = 0; i103 <= NM; i103++) {
                                                        for (int j102 = 0; j102 <= NM - i103; j102++) {
                                                            sum_musL_10 += pi[xH][xL + B][0][B][MAX_K][i103][j102];
                                                        }
                                                    }
                                                    double sum_mufL_10 = 0.0;
                                                    for (int k102 = 0; k102 <= MAX_K - 1; k102++) {
                                                        for (int i104 = 0; i104 <= NM + 1; i104++) {
                                                            sum_mufL_10 += pi[xH][xL + B][0][B][k102][i104][NM + 1 - i104];
                                                        }
                                                    }
                                                    pi[xH][xL][yH][yL][k][i][j] = ( musH * sum_musH_10 +
                                                                                    mufH * sum_mufH_10 +
                                                                                    musL * sum_musL_10 +
                                                                                    mufL * sum_mufL_10 +
                                                                                    lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                                                    lamL * pi[xH][xL + (B - 1)][yH][yL - B][k][i][j] +
                                                                                    (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                                                    (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufL);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufL);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + musL);
                                                }
                                                // 12 end
                                            }
                                        }
                                        else if (yH == B) {
                                            if (yL == 0) {
                                                // 13
                                                if (k == 0 && i == 0 && j == 0) {
                                                    double sum_musH_10 = 0.0;
                                                    for (int i101 = 0; i101 <= NM; i101++) {
                                                        for (int j101 = 0; j101 <= NM - i101; j101++) {
                                                            sum_musH_10 += pi[xH + B][xL][B][0][MAX_K][i101][j101];
                                                        }
                                                    }
                                                    double sum_mufH_10 = 0.0;
                                                    for (int k101 = 0; k101 <= MAX_K - 1; k101++) {
                                                        for (int i102 = 0; i102 <= NM + 1; i102++) {
                                                            sum_mufH_10 += pi[xH + B][xL][B][0][k101][i102][NM + 1 - i102];
                                                        }
                                                    }
                                                    double sum_musL_10 = 0.0;
                                                    for (int i103 = 0; i103 <= NM; i103++) {
                                                        for (int j102 = 0; j102 <= NM - i103; j102++) {
                                                            sum_musL_10 += pi[xH + B][xL][0][B][MAX_K][i103][j102];
                                                        }
                                                    }
                                                    double sum_mufL_10 = 0.0;
                                                    for (int k102 = 0; k102 <= MAX_K - 1; k102++) {
                                                        for (int i104 = 0; i104 <= NM + 1; i104++) {
                                                            sum_mufL_10 += pi[xH + B][xL][0][B][k102][i104][NM + 1 - i104];
                                                        }
                                                    }
                                                    pi[xH][xL][yH][yL][k][i][j] = ( musH * sum_musH_10 +
                                                                                    mufH * sum_mufH_10 +
                                                                                    musL * sum_musL_10 +
                                                                                    mufL * sum_mufL_10 +
                                                                                    lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                                                    (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                                                    (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufH);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufH);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + musH);
                                                }
                                                // 13 emd
                                            }
                                        }
                                    }
                                    else if (xL >= 1 && xL <= B - 1) {
                                        if (yH == 0) {
                                            if (yL == 0) {
                                                // 14
                                                double sum_musH_1 = 0.0;
                                                for (int i11 = 0; i11 <= NM; i11++) {
                                                    for (int j11 = 0; j11 <= NM - i11; j11++) {
                                                        sum_musH_1 += pi[xH][xL][B][0][MAX_K][i11][j11];
                                                    }
                                                }
                                                double sum_mufH_1 = 0.0;
                                                for (int k11 = 0; k11 <= MAX_K - 1; k11++) {
                                                    for (int i12 = 0; i12 <= NM + 1; i12++) {
                                                        sum_mufH_1 += pi[xH][xL][B][0][k11][i12][NM + 1 - i12];
                                                    }
                                                }
                                                double sum_musL_1 = 0.0;
                                                for (int i13 = 0; i13 <= NM; i13++) {
                                                    for (int j12 = 0; j12 <= NM - i13; j12++) {
                                                        sum_musL_1 += pi[xH][xL][0][B][MAX_K][i13][j12];
                                                    }
                                                }
                                                double sum_mufL_1 = 0.0;
                                                for (int k12 = 0; k12 <= MAX_K - 1; k12++) {
                                                    for (int i14 = 0; i14 <= NM + 1; i14++) {
                                                        sum_mufL_1 += pi[xH][xL][0][B][k12][i14][NM + 1 - i14];
                                                    }
                                                }
                                                pi[xH][xL][yH][yL][k][i][j] = ( musH * sum_musH_1 +
                                                                                mufH * sum_mufH_1 +
                                                                                musL * sum_musL_1 +
                                                                                mufL * sum_mufL_1 +
                                                                                lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                                                lamL * pi[xH][xL - 1][yH][yL][k][i][j]) /
                                                                                (lamH + lamL);
                                                // 14 end
                                            }
                                            else if (yL == B) {
                                                if (xH + xL < B) {
                                                    // 15
                                                    // k = 0
                                                    if (k == 0 && i == 0 && j == 0) {
                                                        double sum_musH_10 = 0.0;
                                                        for (int i101 = 0; i101 <= NM; i101++) {
                                                            for (int j101 = 0; j101 <= NM - i101; j101++) {
                                                                sum_musH_10 += pi[xH][xL + B][B][0][MAX_K][i101][j101];
                                                            }
                                                        }
                                                        double sum_mufH_10 = 0.0;
                                                        for (int k101 = 0; k101 <= MAX_K - 1; k101++) {
                                                            for (int i102 = 0; i102 <= NM + 1; i102++) {
                                                                sum_mufH_10 += pi[xH][xL + B][B][0][k101][i102][NM + 1 - i102];
                                                            }
                                                        }
                                                        double sum_musL_10 = 0.0;
                                                        for (int i103 = 0; i103 <= NM; i103++) {
                                                            for (int j102 = 0; j102 <= NM - i103; j102++) {
                                                                sum_musL_10 += pi[xH][xL + B][0][B][MAX_K][i103][j102];
                                                            }
                                                        }
                                                        double sum_mufL_10 = 0.0;
                                                        for (int k102 = 0; k102 <= MAX_K - 1; k102++) {
                                                            for (int i104 = 0; i104 <= NM + 1; i104++) {
                                                                sum_mufL_10 += pi[xH][xL + B][0][B][k102][i104][NM + 1 - i104];
                                                            }
                                                        }
                                                        pi[xH][xL][yH][yL][k][i][j] = ( musH * sum_musH_10 +
                                                                                        mufH * sum_mufH_10 +
                                                                                        musL * sum_musL_10 +
                                                                                        mufL * sum_mufL_10 +
                                                                                        lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufL);
                                                    }
                                                    // k = 1 ~ 2NM
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufL);
                                                    }
                                                    // k = MAX_K
                                                    else if (k == MAX_K && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j == NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musL);
                                                    }
                                                    else if (k == MAX_K && i == NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + musL);
                                                    }
                                                    // 15 end
                                                }
                                                else if (xH + xL >= B && xH + xL <= Q_MAX - 1) {
                                                    // 16
                                                    // k = 0
                                                    if (k == 0 && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufL);
                                                    }
                                                    // k = 1 ~ 2NM
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufL);
                                                    }
                                                    // k = MAX_K
                                                    else if (k == MAX_K && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j == NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musL);
                                                    }
                                                    else if (k == MAX_K && i == NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + musL);
                                                    }
                                                    // 16 end
                                                }
                                            }
                                        }
                                        else if (yH == B) {
                                            if (yL == 0) {
                                                if (xH + xL < B) {
                                                    // 17
                                                    // k = 0
                                                    if (k == 0 && i == 0 && j == 0) {
                                                        double sum_musH_10 = 0.0;
                                                        for (int i101 = 0; i101 <= NM; i101++) {
                                                            for (int j101 = 0; j101 <= NM - i101; j101++) {
                                                                sum_musH_10 += pi[xH + B][xL][B][0][MAX_K][i101][j101];
                                                            }
                                                        }
                                                        double sum_mufH_10 = 0.0;
                                                        for (int k101 = 0; k101 <= MAX_K - 1; k101++) {
                                                            for (int i102 = 0; i102 <= NM + 1; i102++) {
                                                                sum_mufH_10 += pi[xH + B][xL][B][0][k101][i102][NM + 1 - i102];
                                                            }
                                                        }
                                                        double sum_musL_10 = 0.0;
                                                        for (int i103 = 0; i103 <= NM; i103++) {
                                                            for (int j102 = 0; j102 <= NM - i103; j102++) {
                                                                sum_musL_10 += pi[xH + B][xL][0][B][MAX_K][i103][j102];
                                                            }
                                                        }
                                                        double sum_mufL_10 = 0.0;
                                                        for (int k102 = 0; k102 <= MAX_K - 1; k102++) {
                                                            for (int i104 = 0; i104 <= NM + 1; i104++) {
                                                                sum_mufL_10 += pi[xH + B][xL][0][B][k102][i104][NM + 1 - i104];
                                                            }
                                                        }
                                                        pi[xH][xL][yH][yL][k][i][j] = ( musH * sum_musH_10 +
                                                                                        mufH * sum_mufH_10 +
                                                                                        musL * sum_musL_10 +
                                                                                        mufL * sum_mufL_10 +
                                                                                        lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                                                        lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufH);
                                                    }
                                                    // k = 1 ~ 2NM
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufH);
                                                    }
                                                    // k = MAX_K
                                                    else if (k == MAX_K && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j == NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musH);
                                                    }
                                                    else if (k == MAX_K && i == NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + musH);
                                                    }
                                                    // 17 end
                                                }
                                                else if (xH + xL >= B && xH + xL <= Q_MAX - 1) {
                                                    // 18
                                                    // k = 0
                                                    if (k == 0 && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufH);
                                                    }
                                                    // k = 1 ~ 2NM
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufH);
                                                    }
                                                    // k = MAX_K
                                                    else if (k == MAX_K && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j == NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musH);
                                                    }
                                                    else if (k == MAX_K && i == NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + musH);
                                                    }
                                                    // 18 end
                                                }
                                            }
                                        }
                                    }
                                    else if (xL >= B && xL <= Q_MAX - 1) {
                                        if (yH == 0) {
                                            if (yL == B) {
                                                if (xH + xL >= B && xH + xL <= Q_MAX - 1) {
                                                    // 19
                                                    // k = 0
                                                    if (k == 0 && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufL);
                                                    }
                                                    // k = 1 ~ 2NM
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufL);
                                                    }
                                                    // k = MAX_K
                                                    else if (k == MAX_K && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j == NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musL);
                                                    }
                                                    else if (k == MAX_K && i == NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + musL);
                                                    }
                                                    // 19 end
                                                }
                                                else if (xH + xL == Q_MAX) {
                                                    // 20
                                                    // k = 0
                                                    if (k == 0 && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            ((j)*epsilon + mufL);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((j)*epsilon + mufL);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (mufL);
                                                    }
                                                    // k = 1 ~ 2NM
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            ((j)*epsilon + mufL);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((j)*epsilon + mufL);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (mufL);
                                                    }
                                                    // k = MAX_K
                                                    else if (k == MAX_K && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j == NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + musL);
                                                    }
                                                    else if (k == MAX_K && i == NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (musL);
                                                    }
                                                    // 20 end
                                                }
                                            }
                                        }
                                        else if (yH == B) {
                                            if (yL == 0) {
                                                if (xH + xL >= B && xH + xL <= Q_MAX - 1) {
                                                    // 21
                                                    // k = 0
                                                    if (k == 0 && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufH);
                                                    }
                                                    // k = 1 ~ 2NM
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufH);
                                                    }
                                                    // k = MAX_K
                                                    else if (k == MAX_K && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j == NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musH);
                                                    }
                                                    else if (k == MAX_K && i == NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + musH);
                                                    }
                                                    // 21 end
                                                }
                                                else if (xH + xL == Q_MAX) {
                                                    // 22
                                                    // k = 0
                                                    if (k == 0 && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            ((j)*epsilon + mufH);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((j)*epsilon + mufH);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (mufH);
                                                    }
                                                    // k = 1 ~ 2NM
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            ((j)*epsilon + mufH);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((j)*epsilon + mufH);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (mufH);
                                                    }
                                                    // k = MAX_K
                                                    else if (k == MAX_K && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j == NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + musH);
                                                    }
                                                    else if (k == MAX_K && i == NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (musH);
                                                    }
                                                    // 22 end
                                                }
                                            }
                                        }
                                    }
                                }                                
                                else if (xH >= B && xH <= Q_MAX - 1) {
                                    if (xL == 0) {
                                        if (yH == 0) {
                                            if (yL == B) {
                                                // 23
                                                // k = 0
                                                if (k == 0 && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufL);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufL);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + musL);
                                                }
                                                // 23 end
                                            }
                                        }
                                        else if (yH == B) {
                                            if (yL == 0) {
                                                // 24
                                                // k = 0
                                                if (k == 0 && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufH);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + (j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (lamH + lamL + mufH);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + (j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (lamH + lamL + musH);
                                                }
                                                // 24 end
                                            }
                                        }
                                    }
                                    else if (xL >= 1 && xL <= B - 1) {
                                        if (yH == 0) {
                                            if (yL == B) {
                                                if (xH + xL >= B && xH + xL <= Q_MAX - 1) {
                                                    // 25
                                                    // k = 0
                                                    if (k == 0 && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufL);
                                                    }
                                                    // k = 1 ~ 2NM
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufL);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufL);
                                                    }
                                                    // k = MAX_K
                                                    else if (k == MAX_K && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j == NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musL);
                                                    }
                                                    else if (k == MAX_K && i == NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + musL);
                                                    }
                                                    // 25 end
                                                }
                                                else if (xH + xL == Q_MAX) {
                                                    // 26
                                                    // k = 0
                                                    if (k == 0 && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            ((j)*epsilon + mufL);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((j)*epsilon + mufL);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (mufL);
                                                    }
                                                    // k = 1 ~ 2NM
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            ((j)*epsilon + mufL);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((j)*epsilon + mufL);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (mufL);
                                                    }
                                                    // k = MAX_K
                                                    else if (k == MAX_K && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j == NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + musL);
                                                    }
                                                    else if (k == MAX_K && i == NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (musL);
                                                    }
                                                    // 26 end
                                                }
                                            }
                                        }
                                        else if (yH == B) {
                                            if (yL == 0) {
                                                if (xH + xL >= B && xH + xL <= Q_MAX - 1) {
                                                    // 27
                                                    // k = 0
                                                    if (k == 0 && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufH);
                                                    }
                                                    // k = 1 ~ 2NM
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + (j)*epsilon + mufH);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (lamH + lamL + mufH);
                                                    }
                                                    // k = MAX_K
                                                    else if (k == MAX_K && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j == NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + (j)*epsilon + musH);
                                                    }
                                                    else if (k == MAX_K && i == NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (lamH + lamL + musH);
                                                    }
                                                    // 27 end
                                                }
                                                else if (xH + xL == Q_MAX) {
                                                    // 28
                                                    // k = 0
                                                    if (k == 0 && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            ((j)*epsilon + mufH);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((j)*epsilon + mufH);
                                                    }
                                                    else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k == 0 && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (mufH);
                                                    }
                                                    // k = 1 ~ 2NM
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                            ((j)*epsilon + mufH);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            ((j)*epsilon + mufH);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                    }
                                                    else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                            (mufH);
                                                    }
                                                    // k = MAX_K
                                                    else if (k == MAX_K && i == 0 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i == 0 && j == NM) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                    }
                                                    else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            ((j)*epsilon + musH);
                                                    }
                                                    else if (k == MAX_K && i == NM && j == 0) {
                                                        pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                            lamL * pi[xH][xL - 1][yH][yL][k][i][j] +
                                                            (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                            (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                            (musH);
                                                    }
                                                    // 28 end
                                                }
                                            }
                                        }
                                    }
                                }
                                else if (xH == Q_MAX) {
                                    if (xL == 0) {
                                        if (yH == 0) {
                                            if (yL == B) {
                                                // 29
                                                // k = 0
                                                if (k == 0 && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        ((j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((j)*epsilon + mufL);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (mufL);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        ((j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((j)*epsilon + mufL);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (mufL);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musL);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + musL);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (musL);
                                                }
                                                // 29 end
                                            }
                                        }
                                        else if (yH == B) {
                                            if (yL == 0) {
                                                // 30
                                                // k = 0
                                                if (k == 0 && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        ((j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((j)*epsilon + mufH);
                                                }
                                                else if (k == 0 && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k == 0 && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (mufH);
                                                }
                                                // k = 1 ~ 2NM
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j >= 1 && j <= NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == 0 && j == NM + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1]) /
                                                        ((j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM && j == NM - i + 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        ((j)*epsilon + mufH);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i >= 1 && i <= NM - 1 && j >= 1 && j <= NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + (N - k - i - j) * r * p);
                                                }
                                                else if (k >= 1 && k <= 2 * NM && i == NM + 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j]) /
                                                        (mufH);
                                                }
                                                // k = MAX_K
                                                else if (k == MAX_K && i == 0 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j >= 1 && j <= NM - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i == 0 && j == NM) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 2 && j >= 1 && j <= NM - i - 1) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (j + 1) * epsilon * pi[xH][xL][yH][yL][k][i][j + 1] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + (N - k - i - j) * theta + (N - k - i - j) * r * (1.0 - p) + musH);
                                                }
                                                else if (k == MAX_K && i >= 1 && i <= NM - 1 && j == NM - i) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - i - (j - 1)) * theta * pi[xH][xL][yH][yL][k][i][j - 1] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        ((j)*epsilon + musH);
                                                }
                                                else if (k == MAX_K && i == NM && j == 0) {
                                                    pi[xH][xL][yH][yL][k][i][j] = (lamH * pi[xH - 1][xL][yH][yL][k][i][j] +
                                                        (N - k - (i - 1) - j) * r * (1.0 - p) * pi[xH][xL][yH][yL][k][i - 1][j] +
                                                        (N - (k - 1) - i - j) * r * p * pi[xH][xL][yH][yL][k - 1][i][j]) /
                                                        (musH);
                                                }
                                                // 30 end
                                            }
                                        }
                                    }
                                }

                            }
                        }
                    }
                }
            }
        }
    }

    return;
}

void Normalize() {
    sum = 0;
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int yH = 0; yH <= B; yH = yH + B) {
                for (int yL = 0; yL <= B; yL = yL + B) {

                    for (int k = 0; k <= MAX_K; k++) {
                        for (int i = 0; i <= NM + 1; i++) {
                            for (int j = 0; j <= NM + 1 - i; j++) {
                                sum += pi[xH][xL][yH][yL][k][i][j];
                            }
                        }
                    }
                }
            }
        }
    }

    normalize = 0;
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int yH = 0; yH <= B; yH = yH + B) {
                for (int yL = 0; yL <= B; yL = yL + B) {

                    for (int k = 0; k <= MAX_K; k++) {
                        for (int i = 0; i <= NM + 1; i++) {
                            for (int j = 0; j <= NM + 1 - i; j++) {
                                normalize += pi[xH][xL][yH][yL][k][i][j];
                            }
                        }
                    }
                }
            }
        }
    }
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int yH = 0; yH <= B; yH = yH + B) {
                for (int yL = 0; yL <= B; yL = yL + B) {

                    for (int k = 0; k <= MAX_K; k++) {
                        for (int i = 0; i <= NM + 1; i++) {
                            for (int j = 0; j <= NM + 1 - i; j++) {
                                pi[xH][xL][yH][yL][k][i][j] = pi[xH][xL][yH][yL][k][i][j] / normalize;
                            }
                        }
                    }
                }
            }
        }
    }
    // step 4: calc. error;
    error = 0;
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int yH = 0; yH <= B; yH = yH + B) {
                for (int yL = 0; yL <= B; yL = yL + B) {

                    for (int k = 0; k <= MAX_K; k++) {
                        for (int i = 0; i <= NM + 1; i++) {
                            for (int j = 0; j <= NM + 1 - i; j++) {
                                error += pow(old_pi[xH][xL][yH][yL][k][i][j] - pi[xH][xL][yH][yL][k][i][j], 2);
                                old_pi[xH][xL][yH][yL][k][i][j] = pi[xH][xL][yH][yL][k][i][j];
                            }
                        }
                    }
                }
            }
        }
    }
    convergence = sqrt(error);
    return;
}

void Calculation() {
    // Performance measure

    /////////// Lq //////////
    // Queue Len all
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int yH = 0; yH <= B; yH = yH + B) {
                for (int yL = 0; yL <= B; yL = yL + B) {

                    for (int k = 0; k <= MAX_K; k++) {
                        for (int i = 0; i <= NM + 1; i++) {
                            for (int j = 0; j <= NM + 1 - i; j++) {
                                Lq_T += (xH + xL) * pi[xH][xL][yH][yL][k][i][j];
                            }
                        }
                    }
                }
            }
        }
    }
    // Queue Len High
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int yH = 0; yH <= B; yH = yH + B) {
                for (int yL = 0; yL <= B; yL = yL + B) {

                    for (int k = 0; k <= MAX_K; k++) {
                        for (int i = 0; i <= NM + 1; i++) {
                            for (int j = 0; j <= NM + 1 - i; j++) {
                                Lq_H += (xH)*pi[xH][xL][yH][yL][k][i][j];
                            }
                        }
                    }
                }
            }
        }
    }
    // Queue Len Low
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int yH = 0; yH <= B; yH = yH + B) {
                for (int yL = 0; yL <= B; yL = yL + B) {

                    for (int k = 0; k <= MAX_K; k++) {
                        for (int i = 0; i <= NM + 1; i++) {
                            for (int j = 0; j <= NM + 1 - i; j++) {
                                Lq_L += (xL)*pi[xH][xL][yH][yL][k][i][j];
                            }
                        }
                    }
                }
            }
        }
    }

    /////////// Lbl //////////
    // block len all
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int k = 0; k <= MAX_K; k++) {
                for (int i = 0; i <= NM + 1; i++) {
                    for (int j = 0; j <= NM + 1 - i; j++) {
                        Lb_T += (B)*pi[xH][xL][B][0][k][i][j];
                    }
                }
            }
        }
    }
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int k = 0; k <= MAX_K; k++) {
                for (int i = 0; i <= NM + 1; i++) {
                    for (int j = 0; j <= NM + 1 - i; j++) {
                        Lb_T += (B)*pi[xH][xL][0][B][k][i][j];
                    }
                }
            }

        }
    }
    // block len High
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int k = 0; k <= MAX_K; k++) {
                for (int i = 0; i <= NM + 1; i++) {
                    for (int j = 0; j <= NM + 1 - i; j++) {
                        Lb_H += (B)*pi[xH][xL][B][0][k][i][j];
                    }
                }
            }
        }
    }
    // block len Low
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int k = 0; k <= MAX_K; k++) {
                for (int i = 0; i <= NM + 1; i++) {
                    for (int j = 0; j <= NM + 1 - i; j++) {
                        Lb_L += (B)*pi[xH][xL][0][B][k][i][j];
                    }
                }
            }

        }
    }

    /////////// Le //////////
    // len Total
    Le_T = Lq_T + Lb_T;
    // len High
    Le_H = Lq_H + Lb_H;
    // len Low
    Le_L = Lq_L + Lb_L;

    /////////// Nb //////////
    // block num all
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int k = 0; k <= MAX_K; k++) {
                for (int i = 0; i <= NM + 1; i++) {
                    for (int j = 0; j <= NM + 1 - i; j++) {
                        Nb_T += pi[xH][xL][B][0][k][i][j];
                    }
                }
            }
        }
    }
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int k = 0; k <= MAX_K; k++) {
                for (int i = 0; i <= NM + 1; i++) {
                    for (int j = 0; j <= NM + 1 - i; j++) {
                        Nb_T += pi[xH][xL][0][B][k][i][j];
                    }
                }
            }

        }
    }
    // block num High
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int k = 0; k <= MAX_K; k++) {
                for (int i = 0; i <= NM + 1; i++) {
                    for (int j = 0; j <= NM + 1 - i; j++) {
                        Nb_H += pi[xH][xL][B][0][k][i][j];
                    }
                }
            }
        }
    }
    // block num Low
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int k = 0; k <= MAX_K; k++) {
                for (int i = 0; i <= NM + 1; i++) {
                    for (int j = 0; j <= NM + 1 - i; j++) {
                        Nb_L += pi[xH][xL][0][B][k][i][j];
                    }
                }
            }

        }
    }

    /////////// Pb //////////
    // Blocking prob. All = High = Low
    for (int xH = 0; xH <= Q_MAX; xH++) {
        int xL_MAX1 = Q_MAX - xH;
        for (int k = 0; k <= MAX_K; k++) {
            for (int i = 0; i <= NM + 1; i++) {
                for (int j = 0; j <= NM + 1 - i; j++) {
                    //if (pi[xH][xL_MAX1][B][0][k][i][j] != 0) {
                    //    printf("(%d, %d, B, 0, %d, %d, %d) = %f\n", xH, xL_MAX1, k, i, j, pi[xH][xL_MAX1][B][0][k][i][j]);
                    //}
                    Pb_T += pi[xH][xL_MAX1][B][0][k][i][j];
                }
            }
        }
    }
    for (int xH = 0; xH <= Q_MAX; xH++) {
        int xL_MAX2 = Q_MAX - xH;
        for (int k = 0; k <= MAX_K; k++) {
            for (int i = 0; i <= NM + 1; i++) {
                for (int j = 0; j <= NM + 1 - i; j++) {
                    Pb_T += pi[xH][xL_MAX2][0][B][k][i][j];
                }
            }
        }
    }
    Pb_H = Pb_T;
    Pb_L = Pb_T;

    /////////// Lvs //////////
    // Successful len Total
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int i = 0; i <= NM; i++) {
                for (int j = 0; j <= NM - i; j++) {
                    Lvs_T += (B)*pi[xH][xL][B][0][MAX_K][i][j];
                }
            }
        }
    }
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int i = 0; i <= NM; i++) {
                for (int j = 0; j <= NM - i; j++) {
                    Lvs_T += (B)*pi[xH][xL][0][B][MAX_K][i][j];

                }
            }

        }
    }
    // Successful len High
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int i = 0; i <= NM; i++) {
                for (int j = 0; j <= NM - i; j++) {
                    Lvs_H += (B)*pi[xH][xL][B][0][MAX_K][i][j];
                }
            }
        }
    }
    // Successful len Low
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int i = 0; i <= NM; i++) {
                for (int j = 0; j <= NM - i; j++) {
                    Lvs_L += (B)*pi[xH][xL][0][B][MAX_K][i][j];

                }
            }

        }
    }

    /////////// Lvf //////////
    // Failed len Total
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int k = 0; k <= MAX_K - 1; k++) {
                for (int i = 0; i <= NM + 1; i++) {
                    Lvf_T += (B)*pi[xH][xL][B][0][k][i][NM + 1 - i];
                }
            }
        }
    }
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int k = 0; k <= MAX_K - 1; k++) {
                for (int i = 0; i <= NM + 1; i++) {
                    Lvf_T += (B)*pi[xH][xL][0][B][k][i][NM + 1 - i];
                }
            }

        }
    }
    // Failed len High
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int k = 0; k <= MAX_K - 1; k++) {
                for (int i = 0; i <= NM + 1; i++) {
                    Lvf_H += (B)*pi[xH][xL][B][0][k][i][NM + 1 - i];
                }
            }
        }
    }
    // Failed len Low
    for (int xH = 0; xH <= Q_MAX; xH++) {
        for (int xL = 0; xL <= Q_MAX; xL++) {

            for (int k = 0; k <= MAX_K - 1; k++) {
                for (int i = 0; i <= NM + 1; i++) {
                    Lvf_L += (B)*pi[xH][xL][0][B][k][i][NM + 1 - i];
                }
            }

        }
    }

    /////////// Th //////////
    // Throughput Total
    Th_T = (lamH + lamL) * (1 - Pb_T);
    // Throughput High
    Th_H = (lamH) * (1 - Pb_H);
    // Throughput Low
    Th_L = (lamL) * (1 - Pb_L);

    /////////// Ths //////////
    // Successful Throughput High
    Ths_H = (musH)*Lvs_H;
    // Successful Throughput Low
    Ths_L = (musL)*Lvs_L;
    // Successful Throughput Total
    Ths_T = Ths_H + Ths_L;

    /////////// Thf //////////
    // Failed Throughput High
    Thf_H = (mufH)*Lvf_H;
    // Failed Throughput Low
    Thf_L = (mufL)*Lvf_L;
    // Failed Throughput Total
    Thf_T = Thf_H + Thf_L;

    /////////// W //////////
    W_T = Le_T / Th_T;
    W_H = Le_H / Th_H;
    W_L = Le_L / Th_L;

    /////////// Wq //////////
    Wq_T = Lq_T / Th_T;
    Wq_H = Lq_H / Th_H;
    Wq_L = Lq_L / Th_L;

    /////////// Wb //////////
    Wb_T = Lb_T / Th_T;
    Wb_H = Lb_H / Th_H;
    Wb_L = Lb_L / Th_L;

    /////////// Rs //////////
    Rs_T = Ths_T / (Ths_T + Thf_T);
    Rs_H = Ths_H / (Ths_H + Thf_H);
    Rs_L = Ths_L / (Ths_L + Thf_L);

    return;
}

int main() {
    csv1 csv1("Model 4 analytical report.csv");

    csv1.clear_file();

    vector<string> analytical_columns = {
        "lamH", "lamL", "musH", "musL","mufH", "mufL", "theta", "epsilon", "r", "p", "B", "NM",
        "Le.T_ana", "Lq.T_ana", "Lb.T_ana",  "W.T_ana",   "Wq.T_ana", "Wb.T_ana",
        "Pb.T_ana", "Th.T_ana", "Ths.T_ana", "Thf.T_ana", "Rs.T_ana", "Nb.T_ana",
        "Le.H_ana", "Lq.H_ana", "Lb.H_ana",  "W.H_ana",   "Wq.H_ana", "Wb.H_ana",
        "Pb.H_ana", "Th.H_ana", "Ths.H_ana", "Thf.H_ana", "Rs.H_ana", "Nb.H_ana",
        "Le.L_ana", "Lq.L_ana", "Lb.L_ana",  "W.L_ana",   "Wq.L_ana", "Wb.L_ana",
        "Pb.L_ana", "Th.L_ana", "Ths.L_ana", "Thf.L_ana", "Rs.L_ana", "Nb.L_ana" };

    int pointer = 1;
    while (pointer <= 1) {
        cout << "Starting case " << pointer << endl;
        csv1.print_title(pointer, analytical_columns);

        for (int aa = 0; aa < 5; aa++) {
            Initialization(pointer, aa);
            while (convergence > 1.0E-8) {
                Iteration();
                Normalize();
            }

            Calculation();
            vector<double> analytical_data = {
                lamH,  lamL,    musH,  musL,  mufH,  mufL,
                theta, epsilon, r,     p,     B,     NM,
                Le_T,  Lq_T,    Lb_T,  W_T,   Wq_T,  Wb_T,
                Pb_T,  Th_T,    Ths_T, Thf_T, Rs_T,  Nb_T,
                Le_H,  Lq_H,    Lb_H,  W_H,   Wq_H,  Wb_H,
                Pb_H,  Th_H,    Ths_H, Thf_H, Rs_H,  Nb_H,
                Le_L,  Lq_L,    Lb_L,  W_L,   Wq_L,  Wb_L,
                Pb_L,  Th_L,    Ths_L, Thf_L, Rs_L,  Nb_L };
            csv1.print_data(analytical_data);
            cout << endl;
            printf("%f\n", sum); // 正常要是1.0
        }
        pointer++;
    }

    return 0;
}