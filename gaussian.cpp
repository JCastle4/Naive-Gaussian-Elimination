#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>

void printMatrix(std::vector<std::vector<long double>> matrix, std::vector<long double> vec);

void FwdElimination(std::vector<std::vector<long double>>& matrix, std::vector<long double>& vec);
void BackSubst(std::vector<std::vector<long double>>& matrix, std::vector<long double>& vec, std::vector<long double>& sol);
void NaiveGaussian(std::vector<std::vector<long double>>& matrix, std::vector<long double>& vec, std::string fileName);

void SPPfwdElimination(std::vector<std::vector<long double>>& matrix, std::vector<long double>& vec, std::vector<long double>& ind);
void SPPbackSubst(std::vector<std::vector<long double>>& matrix, std::vector<long double>& vec, std::vector<long double>& sol, std::vector<long double>& ind);
void SPPGaussian(std::vector<std::vector<long double>>& matrix, std::vector<long double>& vec, std::string fileName);

int main(int argc, char* argv[])
{
    bool flag = false;
    std::string input = "";
    if (argc <= 1)
    {
        printf("usage: %s ", argv[0]);
    }
    else if (argc == 2) {
        input = argv[1];
    }
    else if (argc >= 3) {
        std::string str = argv[1];
        if (str == "--spp") { flag = true; input = argv[2]; }
        else std::cout << "Invalid flag" << std::endl;
    }

    std::ifstream MyFile(input);
    if (!MyFile.is_open()) {
        std::cout << "File Not Found." << std::endl;
    }
    else {
        input = input.substr(0, input.size() - 4);
        //std::cout << input << std::endl;
        int size = 0;
        MyFile >> size;
        //std::cout << size << std::endl;

        //2D vector with rows and columns of size 'size' and value of 0
        std::vector<std::vector<long double>> matrix(size, std::vector<long double>(size, 0));;
        std::vector<long double> vec(size);

        for (int i = 0; i < vec.size(); i++) {
            for (int j = 0; j < vec.size(); j++) {
                MyFile >> matrix[i][j];
            }
        }
        for (int k = 0; k < size; k++) { MyFile >> vec[k]; }

        if (flag) { SPPGaussian(matrix, vec, input); }
        else {
            NaiveGaussian(matrix, vec, input);
        }
        //printMatrix(matrix, vec);
        MyFile.close();
    }

}

void printMatrix(std::vector<std::vector<long double>> matrix, std::vector<long double> vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << "|";
        for (int j = 0; j < vec.size(); j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << "|" << vec[i] << "\n";
    }
}

void FwdElimination(std::vector<std::vector<long double>>& matrix, std::vector<long double>& vec) {
    for (int k = 0; k < vec.size() - 1; k++) {
        for (int i = k + 1; i < vec.size(); i++) {
            long double mult = matrix[i][k] / matrix[k][k];
            for (int j = k; j < vec.size(); j++) {
                matrix[i][j] = matrix[i][j] - (mult * matrix[k][j]);
            }
            vec[i] = vec[i] - (mult * vec[k]);
        }
    }
}

void BackSubst(std::vector<std::vector<long double>>& matrix, std::vector<long double>& vec, std::vector<long double>& sol) {
    sol[sol.size() - 1] = vec[sol.size() - 1] / matrix[sol.size() - 1][sol.size() - 1];
    for (int i = sol.size() - 2; i >= 0; i--) {
        long double sum = vec[i];
        for (int j = i + 1; j < sol.size(); j++) {
            sum = sum - matrix[i][j] * sol[j];
        }
        sol[i] = sum / matrix[i][i];
    }
}

void NaiveGaussian(std::vector<std::vector<long double>>& matrix, std::vector<long double>& vec, std::string fileName) {
    std::vector<long double> sol(vec.size());
    FwdElimination(matrix, vec);
    BackSubst(matrix, vec, sol);

    std::ofstream newFile;
    newFile.open(fileName + ".sol");
    for (int i = 0; i < sol.size(); i++) { newFile << std::setprecision(20) << sol[i] << " "; }
    newFile.close();
}

void SPPfwdElimination(std::vector<std::vector<long double>>& matrix, std::vector<long double>& vec, std::vector<long double>& ind) {
    std::vector<long double> scaling(vec.size());

    //initialize index and scaling vectors
    for (int i = 0; i < vec.size(); i++) {
        long double smax = 0;

        for (int j = 0; j < vec.size(); j++) {
            if (std::abs(matrix[i][j]) > smax) smax = matrix[i][j];
        }
        scaling[i] = smax;
    }

    for (int k = 0; k < vec.size() - 1; k++) {
        long double rmax = 0;
        int maxInd = k;

        for (int i = k; i < vec.size(); i++) {
            long double r = std::abs(matrix[ind[i]][k] / scaling[ind[i]]);

            if (r > rmax) {
                rmax = r;
                maxInd = i;
            }
        }
        long double temp = ind[maxInd];
        ind[maxInd] = ind[k];
        ind[k] = temp;

        for (int i = k + 1; i < vec.size(); i++) {
            long double mult = matrix[ind[i]][k] / matrix[ind[k]][k];
            for (int j = k; j < vec.size(); j++) {
                matrix[ind[i]][j] = matrix[ind[i]][j] - (mult * matrix[ind[k]][j]);
            }
            vec[ind[i]] = vec[ind[i]] - (mult * vec[ind[k]]);
        }
    }
}

void SPPbackSubst(std::vector<std::vector<long double>>& matrix, std::vector<long double>& vec, std::vector<long double>& sol, std::vector<long double>& ind) {
    sol[vec.size() - 1] = vec[ind[vec.size() - 1]] / matrix[ind[vec.size() - 1]][vec.size() - 1];

    for (int i = vec.size() - 2; i >= 0; i--) {
        long double sum = vec[ind[i]];
        for (int j = i + 1; j < vec.size(); j++) {
            sum = sum - matrix[ind[i]][j] * sol[j];
        }
        sol[i] = sum / matrix[ind[i]][i];
    }
}

void SPPGaussian(std::vector<std::vector<long double>>& matrix, std::vector<long double>& vec, std::string fileName) {
    std::vector<long double> sol(vec.size());
    std::vector<long double> ind(vec.size());

    for (int i = 0; i < vec.size(); i++) {
        ind[i] = i;
    }

    SPPfwdElimination(matrix, vec, ind);
    SPPbackSubst(matrix, vec, sol, ind);

    std::ofstream newFile;
    newFile.open(fileName + ".sol");
    for (int i = 0; i < sol.size(); i++) { newFile << std::setprecision(20) << sol[i] << " "; }
    newFile.close();
}
