#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>

#define matrix std::vector<std::vector<double>>

//*************************************************************************
double polinom(int n, double x, double* k)
//полином вида x^n+k[n-1]*x^(n-1)+k[n-2]*x^(n-2)+...+k[1]*x+k[0]
//со старшим коэффициентом, равным единице
{
  double s = 1;
  for(int i = n - 1; i >= 0; i--)
    s = s * x + k[i];
  return s;
}//polinom

//расчёт корня полинома методом деления пополам отрезка, содержащего этот корень
double dihot(int degree, double edgeNegativ, double edgePositiv, double* kf) {
  for(;;) {//цикл деления отрезка пополам
    double x = 0.5 * (edgeNegativ + edgePositiv);
    if (x == edgeNegativ || x == edgePositiv)return x;
    if (polinom(degree, x, kf) < 0)edgeNegativ = x;
    else edgePositiv = x;
  }//цикл деления отрезка пополам
}//dihot

//один шаг подъёма по лестнице из производных полиномов
void stepUp(
    int level, //индекс достигаемой ступеньки лестницы
    double** A, //вспомогательные массивы
    double** B, //параметров производных полиномов
    int* currentRootsCount //сформированные в вызывающей процедуре
) {
//аргумент полинома, равносильный его бесконечно большому значению
  double major = 0;
  for(int i = 0; i < level; i++) {//формирование major
    double s = fabs(A[level][i]);
    if (s > major)major = s;
  }//формирование major
  major += 1.0;

  currentRootsCount[level] = 0; //рабочая инициализация

//основной цикл поиска корня уравнения
//A[level][0]+A[level][1]*x+...+A[level][level-1]*x^(level-1)+x^level=0
//на очередном интервале монотонности
  for(int i = 0; i <= currentRootsCount[level - 1]; i++) {//очередной интервал монотонности
//signLeft signRight - знаки текущего A-полинома на левой и правой границе интервала монотонности
    int signLeft, signRight;

//предварительная левая и правая границы интервала поиска
    double edgeLeft, edgeRight;

//границы интервала монотонности, несущие информацию о знаке полинома на них
    double edgeNegativ, edgePositiv;

//формирование левой границы поиска
    if (i == 0)edgeLeft = -major;
    else edgeLeft = B[level - 1][i - 1];

//значение текущего A-полинома на левой границе
    double rb = polinom(level, edgeLeft, A[level]);

    if (rb == 0) {//маловероятный случай попадания в корень
      B[level][currentRootsCount[level]] = edgeLeft;
      currentRootsCount[level]++;
      continue;
    }//маловероятный случай попадания в корень

//запомнить знак текущего A-полинома на левой границе
    if (rb > 0)signLeft = 1; else signLeft = -1;

//формирование правой границы поиска
    if (i == currentRootsCount[level - 1])edgeRight = major;
    else edgeRight = B[level - 1][i];

//значение текущего A-полинома на правой границе
    rb = polinom(level, edgeRight, A[level]);

    if (rb == 0) {//маловероятный случай попадания в корень
      B[level][currentRootsCount[level]] = edgeRight;
      currentRootsCount[level]++;
      continue;
    }//маловероятный случай попадания в корень

//запомнить знак текущего A-полинома на правой границе
    if (rb > 0)signRight = 1; else signRight = -1;

//если знаки полинома на границах интервала монотонности совпадают,
//то корня нет
    if (signLeft == signRight)continue;

//теперь можно определить плюс границу и минус границу поиска корня
    if (signLeft < 0) {
      edgeNegativ = edgeLeft;
      edgePositiv = edgeRight;
    }
    else {
      edgeNegativ = edgeRight;
      edgePositiv = edgeLeft;
    }

//всё готово для локализации корня методом деления пополам интервала поиска
    B[level][currentRootsCount[level]] = dihot(level, edgeNegativ, edgePositiv, A[level]);
    currentRootsCount[level]++;
  }//очередной интервал монотонности
  return;
}//stepUp

//процедура находит все вещественные корни полинома любой степени
void polynomRealRoots(
    int n, //степень исходного полинома
    double* kf, //массив коэффициентов исходного полинома
    double* rootsArray, //выходной массив вычисленных корней
    int& rootsCount //количество найденных корней
) {
//используются вспомогательные массивы A и B, имеющие следующее содержание
//A это коэффициенты а B корни производных полиномов
//все A-полиномы нормируются так,
//чтобы коэффициент при старшей степени был равен единице
//A[n] - это массив нормированных коэффициентов исходного полинома
//B[n] - это массив корней исходного полинома
//A[n-1] - это массив нормированных коэффициентов производного полинома
//B[n-1] - это массив корней производного полинома
//аналогичным образом
//A[n-2] и B[n-2] - это коэффициенты и корни дважды производного полинома
//наконец A[1] - это массив коэффициентов последнего полинома
//в цепочке производных полиномов
//это линейный полином и B[1] - это массив его корней,
//представленный единственным значимым элементом

//выделение памяти для вспомогательных массивов
  double** A = new double* [n + 1];
  double** B = new double* [n + 1];

//количество вещественных корней для каждого из A-полиномов
  int* currentRootsCount = new int[n + 1];

  for(int i = 1; i <= n; i++) {
    A[i] = new double[i];
    B[i] = new double[i];
  }

//нормировка исходного полинома
  for(int i = 0; i < n; i++)A[n][i] = kf[i] / kf[n];

//расчёт производных A-полиномов
  for(int i1 = n, i = n - 1; i > 0; i1 = i, i--) {
    for(int j1 = i, j = i - 1; j >= 0; j1 = j, j--) {
      A[i][j] = A[i1][j1] * j1 / i1;
    }
  }

//формирование исходного корня последнего производного полинома
  currentRootsCount[1] = 1;
  B[1][0] = -A[1][0];

//подъём по лестнице производных полиномов
  for(int i = 2; i <= n; i++)stepUp(i, A, B, currentRootsCount);

//формирование результата
  rootsCount = currentRootsCount[n];
  for(int i = 0; i < rootsCount; i++)rootsArray[i] = B[n][i];

//очистка памяти
  for(int i = 1; i <= n; i++) {
    delete[]A[i];
    delete[]B[i];
  }
  delete[]A;
  delete[]B;
  delete[]currentRootsCount;
  return;
}//polynomRealRoots

//*************************************************************************



void PrintMatrix(const matrix m) {
  for(int i = 0; i < m.size(); ++i) {
    for(int j = 0; j < m.size(); ++j) {
      std::cout << m[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

int main() {
  int n;
  std::cin >> n;
  matrix A(n);
  for(int i = 0; i < n; ++i) {
    A[i].resize(n);
    for(int j = 0; j < n; ++j) {
      std::cin >> A[i][j];
    }
  }
  std::vector<int> coeffs;
  coeffs.push_back(n);

///// Danilevski
  for(int i = n - 2; i >= 0; --i) {

    // Swap
    if(std::abs(A[i+1][i]) < 1e-8 ) {
      bool swaped = false;
      for(int j = i-1; j >= 0; --j) {
        if(std::abs(A[i+1][j]) > 1e-8 ) {
          std::swap(A[i], A[j]);
          for(int t = 0; t < n; ++t) {
            std::swap(A[t][i], A[t][j]);
          }
          swaped = true;
          break;
        }
      }
      if(!swaped) {
        coeffs.push_back(i+1);
        --i;
        if(i == -1) continue;
      }
    }

    // Dividing & Multiplication
    double div_mul_value = A[i + 1][i];
    for(int j = 0; j < n; ++j) {
      A[j][i] /= div_mul_value;
      A[i][j] *= div_mul_value;
    }
    // Subtracting
    for(int j = 0; j < n; ++j) {
      if (j != i) {
        double pm_value = A[i + 1][j];
        for(int t = 0; t < n; ++t)
          A[t][j] -= pm_value * A[t][i];
        for(int t = 0; t < n; ++t) {
          A[i][t] += pm_value * A[j][t];
        }
      }
    }
  }
  PrintMatrix(A);
  coeffs.push_back(0);
  std::reverse(coeffs.begin(), coeffs.end());

  std::vector<std::vector<double>> polynoms(coeffs.size()-1);


  for(int i = 0; i < coeffs.size()-1; ++i) {
    for(int j = coeffs[i]; j < coeffs[i+1]; ++j) {
      polynoms[i].push_back(A[coeffs[i]][j]);
    }
  }

///// SOLVING EQUATION
  std::map<double, int> roots;
  //Computing roots
  for(auto & polynom : polynoms) {
    int end_index = polynom.size()-1;
    while(fabs(polynom[end_index]) < 1e-8) {
      ++roots[0];
      --end_index;
    }

    int polynom_order = end_index+1;
    double polynom_coeffs[polynom_order+1];
    polynom_coeffs[polynom_order] = 1;
    for(int j = 0; j < polynom_order; ++j) {
      polynom_coeffs[j] = -polynom[polynom_order - j - 1];
    }

    double f_roots[polynom_order];
    int roots_number = 0;
    polynomRealRoots(polynom_order, polynom_coeffs, f_roots, roots_number);

    for(int i = 0; i < roots_number; ++i) {
      ++roots[f_roots[i]];
    }
  }
  for(const auto& elem : roots) {
    std::cout << "Root - " << elem.first << ", d = " << elem.second << std::endl;
  }

  return 0;
}