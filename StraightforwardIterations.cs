using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ЧМ_лаба_3
{
    internal class StraightforwardIterations
    {
        private int N;
        private double[,] A;
        private double L_n;
        private double[] x_n;
        private double L_n1;
        private double[] x_n1;
        private const int IterationsCount = 1000;

        //private double L;
        private double[] x_prev;
        private double[] x;
        private int K;

        private double[,] A_1;
        private double[] v_prev;
        private double[] v;
        private double G_prev;
        private double G;
        private double[,] H;
        private double[,] Y;

        //private double E = 1e-5;
        private double E;
        private double AccuracyVectors;
        private double AccuracyEigens;
        private double r; //мера точности

        public double getAccuracyVectors()
        {
            return AccuracyVectors;
        }
        public double getAccuracyEigens()
        {
            return AccuracyEigens;
        }
        public double getMeasureOfAccuracy()
        {
            return r;
        }
        public int getIterationsCount()
        {
            return K;
        }
        public StraightforwardIterations(String fileName, double minValue, double maxValue)
        {
            StreamReader sr = new StreamReader(fileName);
            string? readNumber = sr.ReadLine();
            N = Convert.ToInt32(readNumber);
            A = new double[N, N];
            x_n = new double[N];
            x_n1 = new double[N];
            A_1 = new double[N, N];
            x_prev = new double[N];
            x = new double[N];
            v_prev = new double[N];
            v = new double[N];
            ReadMatrixFromFile(sr);
            readNumber = sr.ReadLine();
            L_n = Convert.ToDouble(readNumber);
            ReadArrayFromFile(sr, x_n);
            readNumber = sr.ReadLine();
            L_n1 = Convert.ToDouble(readNumber);
            ReadArrayFromFile(sr, x_n1);
            x = GenerateRandomVector(minValue, maxValue);
        }
        public StraightforwardIterations(int N, double minValue, double MaxValue, double eps)
        {
            this.N = N;
            A = new double[N, N];
            A_1 = new double[N, N];
            H = new double[N, N];
            x_n = new double[N];
            x_n1 = new double[N];
            x_prev = new double[N];
            x = new double[N];
            v_prev = new double[N];
            v = new double[N];
            E = eps;
            CalculateA(minValue, MaxValue);
            GetFirstAndSecondEigenValuesAndVectors();
            x = GenerateRandomVector(minValue, MaxValue);
        }

        public void Solution()
        {
            CalculateA_1();
            //Console.WriteLine();
            //Console.WriteLine("A_1");
            //PrintMatrix(A_1);
            v = VectorNormalazing(x);
            x = MultiplyMatrixOnVector(A_1, v);
            G = MultiplayRowOnColumn(x, v);
            ++K;
            do
            {
                AccuracyVectors = AccuracyEigenVectors();

                v_prev = v;
                x_prev = x;
                G_prev = G;

                v = VectorNormalazing(x);
                x = MultiplyMatrixOnVector(A_1, v);
                G = MultiplayRowOnColumn(x, v);

                ++K;
                AccuracyEigens = AccuracyEigenValues();
            } while ((AccuracyEigenVectors() >= E || AccuracyEigenValues() >= E) && K < IterationsCount);

            r = Calculate_r();

            //Console.WriteLine();
            //Console.WriteLine("Точность собственных значений:");
            //Console.WriteLine(AccuracyEigens);
            //Console.WriteLine("Точность собственных векторов:");
            //Console.WriteLine(AccuracyVectors);
            Console.WriteLine("v");
            PrintVector(v);
            Console.WriteLine("Третье максимальное по модулю собственное значение");
            Console.WriteLine(G);
            Console.WriteLine("Третий собственный вектор:");
            PrintVector(x);
            Console.WriteLine("Мера точности полученной пары - r:");
            Console.WriteLine(r);
            Console.WriteLine("Число выполненных итераций:");
            Console.WriteLine(K);
        }

        public void SolutionTest()
        {
            CalculateA_1();
            v = VectorNormalazing(x);
            x = MultiplyMatrixOnVector(A_1, v);
            G = MultiplayRowOnColumn(x, v);
            ++K;
            do
            {
                AccuracyVectors = AccuracyEigenVectors();

                v_prev = v;
                x_prev = x;
                G_prev = G;

                v = VectorNormalazing(x);
                x = MultiplyMatrixOnVector(A_1, v);
                G = MultiplayRowOnColumn(x, v);

                ++K;
                AccuracyEigens = AccuracyEigenValues();
            } while ((AccuracyEigenVectors() >= E || AccuracyEigenValues() >= E) && K < IterationsCount);

            r = Calculate_r();
        }
        private void CalculateA(double minValue, double maxValue) //
        {
            H = CreateHouseholderMatrix(minValue, maxValue);
            Y = GenerateEigenValuesMatrix(minValue, maxValue);
            Console.WriteLine("Householder matrix:");
            PrintMatrix(H);
            Console.WriteLine();
            Console.WriteLine("Matrix with eigen values on main diagonal:");
            PrintMatrix(Y);
            Console.WriteLine();
            A = MultiplyMatrixOnMatrix(H, Y);
            A = MultiplyMatrixOnMatrix(A, H);
            Console.WriteLine("A");
            PrintMatrix(A);
            Console.WriteLine();
        }
        private void CalculateA_1() //
        {
            A_1 = MatrixDifference(A, MultiplyColumnOnRowAndCoeff(x_n, x_n, L_n));
            A_1 = MatrixDifference(A_1, MultiplyColumnOnRowAndCoeff(x_n1, x_n1, L_n1));
        }
        private double AccuracyEigenValues()
        {
            return Math.Abs(G - G_prev);
        }
        private double AccuracyEigenVectors()
        {
            double numerator = MultiplayRowOnColumn(v, v_prev);
            double denumerator = CalculateVectorLen(v) * CalculateVectorLen(v_prev);
            double cosAngle = numerator / denumerator;
            return Math.Acos(cosAngle);
        }

        //расчет меры точности
        private double Calculate_r()
        {
            double[] tmp_x = new double[N];
            tmp_x = VectorDifference(MultiplyMatrixOnVector(A, x), MultiplyVectorOnCoeff(x, G));
            //return CalculateFirstNormOfVector(tmp_x);
            return FindMaxElem(tmp_x);
        }

        //private double CalculateFirstNormOfVector(double[] vector)
        //{
        //    double result = 0;
        //    for (int i = 0; i < N; ++i)
        //        result += Math.Abs(vector[i]);
        //    return result;
        //}

        private double FindMaxElem(double[] tmp_x)
        {
            double max = 0;
            for (int i = 0; i < N; ++i)
            {
                if (tmp_x[i] > max)
                    max = tmp_x[i];
            }
            return max;
        }
        private double[,] CreateHouseholderMatrix(double minValue, double maxValue) //
        {
            double[] w = GenerateRandomVector(minValue, maxValue);
            w = VectorNormalazing(w);
            double[,] result = MatrixDifference(CreateIdentityMatrix(), MultiplyColumnOnRowAndCoeff(w, w, 2));
            return result;
        }
        private double[,] CreateIdentityMatrix() //
        {
            double[,] matrix = new double[N, N];
            for (int i = 0; i < N; ++i)
            {
                matrix[i, i] = 1;
            }
            return matrix;
        }
        private double[,] GenerateEigenValuesMatrix(double minValue, double maxValue)//
        {
            double[,] result = new double[N, N];
            Random rnd = new Random();
            double d = maxValue - minValue;
            for (int i = 0; i < N; ++i)
            {
                result[i, i] = minValue + rnd.NextDouble() * d;
            }
            return result;

        }
        private double CalculateVectorLen(double[] vector) //
        {
            double result = 0;
            for (int i = 0; i < N; ++i)
            {
                result += Math.Pow(vector[i], 2);
            }
            return Math.Sqrt(result);
        }
        private double[] VectorNormalazing(double[] vector) //
        {
            double[] result = new double[N];
            double vector_len = CalculateVectorLen(vector);
            for (int i = 0; i < N; ++i)
            {
                result[i] = vector[i] / vector_len;
            }
            return result;
        }
        private double[,] MatrixDifference(double[,] matrix1, double[,] matrix2)//
        {
            double[,] result = new double[N, N];
            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < N; ++j)
                {
                    result[i, j] = matrix1[i, j] - matrix2[i, j];
                }
            }
            return result;
        }
        private double[] VectorDifference(double[] vector1, double[] vector2)
        {
            double[] result = new double[N];
            for (int i = 0; i < N; ++i)
            {
                result[i] = vector1[i] - vector2[i];
            }
            return result;
        }
        private double[,] MultiplyColumnOnRowAndCoeff(double[] column, double[] row, double coeff) //
        {
            double[,] result = new double[N, N];
            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < N; ++j)
                {
                    result[i, j] = column[i] * row[j] * coeff;
                }
            }
            return result;
        }
        private double[] MultiplyVectorOnCoeff(double[] vector, double coeff)
        {
            double[] result = new double[N];
            for (int i = 0; i < N; ++i)
            {
                result[i] = vector[i] * coeff;
            }
            return result;
        }
        private double MultiplayRowOnColumn(double[] vector, double[] vectorT) //
        {
            double result = 0;
            for (int i = 0; i < N; ++i)
            {
                result += vector[i] * vectorT[i];
            }
            return result;
        }
        private double[] MultiplyMatrixOnVector(double[,] matrix, double[] vector) //
        {
            double[] result = new double[N];
            for (int i = 0; i < N; ++i)
            {
                result[i] = 0;
                for (int j = 0; j < N; ++j)
                {
                    result[i] += matrix[i, j] * vector[j];
                }
            }
            return result;
        }
        private double[,] MultiplyMatrixOnMatrix(double[,] matrix1, double[,] matrix2) //
        {
            double[,] result = new double[N, N];
            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < N; ++j)
                {
                    for (int k = 0; k < N; ++k)
                    {
                        result[i, j] += matrix1[i, k] * matrix2[k, j];
                    }
                }
            }
            return result;
        }
        private void GetFirstAndSecondEigenValuesAndVectors() //
        {
            int first_max_index = 0, second_max_index = 0;
            double first_max_value = 0, second_max_value = 0;
            for (int i = 0; i < N; ++i)
            {
                if (Math.Abs(Y[i, i]) > Math.Abs(first_max_value))
                {
                    first_max_value = Y[i, i];
                    first_max_index = i;
                }
            }
            for (int i = 0; i < N; ++i)
            {
                if (Math.Abs(Y[i, i]) > Math.Abs(second_max_value) && Math.Abs(Y[i, i]) < Math.Abs(first_max_value))
                {
                    second_max_value = Y[i, i];
                    second_max_index = i;
                }
            }
            L_n = first_max_value;
            L_n1 = second_max_value;
            for (int i = 0; i < N; ++i)
            {
                x_n[i] = H[i, first_max_index];
                x_n1[i] = H[i, second_max_index];
            }
            //Console.WriteLine("L_n");
            //Console.WriteLine(L_n);
            //Console.WriteLine("x_n");
            //PrintVector(x_n);
            //Console.WriteLine("L_n1");
            //Console.WriteLine(L_n1);
            //Console.WriteLine("x_n1");
            //PrintVector(x_n1);
        }
        private double[] GenerateRandomVector(double minValue, double maxValue) //
        {
            double[] result = new double[N];
            Random rnd = new Random();
            double d = maxValue - minValue;
            for (var index = 0; index < N; ++index)
            {
                result[index] = minValue + rnd.NextDouble() * d;
            }
            return result;
        }
        private void ReadArrayFromFile(StreamReader sr, double[] array) //
        {
            string line;
            if ((line = sr.ReadLine()) != null)
            {
                string[] text = line.Split(' ');
                for (int ind = 0; ind < text.Length; ++ind)
                    array[ind] = Convert.ToDouble(text[ind]);
            }
            else
                Console.WriteLine("Can't read array from file");
        }
        private void ReadMatrixFromFile(StreamReader sr) //
        {
            for (int i = 0; i < N; i++)
            {
                string line = sr.ReadLine();
                string[] text = line.Split(' ');
                for (int j = 0; j < N; j++)
                {
                    A[i, j] = Convert.ToDouble(text[j]);
                }
            }
        }
        public void PrintMatrix(double[,] matrix) //
        {
            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < N; ++j)
                {
                    //Console.Write(matrix[i, j]);
                    Console.Write("{0,7} ", string.Format("{0:F5}", matrix[i, j]));
                    Console.Write(' ');
                }
                Console.WriteLine();
            }
        }

        public void PrintVector(double[] vector) //
        {
            for (int i = 0; i < N; ++i)
            {
                Console.Write("{0,7} ", string.Format("{0:F5}", vector[i]));
                Console.Write(' ');
            }
            Console.WriteLine();
        }
    }
}
