Необходимо реализовать метод прямых итераций c исчерпыванием определения пары с третьим максимальным по модулю собственным значением симметричной матрицы простой структуры. 

	Входные параметры основной процедуры: 

	N – размерность матрицы; 
 
	A – двумерный массив размерности N × N; 
 
	λ_n – максимальное по модулю собственное значение матрицы A; 
 
	x_n – собственный вектор, соответствующий собственному значению λ_n; 
 
	λ_n1 – второе максимальное по модулю собственное значение; 
 
	x_n1 – собственный вектор, соответствующий собственному значению λ_n1;
 
	Выходные параметры основной процедуры: 
 
	λ – третье максимальное по модулю собственное значение; 
 
	x – третий собственный вектор; 
 
	K – число выполненных итераций; 
 
	r – мера точности полученной пары (λ, x).
