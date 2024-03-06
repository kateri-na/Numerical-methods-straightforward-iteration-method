namespace ЧМ_лаба_3
{
	public class Program
	{
		static void Main(string[] args)
		{
			//StraightforwardIterations sfI = new StraightforwardIterations("C:\\Users\\user\\source\\repos\\ЧМ лаба 3\\data.txt");
			//sfI.Solution();
			Console.WriteLine("Введите требуемую точность:");
			double eps = Convert.ToDouble(Console.ReadLine());
			StraightforwardIterations sfI = new StraightforwardIterations(5, -10, 10, eps);
			sfI.Solution();
			Console.ReadLine();
		}
	}
}