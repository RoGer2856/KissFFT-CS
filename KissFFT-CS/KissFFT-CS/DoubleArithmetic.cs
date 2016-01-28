using System;

namespace KissFFT
{
	// TODO(@eldor): implement IntegerArithmetic
	public class DoubleArithmetic : IArithmetic<double>
	{
		public double Negate(double a)
		{
			return -a;
		}

		public double Add(double a, double b)
		{
			return a + b;
		}

		public double Subtract(double a, double b)
		{
			return a - b;
		}

		public double Multiply(double a, double b)
		{
			return a * b;
		}

		public double Divide(double a, double b)
		{
			return a / b;
		}

		public double Half(double a)
		{
			return a * 0.5f;
		}

		public kiss_fft_cpx<double> Add(kiss_fft_cpx<double> a, kiss_fft_cpx<double> b)
		{
			return new kiss_fft_cpx<double>(a.r + b.r, a.i + b.i);
		}

		public kiss_fft_cpx<double> Subtract(kiss_fft_cpx<double> a, kiss_fft_cpx<double> b)
		{
			return new kiss_fft_cpx<double>(a.r - b.r, a.i - b.i);
		}

		public kiss_fft_cpx<double> Multiply(kiss_fft_cpx<double> a, kiss_fft_cpx<double> b)
		{
			return new kiss_fft_cpx<double>(
				a.r * b.r - a.i * b.i,
				a.r * b.i + a.i * b.r
			);
		}

		public kiss_fft_cpx<double> Multiply(kiss_fft_cpx<double> a, double b)
		{
			return new kiss_fft_cpx<double>(a.r * b, a.i * b);
		}

		public kiss_fft_cpx<double> FixDivide(kiss_fft_cpx<double> a, int b)
		{
			return a;
		}

		public kiss_fft_cpx<double> Exp(double phase)
		{
			return new kiss_fft_cpx<double>(Math.Cos(phase), Math.Sin(phase));
		}
	}
}
