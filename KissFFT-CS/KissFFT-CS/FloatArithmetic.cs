using System;

namespace KissFFT
{
	// TODO(@eldor): implement IntegerArithmetic
	public class FloatArithmetic : IArithmetic<float>
	{
		public float Negate(float a)
		{
			return -a;
		}

		public float Add(float a, float b)
		{
			return a + b;
		}

		public float Subtract(float a, float b)
		{
			return a - b;
		}

		public float Multiply(float a, float b)
		{
			return a * b;
		}

		public float Divide(float a, float b)
		{
			return a / b;
		}

		public float Half(float a)
		{
			return a * 0.5f;
		}

		public kiss_fft_cpx<float> Add(kiss_fft_cpx<float> a, kiss_fft_cpx<float> b)
		{
			return new kiss_fft_cpx<float>(a.r + b.r, a.i + b.i);
		}

		public kiss_fft_cpx<float> Subtract(kiss_fft_cpx<float> a, kiss_fft_cpx<float> b)
		{
			return new kiss_fft_cpx<float>(a.r - b.r, a.i - b.i);
		}

		public kiss_fft_cpx<float> Multiply(kiss_fft_cpx<float> a, kiss_fft_cpx<float> b)
		{
			return new kiss_fft_cpx<float>(
				a.r * b.r - a.i * b.i,
				a.r * b.i + a.i * b.r
			);
		}

		public kiss_fft_cpx<float> Multiply(kiss_fft_cpx<float> a, float b)
		{
			return new kiss_fft_cpx<float>(a.r * b, a.i * b);
		}

		public kiss_fft_cpx<float> FixDivide(kiss_fft_cpx<float> a, int b)
		{
			return a;
		}

		public kiss_fft_cpx<float> Exp(double phase)
		{
			return new kiss_fft_cpx<float>((float) Math.Cos(phase), (float) Math.Sin(phase));
		}
	}
}
