namespace KissFFT
{
	public interface IArithmetic<kiss_fft_scalar>
	{
		kiss_fft_scalar Negate(kiss_fft_scalar a);
		kiss_fft_scalar Add(kiss_fft_scalar a, kiss_fft_scalar b);
		kiss_fft_scalar Subtract(kiss_fft_scalar a, kiss_fft_scalar b);
		kiss_fft_scalar Multiply(kiss_fft_scalar a, kiss_fft_scalar b);
		kiss_fft_scalar Divide(kiss_fft_scalar a, kiss_fft_scalar b);
		kiss_fft_scalar Half(kiss_fft_scalar a);

		/*
		Requirements of complex arithmetic functions:

		Add: result = a + b
		Subtract: result = a - b
		Multiplication: result = a * b
		fixDivide: if a fixed point impl, result = c / div, noop otherwise
		Exp: result = (Cos(p), Sin(i))
		*/
		kiss_fft_cpx<kiss_fft_scalar> Add(kiss_fft_cpx<kiss_fft_scalar> a, kiss_fft_cpx<kiss_fft_scalar> b);
		kiss_fft_cpx<kiss_fft_scalar> Subtract(kiss_fft_cpx<kiss_fft_scalar> a, kiss_fft_cpx<kiss_fft_scalar> b);
		kiss_fft_cpx<kiss_fft_scalar> Multiply(kiss_fft_cpx<kiss_fft_scalar> a, kiss_fft_cpx<kiss_fft_scalar> b);
		kiss_fft_cpx<kiss_fft_scalar> Multiply(kiss_fft_cpx<kiss_fft_scalar> a, kiss_fft_scalar b);
		kiss_fft_cpx<kiss_fft_scalar> FixDivide(kiss_fft_cpx<kiss_fft_scalar> a, int b);
		kiss_fft_cpx<kiss_fft_scalar> Exp(double phase);
	}
}
