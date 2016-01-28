using System;
using System.Collections.Generic;
using System.Text;

namespace KissFFT
{
	public class kiss_fft_cpx<kiss_fft_scalar>
	{
		public kiss_fft_cpx()
		{
		}

		public kiss_fft_cpx(kiss_fft_cpx<kiss_fft_scalar> other)
		{
			r = other.r;
			i = other.i;
		}

		public kiss_fft_cpx(kiss_fft_scalar r, kiss_fft_scalar i)
		{
			this.r = r;
			this.i = i;
		}

		public kiss_fft_scalar r;
		public kiss_fft_scalar i;
	}
}
