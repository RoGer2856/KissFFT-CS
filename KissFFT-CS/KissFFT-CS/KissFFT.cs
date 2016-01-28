using System;
using System.Collections.Generic;
using System.Text;

namespace KissFFT
{
	public class KissFFT<kiss_fft_scalar>
	{
		private const int MaxFactors = 32;
		private int Nfft;
		private bool Inverse;
		private int[] Factors = new int[2 * MaxFactors];
		private kiss_fft_cpx<kiss_fft_scalar>[] Twiddles;

		private IArithmetic<kiss_fft_scalar> A;

		public KissFFT(int nfft, bool inverse, IArithmetic<kiss_fft_scalar> arithmetic)
		{
			this.A = arithmetic;
			kiss_fft_alloc(nfft, inverse);
		}		

		protected void kf_bfly2(Array<kiss_fft_cpx<kiss_fft_scalar>> Fout, int fstride, int m)
		{
			Array<kiss_fft_cpx<kiss_fft_scalar>> Fout2;
			Array<kiss_fft_cpx<kiss_fft_scalar>> tw1 = new Array<kiss_fft_cpx<kiss_fft_scalar>>(Twiddles);
			kiss_fft_cpx<kiss_fft_scalar> t;
			Fout2 = Fout + m;
			do
			{
				Fout[0] = A.FixDivide(Fout[0], 2);
				Fout2[0] = A.FixDivide(Fout2[0], 2);

				t = A.Multiply(Fout2[0], tw1[0]);
				tw1 += fstride;
				Fout2[0] = A.Subtract(Fout[0], t);
				Fout[0] = A.Add(Fout[0], t);
				++Fout2;
				++Fout;
			} while (--m != 0);
		}

		protected void kf_bfly4(Array<kiss_fft_cpx<kiss_fft_scalar>> Fout, int fstride, int m)
		{
			Array<kiss_fft_cpx<kiss_fft_scalar>> tw1, tw2, tw3;
			kiss_fft_cpx<kiss_fft_scalar>[] scratch = new kiss_fft_cpx<kiss_fft_scalar>[6];
			int k = m;
			int m2 = 2 * m;
			int m3 = 3 * m;

			tw1 = new Array<kiss_fft_cpx<kiss_fft_scalar>>(Twiddles);
			tw2 = new Array<kiss_fft_cpx<kiss_fft_scalar>>(Twiddles);
			tw3 = new Array<kiss_fft_cpx<kiss_fft_scalar>>(Twiddles);

			do
			{
				Fout[0] = A.FixDivide(Fout[0], 4);
				Fout[m] = A.FixDivide(Fout[m], 4);
				Fout[m2] = A.FixDivide(Fout[m2], 4);
				Fout[m3] = A.FixDivide(Fout[m3], 4);

				scratch[0] = A.Multiply(Fout[m], tw1[0]);
				scratch[1] = A.Multiply(Fout[m2], tw2[0]);
				scratch[2] = A.Multiply(Fout[m3], tw3[0]);

				scratch[5] = A.Subtract(Fout[0], scratch[1]);
				Fout[0] = A.Add(Fout[0], scratch[1]);
				scratch[3] = A.Add(scratch[0], scratch[2]);
				scratch[4] = A.Subtract(scratch[0], scratch[2]);
				Fout[m2] = A.Subtract(Fout[0], scratch[3]);
				tw1 += fstride;
				tw2 += fstride * 2;
				tw3 += fstride * 3;
				Fout[0] = A.Add(Fout[0], scratch[3]);

				if (Inverse)
				{
					Fout[m].r = A.Subtract(scratch[5].r, scratch[4].i);
					Fout[m].i = A.Add(scratch[5].i, scratch[4].r);
					Fout[m3].r = A.Add(scratch[5].r, scratch[4].i);
					Fout[m3].i = A.Subtract(scratch[5].i, scratch[4].r);
				}
				else
				{
					Fout[m].r = A.Add(scratch[5].r, scratch[4].i);
					Fout[m].i = A.Subtract(scratch[5].i, scratch[4].r);
					Fout[m3].r = A.Subtract(scratch[5].r, scratch[4].i);
					Fout[m3].i = A.Add(scratch[5].i, scratch[4].r);
				}
				++Fout;
			} while (--k != 0);
		}

		protected void kf_bfly3(Array<kiss_fft_cpx<kiss_fft_scalar>> Fout, int fstride, int m)
		{
			int k = m;
			int m2 = 2 * m;
			Array<kiss_fft_cpx<kiss_fft_scalar>> tw1, tw2;
			kiss_fft_cpx<kiss_fft_scalar>[] scratch = new kiss_fft_cpx<kiss_fft_scalar>[5];
			kiss_fft_cpx<kiss_fft_scalar> epi3;
			epi3 = Twiddles[fstride * m];

			tw1 = new Array<kiss_fft_cpx<kiss_fft_scalar>>(Twiddles);
			tw2 = new Array<kiss_fft_cpx<kiss_fft_scalar>>(Twiddles);

			do
			{
				Fout[0] = A.FixDivide(Fout[0], 3);
				Fout[m] = A.FixDivide(Fout[m], 3);
				Fout[m2] = A.FixDivide(Fout[m2], 3);

				scratch[1] = A.Multiply(Fout[m], tw1[0]);
				scratch[2] = A.Multiply(Fout[m2], tw2[0]);

				scratch[3] = A.Add(scratch[1], scratch[2]);
				scratch[0] = A.Subtract(scratch[1], scratch[2]);
				tw1 += fstride;
				tw2 += fstride * 2;

				Fout[m].r = A.Subtract(Fout[0].r, A.Half(scratch[3].r));
				Fout[m].i = A.Subtract(Fout[0].i, A.Half(scratch[3].i));

				scratch[0] = A.Multiply(scratch[0], epi3.i);

				Fout[0] = A.Add(Fout[0], scratch[3]);

				Fout[m2].r = A.Add(Fout[m].r, scratch[0].i);
				Fout[m2].i = A.Subtract(Fout[m].i, scratch[0].r);

				Fout[m].r = A.Subtract(Fout[m].r, scratch[0].i);
				Fout[m].i = A.Add(Fout[m].i, scratch[0].r);

				++Fout;
			} while (--k != 0);
		}

		protected void kf_bfly5(Array<kiss_fft_cpx<kiss_fft_scalar>> Fout, int fstride, int m)
		{
			Array<kiss_fft_cpx<kiss_fft_scalar>> Fout0, Fout1, Fout2, Fout3, Fout4;
			int u;
			kiss_fft_cpx<kiss_fft_scalar>[] scratch = new kiss_fft_cpx<kiss_fft_scalar>[13];
			Array<kiss_fft_cpx<kiss_fft_scalar>> twiddles = new Array<kiss_fft_cpx<kiss_fft_scalar>>(Twiddles);
			Array<kiss_fft_cpx<kiss_fft_scalar>> tw;
			kiss_fft_cpx<kiss_fft_scalar> ya, yb;
			ya = twiddles[fstride * m];
			yb = twiddles[fstride * 2 * m];

			Fout0 = Fout;
			Fout1 = Fout0 + m;
			Fout2 = Fout0 + 2 * m;
			Fout3 = Fout0 + 3 * m;
			Fout4 = Fout0 + 4 * m;

			tw = new Array<kiss_fft_cpx<kiss_fft_scalar>>(twiddles);
			for (u = 0; u < m; ++u)
			{
				Fout[0] = A.FixDivide(Fout[0], 5);
				Fout1[0] = A.FixDivide(Fout1[0], 5);
				Fout2[0] = A.FixDivide(Fout2[0], 5);
				Fout3[0] = A.FixDivide(Fout3[0], 5);
				Fout4[0] = A.FixDivide(Fout4[0], 5);
				scratch[0] = new kiss_fft_cpx<kiss_fft_scalar>(Fout[0]);

				scratch[1] = A.Multiply(Fout1[0], tw[u * fstride]);
				scratch[2] = A.Multiply(Fout2[0], tw[2 * u * fstride]);
				scratch[3] = A.Multiply(Fout3[0], tw[3 * u * fstride]);
				scratch[4] = A.Multiply(Fout4[0], tw[4 * u * fstride]);

				scratch[7] = A.Add(scratch[1], scratch[4]);
				scratch[10] = A.Subtract(scratch[1], scratch[4]);
				scratch[8] = A.Add(scratch[2], scratch[3]);
				scratch[9] = A.Subtract(scratch[2], scratch[3]);

				Fout[0].r = A.Add(Fout[0].r, A.Add(scratch[7].r, scratch[8].r));
				Fout[0].i = A.Add(Fout[0].i, A.Add(scratch[7].i, scratch[8].i));

				scratch[5].r = A.Add(scratch[0].r, A.Add(A.Multiply(scratch[7].r, ya.r), A.Multiply(scratch[8].r, yb.r)));
				scratch[5].i = A.Add(scratch[0].i, A.Add(A.Multiply(scratch[7].i, ya.r), A.Multiply(scratch[8].i, yb.r)));

				scratch[6].r = A.Add(A.Multiply(scratch[10].i, ya.i), A.Multiply(scratch[9].i, yb.i));
				scratch[6].r = A.Negate(A.Add(A.Multiply(scratch[10].r, ya.i), A.Multiply(scratch[9].r, yb.i)));

				Fout1[0] = A.Subtract(scratch[5], scratch[6]);
				Fout4[0] = A.Add(scratch[5], scratch[6]);

				scratch[11].r = A.Add(scratch[0].r, A.Add(A.Multiply(scratch[7].r, yb.r), A.Multiply(scratch[8].r, ya.r)));
				scratch[11].i = A.Add(scratch[0].i, A.Add(A.Multiply(scratch[7].i, yb.r), A.Multiply(scratch[8].i, ya.r)));
				scratch[12].r = A.Subtract(A.Multiply(scratch[9].i, ya.i), A.Multiply(scratch[10].i, yb.i));
				scratch[12].i = A.Subtract(A.Multiply(scratch[10].r, yb.i), A.Multiply(scratch[9].r, ya.i));

				Fout2[0] = A.Add(scratch[11], scratch[12]);
				Fout3[0] = A.Subtract(scratch[11], scratch[12]);

				++Fout0;
				++Fout1;
				++Fout2;
				++Fout3;
				++Fout4;
			}
		}

		/* perform the butterfly for one stage of a mixed radix FFT */
		protected void kf_bfly_generic(Array<kiss_fft_cpx<kiss_fft_scalar>> Fout, int fstride, int m, int p)
		{
			int u, k, q1, q;
			Array<kiss_fft_cpx<kiss_fft_scalar>> twiddles = new Array<kiss_fft_cpx<kiss_fft_scalar>>(Twiddles);
			kiss_fft_cpx<kiss_fft_scalar> t;
			int Norig = Nfft;

			kiss_fft_cpx<kiss_fft_scalar>[] scratch = new kiss_fft_cpx<kiss_fft_scalar>[p];

			for (u = 0; u < m; ++u)
			{
				k = u;
				for (q1 = 0; q1 < p; ++q1)
				{
					scratch[q1] = new kiss_fft_cpx<kiss_fft_scalar>(Fout[k]);
					scratch[q1] = A.FixDivide(scratch[q1], p);
					k += m;
				}

				k = u;
				for (q1 = 0; q1 < p; ++q1)
				{
					int twidx = 0;
					Fout[k] = new kiss_fft_cpx<kiss_fft_scalar>(scratch[0]);
					for (q = 1; q < p; ++q)
					{
						twidx += fstride * k;
						if (twidx >= Norig)
							twidx -= Norig;
						t = A.Multiply(scratch[q], twiddles[twidx]);
						Fout[k] = A.Add(Fout[k], t);
					}
					k += m;
				}
			}
		}

		protected void kf_work(Array<kiss_fft_cpx<kiss_fft_scalar>> Fout, Array<kiss_fft_cpx<kiss_fft_scalar>> f, int fstride, int in_stride, Array<int> factors)
		{
			Array<kiss_fft_cpx<kiss_fft_scalar>> Fout_beg = Fout;
			int p = factors++[0]; /* the radix  */
			int m = factors++[0]; /* stage's fft length/p */
			Array<kiss_fft_cpx<kiss_fft_scalar>> Fout_end = Fout + p * m;

			// OpenMP
			/*
			// use openmp extensions at the 
			// top-level (not recursive)
			if (fstride == 1 && p <= 5)
			{
				int k;

				// execute the p different work units in different threads
#pragma omp parallel for
				for (k = 0; k < p; ++k)
					kf_work(Fout + k * m, f + fstride * in_stride * k, fstride * p, in_stride, factors, st);
				// all threads have joined by this point

				switch (p)
				{
					case 2:
						kf_bfly2(Fout, fstride, st, m);
						break;
					case 3:
						kf_bfly3(Fout, fstride, st, m);
						break;
					case 4:
						kf_bfly4(Fout, fstride, st, m);
						break;
					case 5:
						kf_bfly5(Fout, fstride, st, m);
						break;
					default:
						kf_bfly_generic(Fout, fstride, st, m, p);
						break;
				}
				return;
			}
			*/

			if (m == 1)
			{
				do
				{
					Fout[0] = new kiss_fft_cpx<kiss_fft_scalar>(f[0]);
					f += fstride * in_stride;
				} while (++Fout != Fout_end);
			}
			else
			{
				do
				{
					// recursive call:
					// DFT of size m*p performed by doing
					// p instances of smaller DFTs of size m, 
					// each one takes a decimated version of the input
					kf_work(Fout, f, fstride * p, in_stride, factors);
					f += fstride * in_stride;
				} while ((Fout += m) != Fout_end);
			}

			Fout = Fout_beg;

			// recombine the p smaller DFTs 
			switch (p)
			{
				case 2:
					kf_bfly2(Fout, fstride, m);
					break;
				case 3:
					kf_bfly3(Fout, fstride, m);
					break;
				case 4:
					kf_bfly4(Fout, fstride, m);
					break;
				case 5:
					kf_bfly5(Fout, fstride, m);
					break;
				default:
					kf_bfly_generic(Fout, fstride, m, p);
					break;
			}
		}

		/*  facbuf is populated by p1,m1,p2,m2, ...
			where 
			p[i] * m[i] = m[i-1]
			m0 = n                  */
		protected void kf_factor(int n, Array<int> facbuf)
		{
			int p = 4;
			double floor_sqrt;
			floor_sqrt = Math.Floor(Math.Sqrt((double) n));

			/*factor out powers of 4, powers of 2, then any remaining primes */
			do
			{
				while (n % p != 0)
				{
					switch (p)
					{
						case 4:
							p = 2;
							break;
						case 2:
							p = 3;
							break;
						default:
							p += 2;
							break;
					}
					if (p > floor_sqrt)
						p = n;          /* no more factors, skip to end */
				}
				n /= p;
				facbuf++[0] = p;
				facbuf++[0] = n;
			} while (n > 1);
		}

		/*
		 *
		 * User-callable function to allocate all necessary storage space for the fft.
		 *
		 * The return value is a contiguous block of memory, allocated with malloc.  As such,
		 * It can be freed with free(), rather than a kiss_fft-specific function.
		 * */
		protected void kiss_fft_alloc(int nfft, bool inverse_fft)
		{
			Twiddles = new kiss_fft_cpx<kiss_fft_scalar>[nfft];
			
			int i;
			Nfft=nfft;
			Inverse = inverse_fft;

			for (i=0;i<nfft;++i) {
				const double pi = 3.141592653589793238462643383279502884197169399375105820974944;
	double phase = -2 * pi * i / nfft;
				if (Inverse)
					phase *= -1;

				Twiddles[i] = A.Exp(phase);
			}


			kf_factor(nfft, new Array<int>(Factors));
		}


		public void kiss_fft_stride(Array<kiss_fft_cpx<kiss_fft_scalar>> fin, Array<kiss_fft_cpx<kiss_fft_scalar>> fout, int in_stride)
		{
			if (fin == fout)
			{
				//NOTE: this is not really an in-place FFT algorithm.
				//It just performs an out-of-place FFT into a temp buffer
				kiss_fft_cpx<kiss_fft_scalar>[] tmpbuf = new kiss_fft_cpx<kiss_fft_scalar>[Nfft];
				kf_work(new Array<kiss_fft_cpx<kiss_fft_scalar>>(tmpbuf), fin, 1, in_stride, new Array<int>(Factors));
				for (int i = 0; i < Nfft; ++i)
				{
					fout[i] = new kiss_fft_cpx<kiss_fft_scalar>(tmpbuf[i]);
				}
			}
			else
			{
				kf_work(fout, fin, 1, in_stride, new Array<int>(Factors));
			}
		}

		public void kiss_fft(Array<kiss_fft_cpx<kiss_fft_scalar>> fin, Array<kiss_fft_cpx<kiss_fft_scalar>> fout)
		{
			kiss_fft_stride(fin, fout, 1);
		}

		public static int kiss_fft_next_fast_size(int n)
		{
			while (true)
			{
				int m = n;
				while ((m % 2) == 0)
					m /= 2;
				while ((m % 3) == 0)
					m /= 3;
				while ((m % 5) == 0)
					m /= 5;
				if (m <= 1)
					break; /* n is completely factorable by twos, threes, and fives */
				n++;
			}
			return n;
		}
	}
}
