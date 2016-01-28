using System;

using KissFFT;
using System.Drawing;

namespace FFT
{
	class Program
	{
		private static void PrintEuler(float r, float i)
		{
			float magnitude = r * r + i * i;
			float phase;
			const float treshold = 0.00001f;
			if (Math.Abs(r) < treshold)
			{
				if (i < 0.0f)
				{
					phase = -(float) Math.PI / 2.0f;
				}
				else
				{
					phase = (float) Math.PI / 2.0f;
				}
			}
			else
			{
				phase = (float) Math.Tan(i / r);
			}
			System.Console.Write("magnitude = {0}, phase = {0}", magnitude, phase);
		}

		private static void TestFFT(string title, kiss_fft_cpx<float>[] timeDomain, kiss_fft_cpx<float>[] frequencyDomain)
		{
			System.Console.WriteLine(title);
			{ // FFT
				KissFFT<float> kissFft = new KissFFT<float>(timeDomain.Length, false, new FloatArithmetic());

				kissFft.kiss_fft(new Array<kiss_fft_cpx<float>>(timeDomain), new Array<kiss_fft_cpx<float>>(frequencyDomain));

				for (int i = 0; i < frequencyDomain.Length; ++i)
				{
					frequencyDomain[i].r /= (float) frequencyDomain.Length;
					frequencyDomain[i].i /= (float) frequencyDomain.Length;
				}

				System.Console.Write("Offset = {0} + {1} * i\n", frequencyDomain[0].r, frequencyDomain[0].i);
				for (int i = 1; i < frequencyDomain.Length / 2 + 1; ++i)
				{
					System.Console.Write("f({0}): {1} + {2} * i, ", i, frequencyDomain[i].r * 2.0f, frequencyDomain[i].i * 2.0f);
					PrintEuler(frequencyDomain[i].r * 2.0f, frequencyDomain[i].i * 2.0f);
					System.Console.WriteLine();
				}
			}
			{ // inverse FFT
				kiss_fft_cpx<float>[] inverted = new kiss_fft_cpx<float>[frequencyDomain.Length];

				KissFFT<float> kissFft = new KissFFT<float>(frequencyDomain.Length, true, new FloatArithmetic());

				kissFft.kiss_fft(new Array<kiss_fft_cpx<float>>(frequencyDomain), new Array<kiss_fft_cpx<float>>(inverted));

				for (int i = 0; i < frequencyDomain.Length; ++i)
				{
					const float treshold = 0.00001f;
					if (!(Math.Abs(timeDomain[i].r - inverted[i].r) < treshold && Math.Abs(timeDomain[i].i - inverted[i].i) < treshold))
					{
						System.Console.Write("{0} ({1}, {2}) != ({3}, {4})\n", i, timeDomain[i].r, timeDomain[i].i, inverted[i].r, inverted[i].i);
					}
				}
				System.Console.WriteLine();
			}
		}

		private static void Test1D()
		{
			// scaling factor is N

			const int N = 16;
			kiss_fft_cpx<float>[] timeDomain = new kiss_fft_cpx<float>[N];
			kiss_fft_cpx<float>[] frequencyDomain = new kiss_fft_cpx<float>[N];

			for (int i = 0; i < N; ++i)
			{
				timeDomain[i] = new kiss_fft_cpx<float>();
				frequencyDomain[i] = new kiss_fft_cpx<float>();
			}

			for (int i = 0; i < N; i++)
			{
				timeDomain[i].r = 0.0f;
				timeDomain[i].i = 0.0f;
			}
			TestFFT("Zeroes", timeDomain, frequencyDomain);

			for (int i = 0; i < N; i++)
			{
				timeDomain[i].r = 1.0f;
				timeDomain[i].i = 0.0f;
			}
			TestFFT("Ones", timeDomain, frequencyDomain);

			for (int i = 0; i < N; i++)
			{
				timeDomain[i].r = (float) Math.Sin(2.0f * Math.PI * 4.0f * i / N);
				timeDomain[i].i = 0.0f;
			}
			TestFFT("SineWave", timeDomain, frequencyDomain);

			for (int i = 0; i < N; i++)
			{
				timeDomain[i].r = (float) Math.Cos(2.0f * Math.PI * 4.0f * i / N);
				timeDomain[i].i = 0.0f;
			}
			TestFFT("CosineWave", timeDomain, frequencyDomain);

			for (int i = 0; i < N; i++)
			{
				timeDomain[i].r = 1.0f + (float) Math.Sin(4.0f * 2.0f * Math.PI * i / N) + (float) Math.Sin(6.0f * 2.0f * Math.PI * i / N);
				timeDomain[i].i = 0;
			}
			TestFFT("sin(4x)+sin(6x)", timeDomain, frequencyDomain);
		}

		private delegate float FilterFunction(int x, int y, int width, int height);
		private static void Filter2D(Bitmap input, ref Bitmap output, FilterFunction filterFunction)
		{
			int[] dims = new int[2] { input.Width, input.Height };

			kiss_fft_cpx<float>[] timeDomain = new kiss_fft_cpx<float>[input.Width * input.Height];
			kiss_fft_cpx<float>[] frequencyDomain = new kiss_fft_cpx<float>[input.Width * input.Height];
			for (int i = 0; i < input.Width * input.Height; ++i)
			{
				timeDomain[i] = new kiss_fft_cpx<float>();
				frequencyDomain[i] = new kiss_fft_cpx<float>();
			}

			{ // FFT
				for (int x = 0; x < input.Width; ++x)
				{
					for (int y = 0; y < input.Height; ++y)
					{
						timeDomain[y * input.Width + x].r = input.GetPixel(x, y).R;
						timeDomain[y * input.Width + x].i = 0.0f;
					}
				}

				KissFFTnd<float> kissFftnd = new KissFFTnd<float>(dims, 2, false, new FloatArithmetic());
				kissFftnd.kiss_fftnd(new Array<kiss_fft_cpx<float>>(timeDomain), new Array<kiss_fft_cpx<float>>(frequencyDomain));
				for (int i = 0; i < input.Width * input.Height; ++i)
				{
					frequencyDomain[i].r /= input.Width * input.Height;
					frequencyDomain[i].i /= input.Width * input.Height;
				}
			}

			{ // filter
				for (int x = 0; x < input.Width; ++x)
				{
					for (int y = 0; y < input.Height; ++y)
					{
						float q = filterFunction(x, y, input.Width, input.Height);
						frequencyDomain[y * input.Width + x].r *= q;
						frequencyDomain[y * input.Width + x].i *= q;
					}
				}
			}

			{ // IFFT
				output = new Bitmap(input.Width, input.Height);

				KissFFTnd<float> kissFftnd = new KissFFTnd<float>(dims, 2, true, new FloatArithmetic());
				kissFftnd.kiss_fftnd(new Array<kiss_fft_cpx<float>>(frequencyDomain), new Array<kiss_fft_cpx<float>>(timeDomain));
				for (int x = 0; x < input.Width; ++x)
				{
					for (int y = 0; y < input.Height; ++y)
					{
						int intensity = (int) timeDomain[y * input.Width + x].r;
						intensity = Math.Min(255, Math.Max(0, intensity));
						output.SetPixel(x, y, Color.FromArgb(intensity, intensity, intensity));
					}
				}
			}
		}

		private static void Test2D()
		{
			Bitmap input = Image.FromFile("lena.png") as Bitmap;

			Bitmap output = input.Clone() as Bitmap;

			// copy
			//Filter2D(input, ref output, (int x, int y, int width, int height) =>
			//{
			//	return 1.0f;
			//});

			// low - pass
			//Filter2D(input, ref output, (int x, int y, int width, int height) =>
			//{
			//	int freq = 16;
			//	if (x > freq && x < width - freq && y > freq && y < height - freq && (x != 0 || y != 0))
			//	{
			//		return 0.0f;
			//	}
			//	return 1.0f;
			//});

			// high-pass
			Filter2D(input, ref output, (int x, int y, int width, int height) =>
			{
				int freq = 10;
				if (!(x > freq && x < width - freq && y > freq && y < height - freq) && (x != 0 || y != 0))
				{
					return 0.0f;
				}
				return 1.0f;
			});

			output.Save("output.png");
		}

		static void Main(string[] args)
		{
			Test1D();
			Test2D();
			System.Console.ReadKey();
		}
	}
}
