using System;
using System.Collections.Generic;
using System.Text;

namespace KissFFT
{
	public class KissFFTnd<kiss_fft_scalar>
	{
		private int Dimprod; /* dimsum would be mighty tasty right now */
		private int Ndims;
		private int[] Dims;
		private KissFFT<kiss_fft_scalar>[] States; /* cfg states for each dimension */
		private kiss_fft_cpx<kiss_fft_scalar>[] Tmpbuf; /*buffer capable of hold the entire input */

		private IArithmetic<kiss_fft_scalar> A;

		public KissFFTnd(int[] dims, int ndims, bool inverse, IArithmetic<kiss_fft_scalar> arithmetic)
		{
			this.A = arithmetic;
			kiss_fftnd_alloc(dims, ndims, inverse);
		}

		private void kiss_fftnd_alloc(int[] dims, int ndims, bool inverse_fft)
		{
			Dimprod = 1;
			for (int i = 0; i < ndims; ++i)
			{
				Dimprod *= dims[i];
			}
			Ndims = ndims;
			Dims = new int[dims.Length];
			for (int i = 0; i < dims.Length; ++i)
			{
				Dims[i] = dims[i];
			}

			States = new KissFFT<kiss_fft_scalar>[ndims];
			for (int i = 0; i < States.Length; ++i)
			{
				States[i] = new KissFFT<kiss_fft_scalar>(dims[i], inverse_fft, A);
			}

			Tmpbuf = new kiss_fft_cpx<kiss_fft_scalar>[Dimprod];
		}

		/*
		 This works by tackling one dimension at a time.

		 In effect,
		 Each stage starts out by reshaping the matrix into a DixSi 2d matrix.
		 A Di-sized fft is taken of each column, transposing the matrix as it goes.

		Here's a 3-d example:
		Take a 2x3x4 matrix, laid out in memory as a contiguous buffer
		 [ [ [ a b c d ] [ e f g h ] [ i j k l ] ]
		   [ [ m n o p ] [ q r s t ] [ u v w x ] ] ]

		Stage 0 ( D=2): treat the buffer as a 2x12 matrix
		   [ [a b ... k l]
			 [m n ... w x] ]

		   FFT each column with size 2.
		   Transpose the matrix at the same time using kiss_fft_stride.

		   [ [ a+m a-m ]
			 [ b+n b-n]
			 ...
			 [ k+w k-w ]
			 [ l+x l-x ] ]

		   Note fft([x y]) == [x+y x-y]

		Stage 1 ( D=3) treats the buffer (the output of stage D=2) as an 3x8 matrix,
		   [ [ a+m a-m b+n b-n c+o c-o d+p d-p ] 
			 [ e+q e-q f+r f-r g+s g-s h+t h-t ]
			 [ i+u i-u j+v j-v k+w k-w l+x l-x ] ]

		   And perform FFTs (size=3) on each of the columns as above, transposing 
		   the matrix as it goes.  The output of stage 1 is 
			   (Legend: ap = [ a+m e+q i+u ]
						am = [ a-m e-q i-u ] )
   
		   [ [ sum(ap) fft(ap)[0] fft(ap)[1] ]
			 [ sum(am) fft(am)[0] fft(am)[1] ]
			 [ sum(bp) fft(bp)[0] fft(bp)[1] ]
			 [ sum(bm) fft(bm)[0] fft(bm)[1] ]
			 [ sum(cp) fft(cp)[0] fft(cp)[1] ]
			 [ sum(cm) fft(cm)[0] fft(cm)[1] ]
			 [ sum(dp) fft(dp)[0] fft(dp)[1] ]
			 [ sum(dm) fft(dm)[0] fft(dm)[1] ]  ]

		Stage 2 ( D=4) treats this buffer as a 4*6 matrix,
		   [ [ sum(ap) fft(ap)[0] fft(ap)[1] sum(am) fft(am)[0] fft(am)[1] ]
			 [ sum(bp) fft(bp)[0] fft(bp)[1] sum(bm) fft(bm)[0] fft(bm)[1] ]
			 [ sum(cp) fft(cp)[0] fft(cp)[1] sum(cm) fft(cm)[0] fft(cm)[1] ]
			 [ sum(dp) fft(dp)[0] fft(dp)[1] sum(dm) fft(dm)[0] fft(dm)[1] ]  ]

		   Then FFTs each column, transposing as it goes.

		   The resulting matrix is the 3d FFT of the 2x3x4 input matrix.

		   Note as a sanity check that the first element of the final 
		   stage's output (DC term) is 
		   sum( [ sum(ap) sum(bp) sum(cp) sum(dp) ] )
		   , i.e. the summation of all 24 input elements. 

		*/
		public void kiss_fftnd(Array<kiss_fft_cpx<kiss_fft_scalar>> fin, Array<kiss_fft_cpx<kiss_fft_scalar>> fout)
		{
			Array<kiss_fft_cpx<kiss_fft_scalar>> bufin = new Array<kiss_fft_cpx<kiss_fft_scalar>>(fin);
			Array<kiss_fft_cpx<kiss_fft_scalar>> bufout;

			/*arrange it so the last bufout == fout*/
			if ((Ndims & 1) != 0)
			{
				bufout = new Array<kiss_fft_cpx<kiss_fft_scalar>>(fout);
				if (fin == fout)
				{
					for (int i = 0; i < Tmpbuf.Length; ++i)
					{
						Tmpbuf[i] = new kiss_fft_cpx<kiss_fft_scalar>(fin[i]);
					}
					bufin = new Array<kiss_fft_cpx<kiss_fft_scalar>>(Tmpbuf);
				}
			}
			else
				bufout = new Array<kiss_fft_cpx<kiss_fft_scalar>>(Tmpbuf);

			for (int k = 0; k < Ndims; ++k)
			{
				int curdim = Dims[k];
				int stride = Dimprod / curdim;

				for (int i = 0; i < stride; ++i)
					States[k].kiss_fft_stride(bufin + i, bufout + i * curdim, stride);

				/*toggle back and forth between the two buffers*/
				if (bufout == Tmpbuf)
				{
					bufout = new Array<kiss_fft_cpx<kiss_fft_scalar>>(fout);
					bufin = new Array<kiss_fft_cpx<kiss_fft_scalar>>(Tmpbuf);
				}
				else
				{
					bufout = new Array<kiss_fft_cpx<kiss_fft_scalar>>(Tmpbuf);
					bufin = new Array<kiss_fft_cpx<kiss_fft_scalar>>(fout);
				}
			}
		}
	}
}
