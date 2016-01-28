using System;
using System.Collections.Generic;
using System.Text;

namespace KissFFT
{
	public class Array<T>
	{
		private T[] FullArray;
		private int MaxSize = 0;

		private int Offset = 0;

		public Array(T[] fullArray)
		{
			FullArray = fullArray;
			MaxSize = fullArray.Length;
		}

		public Array(Array<T> array)
		{
			FullArray = array.FullArray;
			MaxSize = array.MaxSize;

			Offset = array.Offset;
		}

		public Array(Array<T> array, int offset)
		{
			FullArray = array.FullArray;
			MaxSize = array.MaxSize;

			Offset = array.Offset + offset;
		}

		public T this[int index]
		{
			get
			{
				return FullArray[Offset + index];
			}

			set
			{
				FullArray[Offset + index] = value;
			}
		}

		public override bool Equals(Object o)
		{
			Array<T> a = o as Array<T>;
			if (!ReferenceEquals(a, null))
			{
				return this == a;
			}
			return (Object) this == o;
		}

		public override int GetHashCode()
		{
			return base.GetHashCode();
		}

		public static bool operator ==(Array<T> a, Array<T> b)
		{
			if (!ReferenceEquals(a, null) && !ReferenceEquals(b, null))
			{
				return (a.FullArray == b.FullArray) && (a.Offset == b.Offset) && (a.MaxSize == b.MaxSize);
			}
			return false;
		}

		public static bool operator !=(Array<T> a, Array<T> b)
		{
			return !(a == b);
		}

		public static bool operator ==(Array<T> a, T[] b)
		{
			if (!ReferenceEquals(a, null) && !ReferenceEquals(b, null))
			{
				return a.FullArray == b;
			}
			return ReferenceEquals(a, b);
		}

		public static bool operator !=(Array<T> a, T[] b)
		{
			return !(a == b);
		}

		public static Array<T> operator ++(Array<T> a)
		{
			return new Array<T>(a, 1);
		}

		public static Array<T> operator +(Array<T> a, int offset)
		{
			return new Array<T>(a, offset);
		}
	}
}
