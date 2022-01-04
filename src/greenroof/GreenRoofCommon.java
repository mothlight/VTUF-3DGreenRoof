package greenroof;

public class GreenRoofCommon
{
	
	public static final double threeDigitDelta = 1e3;
	public static final double twoDigitDelta = 1e2;
	public static final double oneDigitDelta = 1e1;
	
	public static final int ONE = 0;
	public static final int TWO = 1;
	public static final int THREE = 2;
	public static final int FOUR = 3;
	public static final int FIVE = 4;
	public static final int SIX = 5;

	public static double[] fillArray(int length, double value)
	{
		double[] newArray = new double[length];
		for (int i=0;i<length;i++)
		{
			newArray[i]=value;
		}
		return newArray;
	}
	public static double[] multiplyArrays(double[] array, int[] array2)
	{
		if (array.length == 1)
		{
			double[] returnArray = new double[array2.length];
			
			for (int i=0;i<array.length;i++)
			{
				returnArray[i] = array[0] * array2[i];
			}
			
			return returnArray;
		}
		if (array2.length == 1)
		{
			double[] returnArray = new double[array.length];
			
			for (int i=0;i<array.length;i++)
			{
				returnArray[i] = array[i] * array2[0];
			}
			
			return returnArray;
		}
		
		double[] returnArray = new double[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = array[i] * array2[i];
		}
		
		return returnArray;
	}
	
	public static double[] multiplyArrays(int[] array, double[] array2)
	{
		double[] returnArray = new double[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = array[i] * array2[i];
		}
		
		return returnArray;
	}
	
	public static double[] multiplyArrays(double[] array, double[] array2)
	{
		if (array.length == 1 && array2.length != 1)
		{
			double[] returnArray = new double[array2.length];
			
			for (int i=0;i<array2.length;i++)
			{
				returnArray[i] = array[0] * array2[i];
			}
			
			return returnArray;
		}
		
		double[] returnArray = new double[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = array[i] * array2[i];
		}
		
		return returnArray;
	}
	
	public static double[] addArrays(double[] array, double[] array2)
	{
		double[] returnArray = new double[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = array[i] + array2[i];
		}
		
		return returnArray;
	}
	
	public static double[] cosArray(double[] array)
	{
		double[] returnArray = new double[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = Math.cos(array[i]) ;
		}
		
		return returnArray;
	}
	public static double[] sinArray(double[] array)
	{
		double[] returnArray = new double[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = Math.sin(array[i]) ;
		}
		
		return returnArray;
	}
	
	public static double[] minAndPosition(double[] array)
	{
		
		double min = Double.MAX_VALUE;
		int minPosition = -9999;
		for (int i=0;i<array.length;i++)
		{
			if (array[i] == Double.NaN)
			{
				continue;
			}
			if (min > array[i])
			{
				min = array[i];
				minPosition = i;
			}
		}
		double[] minValue = new double[]{min,minPosition};
		return minValue;
	}
	
	public static double[] max(double[] array)
	{
		double[] maxReturn = new double[1];
		double max = Double.MIN_VALUE;
		for (int i=0;i<array.length;i++)
		{
			if (max < array[i])
			{
				max = array[i];
			}
		}
		maxReturn[0]=max;
		return maxReturn;
	}
	
	public static int max(int[] array)
	{
		int maxValue = Integer.MIN_VALUE;
		
		for (int item : array)
		{
			if (item > maxValue)
			{
				maxValue = item;
			}
		}
		
		return maxValue;
	}
	
	public static int min(int[] array)
	{
		int minValue = Integer.MAX_VALUE;
		
		for (int item : array)
		{
			if (item < minValue)
			{
				minValue = item;
			}
		}
		
		return minValue;
	}
	
	public static double[] min(double[] array)
	{
		double[] minReturn = new double[1];
		double min = Double.MIN_VALUE;
		for (int i=0;i<array.length;i++)
		{
			if (min > array[i])
			{
				min = array[i];
			}
		}
		minReturn[0]=min;
		return minReturn;
	}
	
	public static double[] min(double[] array, double[] array2)
	{
		double[] minReturn = new double[array.length];

		for (int i=0;i<array.length;i++)
		{
			minReturn[i] = min(array[i],array2[i]);
		}
		return minReturn;
	}
	
	public static double max(double array)
	{
		return array;
	}
	public static double max(double value, double[] array)
	{
		double[] maxReturn = new double[1];
		double max = Double.MIN_VALUE;
		for (int i=0;i<array.length;i++)
		{
			if (max < array[i])
			{
				max = array[i];
			}
		}
		if (max < value)
		{
			max = value;
		}
		maxReturn[0]=max;
//		return maxReturn;
		return max;
	}
	
	public static double[] get1DSliceOf2D(double[][] array, int dimensionToSlice, int sliceItem)
	{
		int length;
		if (dimensionToSlice == 1)
		{
			length = array.length;
		}
		else
		{
			length = array[0].length;
		}
		double[] returnArray = new double[length];
		

		for (int i=0;i<length;i++)
		{
			if (dimensionToSlice == 1)
			{
				returnArray[i] = array[i][sliceItem];
				
			}
			else
			{
				returnArray[i] = array[sliceItem][i];
			}
			
		}
		return returnArray;
	}
	
	public static double[] nans(int length)
	{
		double[] returnArray = new double[length];
		for (int i=0;i<length;i++)
		{
			returnArray[i] = Double.NaN;
		}
		return returnArray;
	}
	
	public static double[] ones(int length)
	{
		double[] returnArray = new double[length];
		for (int i=0;i<length;i++)
		{
			returnArray[i] = 1.0;
		}
		return returnArray;
	}
	
	public static double[] multiply(double[] array, double otherTerm)
	{
		double[] returnArray = new double[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = array[i] * otherTerm;
		}
		
		return returnArray;
	}
	public static double[] multiply(double otherTerm, double[] array)
	{
		double[] returnArray = new double[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] =otherTerm * array[i] ;
		}
		
		return returnArray;
	}
	
	public static double min(double item1, double item2)
	{
		if (Double.isNaN(item1))
		{
			return item2;
		}
		else if (Double.isNaN(item2))
		{
			return item1;
		}
		else
		{
			return Math.min(item1, item2);
		}
	}
	
	
	public static double[] subtract(double[] array, double otherTerm)
	{
		double[] returnArray = new double[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = array[i] - otherTerm;
		}
		
		return returnArray;
	}
	
	public static int[] subtract(int[] array, int otherTerm)
	{
		int[] returnArray = new int[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = array[i] - otherTerm;
		}
		
		return returnArray;
	}
	
	public static double[] subtract(double[] array, double[] array2)
	{
		double[] returnArray = new double[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = array[i] - array2[i];
		}
		
		return returnArray;
	}
	
	public static double[] subtract(double item, double[] array)
	{
		double[] returnArray = new double[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = item - array[i] ;
		}
		
		return returnArray;
	}
	
	public static int[] subtract(int item, int[] array)
	{
		int[] returnArray = new int[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = item - array[i] ;
		}
		
		return returnArray;
	}
	
	public static double[] getSumOfLayersOf2D(double[][] array, int dimensionToSlice)
	{
		double[] returnArray = new double[0];
		if (dimensionToSlice == 1)
		{
			int length,width;
			length = array.length;
			width = array[0].length;
			returnArray = new double[length];
			for (int i=0;i<length;i++)
			{		
				for (int sliceItem=0;sliceItem<width;sliceItem++)
				{
					returnArray[i] = array[i][sliceItem] + returnArray[i];
				}				
			}
		}
		
		if (dimensionToSlice == 2)
		{
			int length,width;
			length = array.length;
			width = array[0].length;
			returnArray = new double[width];
			for (int i=0;i<width;i++)
			{		
				for (int sliceItem=0;sliceItem<length;sliceItem++)
				{
					returnArray[i] = array[sliceItem][i] + returnArray[i];
				}				
			}
		}
		

		return returnArray;
	}
	
	


	public static double[][] update2DWith1DSlice(double[][] array, double[] arraySlice, int dimensionToSlice, int sliceItem)
	{
		int length;
		if (dimensionToSlice == 1)
		{
			length = array.length;
		}
		else
		{
			length = array[0].length;
		}
//		double[] returnArray = new double[length];
		

		for (int i=0;i<length;i++)
		{
			if (dimensionToSlice == 1)
			{
				array[i][sliceItem] = arraySlice[i]  ;
				
			}
			else
			{
				array[sliceItem][i] = arraySlice[i];
			}
			
		}
		return array;
	}
	
	
	public static double[] add(double[] array, double otherTerm)
	{
		double[] returnArray = new double[array.length];
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = array[i] + otherTerm;
		}
		
		return returnArray;
	}
	
	public static double[] add(double[] array, double[] array2)
	{
		double[] returnArray = new double[array.length];
		if (array2.length == 1)
		{
			return add(array,array2[0]);
		}
		
		for (int i=0;i<array.length;i++)
		{
			returnArray[i] = array[i] + array2[i];
			if (returnArray[i] < 0.00000001 && returnArray[i]>= 0)
			{
				returnArray[i]  = 0.0;
			}
		}
		
		return returnArray;
	}
	public static double[] subsetArray(double[] array, int start, int end)
	{
		int length = end - start+1;
		double[] returnArray = new double[length];
		int count = 0;
		for (int i=start;i<=end;i++)
		{
			returnArray[count] = array[i];
			count++;
		}
		return returnArray;
	}
	
	public static int[] subsetArray(int[] array, int start, int end)
	{
		int length = end - start+1;
		int[] returnArray = new int[length];
		int count = 0;
		for (int i=start;i<=end;i++)
		{
			returnArray[count] = array[i];
			count++;
		}
		return returnArray;
	}
}
