package greenroof;

public class GreenRoofCommon
{

	public static double[] fillArray(int length, double value)
	{
		double[] newArray = new double[length];
		for (int i=0;i<length;i++)
		{
			newArray[i]=value;
		}
		return newArray;
	}
}
