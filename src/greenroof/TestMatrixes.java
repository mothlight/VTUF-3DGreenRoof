package greenroof;

//import org.apache.commons.math.linear.BlockRealMatrix;
//import org.apache.commons.math3.linear.Array2DRowRealMatrix;
//import org.apache.commons.math3.linear.MatrixUtils;
//import org.apache.commons.math3.linear.RealMatrix;

import jml.matlab.Matlab;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.BlockRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

public class TestMatrixes
{

	public static void main(String[] args)
	{
		TestMatrixes t = new TestMatrixes();
		t.test();

	}
	
	public void test()
	{
//		double[][] double1 = new double[][]{{1,2,3},{4,5,6},{7,8,9}};
//		double[][] double2 = new double[][]{{7,8,9},{4,5,6},{1,2,3}};
//		RealMatrix matrix1 = new Array2DRowRealMatrix(double1);
//		RealMatrix matrix2 = new Array2DRowRealMatrix(double2);
//		RealMatrix result = matrix1.multiply(matrix2);
//		System.out.println(result.toString());
//		RealMatrix result2 = matrix2.multiply(matrix1);
//		System.out.println(result2.toString());
//		
//		RealMatrix result3 = matrix1.add(matrix2);
//		System.out.println(result3.toString());
//		
//		double[][] returnC = GreenRoofCommon.addArrays(double1, double2);//yep, this works.
//		System.out.println(returnC);
//		
//		double[][] eye2 = GreenRoofCommon.eye(2,2);
//		double[][] double4 = new double[][]{{1,2},{3,4}};
//		RealMatrix matrix3 = new Array2DRowRealMatrix(eye2);
//		RealMatrix matrix4 = new Array2DRowRealMatrix(double4);
//		
//		RealMatrix result4 = matrix3.multiply(MatrixUtils.inverse(matrix4));
//		System.out.println(result4);
//		
//		
//		double[][] doubleA = new double[][]{{1,2,3},{4,5,6},{7,8,9}};
////		double[][] doubleD = new double[][]{{3,2,1},{4,5,6},{7,8,9}};
//		RealMatrix matrixA = new Array2DRowRealMatrix(doubleA);
////		RealMatrix matrixD = new Array2DRowRealMatrix(doubleD);
//		RealMatrix matrixAInv = MatrixUtils.inverse(matrixA);
////		RealMatrix matrixDInv = MatrixUtils.inverse(matrixD);
//		RealMatrix resultAD = matrixA.multiply(matrixAInv);
////		RealMatrix resultDA = matrixD.multiply(matrixAInv);
//		System.out.println(resultAD);
////		System.out.println(resultDA);
//		
//		
//		doubleA = new double[][]{{0,2,0,1},{2,2,3,2},{4,-3,0,1.},{6,1,-6,-5}};
//		double[][] doubleB = new double[][]{{0.},{-2.},{-7.},{6.}};
//		matrixA = new Array2DRowRealMatrix(doubleA);
////		RealMatrix matrixB = new Array2DRowRealMatrix(doubleB);
//		matrixAInv = MatrixUtils.inverse(matrixA);
////		matrixAInv = MatrixUtils.inverse(matrixA);
//		RealMatrix resultAB = matrixA.multiply(matrixAInv);
//		System.out.println(resultAB);
		
		
		
		double[][] aData = { {1d,2d}, {3d,4d}};
		double[][] bData = { {3d,6d}, {4d,7d}};
		double[] cData =  {3d,6d};
		RealMatrix A = new BlockRealMatrix(aData);
		RealMatrix B = new BlockRealMatrix(bData);
		RealMatrix C = new Array2DRowRealMatrix(cData);
		
		Matlab.display(A);
//		Matlab.display(B);
		Matlab.display(C);
//		Matlab.display(Matlab.ldivide(A, B));
//		Matlab.display(Matlab.rdivide(A, B));
//		Matlab.display(Matlab.mldivide(A, B));
//		Matlab.display(Matlab.mrdivide(A, B));
		
		Matlab.disp(Matlab.mldivide(A, C));

		
	}

}
