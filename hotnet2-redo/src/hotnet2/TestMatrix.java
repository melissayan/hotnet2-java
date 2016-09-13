package hotnet2;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;

import java.io.IOException;

import org.junit.Test;
//import org.netlib.blas.*;
import com.github.fommil.netlib.*;




public class TestMatrix{

	@Test
	public void testCreateMatrix() throws IOException{
	//	-Dcom.github.fommil.netlib.BLAS=com.github.fommil.netlib.NativeRefBLAS

		System.out.println(BLAS.getInstance().getClass().getName());
		int n = 5;
		double[][] matValues = new double[n][n];
		DenseMatrix matA = new DenseMatrix(matValues);
		DenseMatrix matB = new DenseMatrix(matValues);
		DenseMatrix result = new DenseMatrix(matA.numRows(),matB.numColumns());
		matA.mult(matB,result);
	}
	
}