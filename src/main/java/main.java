import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class main {


	public static void main(String[] args) {

		double[][] matrix;
		BlahutArimoto ba;

		int size_of_alphabet = 3;

		matrix = new double[size_of_alphabet][size_of_alphabet];

		matrix[0][0] = 2d/3d;
		matrix[0][1] = 1d/3d;
		matrix[0][2] = 0d/3d;

		matrix[1][0] = 1d/3d;
		matrix[1][1] = 1d/3d;
		matrix[1][2] = 1d/3d;

		matrix[2][0] = 0d/3d;
		matrix[2][1] = 1d/3d;
		matrix[2][2] = 2d/3d;

		ba = new BlahutArimoto(matrix);

		double capacity = ba.computeCapacity();

		System.out.println("Computed capacity is " + capacity);

		return;


	}
}
