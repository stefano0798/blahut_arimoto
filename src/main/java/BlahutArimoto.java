import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class BlahutArimoto {

	private final int BASE_LOG = 2;
	private final int MAX_ITER = 10;
	private final double THRESHOLD = 1.0E-12;


	private RealMatrix p_y_x;

	private double[] marginal_x, marginal_y;

	public BlahutArimoto(double[][] matrix) {
		if (matrix.length == 0 || matrix[0].length == 0) {
			throw new IllegalArgumentException("Matrix size cannot be 0");
		}

		p_y_x = MatrixUtils.createRealMatrix(matrix);

		// start: create and initialize to 0 the marginal probabilities

		marginal_x = new double[matrix.length];
		for (int i=0; i<matrix.length; i++) {
			marginal_x[i] = 0;
		}
		marginal_y = new double[matrix[0].length];
		for (int i=0; i<matrix[0].length; i++) {
			marginal_y[i] = 0;
		}

		// end: create and initialize to 0 the marginal probabilities
		// start: compute marginal distributions

		for (int i = 0; i< p_y_x.getRowDimension(); i++) {
			for (int j = 0; j< p_y_x.getColumnDimension(); j++) {
				marginal_x[j] += p_y_x.getEntry(i, j);
				marginal_y[i] += p_y_x.getEntry(i, j);
			}
		}

		// end: compute marginal distributions

	}

	public double computeCapacity() {
		boolean running = true;
		double[][] r = new double[1][p_y_x.getColumnDimension()];

		int m = p_y_x.getRowDimension();

		for (int i=0; i<r[0].length; i++) {
			r[0][i] = 1d / r[0].length;
		}

		int iteration = 0;
		RealMatrix q = null;
		while (iteration < MAX_ITER && running == true) {

			RealMatrix rT = MatrixUtils.createRealMatrix(r).transpose();
			q = multiplyByColumn(rT, p_y_x);

			q = normalizeByColumn(q);

			RealMatrix powerMatrix = powerMatrix(q, p_y_x);

			RealMatrix r1 =  MatrixUtils.createRealMatrix(powerMatrix.getRowDimension(), 1);

			for (int i=0; i<r1.getRowDimension(); i++) {
				r1.setEntry(i, 0, 1); // initialize with neutral element of the multiplication
			}

			for (int i=0; i<powerMatrix.getRowDimension(); i++) {
				for (int j=0; j<powerMatrix.getColumnDimension(); j++) {
					r1.setEntry(i, 0, r1.getEntry(i, 0) * powerMatrix.getEntry(i, j));
				}
			}

			double sumOfR1 = 0d;
			for(int i=0; i<r1.getRowDimension(); i++) {
				sumOfR1 += r1.getEntry(i, 0);
			}
			for(int i=0; i<r1.getRowDimension(); i++) {
				r1.setEntry(i, 0, r1.getEntry(i, 0)/sumOfR1);
			}

			RealMatrix diff = MatrixUtils.createRealMatrix(r1.getRowDimension(), 1);

			for (int i=0; i<diff.getRowDimension(); i++) {
				diff.setEntry(i, 0, r1.getEntry(i, 0) - r[0][i]);
			}

			double tolerance = diff.getFrobeniusNorm();

			for (int i=0; i<r[0].length; i++) {
				System.out.println("Updating r[0][" + i + "] from " + r[0][i] + " to " + r1.getEntry(i, 0));
				r[0][i] = r1.getEntry(i, 0);
			}

			if (tolerance < THRESHOLD)
				running = false;

			iteration = iteration +1;
		}

		double c = 0;

		for(int i=0; i<m; i++) {
			if (r[0][i] > 0) {
				double r_i = r[0][i];
				double[] rowP = p_y_x.getRow(i);
				double[] rowQ = q.getRow(i);

				for (int j=0; j<rowQ.length; j++) {
					rowQ[j] = rowQ[j] / r_i + 1e-16;
					rowQ[j] = Math.log(rowQ[j]);
				}
				for (int j=0; j<rowQ.length; j++) {
					double tmp = r_i * rowP[j] * rowQ[j];
					c += tmp;
				}
			}
		}
		return c / Math.log(BASE_LOG);
	}

	private RealMatrix multiplyByColumn(RealMatrix m1, RealMatrix m2) {

		RealMatrix result = MatrixUtils.createRealMatrix(m1.getRowDimension(), m2.getColumnDimension());

		for (int i=0; i<m1.getRowDimension(); i++) {
			for (int j=0; j<m2.getRowDimension(); j++) {
				result.setEntry(i, j, m1.getEntry(i, 0) * m2.getEntry(i, j));
			}
		}

		return result;

	}

	private RealMatrix normalizeByColumn(RealMatrix m) {
		RealMatrix sums = MatrixUtils.createRealMatrix(1, m.getColumnDimension());
		for (int i=0; i<m.getRowDimension(); i++) {
			for (int j=0; j<m.getColumnDimension(); j++) {
				sums.setEntry(0, j, sums.getEntry(0, j) + m.getEntry(i, j));
			}
		}

		for (int i=0; i<m.getRowDimension(); i++) {
			for (int j=0; j<m.getColumnDimension(); j++) {
				m.setEntry(i, j, m.getEntry(i, j)/sums.getEntry(0, j));
			}
		}
		return m;
	}

	private RealMatrix powerMatrix(RealMatrix base, RealMatrix esp) {
		RealMatrix powerMatrix = MatrixUtils.createRealMatrix(esp.getRowDimension(), esp.getColumnDimension());

		for (int i=0; i<powerMatrix.getRowDimension(); i++) {
			for (int j=0; j<powerMatrix.getColumnDimension(); j++) {
				powerMatrix.setEntry(i, j, Math.pow(base.getEntry(i, j), esp.getEntry(i, j)));
			}
		}
		return powerMatrix;
	}

	private void printMatrix(RealMatrix m) {

		double[][] mt = m.getData();

		for (int i=0; i<m.getRowDimension(); i++) {
			System.out.println("Column " + i);
			for (int j=0; j<m.getColumnDimension(); j++) {
				System.out.println("Row " + j + ": " + m.getEntry(i, j));
			}
			System.out.print("\n");
		}

	}

}
