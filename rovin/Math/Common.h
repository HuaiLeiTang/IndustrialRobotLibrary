/*!
*	\file	Common.h
*	\date	2016.01.20
*	\author	Keunjun Choi(ckj.monikaru@gmail.com)
*	\brief	General mathematics functions
*/

#pragma once

#include <cmath>
#include <iostream>

#include "Constant.h"

namespace rovin
{
	/// Calculate sine and cosine simultaneously
	static void fsincos(Real theta, ///< Angle in radians 
		Real& sine, ///< Variable for storing a sine value
		Real& cosine ///< Variable for storing a sine value
		)
	{
		theta -= (int)(theta*Inv_PI_DOUBLE)*PI_DOUBLE;
		if (theta < 0) theta += PI_DOUBLE;

		sine = std::sin(theta);
		if (theta < PI_HALF)
		{
			cosine = std::sqrt(1 - sine*sine);
			return;
		}
		else if (theta < PI + PI_HALF)
		{
			cosine = -std::sqrt(1 - sine*sine);
			return;
		}
		cosine = std::sqrt(1 - sine*sine);
	}

	/// Pseudo-inverse of the matrix
	static MatrixX pInv(const MatrixX& matrix ///< Matrix for inversion
		)
	{
		Eigen::JacobiSVD<MatrixX> svd(matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Real tolerance = RealEps;
		VectorX singular_values = svd.singularValues();
		MatrixX S(matrix.rows(), matrix.cols());
		S.setZero();
		for (int i = 0; i < singular_values.size(); i++)
		{
			if (singular_values(i) > tolerance)
			{
				S(i, i) = 1.0 / singular_values(i);
			}
			else
			{
				S(i, i) = 0;
			}
		}
		return svd.matrixV() * S.transpose() * svd.matrixU().transpose();
	}

	/*!
	 *	\brief Absolute value conversion
	 *	\return f(x)=\|x\|
	*/
	static Real Abs(const Real& op)
	{
		if (op < 0.0) return -op;
		return op;
	}

	/// Return min value
	static Real Min(const Real& op1, const Real& op2)
	{
		if (op1 < op2)
		{
			return op1;
		}
		return op2;
	}

	/// Return max value
	static Real Max(const Real& op1, const Real& op2)
	{
		if (op1 > op2)
		{
			return op1;
		}
		return op2;
	}

	static bool RealEqual(const Real& op1, const Real& op2)
	{
		if (std::abs(op1 - op2) < RealEps + RealEps*Abs(op1)) return true;
		return false;
	}

	template< typename T >
	static bool RealEqual(const Eigen::MatrixBase<T>& op1, const Real& op2)
	{
		int rows = op1.rows();
		int cols = op1.cols();

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				if (!RealEqual(op1(i, j), op2))
				{
					return false;
				}
			}
		}
		return true;
	}

	template< typename T1, typename T2 >
	static bool RealEqual(const Eigen::MatrixBase<T1>& op1, const Eigen::MatrixBase<T2>& op2)
	{
		int rows1 = op1.rows();
		int cols1 = op1.cols();

		int rows2 = op2.rows();
		int cols2 = op2.cols();

		if (!(row1 == row2 && cols1 == cols2)) return false;

		for (int i = 0; i < rows1; i++)
		{
			for (int j = 0; j < cols1; j++)
			{
				if (!RealEqual(op1(i, j), op2(i, j)))
				{
					return false;
				}
			}
		}
		return true;
	}

	static bool RealLess(const Real& op1, const Real& op2)
	{
		if (op1 < op2 - RealEps - RealEps*Abs(op1)) return true;
		return false;
	}

	template< typename T >
	static bool RealLess(const Eigen::MatrixBase<T>& op1, const Real& op2)
	{
		int rows = op1.rows();
		int cols = op1.cols();

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				if (!RealLess(op1(i, j), op2))
				{
					return false;
				}
			}
		}
		return true;
	}

	static bool RealLessEqual(const Real& op1, const Real& op2)
	{
		if (op1 < op2 + RealEps + RealEps*Abs(op1)) return true;
		return false;
	}

	template< typename T >
	static bool RealLessEqual(const Eigen::MatrixBase<T>& op1, const Real& op2)
	{
		int rows = op1.rows();
		int cols = op1.cols();

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				if (!RealLessEqual(op1(i, j), op2))
				{
					return false;
				}
			}
		}
		return true;
	}

	static bool RealBigger(const Real& op1, const Real& op2)
	{
		if (op1 > op2 + RealEps + RealEps*Abs(op1)) return true;
		return false;
	}

	template< typename T >
	static bool RealBigger(const Eigen::MatrixBase<T>& op1, const Real& op2)
	{
		int rows = op1.rows();
		int cols = op1.cols();

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				if (!RealBigger(op1(i, j), op2))
				{
					return false;
				}
			}
		}
		return true;
	}

	static bool RealBiggerEqual(const Real& op1, const Real& op2)
	{
		if (op1 > op2 - RealEps - RealEps*Abs(op1)) return true;
		return false;
	}

	template< typename T >
	static bool RealBiggerEqual(const Eigen::MatrixBase<T>& op1, const Real& op2)
	{
		int rows = op1.rows();
		int cols = op1.cols();

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				if (!RealBiggerEqual(op1(i, j), op2))
				{
					return false;
				}
			}
		}
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////// Hong young suk addition /////////////////////////////////////

	static Real VectorInner(const VectorX& vec1, const VectorX& vec2)
	{
		Real value = 0;
		int size = vec1.size();
		for (int i = 0; i < size; i++)
			value += vec1(i) * vec2(i);
		return value;
	}

	static Real VectorInner(const VectorX& vec1, const VectorX& vec2, const int size)
	{
		Real value = 0;
		for (int i = 0; i < size; i++)
			value += vec1(i) * vec2(i);
		return value;
	}

	static void MatrixVectorMul(const MatrixX& A, const VectorX& b, const int row, const int col, VectorX& x)
	{
		x.setZero();
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < col; j++)
			{
				x(i) += A(i, j) * b(j);
			}
		}
	}

	static void VectorSubtraction(const VectorX& a, const VectorX& b, const int s, VectorX& x)
	{
		for (int i = 0; i < s; i++)
			x(i) = a(i) - b(i);
	}
	
	static void VectorSubtraction(VectorX& a, const VectorX& b, const int s, VectorX& x)
	{
		for (int i = 0; i < s; i++)
			x(i) = a(i) - b(i);
	}

	static void VectorAddition(const VectorX& a, const VectorX& b, const int s, VectorX& x)
	{
		for (int i = 0; i < s; i++)
			x(i) = a(i) + b(i);
	}

	static void VectorAddition(VectorX& a, const VectorX& b, const int s, VectorX& x)
	{
		for (int i = 0; i < s; i++)
			x(i) = a(i) + b(i);
	}

	static void CholeskyDecomposition(const MatrixX& inmat, int dim, MatrixX& outmat)
	{
		outmat.setZero();
		int i, j, k;

		for (k = 0; k < dim; k++) // Cholesky decomposition
		{
			outmat(k, k) = inmat(k, k);
			for (j = 0; j < k; j++)
				outmat(k, k) -= outmat(k, j) * outmat(k, j);
			outmat(k, k) = sqrt(outmat(k, k));

			
			for (i = k + 1; i < dim; i++)
			{
				outmat(i, k) = inmat(i, k);
				for (j = 0; j < k; j++)
					outmat(i, k) -= outmat(i, j) * outmat(k, j);
				outmat(i, k) /= outmat(k, k);
			}
		}
	}

	static void IncompleteCholeskyDecomposition(const MatrixX& inmat, int dim, MatrixX& outmat)
	{
		outmat.setZero();
		int i, j, k;

		for (k = 0; k < dim; k++) // Cholesky decomposition
		{
			outmat(k, k) = inmat(k, k);
			for (j = 0; j < k; j++)
				outmat(k, k) -= outmat(k, j) * outmat(k, j);
			outmat(k, k) = sqrt(outmat(k, k));

			for (i = k + 1; i < dim; i++)
			{
				outmat(i, k) = inmat(i, k);
				for (j = 0; j < k; j++)
					outmat(i, k) -= outmat(i, j) * outmat(k, j);
				outmat(i, k) /= outmat(k, k);

				if (inmat(i, k) == 0)
					outmat(i, k) = 0;
			}
		}

		//for (int k = 0; k < num; k++)
		//{
		//	A(k, k) = std::sqrt(A(k, k));
		//	for (int i = (k + 1); i < num; i++)
		//	{
		//		if (A(i, k) != 0)
		//			A(i, k) = A(i, k) / A(k, k);
		//	}

		//	for (int j = (k + 1); j < num; j++)
		//	{
		//		for (int i = j; i < num; i++)
		//		{
		//			if (A(i, j) != 0)
		//				A(i, j) = A(i, j) - A(i, k)*A(j, k);
		//		}
		//	}
		//}
	}

	static void LowerMatrixInverse(const MatrixX& inmat, int dim, MatrixX& outmat)
	{
		outmat.setZero();
		int i, j, k;
		for (i = 0; i < dim; i++) // Inverse of lower matrix
		{
			outmat(i, i) = 1 / inmat(i, i);
			for (j = i + 1; j < dim; j++)
			{
				for (k = i; k < j; k++)
				{
					outmat(j, i) += inmat(j, k) * outmat(k, i);
				}
				outmat(j, i) /= -inmat(j, j);
			}
		}
	}

	static void MultLowerUpperMatrix(const MatrixX& lower, const MatrixX& Upper, int dim, MatrixX& outmat)
	{
		outmat.setZero();
		int i, j, k;
		for (i = 0; i < dim; i++) // multiple lower matrax & upper matrix
		{
			for (j = i; j < dim; j++)
			{
				for (k = j; k < dim; k++)
				{
					outmat(i, j) += Upper(i, k)*lower(k, j);
				}
				outmat(j, i) = outmat(i, j);
			}
		}
	}

	static void PosDefMatrixInverse(const MatrixX& inmat, int dim, MatrixX& outmat)
	{
		// inverse of positive definite matrix('inmat') by using Cholesky decomposition

		MatrixX Lower(dim, dim), iLower(dim, dim), iLowert(dim, dim);
		Lower.setZero(); iLower.setZero(); iLowert.setZero(); outmat.setZero();

		int i, j, k;

		for (k = 0; k < dim; k++) // Cholesky decomposition
		{
			Lower(k, k) = inmat(k, k);
			for (j = 0; j < k; j++)
				Lower(k, k) -= Lower(k, j) * Lower(k, j);
			Lower(k, k) = sqrt(Lower(k, k));

			for (i = k + 1; i < dim; i++)
			{
				Lower(i, k) = inmat(i, k);
				for (j = 0; j < k; j++)
					Lower(i, k) -= Lower(i, j) * Lower(k, j);
				Lower(i, k) /= Lower(k, k);
			}
		}

		for (i = 0; i < dim; i++) // Inverse of lower matrix
		{
			iLower(i, i) = 1 / Lower(i, i);
			for (j = i + 1; j < dim; j++)
			{
				for (k = i; k < j; k++)
				{
					iLower(j, i) += Lower(j, k) * iLower(k, i);
				}
				iLower(j, i) /= -Lower(j, j);
			}
		}

		for (i = 0; i < dim; i++) // transpose
		{
			iLowert(i, i) = iLower(i, i);
			for (j = 0; j < i; j++)
				iLowert(j, i) = iLower(i, j);
		}

		for (i = 0; i < dim; i++) // multiple lower matrax & upper matrix
		{
			for (j = i; j < dim; j++)
			{
				for (k = j; k < dim; k++)
				{
					outmat(i, j) += iLowert(i, k)*iLower(k, j);
				}
				outmat(j, i) = outmat(i, j);
			}
		}
	}



}