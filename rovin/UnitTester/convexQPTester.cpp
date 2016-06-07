#include <iostream>
#include <conio.h>
#include <time.h>
#include <fstream>

#include "rovin\Math\Common.h"
#include "rovin\Optimizer\convexQP.h"

using namespace std;
using namespace rovin;



int main()
{
	MatrixX G(2, 2);
	G.setIdentity();
	VectorX c(2);
	c.setZero();
	MatrixX A(4, 2);
	A << 1, -1,
		0.5, 1,
		0, 1,
		3, -1;
	VectorX b(4);
	b << -1, 2, 2.5, 3;

	irOpt::convexQP cvQP(G, c, A, b);
	cvQP.solveInequalityConstrainedQP();

	//G << 3, 0, 0, 3;
	//c << 0.3, 0;
	//b << -1, 2, 2, 1;
	//cvQP.setProblem(G, c, A, b);
	//cvQP.solveInequalityConstrainedQP();

	cout << cvQP._resultX << endl << endl;
	cout << cvQP._resultmu << endl << endl;

	cout << "Program Complete" << endl;
	_getch();
	return 0;
}
