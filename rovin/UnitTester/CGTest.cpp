#include <iostream>
#include <conio.h>
#include <time.h>
#include <fstream>

#include "rovin\Math\Common.h"
#include "rovin\Optimizer\UnconstraintedOptimization.h"

using namespace std;
using namespace rovin;

void loadData(MatrixX& data, std::string filename);

int num = 20;

int main()
{
	MatrixX A(num, num);
	VectorX b(num), x0(num), x(num);

	std::string filename = "C:/Users/crazy/Desktop/A.txt";
	loadData(A, filename);

	b << 9.0327,
		5.4126,
		0.9868,
		8.5708,
		0.4274,
		1.4204,
		2.1750,
		7.7746,
		3.5761,
		0.2532,
		8.3378,
		3.6320,
		7.8856,
		8.7966,
		6.7651,
		5.4524,
		8.9889,
		0.9949,
		0.9094,
		0.3424;

	x0.setOnes();

	clock_t time = clock();
	for (int i = 0; i < 1000; i++)
		x = A.inverse() * b;
	cout << clock() - time << endl;
	cout << x << endl << endl;

	time = clock();
	for (int i = 0; i < 1000; i++)
		irLib::Opt::LinearConjugateGradient(A, b, x0, x);
	cout << clock() - time << endl;
	cout << x << endl;

	cout << "Program Complete" << endl;
	_getch();
	return 0;
}


void loadData(MatrixX& data, std::string filename)
{
	ifstream input(filename);
	if (input.fail()) cout << "파일 열기 실패" << endl;
	else cout << "파일 열기 성공" << endl;

	Real tmp;
	unsigned int cnt = 0;
	while (!input.eof()) {
		input >> tmp;
		cnt++;
	}
	input.close();

	ifstream trajectory(filename);
	if (trajectory.fail()) cout << "파일 열기 실패" << endl;
	else cout << "파일 열기 성공" << endl;

	unsigned int dof = num;
	unsigned int data_num = cnt / dof;
	MatrixX q(dof, data_num);

	for (unsigned int i = 0; i < data_num; i++)
		for (unsigned int j = 0; j < dof; j++)
			trajectory >> q(j, i);

	trajectory.close();

	data = q;
}