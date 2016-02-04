#include <iostream>
#include <conio.h>
#include <rovin/Math/LieGroup.h>
#include <rovin\Utils\Diagnostic.h>

#include <Eigen\Dense>

#include <vector>
#include <rovin\Dynamics\State.h>


using namespace std;
using namespace rovin;

class A
{
public:
	int a;
	int b;
	int add() { return a + b; }
	

	A() { 
		//a = 0;
		//b = 0;
	}

	A(int x, int y = int()) : a(x), b(y) {}
	void print() 
	{
		cout << "a : " << this->a << endl;
		cout << "b : " << this->b << endl;
	}
	int geta() { return a; }
	int getb() { return b; }
};

int main()
{
	//SE3 T(SO3::EulerZYX(PI, 0.0, 0.0), Vector3(10, 20, 30));

	//cout << T << endl;

	//int a = 1;
	//LOGIF(a==0, "Hi");
	
	//cout << State::ALL_INFO << endl;

	Matrix6X a(1, 3);
	Matrix6X b(2, 3);
	Matrix6X c(3, 3);
	Matrix6X d(4, 3);
	Matrix6X e(5, 3);
	Matrix6X f(6, 3);
	Matrix6X g(7, 3);
	MatrixX h(3, 3);

	cout << a << endl;
	cout << endl;
	cout << b << endl;
	cout << endl;
	cout << c << endl;
	cout << endl;
	cout << d << endl;
	cout << endl;
	cout << e << endl;
	cout << endl;
	cout << f << endl;
	cout << endl;
	cout << g << endl;
	cout << endl;
	cout << h << endl;

	//VectorX a(3);
	//cout << a.size() << endl;

	_getch();
	return 0;
}