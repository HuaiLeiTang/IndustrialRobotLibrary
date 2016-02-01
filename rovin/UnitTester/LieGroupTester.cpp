#include <iostream>
#include <conio.h>
#include <rovin/Math/LieGroup.h>

#include <Eigen\Dense>


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
	SE3 T(SO3::EulerZYX(PI, 0.0, 0.0), Vector3(10, 20, 30));

	cout << T << endl;

	Real a = std::numeric_limits<Real>::min();
	Real b = std::numeric_limits<Real>::max();

	cout << a << endl;
	cout << b << endl;

	//A aa(1);
	//aa.print();
	//int result = aa.add();
	//cout << result << endl;
	//aa.print();

	//int geta = aa.geta();
	//int getb = aa.getb();
	//cout << geta << endl;
	//cout << getb << endl;

	//Vector6 a;
	//a << 0, 1, 2, 3, 4, 5;
	//cout << a << endl;
	//cout << endl;
	//cout << a.head(3) << endl;
	//Vector3 b = a.tail(3);
	//cout << endl;
	//cout << b << endl;


	_getch();
	return 0;
}