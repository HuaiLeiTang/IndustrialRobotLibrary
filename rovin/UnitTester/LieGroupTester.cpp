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
	
	cout << State::ALL_INFO << endl;

	_getch();
	return 0;
}