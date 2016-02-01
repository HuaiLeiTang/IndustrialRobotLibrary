#include <iostream>
#include <conio.h>
#include <rovin/Math/LieGroup.h>
#include <rovin\Utils\Diagnostic.h>

#include <Eigen\Dense>

#include <vector>


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

	std::vector<int> a;
	a.push_back(0);
	a.push_back(1);
	a.push_back(2);
	a.push_back(3);
	a.push_back(4);

	a.insert(a.begin() + 5, 10);
	cout << a.size() << endl;
	cout << endl;
	for (int i = 0; i < a.size(); i++)
		cout << a[i] << endl;

	_getch();
	return 0;
}