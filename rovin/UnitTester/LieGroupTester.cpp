#include <iostream>
#include <conio.h>
#include <rovin/Math/LieGroup.h>

using namespace std;
using namespace rovin;

int main()
{
	SE3 T(SO3::EulerZYX(PI, 0.0, 0.0), Vector3(10, 20, 30));

	cout << T << endl;

	_getch();
	return 0;
}