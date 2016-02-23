#include <rovin\Optimizer\NonlinearOptimization.h>
#include <rovin\Math\Function.h>

#include <iostream>
#include <conio.h>
#include <list>

using namespace rovin;
using namespace std;

class TestObjFunction : public Function
{
public:
	TestObjFunction() {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(1);
		val(0) = x(0)*x(0) + x(1)*x(1);
		return val;
	}
};

class TestEqFunction : public Function
{
public:
	TestEqFunction() {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(1);
		val(0) = x(0) + x(1) - 5;
		return val;
	}
};

class TestIneqFunction : public Function
{
public:
	TestIneqFunction() {}

	VectorX func(const VectorX& x) const
	{
		VectorX val(1);
		val(0) = -x(1) + 3;
		return val;
	}
};

int main()
{
	//NonlinearOptimization nop;
	//VectorX initX(2);

	//initX.setRandom();

	//nop.setObjectiveFunction(FunctionPtr(new TestObjFunction));
	//nop.setEqualityConstraint(FunctionPtr(new TestEqFunction));
	//nop.setInequalityConstraint(FunctionPtr(new TestIneqFunction));
	//nop.solve(initX);

	//cout << nop.resultX << endl;
	//cout << "f(x) = " << nop.resultFunc << endl;

	list<int> a;
	a.push_back(1);
	a.push_back(2);
	a.push_back(3);

	list<int> b;
	b.push_front(100);
	b.push_front(200);
	b.push_front(300);

	a.merge(b);

	while (!a.empty())
	{
		cout << a.front() << endl;
		a.pop_front();
	}
	cout << "End" << endl;
	cout << "Size : " << a.size() << endl;
	_getch();

	return 0;
}