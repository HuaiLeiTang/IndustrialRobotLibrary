#include <rovin\Optimizer\NonlinearOptimization.h>
#include <rovin\EnergyOptimization\PTPOptimization.h>
#include <rovin\Math\Function.h>

#include <rovin\Renderer\OSG_simpleRender.h>

#include "efortRobot.h"

#include <iostream>
#include <conio.h>
#include <Windows.h>

using namespace rovin;
using namespace std;

unsigned long __stdcall NET_RvThr(void * pParam);
unsigned long __stdcall apply(void * pParam);

HANDLE hPipe;
BOOL Finished;

SerialOpenChainPtr robot;
PTPOptimization* PTPManager;

int main()
{
	
	LPTSTR lpszPipename = TEXT("\\\\.\\pipe\\NamePipeGoogolTech");

	//Thread Init Data
	HANDLE hThread = NULL;
	BOOL Write_St = TRUE;
	Finished = FALSE;

	hPipe = CreateFile(lpszPipename, GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_FLAG_OVERLAPPED, NULL);

	if (hPipe == NULL || hPipe == INVALID_HANDLE_VALUE)
	{
		printf("Please execute the control panel  - (error %d)\n", GetLastError());

	}
	else
	{
		robot = SerialOpenChainPtr(new efortRobot());
		unsigned int dof = robot->getNumOfJoint();

		StatePtr initState, finalState;
		initState = robot->makeState();
		finalState = robot->makeState();

		VectorX init_q(dof), init_qdot(dof), init_qddot(dof);
		VectorX final_q(dof), final_qdot(dof), final_qddot(dof);

		init_q << 1.0854, 1.02654, 0.798359, 2.97849, 1.50724, 1.45496;
		init_qdot.setZero();
		init_qddot.setZero();

		final_q << -1.31465, 0.128787, -0.546992, 2.83671, -1.8583, 2.95029;
		final_qdot.setZero();
		final_qddot.setZero();

		initState->setJointStatePos(init_q);
		initState->setJointStateVel(init_qdot);
		initState->setJointStateAcc(init_qddot);

		finalState->setJointStatePos(init_q);
		finalState->setJointStateVel(init_qdot);
		finalState->setJointStateAcc(init_qddot);

		vector<bool> optJoint(robot->getNumOfJoint());
		optJoint[0] = optJoint[1] = optJoint[2] = true;

		PTPManager = new PTPOptimization(robot, optJoint, 4, 4, 20, 2.0, initState, finalState);
		

		hThread = CreateThread(NULL, 0, &NET_RvThr, NULL, 0, NULL);

		while (1)
		{
			
		}

		CloseHandle(hPipe);
		Finished = TRUE;
	}

	_getch();

	return 0;
}

unsigned long __stdcall NET_RvThr(void * pParam)
{
	BOOL fSuccess;
	char chBuf[500];
	DWORD dwBytesToWrite = (DWORD)strlen(chBuf);
	DWORD cbRead;

	while (1)
	{
		for (int i = 0; i < 500; i++) chBuf[i] = 0;
		fSuccess = ReadFile(hPipe, chBuf, 500, &cbRead, NULL);
		if (fSuccess)
		{
			int ix;
			string command(chBuf);
			while ((ix = command.find('=')) != std::string::npos)
			{
				string scommand = command.substr(0, ix);
				command.erase(0, ix + 1);

				cout << scommand << endl;
				vector<string> subcommand;
				int count = 0;
				int idx = 0;
				while ((idx = scommand.find(' ')) != std::string::npos)
				{
					subcommand.push_back(scommand.substr(0, idx));
					scommand.erase(0, idx + 1);
				}
				subcommand.push_back(scommand);

				if (subcommand[0].compare("sx") == 0)
				{
					
				}
				else if (subcommand[0].compare("apply") == 0)
				{
					// thread 하나 더 만들어줌..
					HANDLE hThread = NULL;
					hThread = CreateThread(NULL, 0, &apply, NULL, 0, NULL);
				}
				
				else if (subcommand[0].compare("save") == 0)
				{
					HANDLE hThread = NULL;
					//hThread = CreateThread(NULL, 0, &savetraj, NULL, 0, NULL);
				}
			}
		}
		if (!fSuccess && GetLastError() != ERROR_MORE_DATA)
		{
			if (Finished)
				break;
		}
	}

	return 0;

}

unsigned long __stdcall apply(void * pParam)
{
	cout << "APPLY BUTTON CLICK" << endl;

	PTPManager->generateTrajectory();

	return 0;
}