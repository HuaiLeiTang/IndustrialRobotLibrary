/*!
 *	\file	SerialOpenChain.h
 *	\date	2016.01.22
 *	\author	Youngsuk (crazyhys@gmail.com)
 *	\brief	SerialOpenChain class
 *          this class has robot properties including assembly
 *          Robot assembly is implemented in this class
*/

#pragma once

#include <vector>
#include <memory>

#include "Link.h"
#include "MotorJoint.h"
#include "State.h"

namespace rovin
{
	class SerialOpenChain;

	typedef std::shared_ptr<SerialOpenChain> SerialOpenChainPtr;

	class SerialOpenChain
	{

	private:
		/*!
		* \brief SerialOpenChain class member variable
		*/


	public:
		/*!
		* \brief SerialOpenChain class member functions
		*/

		// constructor & destructor
		SerialOpenChain();
		~SerialOpenChain();


		// set-function
		


		// get-function
		


		// deep-copy
		


	};


}