#pragma once

#include <rovin\Dynamics\SerialOpenChain.h>
#include <rovin\Math\Interpolation.h>
#include <rovin\Math\GaussianQuadrature.h>
#include <list>
#include <fstream>

namespace rovin {

	class TOPP;
	class SwitchPoint;

	typedef std::shared_ptr<TOPP> TOPPPtr;

	class SwitchPoint
	{
	public:
		enum SPID
		{
			SINGULAR = 0, TANGENT, DISCONTIUOUS
		};
	public:
		SwitchPoint() {}
		SwitchPoint(Real s, Real sdot, SPID id, Real lambda) : _s(s), _sdot(sdot),
			_id(id), _lambda(lambda) {}
	public:
		Real _s;
		Real _sdot;
		SPID _id;
		Real _lambda;
	};

	const enum CONSTRAINT_TYPE
	{
		TORQUE = 1 << 0,
		VELOCITY = 1 << 1,
		ACCELERATION = 1 << 2,

		TORQUE_VEL = TORQUE | VELOCITY,
		TORQUE_ACC = TORQUE | ACCELERATION,
		TORQUE_VEL_ACC = TORQUE | VELOCITY | ACCELERATION,
	};

	class TOPP
	{
	public:
		TOPP(const MatrixX& q_data, const SerialOpenChainPtr& soc, const Real vi, const Real vf,
			const Real ds = 1e-3, const Real si = 0, const Real sf = 1, CONSTRAINT_TYPE constraintType = TORQUE);
		TOPP(const SerialOpenChainPtr& soc, const Real vi = 0, const Real vf = 0,
			const Real ds = 1e-3, const Real si = 0, const Real sf = 1, CONSTRAINT_TYPE constraintType = TORQUE);
		~TOPP() {}

		bool generateTrajectory();
		void initialization();
		void calculateAllMVCPoint();
		void calculateAllSwitchPoint();
		void settingconstraint();
		void makespline();

		void setJointTrajectory(const MatrixX& q_data);
		void setSerialOpenChain(const SerialOpenChainPtr& soc);
		void setInitialParameter(const Real si) { _si = si; }
		void setFinalParameter(const Real sf) { _sf = sf; }
		void setInitialVel(const Real vi) { _vi = vi; }
		void setFinalVel(const Real vf) { _vf = vf; }
		void setStepSize(const Real ds) { _ds = ds; }
		void setConstraintType(CONSTRAINT_TYPE constraintType) { _constraintType = constraintType; }
		
		const std::list<Real>& gets() const { return _s; }
		const std::list<Real>& getsdot() const { return _sdot; }
		const std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& getAllMVCPoint() const { return _allMVCPoints; }
		const std::vector<unsigned int>& getAllMVCPointFlag() const { return _allMVCPointsFlag; }
		const std::vector<SwitchPoint>& getAllSwitchPoint() const { return _switchPoint; }
		const Real getFinalTime() const { return _tf_result; }
		const Real getStepSize() const { return _ds; }
		const Real getInitialParam() const { return _si; }
		const Real getFinalParam() const { return _sf; }
		const MatrixX& getTorqueTrajectory() const { return _torque_result; }
		const unsigned int getdof() const { return _dof; }

		Vector2 determineAlphaBeta(Real s, Real sdot, VectorX& a = VectorX());
		Real calculateMVCPoint(Real s, int& flag);
		Real calculateMVCPointExclude(Real s, int iExclude, int& flag);
		std::vector<VectorX> calculateBandC(Real s);

		unsigned int forward(Real& s_cur, Real& sdot_cur);
		unsigned int forwardVel(Real& s_cur, Real& sdot_cur);
		bool findNearestSwitchPoint(Real s);
		
		void forwardIntegrate(Real& s, Real& sdot, Real sddot);
		void backwardIntegrate(Real& s, Real& sdot, Real sddot);
		void calculateA(Real s, VectorX& a);
		void determineVelminmax(Real s, Vector2& minmax);
		void backward(Real& s_cur, Real& sdot_cur);
		void switchpointfunction(Real& s_cur, Real& sdot_cur, bool& singularPoint_swi);
		
		void calculateFinalTime();
		void calculateTorqueTrajectory();
		
	public:
		CONSTRAINT_TYPE _constraintType;
		SerialOpenChainPtr _soc;
		StatePtr _state;
		unsigned int _dof;

		MatrixX _q_data;
		MatrixX _torque_result;

		BSpline<-1, -1, -1> _q;
		BSpline<-1, -1, -1> _dqds;
		BSpline<-1, -1, -1> _ddqdds;

		std::list<Real> _s;
		std::list<Real> _sdot;
		std::list<Real> _sddot;
		std::vector<SwitchPoint> _switchPoint;
		std::vector<Vector2, Eigen::aligned_allocator<Vector2>> _allMVCPoints;
		std::vector<unsigned int> _allMVCPointsFlag;

		Real _vi, _vf;
		Real _si, _sf;
		Real _ds;
		Real _tf_result;

		VectorX _t;
		VectorX _torqueConstraint;
		VectorX _velConstraint;
		VectorX _accConstraint;
		
		//std::vector<VectorX> _a;
		//std::vector<VectorX> _b;
		//std::vector<VectorX> _c;

		unsigned int _nconstraints;
		unsigned int _nconstraintsWithoutVel;

	public:
		std::string ssibal;

		std::vector<Real> s_FI;
		std::vector<Real> sdot_FI;

		std::vector<Real> s_BI;
		std::vector<Real> sdot_BI;

		void saveRealVector2txt(std::vector<Real> in, std::string filename)
		{
			std::ofstream fout;
			fout.open(filename);

			for (unsigned int i = 0; i < in.size(); i++)
				fout << in[i] << std::endl;

			fout.close();
		}

		void saveVectorX2txt(VectorX in, std::string filename)
		{
			std::ofstream fout;
			fout.open(filename);

			for (unsigned int i = 0; i < in.size(); i++)
				fout << in(i) << std::endl;

			fout.close();
		}

		void saveMatrixX2txt(MatrixX in, std::string filename)
		{
			std::ofstream fout;
			fout.open(filename);

			for (unsigned int i = 0; i < in.cols(); i++)
			{
				for (unsigned int j = 0; j < in.rows(); j++)
					fout << in(j, i) << '\t';
				fout << std::endl;
			}
				
			fout.close();
		}

		void saveIntVector2txt(std::vector<unsigned int> in, std::string filename)
		{
			std::ofstream fout;
			fout.open(filename);

			for (unsigned int i = 0; i < in.size(); i++)
				fout << in[i] << std::endl;

			fout.close();
		}

		void saveMVC(std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& allMVCPoints)
		{
			std::string name = "C:/Users/crazy/Desktop/Time optimization";
			std::vector<Real> s, sdot;
			for (unsigned int i = 0; i < allMVCPoints.size(); i++)
			{
				s.push_back(allMVCPoints[i](0));
				sdot.push_back(allMVCPoints[i](1));
			}
			saveRealVector2txt(s, "C:/Users/crazy/Desktop/Time optimization/avp test/s_MVC.txt");
			saveRealVector2txt(sdot, "C:/Users/crazy/Desktop/Time optimization/avp test/sdot_MVC.txt");
		}

		void savevectorOfVector2(std::vector<Vector2, Eigen::aligned_allocator<Vector2>>& in, std::string filename_first, std::string filename_second)
		{
			std::vector<Real> s, sdot;
			for (unsigned int i = 0; i < in.size(); i++)
			{
				s.push_back(in[i](0));
				sdot.push_back(in[i](1));
			}
			saveRealVector2txt(s, filename_first);
			saveRealVector2txt(sdot, filename_second);
		}

		void saveRealList2txt(std::list<Real> in, std::string filename)
		{
			std::vector<Real> vec;
			for (std::list<Real>::iterator it = in.begin(); it != in.end(); ++it)
			{
				vec.push_back(*(it));
			}
			saveRealVector2txt(vec, filename);
		}

		void saveLC(std::list<Vector2, Eigen::aligned_allocator<Vector2>>& in, std::string filename_s, std::string filename_sdot)
		{
			std::vector<Real> s;
			std::vector<Real> sdot;

			Vector2 tmp;

			for (std::list<Vector2, Eigen::aligned_allocator<Vector2>>::iterator it = in.begin(); it != in.end(); ++it)
			{
				tmp = *(it);
				s.push_back(tmp[0]);
				sdot.push_back(tmp[1]);
			}
			saveRealVector2txt(s, filename_s);
			saveRealVector2txt(sdot, filename_sdot);
		}

		void saveSwitchPoint(std::vector<SwitchPoint>& allSwitchPoint)
		{
			std::vector<Real> s_sw;
			std::vector<Real> sdot_sw;
			for (unsigned int i = 0; i < allSwitchPoint.size(); i++)
			{
				s_sw.push_back(allSwitchPoint[i]._s);
				sdot_sw.push_back(allSwitchPoint[i]._sdot);
			}
			saveRealVector2txt(s_sw, "C:/Users/crazy/Desktop/Time optimization/avp test/s_sw.txt");
			saveRealVector2txt(sdot_sw, "C:/Users/crazy/Desktop/Time optimization/avp test/sdot_sw.txt");
		}

	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

}