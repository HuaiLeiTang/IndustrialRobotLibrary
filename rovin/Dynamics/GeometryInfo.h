/*!
 *	\file	GeometryInfo.h
 *	\date	2015.11.04
 *	\author	Keunjun (ckj@robotics.snu.ac.kr)
 *	\brief	GeometryInfo class
*/

#pragma once

#include <string>
#include <list>
#include <memory>

#include <rovin/Math/Constant.h>
#include <rovin/Math/LieGroup.h>

namespace rovin
{
	class GeometryInfo;

	typedef std::shared_ptr< GeometryInfo > GeometryInfoPtr;

	class GeometryInfo
	{
	public:
		/// Geometry type�� ����
		enum GEOMETRY_TYPE
		{
			_BOX,		///< Box ����
			_SPHERE,	///< Sphere ��
			_CAPSULE,	///< Capsule ĸ��
			_CYLINDER,	///< Cylinder ����
			_MESH,		///< Mesh mesh�� ������ ����
			_USERMODEL	///< UserModel ����ڰ� ���� ����� ��
		};

		/// ������
		GeometryInfo(const GEOMETRY_TYPE& Type, ///< Type
			const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : _Type(Type), _T(T), _Color(Color) {}

		/// �Ҹ���, ���� Ŭ������ ����� �ݴϴ�.
		virtual ~GeometryInfo() = 0;

		/// Type�� �����մϴ�.
		void setType(const GEOMETRY_TYPE& Type)
		{
			_Type = Type;
		}
		/// Frame�� ��ġ�� �����մϴ�.
		void setFrame(const SE3& T)
		{
			_T = T;
		}
		/// (R, G, B, Alpha) ���� �����մϴ�.
		void setColor(const Real& R, ///< Red
			const Real& G, ///< Green
			const Real& B, ///< Blue
			const Real& Alpha = (1.0) /// Alpha
			)
		{
			_Color << R, G, B, Alpha;
		}
		/// (R, G, B, Alpha) ���� �����մϴ�.
		void setColor(const Vector4& Color // �� [R; G; B; Alpha]
			)
		{
			_Color = Color;
		}

		/// ���� geometry�� type�� ������ �ɴϴ�.
		const GEOMETRY_TYPE& getType() const
		{
			return _Type;
		}
		/**
		*	\return Vector4 [R; G; B; Alpha]
		*	\brief ���� ���� �Ǿ��ִ� ���� �˷��ݴϴ�.
		*/
		const Vector4& getColor() const
		{
			return _Color;
		}
		/**
		*	\return SE3 Frame
		*	\brief ���� ���� �Ǿ� �ִ� Frame�� �˷��ݴϴ�.
		*/
		const SE3& getTransform() const
		{
			return _T;
		}

		/// Geometry ���� ����, �ڽ� Ŭ������ ������ copy�Լ� ������ �ʿ��մϴ�.
		virtual GeometryInfoPtr copy() const = 0;

	private:
		GEOMETRY_TYPE _Type; ///< Geometry type�� �����ϴ� ����
		SE3 _T; ///< Frame�� ��ġ�� �����ϰ� �ִ� ����
		Vector4 _Color; ///< R, G, B, Alpha ���� �����ϴ� ����

	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	/**
	*	\class Box
	*	\brief Geometry box�� ������ �����ϴ� Ŭ����, GEOMETRY_TYPE = _BOX
	*/
	class Box :public GeometryInfo
	{
	public:
		/// �⺻ ������, ����=0, ����=0, ����=0�� �ʱ�ȭ�˴ϴ�.
		Box(const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : GeometryInfo(GeometryInfo::_BOX, T, Color),
			_width(0.0), _height(0.0), _depth(0.0) {}
		/// �Է����� ���� ����, ����, ���̷� �ʱ�ȭ�˴ϴ�.
		Box(const Real& width, ///< ����
			const Real& depth, ///< ����
			const Real& height, ///< ����
			const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : GeometryInfo(GeometryInfo::_BOX, T, Color),
			_width(width), _height(height), _depth(depth) {}
		/// �Է����� ���� config�� �̿��Ͽ� �ʱ�ȭ�� �մϴ�.
		Box(const Vector3& config, ///< [����; ����; ����]
			const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : GeometryInfo(GeometryInfo::_BOX, T, Color),
			_width(config(0)), _height(config(1)), _depth(config(2)) {}

		/// �Ҹ���
		~Box() {}

		/// �Է����� ���� ����, ����, ���̷� ���� �ٲߴϴ�.
		void setDimension(const Real& width, ///< ����
			const Real& depth, ///< ����
			const Real& height ///< ����
			);
		/// �Է����� ���� config�� ����, ����, ���̸� �ٲ��ݴϴ�.
		void setDimension(const Vector3& config ///< [����; ����; ����]
			);
		/// �Է����� ���� ���� �Ѻ��� ���̷� �ϴ� ������ü�� ����ϴ�.
		void setCube(const Real& length ///< ������ü�� �Ѻ��� ����
			);

		/**
		*	\return Vector3 [����; ����; ����]
		*	\brief ���� ���� �Ǿ��ִ� ����, ����, ���̸� �˷��ݴϴ�.
		*/
		Vector3 getDimension() const;

		// ���� ����
		GeometryInfoPtr copy() const;

	private:
		Real _width, _depth, _height;
	};

	/**
	*	\class Sphere
	*	\brief Geometry sphere�� ������ �����ϴ� Ŭ����, GEOMETRY_TYPE = _SPHERE
	*/
	class Sphere :public GeometryInfo
	{
	public:
		/// �⺻ ������, ������=0�� �ʱ�ȭ�˴ϴ�.
		Sphere(const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : GeometryInfo(GeometryInfo::_SPHERE, T, Color),
			_radius(0.0) {}
		/// �Է����� ���� ���������� �ʱ�ȭ�˴ϴ�.
		Sphere(const Real& radius, ///< ������
			const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : GeometryInfo(GeometryInfo::_SPHERE, T, Color),
			_radius(radius) {}

		/// �Ҹ���
		~Sphere() {}

		/// �Է����� ���� radius�� ������ ���� �ٲߴϴ�.
		void setRadius(const Real& radius ///< ������
			);

		/**
		*	\return ������
		*	\brief ���� ���� �Ǿ��ִ� ������ ���� �˷��ݴϴ�.
		*/
		const Real& getRadius() const;

		// ���� ����
		GeometryInfoPtr copy() const;

	private:
		Real _radius;
	};

	/**
	*	\class Capsule
	*	\brief Geometry capsule�� ������ �����ϴ� Ŭ����, GEOMETRY_TYPE = _CAPSULE
	*/
	class Capsule :public GeometryInfo
	{
	public:
		/// �⺻ ������, ������=0, ����=0�� �ʱ�ȭ�˴ϴ�.
		Capsule(const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : GeometryInfo(GeometryInfo::_CAPSULE, T, Color),
			_radius(0.0), _height(0.0) {}
		/// �Է����� ���� �������� ���̷� �ʱ�ȭ�˴ϴ�.
		Capsule(const Real& radius, ///< ������
			const Real& height, ///< ����
			const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : GeometryInfo(GeometryInfo::_CAPSULE, T, Color),
			_radius(radius), _height(height) {}

		/// �Ҹ���
		~Capsule() {}

		/// �Է����� ���� ������, ���̷� ���� �ٲߴϴ�.
		void setDimension(const Real& radius, ///< ����
			const Real& height ///< ����
			);
		/// �Է����� ���� config�� ������, ���̸� �ٲ��ݴϴ�.
		void setDimension(const Vector2& config ///< [������: ����]
			);

		/**
		*	\return Vector2 [������; ����]
		*	\brief ���� ���� �Ǿ��ִ� ������, ���̸� �˷��ݴϴ�.
		*/
		Vector2 getDimension() const;

		// ���� ����
		GeometryInfoPtr copy() const;

	private:
		Real _radius, _height;
	};

	/**
	*	\class Cylinder
	*	\brief Geometry cylinder�� ������ �����ϴ� Ŭ����, GEOMETRY_TYPE = _CYLINDER
	*/
	class Cylinder :public GeometryInfo
	{
	public:
		/// �⺻ ������, ������=0, ����=0�� �ʱ�ȭ�˴ϴ�.
		Cylinder(const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : GeometryInfo(GeometryInfo::_CYLINDER, T, Color),
			_radius(0.0), _height(0.0) {}
		/// �Է����� ���� �������� ���̷� �ʱ�ȭ�˴ϴ�.
		Cylinder(const Real& radius, ///< ������
			const Real& height, ///< ����
			const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : GeometryInfo(GeometryInfo::_CYLINDER, T, Color),
			_radius(radius), _height(height) {}

		/// �Ҹ���
		~Cylinder() {}

		/// �Է����� ���� ������, ���̷� ���� �ٲߴϴ�.
		void setDimension(const Real& radius, ///< ����
			const Real& height ///< ����
			);
		/// �Է����� ���� config�� ������, ���̸� �ٲ��ݴϴ�.
		void setDimension(const Vector2& config ///< [������: ����]
			);

		/**
		*	\return Vector2 [������; ����]
		*	\brief ���� ���� �Ǿ��ִ� ������, ���̸� �˷��ݴϴ�.
		*/
		Vector2 getDimension() const;

		// ���� ����
		GeometryInfoPtr copy() const;

	private:
		Real _radius, _height;
	};

	/**
	*	\class Mesh
	*	\brief Geometry mesh�� ������ �����ϴ� Ŭ����, GEOMETRY_TYPE = _MESH
	*/
	class Mesh :public GeometryInfo
	{
	public:
		/// �⺻ ������, �ּҴ� NULL�� �ʱ�ȭ �˴ϴ�,
		Mesh(const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : GeometryInfo(GeometryInfo::_MESH, T, Color),
			_url("") {}
		/// Mesh ������ ����Ǿ� �ִ� �ּ� ���� �̿��Ͽ� �ʱ�ȭ�մϴ�.
		Mesh(const std::string& url, ///< Mesh ������ ����Ǿ� �ִ� �ּ�
			const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : GeometryInfo(GeometryInfo::_MESH, T, Color),
			_url(url) {}

		/// Mesh ������ ����Ǿ� �ִ� �ּ� ���� �޾� �����մϴ�.
		void setUrl(const std::string& url ///< Mesh ������ ����Ǿ� �ִ� �ּ�
			)
		{
			_url = url;
		}

		/// ���� ���� �Ǿ��ִ� �ּ� ���� �˷��ݴϴ�.
		const std::string& getUrl() const
		{
			return _url;
		}

		void		setDimension(Real scale) { _scale = scale; }
		Real	getDimension() const { return _scale; }

		/**
		*	\return �����ϸ� True, �������� ������ False
		*	\brief ���� ���� �Ǿ��ִ� �ּҿ� ������ �����ϴ��� �Ǵ�
		*/
		bool isValid() const;

		// ���� ����
		GeometryInfoPtr copy() const;

	private:
		std::string _url;
		Real	_scale = 1;
	};

	/**
	*	\class UserModel
	*	\brief Geometry ������ ����ڰ� ������� ���� �� �ִ� Ŭ����, GEOMETRY_TYPE = _USERMODEL
	*	\todo addList ���ѷ��� ���� ���� �ʰ� �ϱ�, �˻�, ���� �Լ� �����
	*/
	class UserModel :public GeometryInfo
	{
	public:
		/// ������
		UserModel(const SE3& T = (SE3()), ///< Frame��ġ
			const Vector4& Color = (Vector4(-1, -1, -1, -1)) ///< (R, G, B, Alpha) ��
			) : GeometryInfo(GeometryInfo::_USERMODEL, T, Color), _GeometryList() {}

		/// srGeometry�ּҿ� SE3���� �޾Ƽ� �����մϴ�.
		void addList(const GeometryInfoPtr& shape, const SE3& T);

		/// List�κ��� ���� ã�ų� ������ �ٲٰ� ���� �� ����մϴ�.
		std::list< std::pair< GeometryInfoPtr, SE3 > >& getList()
		{
			return _GeometryList;
		}

		// ���� ����
		GeometryInfoPtr copy() const;

	private:
		std::list< std::pair< GeometryInfoPtr, SE3 > > _GeometryList; ///< Geometry�� �����ϰ� �ִ� list ����
	};
}