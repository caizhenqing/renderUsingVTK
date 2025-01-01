#pragma once
#include "algFoundation.h"
#include "transform.h"


class RigidBody
{
public:
	/*
	* ����/��ȡ RigidBody ��λ��
	*/
	virtual void getGlobalPosition(ncs::Transform& pos) const = 0;
	virtual void setGlobalPosition(const ncs::Transform& pos) = 0;

	/*
	* RigidBody �Ƿ�ֹͣ����
	*/
	virtual	bool isSleeping() const = 0;
	
	/*
	* ����/��ȡ RigidBody ������
	*/
	virtual	void setMass(const double& m) = 0;
	virtual	void getMass(double& m) const = 0;
	/*
	* ����/��ȡ RigidBody �ռ��xyz���ת������
	*/
	virtual void setMassSpaceInertiaTensor(const ncs::Vec3& v) = 0;
	virtual void getMassSpaceInertiaTensor(ncs::Vec3& v) const = 0;
	/*
	* �����������������Ŀ��ľֲ�����
	*/
	/*virtual ncs::Transform& getCMassLocalPose() = 0;
	virtual void setCMassLocalPose(const ncs::Transform& pose) = 0;*/
	/*
	* ����/��ȡ���˶�����
	*/
	virtual void getLinearDamping(double& linDamp) const = 0;
	virtual void setLinearDamping(const double& linDamp) = 0;
	/*
	* ����/��ȡת������
	*/
	virtual void setAngularDamping(const double& angDamp) = 0;
	virtual void getAngularDamping(double& angDamp) const = 0;
	/*
	* ����/��ȡ���ٶ�
	*/
	virtual void setLinearVelocity(const ncs::Vec3 &lineVel) = 0;
	virtual void getLinearVelocity(ncs::Vec3& lineVel) const = 0;
	/*
	* ����/��ȡ���ٶ�
	*/
	virtual void setAngularVelocity(const ncs::Vec3 &angVel) = 0;
	virtual void getAngularVelocity(ncs::Vec3& angVel) const = 0;
};