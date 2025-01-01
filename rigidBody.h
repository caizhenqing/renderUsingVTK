#pragma once
#include "algFoundation.h"
#include "transform.h"


class RigidBody
{
public:
	/*
	* 设置/获取 RigidBody 的位置
	*/
	virtual void getGlobalPosition(ncs::Transform& pos) const = 0;
	virtual void setGlobalPosition(const ncs::Transform& pos) = 0;

	/*
	* RigidBody 是否停止计算
	*/
	virtual	bool isSleeping() const = 0;
	
	/*
	* 设置/获取 RigidBody 的质量
	*/
	virtual	void setMass(const double& m) = 0;
	virtual	void getMass(double& m) const = 0;
	/*
	* 设置/获取 RigidBody 空间对xyz轴的转动惯量
	*/
	virtual void setMassSpaceInertiaTensor(const ncs::Vec3& v) = 0;
	virtual void getMassSpaceInertiaTensor(ncs::Vec3& v) const = 0;
	/*
	* 设置质量中心相对于目标的局部坐标
	*/
	/*virtual ncs::Transform& getCMassLocalPose() = 0;
	virtual void setCMassLocalPose(const ncs::Transform& pose) = 0;*/
	/*
	* 设置/获取线运动阻尼
	*/
	virtual void getLinearDamping(double& linDamp) const = 0;
	virtual void setLinearDamping(const double& linDamp) = 0;
	/*
	* 设置/获取转动阻尼
	*/
	virtual void setAngularDamping(const double& angDamp) = 0;
	virtual void getAngularDamping(double& angDamp) const = 0;
	/*
	* 设置/获取线速度
	*/
	virtual void setLinearVelocity(const ncs::Vec3 &lineVel) = 0;
	virtual void getLinearVelocity(ncs::Vec3& lineVel) const = 0;
	/*
	* 设置/获取角速度
	*/
	virtual void setAngularVelocity(const ncs::Vec3 &angVel) = 0;
	virtual void getAngularVelocity(ncs::Vec3& angVel) const = 0;
};