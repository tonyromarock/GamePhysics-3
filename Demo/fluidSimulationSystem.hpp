#pragma once
#ifndef FLUIDSIMULATIONSYSTEM_H
#define FLUIDSIMULATIONSYSTEM_H

#include <vector>
#include <DirectXMath.h>
#include <AntTweakBar.h>
#include <iostream>
using namespace DirectX;
using std::cout;

class FluidSimulationSystem {
public:
	//FluidSimulationSystem();
	//virtual ~FluidSimulationSystem();

	struct Particle
	{
		XMFLOAT3 pos;
		XMFLOAT3 vel;
		float pressure;
		float density;
	};

	bool equalParticles (Particle& lhs, Particle& rhs) {
		return (lhs.pos.x == rhs.pos.x) && (lhs.pos.y == rhs.pos.y) && (lhs.pos.z == rhs.pos.z);
	}

	void setupScene(int noOfParticles);
	//void sphAlgorithm();
	//void reset();
	const std::vector<Particle>&  GetParticles() { return m_particles; }

	void SetMass(float mass) { m_mass = mass; }
	void SetKernelSize(float kernelSize) { m_kernelSize = kernelSize; }
	void SetStiffness(float stiffness) { m_stiffnessFactor = stiffness; }
	void SetExponent(float exponent) { m_exponent = exponent; }
	void SetRestDensity(int restDensity) { m_restDensity = restDensity; }
	void SetTestCase(int testCase) { m_testCase = testCase; }
	void SetVelocityDamping(float damping) { m_damping = damping; }

	void DoEulerStep(float deltaTime);
	void BoundingBoxCheck(float times = 1.0f);

private:
	std::vector<Particle>  m_particles;
	Particle** m_particleGrid;
	int* m_particlesPerCell;

	float m_mass;
	float m_kernelSize;
	XMFLOAT3 m_gravity = XMFLOAT3(0, -0.05f, 0);

	float m_exponent;
	float m_stiffnessFactor;
	int m_restDensity;
	float m_damping;
	int m_testCase;
	int m_maxPerCell = 20;
	int m_cellsPerRow;

	void UpdatePositions(float deltaTime);
	void UpdateVelocities(float deltaTime, const std::vector<XMVECTOR>& forces);
	std::vector<XMVECTOR> ComputeForces();
	void UpdateDensitiesNaive();
	void UpdateDensitiesUG();
	void UpdatePressures();
	float KernelFunction(Particle pi, Particle pj);
	float KernelFunctionGradient(Particle pi, Particle pj);

	float Distance(Particle pi, Particle pj);
	float ToNextCell(float f);
	void AddToGrid(Particle p);
	int GetPosInGrid(Particle p);
	void SortParticlesIntoGrid();
	std::vector<Particle*> GetList(int posInGrid);
};

#endif