#include "fluidSimulationSystem.hpp"

#include <GeometricPrimitive.h>
#include <Effects.h>
#include <PrimitiveBatch.h>
#include <VertexTypes.h>
using namespace DirectX;

extern DirectX::BasicEffect* g_pEffectPositionNormal;

void FluidSimulationSystem::setupScene(int noOfParticles) {

	//cellsPerRow = 1 / (2 * kernelSize) rounded up
	if (abs((int)(1 / (2 * m_kernelSize)) - (1 / (2 * m_kernelSize))) < 0.00000001f)
		m_cellsPerRow = (int)(1 / (2 * m_kernelSize));
	else
		m_cellsPerRow = (int)(1 / (2 * m_kernelSize)) + 1;

	m_particlesPerCell = (int*)malloc(sizeof(int) * (int)pow(m_cellsPerRow, 3));
	m_particleGrid = (Particle**)malloc(sizeof(Particle*) * (int)pow(m_cellsPerRow, 3) * m_maxPerCell);

	float distance = m_kernelSize;
	float length = (float) pow(noOfParticles, 1.0f/3);
	float startX = 0.2f, startY = 0.2f, startZ = 0.2f;
	float x = startX, y = startY, z = startZ;
	for (int i = 0; i < noOfParticles; i++) {
		Particle p;
		p.pos = XMFLOAT3(x, y, z);
		p.vel = XMFLOAT3(0, 0, 0);
		p.density = 0;
		p.pressure = 0;
		m_particles.push_back(p);
		if ((x -= distance) <= startX - distance * length) {
			x = startX;
			if ((y -= distance) <= startY - distance * length) {
				y = startY;
				z -= distance;
			}
		}
	}

}

void FluidSimulationSystem::DoEulerStep(float deltaTime) {

	switch (m_testCase) {
	case 1:
		UpdateDensitiesNaive();
		UpdatePressures();
		break;
	case 2:
		SortParticlesIntoGrid();
		UpdateDensitiesUG();
		UpdatePressures();
		break;
	}

	std::vector<XMVECTOR> forces;
	forces = ComputeForces();
	UpdatePositions(deltaTime);
	UpdateVelocities(deltaTime, forces);

}

void FluidSimulationSystem::UpdatePositions(float deltaTime) {

	for (Particle& p : m_particles) {
		XMVECTOR pos = XMLoadFloat3(&p.pos);
		XMVECTOR vel = XMLoadFloat3(&p.vel);

		pos += deltaTime * vel;

		XMStoreFloat3(&p.pos, pos);
	}

}

void FluidSimulationSystem::UpdateVelocities(float deltaTime, const std::vector<XMVECTOR>& forces) {

	for (int i = 0; i < m_particles.size(); i++) {
		XMVECTOR vel = XMLoadFloat3(&m_particles[i].vel);
		vel += deltaTime * (forces[i] / m_mass);
		XMStoreFloat3(&m_particles[i].vel, vel);
	}

}

void FluidSimulationSystem::UpdateDensitiesNaive() {

	for (Particle& pi : m_particles) {
		pi.density = 0;
		
		for (Particle& pj : m_particles) {
			pi.density += KernelFunction(pi, pj);
		}

		pi.density *= m_mass;
	}

}

void FluidSimulationSystem::UpdateDensitiesUG() {

	for (Particle& pi : m_particles) {
		pi.density = 0;
		
		int posInGrid = GetPosInGrid(pi);
		std::vector<Particle*> neighbours = GetList(posInGrid);

		for (Particle* pj : neighbours) {
			pi.density += KernelFunction(pi, *pj);
		}

		pi.density *= m_mass;
	}

}

void FluidSimulationSystem::UpdatePressures(){

	for (Particle& p : m_particles) {
		p.pressure = m_stiffnessFactor * (pow(p.density / m_restDensity, m_exponent) - 1);
	}

}

float FluidSimulationSystem::KernelFunction(Particle pi, Particle pj) {

	float d = m_kernelSize;
	float q = Distance(pi, pj) / d;
	float W = 3 / (2.0f * XM_PI * pow(d, 3));
	if (0 <= q && q < 1) {
		W *= (2.0f / 3) - pow(q, 2) + 0.5f * pow(q, 3);
	}
	else if (1 <= q && q < 2) {
		W *= (1.0f / 6) * pow(2 - q, 3);
	}
	else {
		W *= 0;
	}

	return W;
}

float FluidSimulationSystem::KernelFunctionGradient(Particle pi, Particle pj) {
	
	float d = m_kernelSize;
	float q = Distance(pi, pj) / d;
	float gradientW = 9 / (4 * XM_PI * pow(d, 5));
	if (0 <= q && q < 1) {
		gradientW *= (q - (4.0f / 3)) * q * d;
	}
	else if (1 <= q && q < 2) {
		gradientW *= -pow(2 - q, 2) * (d / 3);
	}
	else {
		gradientW *= 0;
	}

	return gradientW;
}
std::vector<XMVECTOR> FluidSimulationSystem::ComputeForces() {

	std::vector<XMVECTOR> forces;
	XMVECTOR gravity = XMLoadFloat3(&m_gravity);

	for (Particle& pi : m_particles) {
		XMVECTOR forceI = XMLoadFloat3(&XMFLOAT3(0, 0, 0));
		XMVECTOR posI = XMLoadFloat3(&pi.pos);
		XMVECTOR vel = XMLoadFloat3(&pi.vel);

		switch (m_testCase) {
		case 1:
			for (Particle& pj : m_particles) {
				if (!(equalParticles(pi, pj))) {
					XMVECTOR posJ = XMLoadFloat3(&pj.pos);

					//Using the given force function for some reason made the simulation explode, so I used an alternative one I found.
					//forceI += ((pressures[i] + pressures[j]) / 2) * (m_mass / densities[j]) * KernelFunctionGradient(pi, pj) * ((posI - posJ) / Distance(pi, pj));
					forceI += pow(m_mass, 2) * ((pi.pressure / pow(pi.density, 2)) + (pj.pressure / pow(pj.density, 2))) * KernelFunctionGradient(pi, pj) * ((posI - posJ) / Distance(pi, pj));
				}
			}
			break;
		case 2:
			int posInGrid = GetPosInGrid(pi);
			std::vector<Particle*> neighbours = GetList(posInGrid);

			for (Particle* pj : neighbours) {
				if (!(equalParticles(pi, *pj))) {
					XMVECTOR posJ = XMLoadFloat3(&(pj->pos));
					forceI += pow(m_mass, 2) * ((pi.pressure / pow(pi.density, 2)) + (pj->pressure / pow(pj->density, 2))) * KernelFunctionGradient(pi, *pj) * ((posI - posJ) / Distance(pi, *pj));

				}
			}
			break;
		}

		forceI *= -1;
		forceI += gravity;
		forceI += -m_damping * vel;

		forces.push_back(forceI);
	}

	return forces;
}

void FluidSimulationSystem::BoundingBoxCheck(float times) {

	for (int i = 0; i < m_particles.size(); i++) {
		XMVECTOR pos = XMLoadFloat3(&m_particles[i].pos);
		XMVECTOR vel = XMLoadFloat3(&m_particles[i].vel);

		for (int f = 0; f < 6; f++)
		{
			float sign = (f % 2 == 0) ? -1.0f : 1.0f;
			if (sign * XMVectorGetByIndex(pos, f / 2) < -0.5f * times)
			{
				pos = XMVectorSetByIndex(pos, sign * -0.5f * times, f / 2);
				vel = XMVectorSetByIndex(vel, 0, f / 2);
			}
		}

		XMStoreFloat3(&m_particles[i].pos, pos);
		XMStoreFloat3(&m_particles[i].vel, vel);
	}

}

float FluidSimulationSystem::Distance(Particle pi, Particle pj) {

	return sqrtf(pow(pi.pos.x - pj.pos.x, 2) + pow(pi.pos.y - pj.pos.y, 2) + pow(pi.pos.z - pj.pos.z, 2));

}

float FluidSimulationSystem::ToNextCell(float f) {

	return ((int)(f / (2 * m_kernelSize))) * 2 * m_kernelSize;

}

void FluidSimulationSystem::SortParticlesIntoGrid() {

	for (int i = 0; i < pow(m_cellsPerRow, 3); i++) {
		m_particlesPerCell[i] = 0;
	}

	for (Particle& p : m_particles) {
		int posInGrid = GetPosInGrid(p);
		if (m_particlesPerCell[posInGrid] < m_maxPerCell) {
			m_particleGrid[posInGrid * m_maxPerCell + m_particlesPerCell[posInGrid]] = &p;
			m_particlesPerCell[posInGrid]++;
		}
	}

}

int FluidSimulationSystem::GetPosInGrid(Particle p) {

	return (int)(ToNextCell(p.pos.x + 0.5f) / (2 * m_kernelSize) + ToNextCell(p.pos.y + 0.5f) / pow(2 * m_kernelSize, 2) + ToNextCell(p.pos.z + 0.5f) / pow(2 * m_kernelSize, 3));

}

std::vector<FluidSimulationSystem::Particle*> FluidSimulationSystem::GetList(int posInGrid) {

	std::vector<Particle*> list;
	int currentPos = posInGrid;
	/*
	for (int count = 0; count < 7; count++) {
		currentPos = posInGrid;
		switch (count) {
		case 0: break;
		case 1:
			if ((int)(currentPos * pow(2 * m_kernelSize, 2)) != 0)
				currentPos -= (int)(pow(m_cellsPerRow, 2));
			else 
				continue;
			break;
		case 2:
			if ((int)(currentPos * pow(2 * m_kernelSize, 2)) < (m_cellsPerRow - 1))
				currentPos += (int)(pow(m_cellsPerRow, 2));
			else
				continue;
			break;
		case 3:
			if ((int)((currentPos % (int)(pow(m_cellsPerRow, 2))) * 2 * m_kernelSize) != 0)
				currentPos -= (int)(m_cellsPerRow);
			else
				continue;
			break;
		case 4:
			if (((int)(currentPos % (int)(pow(m_cellsPerRow, 2))) * 2 * m_kernelSize) < (m_cellsPerRow - 1))
				currentPos += (int)(m_cellsPerRow);
			else
				continue;
			break;
		case 5:
			if (currentPos % (int)(m_cellsPerRow) != 0)
				currentPos--;
			else
				continue;
			break;
		case 6:
			if (currentPos % (int)(m_cellsPerRow) < (m_cellsPerRow - 1))
				currentPos++;
			else
				continue;
			break;
		}
		int no = m_particlesPerCell[currentPos];
		for (int i = 0; i < no; i++) {
			list.push_back(m_particleGrid[currentPos * m_maxPerCell + i]);
		}
	}
	*/
	/* For all cells adjacent to 'posInGrid' (doesn't recognize, when cell is against wall, so it may loop over a few unnecessary cells): */
	for (int i = -1; i <= 1; i++) {
		if ((i == -1 && (int)(posInGrid * pow(2 * m_kernelSize, 2)) == 0) || (i == 1 && (int)(posInGrid * pow(2 * m_kernelSize, 2)) >= (m_cellsPerRow - 1))) {
			continue;
		}
		else {
			for (int j = -1; j <= 1; j++) {
				if ((j == -1 && (int)((posInGrid % (int)(pow(m_cellsPerRow, 2))) * 2 * m_kernelSize) == 0) || (j == 1 && ((int)(posInGrid % (int)(pow(m_cellsPerRow, 2))) * 2 * m_kernelSize) >= (m_cellsPerRow - 1))) {
					continue;
				}
				else {
					for (int k = -1; k <= 1; k++) {
						if ((k == -1 && posInGrid % (int)(m_cellsPerRow) == 0) || (k == 1 && posInGrid % (int)(m_cellsPerRow) >= (m_cellsPerRow - 1))) {
							continue;
						}
						else {
							currentPos = posInGrid + i * (int)(pow(m_cellsPerRow, 2)) + j * (int)(m_cellsPerRow) + k;
							int no = m_particlesPerCell[currentPos];
							for (int l = 0; l < no; l++) {
								list.push_back(m_particleGrid[currentPos * m_maxPerCell + l]);
							}
						}
					}
				}
			}
		}
	}

	return list;
}