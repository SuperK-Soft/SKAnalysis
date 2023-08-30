#include "Particle.h"

Particle::~Particle() = default;

Particle::Particle(int id, float time, TVector3 vtx, TVector3 mom, int parPID, unsigned int interactionID)
: pid(id), t(time), v(vtx), p(mom), parentPID(parPID), intID(interactionID) {}
