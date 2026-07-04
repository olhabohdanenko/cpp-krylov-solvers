#pragma once
#include <iostream>

#include "Vector.h"
#include "Matrix.h"

class Preconditioner
{
public:
	virtual ~Preconditioner () = default;

	virtual Vector apply (const Vector& x) const = 0;

	virtual Vector matvec (const Vector& v) const = 0;
	virtual Vector rmatvec (const Vector& v) const = 0;
};

class NoPreconditioner : public Preconditioner
{
	public:
	Vector apply (const Vector& x) const override { return x; }
	Vector matvec (const Vector& v) const override { return v; }
	Vector rmatvec (const Vector& v) const override { return v; }
};
