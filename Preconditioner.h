#pragma once
#include <iostream>
#include "Vector.h"
#include "Matrix.h"

using namespace std;

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
};

class JacobiPreconditioner : public Preconditioner
{
    private:
    Vector inv_diag;

    public:
    explicit JacobiPreconditioner (const Matrix& A);

    Vector apply (const Vector& x) const override;

    Vector matvec (const Vector& v) const override;
    Vector rmatvec (const Vector& v) const override;
};