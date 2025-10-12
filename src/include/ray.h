#pragma once
#include "core.h"

struct Ray
{
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_.norm()) {}
};