#include "Sylinder.hpp"
#include "Util/Base64.hpp"

/*****************************************************
 *  Sphero-cylinder
 ******************************************************/

Sylinder::Sylinder(const int &gid_, const double &radius_, const double &radiusCollision_, const double &length_,
                   const double &lengthCollision_, const double pos_[3], const double orientation_[4]) {
    gid = gid_;
    radius = radius_;
    radiusCollision = radiusCollision_;
    length = length_;
    lengthCollision = lengthCollision_;
    if (pos_ == nullptr) {
        Emap3(pos).setZero();
    } else {
        for (int i = 0; i < 3; i++) {
            pos[i] = pos_[i];
        }
    }
    if (orientation_ == nullptr) {
        Emapq(orientation).setIdentity();
    } else {
        for (int i = 0; i < 4; i++) {
            orientation[i] = orientation_[i];
        }
    }

    clear();
    return;
}

void Sylinder::clear() {
    Emap3(vel).setZero();
    Emap3(omega).setZero();
    Emap3(velCol).setZero();
    Emap3(omegaCol).setZero();
    Emap3(velBi).setZero();
    Emap3(omegaBi).setZero();
    Emap3(velNonB).setZero();
    Emap3(omegaNonB).setZero();

    Emap3(force).setZero();
    Emap3(torque).setZero();
    Emap3(forceCol).setZero();
    Emap3(torqueCol).setZero();
    Emap3(forceBi).setZero();
    Emap3(torqueBi).setZero();
    Emap3(forceNonB).setZero();
    Emap3(torqueNonB).setZero();

    Emap3(velBrown).setZero();
    Emap3(omegaBrown).setZero();

    sepmin = std::numeric_limits<double>::max();
    globalIndex = GEO_INVALID_INDEX;
    rank = -1;
}

bool Sylinder::isSphere(bool collision) const {
    if (collision) {
        return lengthCollision < radiusCollision * 2;
    } else {
        return length < radius * 2;
    }
}

void Sylinder::calcDragCoeff(const double viscosity, double &dragPara, double &dragPerp, double &dragRot) const {
    if (isSphere()) { // use drag for sphere
        const double rad = 0.5 * length + radius;
        dragPara = 6 * Pi * rad * viscosity;
        dragPerp = dragPara;
        dragRot = 8 * Pi * rad * rad * rad * viscosity;
    } else { // use spherocylinder drag and slender body theory
        const double b = -(1 + 2 * log(radius / (length)));
        dragPara = 8 * Pi * length * viscosity / (2 * b);
        dragPerp = 8 * Pi * length * viscosity / (b + 2);
        dragRot = 2 * Pi * viscosity * length * length * length / (3 * (b + 2));
    }
    return;
}

void Sylinder::dumpSylinder() const {
    printf("gid %d, R %g, RCol %g, L %g, LCol %g, pos %g, %g, %g\n", gid, radius, radiusCollision, length,
           lengthCollision, pos[0], pos[1], pos[2]);
    printf("vel %g, %g, %g; omega %g, %g, %g\n", vel[0], vel[1], vel[2], omega[0], omega[1], omega[2]);
    printf("orient %g, %g, %g, %g\n", orientation[0], orientation[1], orientation[2], orientation[3]);
}

void Sylinder::stepEuler(double dt) {
    Emap3(pos) += Emap3(vel) * dt;
    Equatn currOrient = Emapq(orientation);
    EquatnHelper::rotateEquatn(currOrient, Emap3(omega), dt);
    Emapq(orientation).x() = currOrient.x();
    Emapq(orientation).y() = currOrient.y();
    Emapq(orientation).z() = currOrient.z();
    Emapq(orientation).w() = currOrient.w();
}

void Sylinder::writeAscii(FILE *fptr) const {
    Evec3 direction = ECmapq(orientation) * Evec3(0, 0, 1);
    Evec3 minus = ECmap3(pos) - 0.5 * length * direction;
    Evec3 plus = ECmap3(pos) + 0.5 * length * direction;

    char typeChar = isImmovable ? 'S' : 'C';
    if (isAnchor) {
        typeChar = 'A';
    }
    fprintf(fptr, "%c %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %d\n", //
            typeChar, gid, radius,                                 //
            minus[0], minus[1], minus[2], plus[0], plus[1], plus[2], group);
}
