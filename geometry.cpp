
#include "geometry.h"
#include "matrix_math.h"

namespace Geom {

   // ================================================================== //
   // ====================== Geometry Declarations ===================== //
   // ================================================================== //

   Rayf::Rayf() : start(Eigen::Vector3f(0,0,0)), direction(Eigen::Vector3f(0,1,0)) {}

   Rayf::Rayf(Eigen::Vector3f pnt, Eigen::Vector3f dir)
   : start(pnt), direction(dir) {}

   Eigen::Vector3f Rayf::getPointByDist(float dist) {
      return start + dist * direction;
   }

   float Rayf::distToPoint(Eigen::Vector3f pnt) {
      return (pnt - start).cross(pnt - start - direction).norm();
   }

   float Rayf::squaredDistToPoint(Eigen::Vector3f pnt) {
      return (pnt - start).cross(pnt - start - direction).squaredNorm();
   }

   Planef::Planef() : point(Eigen::Vector3f(0,0,0)), normal(Eigen::Vector3f(0,1,0)) {}

   Planef::Planef(Eigen::Vector3f pnt, Eigen::Vector3f norm)
   : point(pnt), normal(norm) {}

   Planef::Planef(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c) {
      point = b;
      normal = (c - b).cross(a - b).normalized();
   }

   float Planef::distToPoint(Eigen::Vector3f pnt) {
      return normal.dot(pnt - point);
   }

   Planef operator *(Eigen::Matrix4f transformM, Planef planeOp) {
      return Planef((transformM * Mmath::vec3To4(planeOp.point, 1.0f)).head<3>(),
                    (transformM * Mmath::vec3To4(planeOp.normal, 0.0f)).head<3>());
   }

   Spheref::Spheref() : center(Eigen::Vector3f(0.0f,0.0f,0.0f)), radius(1.0f) {}

   Spheref::Spheref(Eigen::Vector3f cent, float rad)
   : center(cent), radius(rad) {}

   Trianglef::Trianglef(Eigen::Vector3f pA, Eigen::Vector3f pB, Eigen::Vector3f pC)
   : pA(pA), pB(pB), pC(pC) {
      normal = (pC - pB).cross(pA - pB);
   }

   Trianglef::Trianglef() {}

   Trianglef::Trianglef(Eigen::Vector3f pA, Eigen::Vector3f pB, Eigen::Vector3f pC, Eigen::Vector3f normal)
   : pA(pA), pB(pB), pC(pC), normal(normal) {}

   bool Trianglef::isPointInside(Eigen::Vector3f pnt) {
      bool inside12 = (pB - pA).cross(pnt - pA).dot(normal) >= 0;
      bool inside23 = (pC - pB).cross(pnt - pB).dot(normal) >= 0;
      bool inside31 = (pA - pC).cross(pnt - pC).dot(normal) >= 0;

      return inside12 && inside23 && inside31;
   }

   AABBf::AABBf() : lowBound(Eigen::Vector3f(0.0f,0.0f,0.0f)), highBound(Eigen::Vector3f(0.0f,0.0f,0.0f)) {}

   AABBf::AABBf(Eigen::Vector3f lowBound, Eigen::Vector3f highBound)
   : lowBound(lowBound), highBound(highBound) {}

   Frustumf::Frustumf() {}

   Frustumf::Frustumf(Planef left, Planef right, Planef bottom, Planef top, Planef near, Planef far)
   : left(left), right(right), bottom(bottom), top(top), near(near), far(far) {}

   Frustumf::Frustumf(Eigen::Vector3f nbl, Eigen::Vector3f nbr, Eigen::Vector3f ntl, Eigen::Vector3f ntr,
                      Eigen::Vector3f fbl, Eigen::Vector3f fbr, Eigen::Vector3f ftl, Eigen::Vector3f ftr) {
      left     = Planef(ftl, ntl, nbl);
      right    = Planef(ntr, ftr, fbr);
      bottom   = Planef(fbr, fbl, nbl);
      top      = Planef(ntr, ntl, ftl);
      near     = Planef(ntl, ntr, nbr);
      far      = Planef(ftr, ftl, fbl);
   }

   Frustumf operator *(Eigen::Matrix4f transformM, Frustumf frustOp) {
      return Frustumf(transformM * frustOp.left, transformM * frustOp.right,
                      transformM * frustOp.bottom, transformM * frustOp.top,
                      transformM * frustOp.near, transformM * frustOp.far);
   }

   bool Frustumf::contains(Eigen::Vector3f pnt) {
      return left.distToPoint(pnt) > 0 && right.distToPoint(pnt) > 0 &&
             bottom.distToPoint(pnt) > 0 && top.distToPoint(pnt) > 0 &&
             near.distToPoint(pnt) > 0 && far.distToPoint(pnt) > 0;
   }

   // ================================================================== //
   // ======================= Intersection Tests ======================= //
   // ================================================================== //

   bool DoesIntersect(Rayf& ray, Planef& plane) {
      return plane.normal.dot(ray.direction) != 0;
   }

   bool DoesIntersect(Rayf& ray, Trianglef& triangle) {
      float t = - triangle.normal.dot(ray.start - triangle.pA) /
                  triangle.normal.dot(ray.direction);
      Eigen::Vector3f interPnt = ray.getPointByDist(t);
      return triangle.isPointInside(interPnt);
   }

   bool DoesIntersect(Rayf& ray, Spheref& sphere) {
      // float t = ray.direction.dot(sphere.center - ray.start);
      // Eigen::Vector3f pClose = ray.getPointByDist(t);
      // float radiusSq = sphere.radius * sphere.radius;
      // float pDistSq = pClose.squaredNorm();
      // return (pDistSq <= radiusSq);

      Eigen::Vector3f vpc = sphere.center - ray.start;
      Eigen::Vector3f pc = Mmath::ProjectVectorOntoVector(vpc, ray.direction); // projection of c on the line
      // printf("vpc = %f %f %f\n", vpc(0), vpc(1), vpc(2));
      // printf("dir = %f %f %f\n", ray.direction(0), ray.direction(1), ray.direction(2));
      // printf("pc =  %f %f %f\n", pc(0), pc(1), pc(2));
      float radSq = sphere.radius * sphere.radius;

      if (vpc.dot(ray.direction) < 0.0f) { // when the sphere is behind the origin p
         // printf("behind\n");
         float vpcMagSq = vpc.squaredNorm();
         if (vpcMagSq > radSq) {
            return false; // there is no intersection
         } else {
            return true;
         }
      } else { // center of sphere projects on the ray
         Eigen::Vector3f distVec = vpc - pc;
         float distSq = distVec.squaredNorm();
         // printf("distSq = %f, radSq = %f\n", distSq, radSq);
         if (distSq > radSq) {
            return false; // there is no intersection
         } else {
            return true;
         }
      }
   }

   bool DoesIntersect(Rayf& ray, AABBf& box) {
      float tmin = -100000.0;
      float tmax =  100000.0;

      if (ray.direction(0) != 0.0) {
         float tx1 = (box.lowBound(0) - ray.start(0)) / ray.direction(0);
         float tx2 = (box.highBound(0) - ray.start(0)) / ray.direction(0);

         tmin = Mmath::max(tmin, Mmath::min(tx1, tx2));
         tmax = Mmath::min(tmax, Mmath::max(tx1, tx2));
      }

      if (ray.direction(1) != 0.0) {
         float ty1 = (box.lowBound(1) - ray.start(1)) / ray.direction(1);
         float ty2 = (box.highBound(1) - ray.start(1)) / ray.direction(1);

         tmin = Mmath::max(tmin, Mmath::min(ty1, ty2));
         tmax = Mmath::min(tmax, Mmath::max(ty1, ty2));
      }

      if (ray.direction(2) != 0.0) {
         float tz1 = (box.lowBound(2) - ray.start(2)) / ray.direction(2);
         float tz2 = (box.highBound(2) - ray.start(2)) / ray.direction(2);

         tmin = Mmath::max(tmin, Mmath::min(tz1, tz2));
         tmax = Mmath::min(tmax, Mmath::max(tz1, tz2));
      }

      return tmax >= tmin;
   }

   bool DoesIntersect(Spheref& sphere, AABBf& box) {
      // False positive when center of circle lies just above or below the corners,
      // But whatever, we can fix it later
      return sphere.center(0) + sphere.radius >= box.lowBound(0) &&
             sphere.center(0) - sphere.radius <= box.highBound(0) &&
             sphere.center(1) + sphere.radius >= box.lowBound(1) &&
             sphere.center(1) - sphere.radius <= box.highBound(1) &&
             sphere.center(2) + sphere.radius >= box.lowBound(2) &&
             sphere.center(2) - sphere.radius <= box.highBound(2);
   }

   bool DoesIntersect(Trianglef& triangle, AABBf& box) {
      // fill this out
      return false;
   }

   Eigen::Vector3f Intersect(Rayf& ray, Planef& plane) {
      float t = - plane.normal.dot(ray.start - plane.point) /
                  plane.normal.dot(ray.direction);
      return ray.getPointByDist(t);
   }

   Eigen::Vector3f Intersect(Geom::Rayf& ray, Trianglef& triangle) {
      float t = - triangle.normal.dot(ray.start - triangle.pA) /
                  triangle.normal.dot(ray.direction);
      return ray.getPointByDist(t);
   }

   Eigen::Vector3f Intersect(Rayf& ray, Spheref& sphere) {
      float t = ray.direction.dot(sphere.center - ray.start);
      Eigen::Vector3f pClose = ray.getPointByDist(t);

      float radiusSq = sphere.radius * sphere.radius;
      float pDistSq = pClose.squaredNorm();

      if (pDistSq > radiusSq) // no hit
         return Eigen::Vector3f(NAN, NAN, NAN);
      else if (pDistSq == radiusSq) // hit sphere tangentially
         return pClose;
      else // hit through the sphere
         return ray.getPointByDist(t - sqrt(radiusSq - pDistSq));


      // Eigen::Vector3f vpc = sphere.center - ray.start;
      // Eigen::Vector3f pc = Mmath::ProjectVectorOntoVector(startToCenter, ray.direction); // projection of c on the line
      // float radSq = sphere.radius * sphere.radius;

      // if ((vpc.dot(ray.direction) < 0.0f) { // when the sphere is behind the origin p
      //    float vpcMagSq = vpc.normSquared();

      //    if (vpcMagSq > radSq) {
      //       return false; // there is no intersection
      //    } else if (vpcMagSq == r) {
      //       return true; // intersection = start
      //    } else { // occurs when p is inside the sphere
      //       // pc = projection of c on the line
      //       // // distance from pc to i1
      //       // dist = sqrt(radius^2 - |pc - c|^2)
      //       // di1 = dist - |pc - p|
      //       // intersection = p + d * di1
      //       return true;
      //    }
      // } else { // center of sphere projects on the ray
      //    Eigen::Vector3f cpc = sphere.center - pc;
      //    float cpcMagSq = cpc.normSquared();

      //    if (cpcMagSq > radSq) {
      //       return false; // there is no intersection
      //    } else {
      //       // // distance from pc to i1
      //       // dist = sqrt(radius^2 - |pc - c|^2)

      //       // if (|vpc| > r) { // origin is outside sphere
      //       //    di1 = |pc - p| - dist
      //       // } else { // origin is inside sphere
      //       //    di1 = |pc - p| + dist
      //       // }
      //       // intersection = p + d * di1
      //       return true;
      //    }
      // }


   }

   Eigen::Vector3f Intersect(Rayf& ray, AABBf& box) {
      return Eigen::Vector3f(0,0,0);
   }
}