#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <Eigen/Dense>
#define EIGEN_DEFAULT_TO_COLUMN_MAJOR

namespace Geom {

   class Rayf {
   public:
      Eigen::Vector3f start;
      Eigen::Vector3f direction;

      Rayf();
      Rayf(Eigen::Vector3f pnt, Eigen::Vector3f dir);
      Eigen::Vector3f getPointByDist(float dist);
      float distToPoint(Eigen::Vector3f pnt);
      float squaredDistToPoint(Eigen::Vector3f pnt);
   };

   class Planef {
   public:
      Eigen::Vector3f point;
      Eigen::Vector3f normal;

      Planef();
      Planef(Eigen::Vector3f pnt, Eigen::Vector3f norm);
      // Constructs a plane who's normal is defined by (c - b)x(a - b)
      // Wind counterclockwise for normal to point towards you
      Planef(Eigen::Vector3f a, Eigen::Vector3f b, Eigen::Vector3f c);

      // Returns a positive distance if the point is on the same side as the direction of the normal.
      // Returns a negative distance if the point is on the other side of the plane
      float distToPoint(Eigen::Vector3f pnt);
   };

   class Spheref {
   public:
      Eigen::Vector3f center;
      float radius;

      Spheref();
      Spheref(Eigen::Vector3f cent, float rad);
   };

   class Trianglef {
   public:
      Eigen::Vector3f pA, pB, pC;
      Eigen::Vector3f normal;

      Trianglef();
      // Constructs a triangle who's normal is defined by (c - b)x(a - b)
      // Wind counterclockwise for normal to point towards you
      Trianglef(Eigen::Vector3f pA, Eigen::Vector3f pB, Eigen::Vector3f pC);
      Trianglef(Eigen::Vector3f pA, Eigen::Vector3f pB, Eigen::Vector3f pC, Eigen::Vector3f normal);

      bool isPointInside(Eigen::Vector3f pnt);
   };

   class AABBf {
   public:
      Eigen::Vector3f lowBound, highBound;

      AABBf();
      AABBf(Eigen::Vector3f lowBound, Eigen::Vector3f highBound);
   };

   class Frustumf {
   public:
      Planef left, right;
      Planef bottom, top;
      Planef near, far;

      Frustumf();
      Frustumf(Planef left, Planef right, Planef bottom, Planef top, Planef near, Planef far);
      Frustumf(Eigen::Vector3f nbl, Eigen::Vector3f nbr, Eigen::Vector3f ntl, Eigen::Vector3f ntr,
               Eigen::Vector3f fbl, Eigen::Vector3f fbr, Eigen::Vector3f ftl, Eigen::Vector3f ftr);
      bool contains(Eigen::Vector3f pnt);
   };

   // Find the point in which the ray intersects the plane
   bool DoesIntersect(Rayf& ray, Planef& plane);
   bool DoesIntersect(Rayf& ray, Trianglef& triangle);
   bool DoesIntersect(Rayf& ray, Spheref& sphere);
   bool DoesIntersect(Rayf& ray, AABBf& box);
   bool DoesIntersect(Spheref& sphere, AABBf& box);
   bool DoesIntersect(Trianglef& triangle, AABBf& box);

   Eigen::Vector3f Intersect(Rayf& ray, Planef& plane);
   Eigen::Vector3f Intersect(Rayf& ray, Trianglef& triangle); // Aame as the above test (ignores the boundaries of the triangle)
   Eigen::Vector3f Intersect(Rayf& ray, Spheref& sphere);
   Eigen::Vector3f Intersect(Rayf& ray, AABBf& box);
}

#endif // __GEOMETRY_H__
