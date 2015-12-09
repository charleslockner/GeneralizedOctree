
#include "test.h"
#include "geometry.h"
#include "octree.h"

using namespace Eigen;
using namespace Geom;

int main() {
   printf("Testing geometry\n");

   // test Rayf's getPointByDist method
   {
      Rayf ray(Vector3f(0,0,0), Vector3f(1,1,0).normalized());
      Vector3f pnt = ray.getPointByDist(sqrt(2));
      equalityFloatCheck(pnt(0), 1, 1e-5);
      equalityFloatCheck(pnt(1), 1, 1e-5);
      equalityFloatCheck(pnt(2), 0, 1e-5);
   }

   // Test Rayf's squaredDistToPoint method
   {
      Vector3f pnt = Vector3f(1,1,1);
      Rayf ray(Vector3f(0,0,0), Vector3f(0,1,0).normalized());
      float res = ray.squaredDistToPoint(pnt);
      equalityFloatCheck(res, 2, 1e-5);
   }

   // Test Rayf's distToPoint method
   {
      Vector3f pnt = Vector3f(0,0,0);
      Rayf ray(Vector3f(0,-1,0), Vector3f(0,1,0).normalized());
      float res = ray.distToPoint(pnt);
      equalityFloatCheck(res, 0, 1e-5);
   }

   {
      Vector3f pnt = Vector3f(1,1,1);
      Rayf ray(Vector3f(0,0,0), Vector3f(0,1,0).normalized());
      float res = ray.distToPoint(pnt);
      equalityFloatCheck(res, sqrt(2), 1e-5);
   }

   // Test Planef's distToPoint method
   {
      Vector3f pnt = Vector3f(0,0,0);
      Planef plane(Vector3f(0,-1,0), Vector3f(0,1,0).normalized());
      float res = plane.distToPoint(pnt);
      equalityFloatCheck(res, 1, 1e-5);
   }

   {
      Vector3f pnt = Vector3f(0,0,0);
      Planef plane(Vector3f(0,-1,0), Vector3f(0,-1,0).normalized());
      float res = plane.distToPoint(pnt);
      equalityFloatCheck(res, -1, 1e-5);
   }

   {
      Vector3f pnt = Vector3f(0,0,0);
      Planef plane(Vector3f(-2,-2,0), Vector3f(-1,-1,0).normalized());
      float res = plane.distToPoint(pnt);
      equalityFloatCheck(res, -sqrt(8), 1e-5);
   }

   printf("Testing DoesIntersect Functions\n");

   // Test DoesIntersect (ray and plane)
   {
      Rayf ray(Vector3f(2,2,0), Vector3f(0,-1,0));
      Planef plane(Vector3f(-1,-1,-1), Vector3f(0,1,0));
      boolCheck(DoesIntersect(ray, plane), true);
   }

   {
      Rayf ray(Vector3f(2,2,2), Vector3f(1,0,0));
      Planef plane(Vector3f(-1,-1,-2), Vector3f(0,-1,0));
      boolCheck(DoesIntersect(ray, plane), false);
   }

   {
      Rayf ray(Vector3f(1,0,0), Vector3f(1,1,0));
      Planef plane(Vector3f(0,-1,0), Vector3f(-1,-1,0).normalized());
      boolCheck(DoesIntersect(ray, plane), true);
   }

   {
      Rayf ray(Vector3f(1,0,0), Vector3f(1,1,0));
      Planef plane(Vector3f(0,-1,0), Vector3f(1,1,0).normalized());
      boolCheck(DoesIntersect(ray, plane), true);
   }

   // Test DoesIntersect (ray and triangle)
   {
      Rayf ray(Vector3f(0,2,-1), Vector3f(0,-1,0));
      Trianglef tri(Vector3f(-2,0,0), Vector3f(2,0,0), Vector3f(0,0,2));
      boolCheck(DoesIntersect(ray, tri), false);
   }

   {
      Rayf ray(Vector3f(0,2,1), Vector3f(0,-1,0));
      Trianglef tri(Vector3f(-2,0,0), Vector3f(2,0,0), Vector3f(0,0,2));
      boolCheck(DoesIntersect(ray, tri), true);
   }

   // Test DoesIntersect (ray and sphere)
   {
      Rayf ray(Vector3f(1,1,0), Vector3f(1,0,0));
      Spheref sphere(Vector3f(0,0,0), 0.5);
      boolCheck(DoesIntersect(ray, sphere), false);
   }

   {
      Rayf ray(Vector3f(1,0.5,0), Vector3f(-1,0,0));
      Spheref sphere(Vector3f(0,0,0), 0.5);
      boolCheck(DoesIntersect(ray, sphere), true);
   }

   {
      Rayf ray(Vector3f(1,0.5,0), Vector3f(-1,0,0));
      Spheref sphere(Vector3f(0,0,0), 0.6);
      boolCheck(DoesIntersect(ray, sphere), true);
   }

   {
      Rayf ray(Vector3f(1,0,0), Vector3f(-1,0,0));
      Spheref sphere(Vector3f(0,0,0), 0.5);
      boolCheck(DoesIntersect(ray, sphere), true);
   }

   {
      Rayf ray(Vector3f(4,5,0), Vector3f(-1,0,0));
      Spheref sphere(Vector3f(0,0,0), 0.5);
      boolCheck(DoesIntersect(ray, sphere), false);
   }

   {
      Rayf ray(Vector3f(3,4,1), Vector3f(-1,0,0));
      Spheref sphere(Vector3f(-1,-1,1), 0.5);
      boolCheck(DoesIntersect(ray, sphere), false);
   }

   {
      Rayf ray(Vector3f(1,0,0), Vector3f(-1,1,0).normalized());
      Spheref sphere(Vector3f(0,1,0), 0.25);
      boolCheck(DoesIntersect(ray, sphere), true);
   }

   // Test DoesIntersect (ray and box)
   {
      Rayf ray(Vector3f(-1.1,0,0), Vector3f(1,1,0).normalized());
      AABBf box(Vector3f(0,0,0), Vector3f(1,1,1));
      boolCheck(DoesIntersect(ray, box), false);
   }

   {
      Rayf ray(Vector3f(-0.9,0,0), Vector3f(1,1,0).normalized());
      AABBf box(Vector3f(0,0,0), Vector3f(1,1,1));
      boolCheck(DoesIntersect(ray, box), true);
   }

   {
      Rayf ray(Vector3f(0,0,0), Vector3f(1,0,0).normalized());
      AABBf box(Vector3f(0,0,0), Vector3f(1,1,1));
      boolCheck(DoesIntersect(ray, box), true);
   }

   {
      Rayf ray(Vector3f(0,1,0), Vector3f(1,0,0).normalized());
      AABBf box(Vector3f(0,0,0), Vector3f(1,1,1));
      boolCheck(DoesIntersect(ray, box), true);
   }

   // Test DoesIntersect (sphere and box)
   {
      Spheref sphere(Vector3f(0.5,0.5,0.5), 0.5);
      AABBf box(Vector3f(0,0,0), Vector3f(1,1,1));
      boolCheck(DoesIntersect(sphere, box), true);
   }

   {
      Spheref sphere(Vector3f(0,0,0), 0.5);
      AABBf box(Vector3f(0,0,0), Vector3f(1,1,1));
      boolCheck(DoesIntersect(sphere, box), true);
   }

   {
      Spheref sphere(Vector3f(-0.5,0,0), 0.5);
      AABBf box(Vector3f(0,0,0), Vector3f(1,1,1));
      boolCheck(DoesIntersect(sphere, box), true);
   }

   {
      Spheref sphere(Vector3f(-0.475,-0.475,0.0), 0.5); // 0.672 away from corner of box
      AABBf box(Vector3f(0,0,0), Vector3f(1,1,1));
      boolCheck(DoesIntersect(sphere, box), false);
   }

   printf("Testing Intersect Functions\n");

   // Test Intersect (ray and plane)
   {
      Rayf ray(Vector3f(2,2,0), Vector3f(0,-1,0));
      Planef plane(Vector3f(-1,-1,-1), Vector3f(0,1,0));
      Vector3f res = Intersect(ray, plane);
      equalityFloatCheck(res(0), 2, 1e-5);
      equalityFloatCheck(res(1), -1, 1e-5);
      equalityFloatCheck(res(2), 0, 1e-5);
   }

   {
      Rayf ray(Vector3f(2,2,2), Vector3f(1,0,0));
      Planef plane(Vector3f(-1,-1,-2), Vector3f(0,-1,0));
      Vector3f res = Intersect(ray, plane);
      nanCheck(res(0));
      nanCheck(res(1));
      nanCheck(res(2));
   }

   {
      Rayf ray(Vector3f(1,0,0), Vector3f(1,1,0));
      Planef plane(Vector3f(0,-1,0), Vector3f(-1,-1,0).normalized());
      Vector3f res = Intersect(ray, plane);
      equalityFloatCheck(res(0), 0, 1e-5);
      equalityFloatCheck(res(1), -1, 1e-5);
      equalityFloatCheck(res(2), 0, 1e-5);
   }

   {
      Rayf ray(Vector3f(1,0,0), Vector3f(1,1,0));
      Planef plane(Vector3f(0,-1,0), Vector3f(1,1,0).normalized());
      Vector3f res = Intersect(ray, plane);
      equalityFloatCheck(res(0), 0, 1e-5);
      equalityFloatCheck(res(1), -1, 1e-5);
      equalityFloatCheck(res(2), 0, 1e-5);
   }

   // Test Intersect (ray and sphere)
   {
      Rayf ray(Vector3f(1,1,0), Vector3f(1,0,0));
      Spheref sphere(Vector3f(0,0,0), 0.5);
      Vector3f res = Intersect(ray, sphere);
      nanCheck(res(0));
      nanCheck(res(1));
      nanCheck(res(2));
   }

   {
      Rayf ray(Vector3f(1,0.5,0), Vector3f(-1,0,0));
      Spheref sphere(Vector3f(0,0,0), 0.5);
      Vector3f res = Intersect(ray, sphere);
      equalityFloatCheck(res(0), 0, 1e-5);
      equalityFloatCheck(res(1), 0.5, 1e-5);
      equalityFloatCheck(res(2), 0, 1e-5);
   }

   {
      Rayf ray(Vector3f(1,0,0), Vector3f(-1,0,0));
      Spheref sphere(Vector3f(0,0,0), 0.5);
      Vector3f res = Intersect(ray, sphere);
      equalityFloatCheck(res(0), 0.5, 1e-5);
      equalityFloatCheck(res(1), 0, 1e-5);
      equalityFloatCheck(res(2), 0, 1e-5);
   }

   {
      Rayf ray(Vector3f(4,5,0), Vector3f(-1,0,0));
      Spheref sphere(Vector3f(0,0,0), 0.5);
      Vector3f res = Intersect(ray, sphere);
      nanCheck(res(0));
      nanCheck(res(1));
      nanCheck(res(2));
   }

   printf("Testing octree\n");

   return 0;
}