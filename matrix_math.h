
#ifndef __MATRIX_MATH_H__
#define __MATRIX_MATH_H__

#include <Eigen/Dense>
#define EIGEN_DEFAULT_TO_COLUMN_MAJOR

namespace Mmath {
   template <typename T>
   T min(const T a, const T b) {
      return a < b ? a : b;
   }

   template <typename T>
   T max(const T a, const T b) {
      return a > b ? a : b;
   }

   template <typename T>
   T clamp(
      const T low,
      const T high,
      T val
   ) {
      if (val < low)
         val = low;
      else if (val > high)
         val = high;
      return val;
   }

   // onto must be a unit vector
   template <typename T>
   Eigen::Matrix<T,3,1> ProjectVectorOntoVector(
      Eigen::Matrix<T,3,1>& projectee,
      Eigen::Matrix<T,3,1>& onto
   ) {
      return projectee.dot(onto) * onto;
   }

   template <typename T>
   Eigen::Matrix<T,3,1> ProjectVectorOntoPlane(
      Eigen::Matrix<T,3,1>& vector,
      Eigen::Matrix<T,3,1>& normal
   ) {
      Eigen::Matrix<T,3,1> v = normal * vector.dot(normal);
      return vector - v;
   }

   template <typename T>
   Eigen::Matrix<T,4,4> TranslationMatrix(
      const Eigen::Matrix<T,3,1> tns
   ) {
      Eigen::Matrix<T,4,4> m = Eigen::Matrix<T,4,4>::Identity();
      m.block(0,3,3,1) = tns;
      return m;
   }

   template <typename T>
   Eigen::Matrix<T,4,4> RotationMatrix(
      const Eigen::Quaternion<T> rot
   ) {
      Eigen::Matrix<T,4,4> m = Eigen::Matrix<T,4,4>::Identity();
      m.block(0,0,3,3) = rot.toRotationMatrix();
      return m;
   }

   template <typename T>
   Eigen::Matrix<T,4,4> ScaleMatrix(
      const Eigen::Matrix<T,3,1> scl
   ) {
      Eigen::Matrix<T,4,4> m = Eigen::Matrix<T,4,4>::Identity();
      m(0,0) = scl(0);
      m(1,1) = scl(1);
      m(2,2) = scl(2);
      return m;
   }

   template <typename T>
   Eigen::Matrix<T,4,4> TransformationMatrix(
      const Eigen::Matrix<T,3,1> tns,
      const Eigen::Quaternion<T> rot,
      const Eigen::Matrix<T,3,1> scl
   ) {
      return TranslationMatrix(tns) * RotationMatrix(rot) * ScaleMatrix(scl);
   }

   template <typename T>
   Eigen::Quaternion<T> AngleAxisQuat(T angle, Eigen::Matrix<T,3,1> axis) {
      Eigen::Quaternion<T> rotQuat(Eigen::AngleAxis<T>(angle, axis));
      return rotQuat;
   }

   template <typename T>
   Eigen::Matrix<T,4,4> AngleAxisMatrix4(T angle, Eigen::Matrix<T,3,1> axis) {
      return RotationMatrix(AngleAxisQuat(angle, axis));
   }

   template <typename T>
   Eigen::Matrix<T,3,3> AngleAxisMatrix3(T angle, Eigen::Matrix<T,3,1> axis) {
      Eigen::Matrix<T,3,3> mat;
      mat = Eigen::AngleAxis<T>(angle, axis);;
      return mat;
   }

   template <typename T>
   Eigen::Matrix<T,3,1> RotateVec3(
      const Eigen::Matrix<T,3,1> subject,
      const T angle,
      const Eigen::Matrix<T,3,1> axis
   ) {
      Eigen::Matrix<T,3,3> rotM;
      rotM = Eigen::AngleAxis<T>(angle, axis);
      return rotM * subject;
   }

   template <typename T>
   Eigen::Matrix<T,3,1> SlerpVec3(
      const Eigen::Matrix<T,3,1> a,
      const Eigen::Matrix<T,3,1> b,
      T ratio
   ) {
      ratio = clamp(T(0),T(1),ratio);
      return ((T(1) - ratio) * a + ratio * b).normalized();
   }

   template <typename T>
   Eigen::Matrix<T,4,4> PerspectiveMatrix(
      const T fovy,
      const T aspect,
      const T zNear,
      const T zFar
   ) {
      T tanHalfFovy = tan(fovy / static_cast<T>(2));

      Eigen::Matrix<T,4,4> Result;
      Result.setZero();
      Result(0,0) = static_cast<T>(1) / (aspect * tanHalfFovy);
      Result(1,1) = static_cast<T>(1) / (tanHalfFovy);
      Result(2,2) = - (zFar + zNear) / (zFar - zNear);
      Result(3,2) = - static_cast<T>(1);
      Result(2,3) = - (static_cast<T>(2) * zFar * zNear) / (zFar - zNear);

      return Result;
   }

   template <typename T>
   Eigen::Matrix<T,4,4> ViewMatrix(
      const Eigen::Matrix<T,3,1> eye,
      const Eigen::Matrix<T,3,1> direction,
      const Eigen::Matrix<T,3,1> up
   ) {
      Eigen::Matrix<T,3,1> f(direction);
      Eigen::Matrix<T,3,1> s((f.cross(up)).normalized());
      Eigen::Matrix<T,3,1> u(s.cross(f));

      Eigen::Matrix<T,4,4> Result = Eigen::Matrix<T,4,4>::Identity();
      Result(0,0) = s(0);
      Result(0,1) = s(1);
      Result(0,2) = s(2);
      Result(1,0) = u(0);
      Result(1,1) = u(1);
      Result(1,2) = u(2);
      Result(2,0) =-f(0);
      Result(2,1) =-f(1);
      Result(2,2) =-f(2);
      Result(0,3) =-s.dot(eye);
      Result(1,3) =-u.dot(eye);
      Result(2,3) = f.dot(eye);

      return Result;
   }

   template <typename T>
   Eigen::Matrix<T,3,1> QuatToEuler(
      const Eigen::Quaternion<T> quat
   ) {
      T sqw = quat.w() * quat.w();
      T sqx = quat.x() * quat.x();
      T sqy = quat.y() * quat.y();
      T sqz = quat.z() * quat.z();
      T unit = sqx + sqy + sqz + sqw; // if normalised is one, otherwise is correction factor
      T test = quat.w() * quat.x() * quat.y() + quat.z();

      if (test > T(0.499)*unit) { // singularity at north pole
         return Eigen::Matrix<T,3,1>(
            T(2) * atan2(quat.x(), quat.w()),
            T(M_PI/2),
            T(0)
         );
      }
      else if (test < -0.499*unit) { // singularity at south pole
         return Eigen::Matrix<T,3,1>(
            T(-2) * atan2(quat.x(), quat.w()),
            T(-M_PI/2),
            T(0)
         );
      }
      else {
         return Eigen::Matrix<T,3,1>(
            atan2(T(2) * quat.y() * quat.w() - T(2) * quat.x() * quat.z(), sqx - sqy - sqz + sqw),
            atan2(T(2) * quat.x() * quat.w() - T(2) * quat.y() * quat.z(), -sqx + sqy - sqz + sqw),
            asin(T(2) * test/unit)
         );
      }
   }

   template <typename T>
   Eigen::Quaternion<T> EulerToQuat(
      const Eigen::Matrix<T,3,1> eulerAngles
   ) {
      T yaw = eulerAngles.x();
      T pitch = eulerAngles.y();
      T roll = eulerAngles.z();

      T c1 = cos(yaw);
      T c2 = cos(pitch);
      T c3 = cos(roll);
      T s1 = sin(yaw);
      T s2 = sin(pitch);
      T s3 = sin(roll);
      T w = T(0.5) * sqrt(T(1.0) + c1*c3 + c1*c2 - s1*s2*s3 + c2*c3);
      T w1o4 = T(1.0) / (T(4.0) * w);

      return Eigen::Quaternion<T>(
         w,
         w1o4 * (c3*s2 + c1*s2 + s1*s3*c2),
         w1o4 * (s1*c3 + s1*c2 + c1*s3*s2),
         w1o4 * (-s1*s2 + c1*s3 * c2+s3)
      );
   }



   template <typename T>
   Eigen::Matrix<T,3,3> TBN(
      const Eigen::Matrix<T,3,1> tangent,
      const Eigen::Matrix<T,3,1> bitangent,
      const Eigen::Matrix<T,3,1> normal
   ) {
      Eigen::Matrix<T,3,3> tbn;
      tbn.block(0,0,3,1) = tangent;
      tbn.block(0,1,3,1) = bitangent;
      tbn.block(0,2,3,1) = normal;
      return tbn;
   }

   template <typename T>
   Eigen::Matrix<T,3,3> InverseTBN(
      const Eigen::Matrix<T,3,1> tangent,
      const Eigen::Matrix<T,3,1> bitangent,
      const Eigen::Matrix<T,3,1> normal
   ) {
      Eigen::Matrix<T,3,3> iTBN;
      iTBN.block(0,0,1,3) = tangent.transpose();
      iTBN.block(1,0,1,3) = bitangent.transpose();
      iTBN.block(2,0,1,3) = normal.transpose();
      return iTBN;
   }

   template <typename T>
   Eigen::Matrix<T,4,1> vec3To4(
      const Eigen::Matrix<T,3,1> vec3,
      const T elem
   ) {
      return Eigen::Matrix<T,4,1>(vec3(0), vec3(1), vec3(2), elem);
   }
}

#endif // __MATRIX_MATH_H__
