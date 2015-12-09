#ifndef __TEST_H__
#define __TEST_H__

#include <iostream>
#include <math.h>

#define PRINT_ERROR(file, line, exp, got) ( \
   std::cout << "FAILED TEST [" << (file) << ":" << (line) << "]: " \
   "Expected(" << (exp) << ") Got(" << (got) << ")\n" \
)

#define PRINT_SUCCESS(file, line, got) ( \
   std::cout << "Passed test [" << (file) << ":" << (line) << "]\n" \
)

inline void _boolCheck_(
   const char * file,
   int line,
   bool got,
   bool exp
) {
   if (got != exp)
      PRINT_ERROR(file, line, (exp ? "true" : "false"), (got ? "true" : "false"));
   else
      PRINT_SUCCESS(file, line, got);
}

inline void _equalityIntCheck_(
   const char * file,
   int line,
   int got,
   int exp
) {
   if (got != exp)
      PRINT_ERROR(file, line, exp, got);
   else
      PRINT_SUCCESS(file, line, got);
}

inline void _equalityFloatCheck_(
   const char * file,
   int line,
   double got,
   double exp,
   double tol
) {
   if ((isnan(got) && !isnan(exp)) || (isinf(got) && !isinf(exp)) || (got > exp && got - exp > tol) || (exp > got && exp - got > tol))
      PRINT_ERROR(file, line, exp, got);
   else
      PRINT_SUCCESS(file, line, got);
}

inline void _nanCheck_(
   const char * file,
   int line,
   double got
) {
   if(!isnan(got) && !isinf(got))
      PRINT_ERROR(file, line, "nan", got);
   else
      PRINT_SUCCESS(file, line, got);
}

#define boolCheck(got, exp) _boolCheck_(__FILE__, __LINE__, got, exp)
#define equalityIntCheck(got, exp) _equalityIntCheck_(__FILE__, __LINE__, got, exp)
#define equalityFloatCheck(got, exp, tol) _equalityFloatCheck_(__FILE__, __LINE__, got, exp, tol)
#define nanCheck(got) _nanCheck_(__FILE__, __LINE__, got)

#endif // __TEST_H__