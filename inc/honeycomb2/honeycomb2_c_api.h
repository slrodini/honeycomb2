#ifndef HC2_HONEYCOMB2_C_API_H
#define HC2_HONEYCOMB2_C_API_H

/**
 *
 *  \brief expose simplified API to work with Honeycomb2 in C/FORTRAN codebasis
 *  WARNING: This API is under development and is neither stable nor optimized
 *           in any way, shape or form.
 *           The use of this API is STRONGLY discouraged in a pure C++ application
 *           and the user should consider to manage Honeycomb2 explicitly.
 *  WARNING: This API is NOT exposed in the full honeycomb2.hpp header
 *           and should be manually included if needed.
 *
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

void set_up_(const char *config_name, int len);

#ifdef __cplusplus
}
#endif

#endif