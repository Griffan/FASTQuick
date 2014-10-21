# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.65])
AC_INIT([MPUTOOLS], [1.1.0], [hmkang@umich.edu])
AM_INIT_AUTOMAKE(MPUTOOLS, 1.1.0, [-Wall -Werror foreign])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_FILES([
	Makefile 
	src/Makefile
	scripts/Makefile
]) 
AC_PROG_CXX
AC_PROG_CC
AM_PROG_LIBTOOL
AC_PROG_LIBTOOL

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_HEADER_MAJOR
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_FUNC_FSEEKO
AC_TYPE_INT16_T
AC_TYPE_INT8_T
AC_TYPE_OFF_T
AC_TYPE_SSIZE_T
AC_CHECK_FUNCS([floor gethostname memset pow select sqrt strchr strstr strtol gethostbyaddr gethostbyname getpagesize socket strdup strerror])
AC_CHECK_LIB(z, gzopen, [], [AC_MSG_ERROR([libz.{so,a} was not found. Please install zlib at http://www.zlib.net/ first])])
AC_CHECK_HEADERS([limits.h stddef.h stdint.h stdlib.h string.h unistd.h arpa/inet.h fcntl.h netdb.h sys/socket.h])

AC_OUTPUT
